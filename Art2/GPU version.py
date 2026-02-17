import math
import secrets
import numpy as np
import cv2
from numba import cuda

# ============================================================
# WINDOW / RENDER
# ============================================================
W, H = 1280, 720
WINDOW_NAME = "Differential Growth (CUDA) — R reset | SPACE pause | S save | Q quit"

# Trails like p5.js background(255, BG_FADE_ALPHA)
BG_FADE_ALPHA = 2  # 0..255 (2 = long trails)

# Single-color look
LINE_GRAY = 60                 # 0..255 (smaller = darker)
LINE_ALPHA = 140 / 255.0       # like p5 stroke(0,140)
LINE_THICKNESS = 1             # 1–2

# Performance / sim rate
ITER_PER_FRAME = 5
GRID_REBUILD_EVERY = 1         # 1 = most accurate, 2 = faster

# Safety cap so randomization doesn’t produce insane CPU draw cost
MAX_TOTAL_POINTS = 45000

# ============================================================
# CUDA kernels
# ============================================================

@cuda.jit
def reset_counts(counts):
    i = cuda.grid(1)
    if i < counts.size:
        counts[i] = 0

@cuda.jit
def build_grid(posx, posy, counts, cell_points, grid_w, grid_h, cell_size, max_per_cell):
    i = cuda.grid(1)
    n = posx.size
    if i >= n:
        return

    x = posx[i]
    y = posy[i]

    cx = int(x / cell_size)
    cy = int(y / cell_size)

    if cx < 0: cx = 0
    if cy < 0: cy = 0
    if cx >= grid_w: cx = grid_w - 1
    if cy >= grid_h: cy = grid_h - 1

    cell = cx + cy * grid_w
    slot = cuda.atomic.add(counts, cell, 1)
    if slot < max_per_cell:
        cell_points[cell, slot] = i

@cuda.jit(device=True)
def _hash_u32(x):
    x ^= (x >> 16)
    x *= 0x7feb352d
    x ^= (x >> 15)
    x *= 0x846ca68b
    x ^= (x >> 16)
    return x

@cuda.jit(device=True)
def _rand2_cell(ix, iy, tick, salt):
    # deterministic pseudo-random vector in [-1, 1], cheaply normalized
    h = _hash_u32(ix * 73856093 ^ iy * 19349663 ^ tick * 83492791 ^ salt)
    a = (h & 1023) / 511.5 - 1.0
    b = ((h >> 10) & 1023) / 511.5 - 1.0
    inv = 1.0 / (abs(a) + abs(b) + 1e-6)
    return a * inv, b * inv

@cuda.jit(fastmath=True)
def update_step(posx_in, posy_in, velx_in, vely_in,
                posx_out, posy_out, velx_out, vely_out,
                prev_idx, next_idx,
                counts, cell_points,
                grid_w, grid_h,
                cell_size, max_per_cell,
                rep_r, rep_strength,
                smooth_strength, step_size, damping,
                drift_cell_size, drift_strength, drift_rate, t, drift_salt,
                edge_margin, edge_push,
                W, H):
    i = cuda.grid(1)
    n = posx_in.size
    if i >= n:
        return

    x = posx_in[i]
    y = posy_in[i]
    vx = velx_in[i]
    vy = vely_in[i]

    # --- smoothing along loop ---
    pi = prev_idx[i]
    ni = next_idx[i]
    px = posx_in[pi]
    py = posy_in[pi]
    nx = posx_in[ni]
    ny = posy_in[ni]

    avgx = 0.5 * (px + nx)
    avgy = 0.5 * (py + ny)

    fx = (avgx - x) * smooth_strength
    fy = (avgy - y) * smooth_strength

    # --- repulsion grid cell ---
    cx = int(x / cell_size)
    cy = int(y / cell_size)
    if cx < 0: cx = 0
    if cy < 0: cy = 0
    if cx >= grid_w: cx = grid_w - 1
    if cy >= grid_h: cy = grid_h - 1

    rr2 = rep_r * rep_r
    repx = 0.0
    repy = 0.0

    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            ncx = cx + dx
            ncy = cy + dy
            if ncx < 0 or ncy < 0 or ncx >= grid_w or ncy >= grid_h:
                continue

            cell = ncx + ncy * grid_w
            cnt = counts[cell]
            if cnt > max_per_cell:
                cnt = max_per_cell

            for k in range(cnt):
                j = cell_points[cell, k]
                if j == i or j == pi or j == ni:
                    continue

                qx = posx_in[j]
                qy = posy_in[j]
                dxp = x - qx
                dyp = y - qy
                d2 = dxp * dxp + dyp * dyp

                if d2 > 0.0 and d2 < rr2:
                    d = math.sqrt(d2)
                    w = 1.0 - d / rep_r
                    invd = 1.0 / d
                    repx += (dxp * invd) * (w * w)
                    repy += (dyp * invd) * (w * w)

    fx += repx * rep_strength
    fy += repy * rep_strength

    # --- drift field (THIS is what was constant before; now drift_salt changes every reset) ---
    dcx = int(x / drift_cell_size)
    dcy = int(y / drift_cell_size)

    tt = t * drift_rate
    tick0 = int(tt)
    a = tt - tick0

    ux0, uy0 = _rand2_cell(dcx, dcy, tick0, drift_salt)
    ux1, uy1 = _rand2_cell(dcx, dcy, tick0 + 1, drift_salt)

    ux = ux0 + (ux1 - ux0) * a
    uy = uy0 + (uy1 - uy0) * a

    fx += ux * drift_strength
    fy += uy * drift_strength

    # --- edge push ---
    if x < edge_margin:
        fx += (edge_margin - x) * edge_push
    if x > W - edge_margin:
        fx -= (x - (W - edge_margin)) * edge_push
    if y < edge_margin:
        fy += (edge_margin - y) * edge_push
    if y > H - edge_margin:
        fy -= (y - (H - edge_margin)) * edge_push

    # --- integrate ---
    vx = (vx + fx * step_size) * damping
    vy = (vy + fy * step_size) * damping
    x = x + vx
    y = y + vy

    posx_out[i] = x
    posy_out[i] = y
    velx_out[i] = vx
    vely_out[i] = vy

# ============================================================
# Random scene generation
# ============================================================

def pick_style(seed: int):
    """Randomize the whole world, but keep it in safe/pretty ranges."""
    rng = np.random.default_rng(seed)

    # Composition
    target_loops = int(rng.integers(200, 5000))
    skip_prob = float(rng.uniform(0.05, 0.35))
    mode = "grid" if rng.random() < 0.70 else "random"

    rmin = float(rng.uniform(30, 90))
    rmax = float(rng.uniform(max(rmin + 20, 80), 180))

    pmin = int(rng.integers(140, 320))
    pmax = int(rng.integers(max(pmin + 80, 320), 900))

    # Dynamics (keep repulsion radius stable-ish for a consistent “lace” look)
    rep_r = float(rng.uniform(16.0, 26.0))
    rep_strength = float(rng.uniform(0.40, 0.75))
    smooth_strength = float(rng.uniform(0.18, 0.34))
    step_size = float(rng.uniform(0.45, 0.75))
    damping = float(rng.uniform(0.82, 0.90))

    # Drift (this is the big visual difference maker)
    base = 1.0 / 0.0022  # ~454, similar scale to your p5 NOISE_SCALE
    drift_cell_size = float(base * rng.uniform(0.55, 1.90))
    drift_strength = float(rng.uniform(0.30, 0.95))
    drift_rate = float(rng.uniform(0.60, 2.40))
    noise_time_speed = float(rng.uniform(0.006, 0.020))

    # Fresh drift seed (OS randomness)
    drift_salt = np.int32(secrets.randbits(31))

    return dict(
        mode=mode,
        target_loops=target_loops,
        skip_prob=skip_prob,
        rmin=rmin,
        rmax=rmax,
        pmin=pmin,
        pmax=pmax,
        rep_r=np.float32(rep_r),
        rep_strength=np.float32(rep_strength),
        smooth_strength=np.float32(smooth_strength),
        step_size=np.float32(step_size),
        damping=np.float32(damping),
        drift_cell_size=np.float32(drift_cell_size),
        drift_strength=np.float32(drift_strength),
        drift_rate=np.float32(drift_rate),
        noise_time_speed=np.float32(noise_time_speed),
        drift_salt=drift_salt,
        edge_margin=np.float32(40.0),
        edge_push=np.float32(0.020),
        max_per_cell=int(96),
    )

def seed_loops(rng: np.random.Generator, style):
    loops = []

    if style["mode"] == "grid":
        area = W * H
        spacing = math.sqrt(area / max(1, style["target_loops"]))
        jitter = spacing * 0.35

        for y in np.arange(spacing * 0.5, H, spacing):
            for x in np.arange(spacing * 0.5, W, spacing):
                if rng.random() < style["skip_prob"]:
                    continue
                cx = float(x + rng.uniform(-jitter, jitter))
                cy = float(y + rng.uniform(-jitter, jitter))
                r = float(rng.uniform(style["rmin"], style["rmax"]))
                n = int(rng.integers(style["pmin"], style["pmax"] + 1))
                loops.append((cx, cy, r, n))

    else:
        # Fully random placement
        n_loops = style["target_loops"] + int(rng.integers(-2, 7))
        n_loops = max(1, n_loops)
        for _ in range(n_loops):
            cx = float(rng.uniform(0, W))
            cy = float(rng.uniform(0, H))
            r = float(rng.uniform(style["rmin"], style["rmax"]))
            n = int(rng.integers(style["pmin"], style["pmax"] + 1))
            loops.append((cx, cy, r, n))

    if not loops:
        loops = [(W * 0.5, H * 0.5, 120.0, style["pmax"])]

    # Cap total points (for live rendering speed)
    total = sum(n for (_, _, _, n) in loops)
    if total > MAX_TOTAL_POINTS:
        f = MAX_TOTAL_POINTS / float(total)
        new_loops = []
        for (cx, cy, r, n) in loops:
            n2 = max(50, int(n * f))
            new_loops.append((cx, cy, r, n2))
        loops = new_loops

    total = sum(n for (_, _, _, n) in loops)

    posx = np.empty(total, dtype=np.float32)
    posy = np.empty(total, dtype=np.float32)
    prev_idx = np.empty(total, dtype=np.int32)
    next_idx = np.empty(total, dtype=np.int32)

    loop_starts = []
    loop_lens = []

    base = 0
    for (cx, cy, r, n) in loops:
        loop_starts.append(base)
        loop_lens.append(n)

        # extra randomness: ellipse + rotation
        rx = r * float(rng.uniform(0.65, 1.35))
        ry = r * float(rng.uniform(0.65, 1.35))
        phi = float(rng.uniform(0.0, 2.0 * math.pi))
        cph = math.cos(phi)
        sph = math.sin(phi)

        # random phase so it doesn't always start at angle 0
        phase = float(rng.uniform(0.0, 2.0 * math.pi))

        for k in range(n):
            t = (k / n) * (2.0 * math.pi) + phase
            wob = float(rng.uniform(0.85, 1.15))

            ux = math.cos(t) * rx * wob
            uy = math.sin(t) * ry * wob

            # rotate ellipse
            x = cx + (ux * cph - uy * sph)
            y = cy + (ux * sph + uy * cph)

            posx[base + k] = x
            posy[base + k] = y

        for k in range(n):
            prev_idx[base + k] = base + (k - 1) % n
            next_idx[base + k] = base + (k + 1) % n

        base += n

    return posx, posy, prev_idx, next_idx, np.array(loop_starts, np.int32), np.array(loop_lens, np.int32)

# ============================================================
# Simulation wrapper
# ============================================================

class CudaSim:
    def __init__(self):
        self.reset()

    def reset(self):
        # Strong randomness (not time-based), so every reset is different
        self.seed = secrets.randbits(32)
        rng = np.random.default_rng(self.seed)

        self.style = pick_style(self.seed)

        posx, posy, prev_i, next_i, self.loop_starts, self.loop_lens = seed_loops(rng, self.style)
        self.n = int(posx.size)

        # Grid uses repulsion radius
        self.cell_size = float(self.style["rep_r"])  # cell size ~= repulsion radius
        self.grid_w = int(math.ceil(W / self.cell_size))
        self.grid_h = int(math.ceil(H / self.cell_size))
        self.n_cells = int(self.grid_w * self.grid_h)

        # Device buffers (double-buffer)
        self.d_posx_a = cuda.to_device(posx)
        self.d_posy_a = cuda.to_device(posy)
        self.d_velx_a = cuda.to_device(np.zeros_like(posx))
        self.d_vely_a = cuda.to_device(np.zeros_like(posy))

        self.d_posx_b = cuda.device_array_like(self.d_posx_a)
        self.d_posy_b = cuda.device_array_like(self.d_posy_a)
        self.d_velx_b = cuda.device_array_like(self.d_velx_a)
        self.d_vely_b = cuda.device_array_like(self.d_vely_a)

        self.d_prev = cuda.to_device(prev_i)
        self.d_next = cuda.to_device(next_i)

        self.d_counts = cuda.device_array(self.n_cells, dtype=np.int32)
        self.d_cell_points = cuda.device_array((self.n_cells, self.style["max_per_cell"]), dtype=np.int32)

        # Pinned host arrays for fast copies
        self.h_posx = cuda.pinned_array(self.n, dtype=np.float32)
        self.h_posy = cuda.pinned_array(self.n, dtype=np.float32)

        self.pts_int = np.empty((self.n, 2), dtype=np.int32)

        self.step_count = 0

        print(
            f"[reset] seed={self.seed} mode={self.style['mode']} "
            f"loops={len(self.loop_lens)} points={self.n} "
            f"rep_r={float(self.style['rep_r']):.2f} "
            f"rep={float(self.style['rep_strength']):.2f} "
            f"smooth={float(self.style['smooth_strength']):.2f} "
            f"drift_cell={float(self.style['drift_cell_size']):.0f} "
            f"drift_str={float(self.style['drift_strength']):.2f} "
            f"drift_rate={float(self.style['drift_rate']):.2f}"
        )

    def step(self, iters: int):
        threads_points = 256
        blocks_points = (self.n + threads_points - 1) // threads_points

        threads_cells = 128
        blocks_cells = (self.n_cells + threads_cells - 1) // threads_cells

        # Pull scalars once (keeps call clean)
        rep_r = self.style["rep_r"]
        rep_strength = self.style["rep_strength"]
        smooth_strength = self.style["smooth_strength"]
        step_size = self.style["step_size"]
        damping = self.style["damping"]
        drift_cell_size = self.style["drift_cell_size"]
        drift_strength = self.style["drift_strength"]
        drift_rate = self.style["drift_rate"]
        noise_time_speed = self.style["noise_time_speed"]
        drift_salt = self.style["drift_salt"]
        edge_margin = self.style["edge_margin"]
        edge_push = self.style["edge_push"]

        for _ in range(iters):
            self.step_count += 1
            t = np.float32(self.step_count) * noise_time_speed

            if (self.step_count % GRID_REBUILD_EVERY) == 0:
                reset_counts[blocks_cells, threads_cells](self.d_counts)
                build_grid[blocks_points, threads_points](
                    self.d_posx_a, self.d_posy_a,
                    self.d_counts, self.d_cell_points,
                    self.grid_w, self.grid_h,
                    np.float32(self.cell_size), int(self.style["max_per_cell"])
                )

            update_step[blocks_points, threads_points](
                self.d_posx_a, self.d_posy_a, self.d_velx_a, self.d_vely_a,
                self.d_posx_b, self.d_posy_b, self.d_velx_b, self.d_vely_b,
                self.d_prev, self.d_next,
                self.d_counts, self.d_cell_points,
                self.grid_w, self.grid_h,
                np.float32(self.cell_size), int(self.style["max_per_cell"]),
                rep_r, rep_strength,
                smooth_strength, step_size, damping,
                drift_cell_size, drift_strength, drift_rate, t, drift_salt,
                edge_margin, edge_push,
                np.float32(W), np.float32(H)
            )

            # swap
            self.d_posx_a, self.d_posx_b = self.d_posx_b, self.d_posx_a
            self.d_posy_a, self.d_posy_b = self.d_posy_b, self.d_posy_a
            self.d_velx_a, self.d_velx_b = self.d_velx_b, self.d_velx_a
            self.d_vely_a, self.d_vely_b = self.d_vely_b, self.d_vely_a

    def fetch_points_int(self):
        self.d_posx_a.copy_to_host(self.h_posx)
        self.d_posy_a.copy_to_host(self.h_posy)
        self.pts_int[:, 0] = self.h_posx
        self.pts_int[:, 1] = self.h_posy
        return self.pts_int

# ============================================================
# Live loop
# ============================================================

def main():
    if not cuda.is_available():
        raise RuntimeError("CUDA not available (numba cuda.is_available() == False).")

    sim = CudaSim()

    # Persistent canvas (float for blending)
    canvas = np.full((H, W, 3), 255.0, dtype=np.float32)

    fade = float(BG_FADE_ALPHA) / 255.0

    mask = np.zeros((H, W), dtype=np.uint8)
    alpha_map = np.zeros((H, W), dtype=np.float32)

    paused = False

    cv2.namedWindow(WINDOW_NAME, cv2.WINDOW_NORMAL)
    cv2.resizeWindow(WINDOW_NAME, W, H)

    while True:
        # Fade to white
        canvas += (255.0 - canvas) * fade

        if not paused:
            sim.step(ITER_PER_FRAME)

        pts_int = sim.fetch_points_int()

        # Draw all loops into a mask
        mask.fill(0)
        for start, length in zip(sim.loop_starts, sim.loop_lens):
            poly = pts_int[start:start + length].reshape(-1, 1, 2)
            cv2.polylines(mask, [poly], True, 255, thickness=LINE_THICKNESS, lineType=cv2.LINE_AA)

        # Alpha blend gray lines onto canvas
        alpha_map[:] = mask
        alpha_map *= (LINE_ALPHA / 255.0)

        a3 = alpha_map[:, :, None]
        canvas *= (1.0 - a3)
        canvas += float(LINE_GRAY) * a3

        frame_u8 = np.clip(canvas, 0, 255).astype(np.uint8)
        cv2.imshow(WINDOW_NAME, frame_u8)

        key = cv2.waitKey(1) & 0xFF
        if key in (27, ord('q'), ord('Q')):
            break
        elif key == ord(' '):
            paused = not paused
        elif key in (ord('r'), ord('R')):
            sim.reset()
            canvas[:] = 255.0
        elif key in (ord('s'), ord('S')):
            fn = f"differential_live_seed{sim.seed}.png"
            cv2.imwrite(fn, frame_u8)
            print("saved:", fn)

    cv2.destroyAllWindows()

if __name__ == "__main__":
    main()
