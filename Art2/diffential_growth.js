// Differential Line Growth — Full Screen + Always Random
// Multiple closed loops seeded across the canvas.
// Repulsion across ALL loops + smoothing within each loop.
// Randomness: randomSeed/noiseSeed on reset + noise drift term.
//
// Controls:
//   R = new random scene
//   SPACE = pause/resume
//   S = save canvas

let lines = [];
let paused = false;

// ===== Coverage knobs =====
const TARGET_LOOPS = 1;          // increase to cover screen more
const LOOP_RADIUS_MIN = 50;
const LOOP_RADIUS_MAX = 120;
const LOOP_POINTS_MIN = 70;
const LOOP_POINTS_MAX = 140;

const ITER_PER_FRAME = 2;         // increase = faster evolution (heavier CPU)
const BG_FADE_ALPHA = 20;         // 0..255 (higher = more fade each frame)

// ===== Growth knobs =====
const REPULSION_RADIUS = 20;
const REPULSION_STRENGTH = 0.55;

const SMOOTHING_STRENGTH = 0.26;
const STEP_SIZE = 0.60;
const DAMPING = 0.86;

const MAX_SEG_LEN = 10.5;

// ===== Random drift (makes it "always random") =====
const NOISE_SCALE = 0.0022;
const NOISE_STRENGTH = 0.55;
const NOISE_TIME_SPEED = 0.010;

// ===== Boundary push (keeps everything inside) =====
const EDGE_MARGIN = 40;
const EDGE_PUSH = 0.020;

// Spatial hash cell size (≈ repulsion radius)
const CELL_SIZE = REPULSION_RADIUS;

// -------- Data structures --------
class Pt {
    constructor(x, y) {
        this.p = createVector(x, y);
        this.v = createVector(0, 0);
    }
}

function setup() {
    createCanvas(windowWidth, windowHeight);
    pixelDensity(1);
    resetScene(true);
}

function windowResized() {
    resizeCanvas(windowWidth, windowHeight);
    resetScene(true);
}

function resetScene(clearBg = false) {
    // Different every time (and different between runs)
    const seed = (Date.now() % 1000000000) | 0;
    randomSeed(seed);
    noiseSeed(seed);

    lines = seedLoopsAcrossScreen();

    if (clearBg) background(255);
}

function seedLoopsAcrossScreen() {
    const out = [];
    // Estimate spacing so we get ~TARGET_LOOPS across the area
    const area = width * height;
    const spacing = sqrt(area / TARGET_LOOPS);
    const jitter = spacing * 0.35;

    for (let y = spacing * 0.5; y < height; y += spacing) {
        for (let x = spacing * 0.5; x < width; x += spacing) {
            // Randomly skip some cells to keep it organic
            if (random() < 0.25) continue;

            const cx = x + random(-jitter, jitter);
            const cy = y + random(-jitter, jitter);

            const r = random(LOOP_RADIUS_MIN, LOOP_RADIUS_MAX);
            const n = floor(random(LOOP_POINTS_MIN, LOOP_POINTS_MAX));

            out.push(makeLoop(cx, cy, r, n));
        }
    }
    return out;
}

function makeLoop(cx, cy, r, nPts) {
    const pts = [];
    for (let i = 0; i < nPts; i++) {
        const t = (i / nPts) * TWO_PI;
        // small random wobble so loops aren't perfect circles
        const wob = random(0.85, 1.15);
        const x = cx + cos(t) * r * wob;
        const y = cy + sin(t) * r * wob;
        pts.push(new Pt(x, y));
    }
    return { pts, closed: true };
}

function draw() {
    background(255, BG_FADE_ALPHA);

    if (!paused) {
        for (let k = 0; k < ITER_PER_FRAME; k++) {
            stepAllLines();
            subdivideAllLines();
        }
    }

    renderAllLines();
}

function stepAllLines() {
    // Build a global spatial hash: cell -> list of {li, pi}
    const grid = new Map();

    for (let li = 0; li < lines.length; li++) {
        const pts = lines[li].pts;
        for (let pi = 0; pi < pts.length; pi++) {
            const p = pts[pi].p;
            const cx = floor(p.x / CELL_SIZE);
            const cy = floor(p.y / CELL_SIZE);
            const key = cx + "," + cy;
            if (!grid.has(key)) grid.set(key, []);
            grid.get(key).push({ li, pi });
        }
    }

    const t = frameCount * NOISE_TIME_SPEED;

    // Compute forces
    for (let li = 0; li < lines.length; li++) {
        const L = lines[li];
        const pts = L.pts;
        const n = pts.length;

        for (let pi = 0; pi < n; pi++) {
            const pt = pts[pi];
            const p = pt.p;

            // Neighbors (within this loop)
            const prev = pts[(pi - 1 + n) % n].p;
            const next = pts[(pi + 1) % n].p;

            // 1) Smoothing (keeps line coherent)
            const neighborAvg = p5.Vector.add(prev, next).mult(0.5);
            const smoothForce = p5.Vector.sub(neighborAvg, p).mult(SMOOTHING_STRENGTH);

            // 2) Repulsion (across ALL loops)
            const repForce = createVector(0, 0);

            const gcx = floor(p.x / CELL_SIZE);
            const gcy = floor(p.y / CELL_SIZE);

            for (let dx = -1; dx <= 1; dx++) {
                for (let dy = -1; dy <= 1; dy++) {
                    const key = (gcx + dx) + "," + (gcy + dy);
                    const bucket = grid.get(key);
                    if (!bucket) continue;

                    for (const ref of bucket) {
                        const li2 = ref.li;
                        const pi2 = ref.pi;

                        // skip self
                        if (li2 === li && pi2 === pi) continue;

                        // skip immediate neighbors in SAME loop (so it doesn't explode)
                        if (li2 === li) {
                            const prevIdx = (pi - 1 + n) % n;
                            const nextIdx = (pi + 1) % n;
                            if (pi2 === prevIdx || pi2 === nextIdx) continue;
                        }

                        const q = lines[li2].pts[pi2].p;
                        const dxp = p.x - q.x;
                        const dyp = p.y - q.y;
                        const d2 = dxp * dxp + dyp * dyp;
                        const r2 = REPULSION_RADIUS * REPULSION_RADIUS;

                        if (d2 > 0 && d2 < r2) {
                            const d = sqrt(d2);
                            const w = 1 - d / REPULSION_RADIUS; // 0..1
                            const ux = dxp / d;
                            const uy = dyp / d;
                            // quadratic falloff looks nice
                            repForce.x += ux * (w * w);
                            repForce.y += uy * (w * w);
                        }
                    }
                }
            }
            repForce.mult(REPULSION_STRENGTH);

            // 3) Noise drift (this makes growth "always random")
            // Use noise to create a small direction field
            const ang = noise(p.x * NOISE_SCALE, p.y * NOISE_SCALE, t) * TWO_PI * 2.0;
            const noiseForce = createVector(cos(ang), sin(ang)).mult(NOISE_STRENGTH);

            // 4) Boundary push (softly keep inside canvas)
            const edgeForce = createVector(0, 0);
            if (p.x < EDGE_MARGIN) edgeForce.x += (EDGE_MARGIN - p.x);
            if (p.x > width - EDGE_MARGIN) edgeForce.x -= (p.x - (width - EDGE_MARGIN));
            if (p.y < EDGE_MARGIN) edgeForce.y += (EDGE_MARGIN - p.y);
            if (p.y > height - EDGE_MARGIN) edgeForce.y -= (p.y - (height - EDGE_MARGIN));
            edgeForce.mult(EDGE_PUSH);

            // Combine
            const force = createVector(0, 0);
            force.add(smoothForce);
            force.add(repForce);
            force.add(noiseForce);
            force.add(edgeForce);

            // Integrate
            pt.v.add(force.mult(STEP_SIZE));
            pt.v.mult(DAMPING);
        }
    }

    // Update positions
    for (const L of lines) {
        for (const pt of L.pts) {
            pt.p.add(pt.v);
        }
    }
}

function subdivideAllLines() {
    for (let li = 0; li < lines.length; li++) {
        const pts = lines[li].pts;
        const n = pts.length;
        const out = [];

        for (let i = 0; i < n; i++) {
            const a = pts[i];
            const b = pts[(i + 1) % n];

            out.push(a);

            const dx = a.p.x - b.p.x;
            const dy = a.p.y - b.p.y;
            const d = sqrt(dx * dx + dy * dy);

            if (d > MAX_SEG_LEN) {
                const mid = p5.Vector.add(a.p, b.p).mult(0.5);
                const m = new Pt(mid.x, mid.y);
                m.v = p5.Vector.add(a.v, b.v).mult(0.5);
                out.push(m);
            }
        }

        lines[li].pts = out;
    }
}

function renderAllLines() {
    noFill();
    stroke(0, 140);
    strokeWeight(1);

    for (const L of lines) {
        beginShape();
        for (const pt of L.pts) vertex(pt.p.x, pt.p.y);
        endShape(CLOSE);
    }
}

function keyPressed() {
    if (key === ' ') paused = !paused;
    if (key === 'r' || key === 'R') resetScene(true);
    if (key === 's' || key === 'S') saveCanvas('differential-growth-fullscreen', 'png');
}
