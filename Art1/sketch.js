// Top-left → bottom-right random mosaic fill (no grid)
// Click or press R to regenerate.

let tiles = [];

const TILE_R = 5;            // tile size (radius)
const GAP = 1;                // grout thickness (space between tiles)
const PADDING = 20;

const BG = "#8f8f8f";
const TILE_COLOR = "#1e5cff";

// Density/quality knobs
const CANDIDATE_MULT = 70;     // more candidates => fuller fill (slower)
const HOLE_SAMPLES = 3500;     // more => better "no big empty space" (slower)
const HOLE_ADDS_MULT = 0.9;    // how aggressively to add tiles during hole fill
const RELAX_ITERS = 80;        // repel-only relaxation (keeps edges filled)

function setup() {
    createCanvas(1500, 800);
    pixelDensity(2);
    noLoop();
    generateMural();
}

function draw() {
    background(BG);
    noStroke();
    fill(TILE_COLOR);

    for (const t of tiles) {
        push();
        translate(t.x, t.y);
        rotate(t.rot);
        beginShape();
        for (const p of t.poly) vertex(p.x, p.y);
        endShape(CLOSE);
        pop();
    }
}

function mousePressed() { regenerate(); }
function keyPressed() { if (key === 'r' || key === 'R') regenerate(); }

function regenerate() {
    generateMural();
    redraw();
}

function generateMural() {
    tiles = [];

    const minX = PADDING + TILE_R;
    const maxX = width - PADDING - TILE_R;
    const minY = PADDING + TILE_R;
    const maxY = height - PADDING - TILE_R;

    const minDist = 2 * TILE_R + GAP;
    const area = (maxX - minX) * (maxY - minY);
    const baseCount = area / (minDist * minDist);

    // 1) Greedy scan-fill from top-left to bottom-right using random candidates
    let centers = scanFillCenters(
        minX, maxX, minY, maxY,
        minDist,
        floor(baseCount * CANDIDATE_MULT)
    );

    // 2) Hole fill: add tiles into the emptiest remaining gaps until none fits
    centers = fillHoles(
        centers,
        minX, maxX, minY, maxY,
        minDist,
        floor(baseCount * HOLE_ADDS_MULT)
    );

    // 3) Repel-only relax (prevents overlaps + tightens) without pulling to center
    centers = relaxRepelOnly(
        centers,
        minX, maxX, minY, maxY,
        minDist,
        RELAX_ITERS
    );

    // 4) Build tiles at those centers
    for (const c of centers) {
        tiles.push({
            x: c.x,
            y: c.y,
            rot: random(TWO_PI),
            poly: generateConvexTile(TILE_R)
        });
    }
}

// ---------- 1) Scan-fill: random points sorted top-left -> bottom-right ----------
function scanFillCenters(minX, maxX, minY, maxY, minDist, candidateCount) {
    const candidates = [];
    for (let i = 0; i < candidateCount; i++) {
        candidates.push({ x: random(minX, maxX), y: random(minY, maxY) });
    }

    // Sort like a “printer”: top to bottom, then left to right
    candidates.sort((a, b) => (a.y === b.y ? a.x - b.x : a.y - b.y));

    const cell = minDist;
    const minDistSq = minDist * minDist;

    const pts = [];
    const hash = new Map();

    // Force the very first tile near top-left corner
    addPoint({ x: minX, y: minY }, pts, hash, cell);

    for (const c of candidates) {
        if (farEnoughHash(c.x, c.y, pts, hash, cell, minDistSq)) {
            addPoint(c, pts, hash, cell);
        }
    }

    return pts;
}

function addPoint(p, pts, hash, cell) {
    const idx = pts.length;
    pts.push(p);

    const gx = floor(p.x / cell);
    const gy = floor(p.y / cell);
    const key = gx + "," + gy;

    if (!hash.has(key)) hash.set(key, []);
    hash.get(key).push(idx);
}

function farEnoughHash(x, y, pts, hash, cell, minDistSq) {
    const gx = floor(x / cell);
    const gy = floor(y / cell);

    for (let ix = gx - 1; ix <= gx + 1; ix++) {
        for (let iy = gy - 1; iy <= gy + 1; iy++) {
            const bucket = hash.get(ix + "," + iy);
            if (!bucket) continue;

            for (const idx of bucket) {
                const p = pts[idx];
                const dx = p.x - x, dy = p.y - y;
                if (dx * dx + dy * dy < minDistSq) return false;
            }
        }
    }
    return true;
}

// ---------- 2) Hole fill: repeatedly insert at “largest empty” until no fit ----------
function fillHoles(pts, minX, maxX, minY, maxY, minDist, maxAdds) {
    const cell = minDist;

    for (let add = 0; add < maxAdds; add++) {
        const hash = buildHash(pts, cell);

        let best = null;
        let bestD2 = -1;

        for (let s = 0; s < HOLE_SAMPLES; s++) {
            const x = random(minX, maxX);
            const y = random(minY, maxY);
            const d2 = nearestDist2(x, y, pts, hash, cell);

            if (d2 > bestD2) {
                bestD2 = d2;
                best = { x, y };
            }
        }

        if (!best) break;
        if (sqrt(bestD2) < minDist) break; // no more space for a new tile center

        pts.push(best);
    }

    return pts;
}

function buildHash(pts, cell) {
    const map = new Map();
    for (let i = 0; i < pts.length; i++) {
        const p = pts[i];
        const gx = floor(p.x / cell);
        const gy = floor(p.y / cell);
        const key = gx + "," + gy;
        if (!map.has(key)) map.set(key, []);
        map.get(key).push(i);
    }
    return map;
}

function nearestDist2(x, y, pts, hash, cell) {
    if (pts.length === 0) return Infinity;

    const cx = floor(x / cell);
    const cy = floor(y / cell);
    let best = Infinity;

    for (let r = 0; r < 18; r++) {
        for (let gx = cx - r; gx <= cx + r; gx++) {
            for (let gy = cy - r; gy <= cy + r; gy++) {
                if (abs(gx - cx) !== r && abs(gy - cy) !== r) continue;

                const bucket = hash.get(gx + "," + gy);
                if (!bucket) continue;

                for (const idx of bucket) {
                    const p = pts[idx];
                    const dx = p.x - x, dy = p.y - y;
                    const d2 = dx * dx + dy * dy;
                    if (d2 < best) best = d2;
                }
            }
        }
        if (best < Infinity && r * cell > sqrt(best) + cell) break;
    }

    return best;
}

// ---------- 3) Repel-only relax: tightens without collapsing to center ----------
function relaxRepelOnly(pts, minX, maxX, minY, maxY, minDist, iters) {
    const cell = minDist;
    const minDistSq = minDist * minDist;
    const pushK = 0.18;

    for (let it = 0; it < iters; it++) {
        const hash = buildHash(pts, cell);
        const next = pts.map(p => ({ x: p.x, y: p.y }));

        for (let i = 0; i < pts.length; i++) {
            const p = pts[i];
            let fx = 0, fy = 0;

            const gx = floor(p.x / cell);
            const gy = floor(p.y / cell);

            for (let ix = gx - 1; ix <= gx + 1; ix++) {
                for (let iy = gy - 1; iy <= gy + 1; iy++) {
                    const bucket = hash.get(ix + "," + iy);
                    if (!bucket) continue;

                    for (const j of bucket) {
                        if (j === i) continue;
                        const q = pts[j];

                        const dx = p.x - q.x;
                        const dy = p.y - q.y;
                        const d2 = dx * dx + dy * dy;
                        if (d2 === 0 || d2 >= minDistSq) continue;

                        const d = sqrt(d2);
                        const overlap = (minDist - d);
                        const ux = dx / d;
                        const uy = dy / d;

                        fx += ux * overlap * pushK;
                        fy += uy * overlap * pushK;
                    }
                }
            }

            next[i].x = constrain(p.x + fx, minX, maxX);
            next[i].y = constrain(p.y + fy, minY, maxY);
        }

        pts = next;
    }

    return pts;
}

// ---------- Tile generation (convex polygon, 3..7 sides, side ratio <= 5) ----------
function generateConvexTile(radius) {
    for (let attempt = 0; attempt < 300; attempt++) {
        const m = floor(random(10, 22));
        const pts = [];

        for (let i = 0; i < m; i++) {
            const a = random(TWO_PI);
            const r = radius * 1.2 * sqrt(random());
            pts.push({ x: r * cos(a), y: r * sin(a) });
        }

        let hull = convexHull(pts);
        if (hull.length < 3 || hull.length > 7) continue;

        hull = normalizeToRadius(hull, radius);

        if (sideRatio(hull) > 5.0) continue;
        if (minSide(hull) < radius * 0.25) continue;

        return hull;
    }

    return regularPoly(floor(random(3, 8)), radius);
}

function regularPoly(n, radius) {
    const pts = [];
    const start = random(TWO_PI);
    for (let i = 0; i < n; i++) {
        const a = start + i * TWO_PI / n;
        pts.push({ x: radius * cos(a), y: radius * sin(a) });
    }
    return pts;
}

function centroid(pts) {
    let cx = 0, cy = 0;
    for (const p of pts) { cx += p.x; cy += p.y; }
    return { x: cx / pts.length, y: cy / pts.length };
}

function normalizeToRadius(pts, radius) {
    const c = centroid(pts);
    let maxD = 0;

    const centered = pts.map(p => {
        const x = p.x - c.x;
        const y = p.y - c.y;
        maxD = max(maxD, sqrt(x * x + y * y));
        return { x, y };
    });

    const s = (maxD > 0) ? (radius / maxD) : 1;
    return centered.map(p => ({ x: p.x * s, y: p.y * s }));
}

function sideRatio(poly) {
    let mn = Infinity, mx = 0;
    for (let i = 0; i < poly.length; i++) {
        const a = poly[i];
        const b = poly[(i + 1) % poly.length];
        const d = dist(a.x, a.y, b.x, b.y);
        mn = min(mn, d);
        mx = max(mx, d);
    }
    return mx / mn;
}

function minSide(poly) {
    let mn = Infinity;
    for (let i = 0; i < poly.length; i++) {
        const a = poly[i];
        const b = poly[(i + 1) % poly.length];
        mn = min(mn, dist(a.x, a.y, b.x, b.y));
    }
    return mn;
}

// Strict convex hull (no collinear), so all interior angles < 180°
function convexHull(points) {
    if (points.length <= 1) return points.slice();

    const pts = points.slice().sort((a, b) => (a.x === b.x) ? a.y - b.y : a.x - b.x);
    const cross = (o, a, b) => (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);

    const lower = [];
    for (const p of pts) {
        while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
            lower.pop();
        }
        lower.push(p);
    }

    const upper = [];
    for (let i = pts.length - 1; i >= 0; i--) {
        const p = pts[i];
        while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
            upper.pop();
        }
        upper.push(p);
    }

    lower.pop();
    upper.pop();
    return lower.concat(upper);
}
