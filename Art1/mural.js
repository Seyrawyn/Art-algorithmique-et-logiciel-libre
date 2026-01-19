// FAST: Painting -> tiny tiles colored from painting pixels
// Click or press R to regenerate.

let pg;
let tiles = [];
let shapes = [];

// ===== Your tile scale =====
const TILE_R = 6;
const GAP = 0.01;
const PADDING = 20;
const GROUT = "#8f8f8f";

// ===== Speed knobs =====
const SHAPE_CACHE = 250; // precomputed tile shapes (reused)
const POISSON_K = 18;    // attempts per active point (higher = denser, slower)

// ===== Painting knobs =====
const STROKE_MIN = 10;
const STROKE_MAX_EXCL = 20; // 6..14
const FINAL_BLUR = 5;

function setup() {
    createCanvas(1500, 800);
    pixelDensity(2);

    pg = createGraphics(width, height);
    pg.pixelDensity(2);

    noLoop();
    regenerate();
}

function draw() {
    background(GROUT);
    noStroke();

    // IMPORTANT: ensure RGB fill for (r,g,b)
    colorMode(RGB, 255);

    for (let i = 0; i < tiles.length; i++) {
        const t = tiles[i];
        fill(t.r, t.g, t.b, 255);

        const poly = shapes[t.si];

        push();
        translate(t.x, t.y);
        rotate(t.rot);

        beginShape();
        for (let k = 0; k < poly.length; k++) vertex(poly[k].x, poly[k].y);
        endShape(CLOSE);

        pop();
    }
}

function mousePressed() { regenerate(); redraw(); }
function keyPressed() { if (key === 'r' || key === 'R') { regenerate(); redraw(); } }

function regenerate() {
    // 1) Paint into pg (HSB correctly)
    generatePainting(pg);

    // 2) Precompute tile shapes once per regen
    shapes = [];
    for (let i = 0; i < SHAPE_CACHE; i++) shapes.push(generateConvexTile(TILE_R));

    // 3) Generate centers quickly (Poisson-disc) seeded from top-left
    const minX = PADDING + TILE_R;
    const maxX = width - PADDING - TILE_R;
    const minY = PADDING + TILE_R;
    const maxY = height - PADDING - TILE_R;
    const minDist = 2 * TILE_R + GAP;

    const centers = poissonDiscFromTopLeft(minX, maxX, minY, maxY, minDist, POISSON_K);

    // 4) Sample color for each center from pg.pixels (FAST)
    pg.loadPixels();
    tiles = [];

    for (let i = 0; i < centers.length; i++) {
        const c = centers[i];
        const { r, g, b } = samplePixelFast(pg, c.x, c.y);

        tiles.push({
            x: c.x,
            y: c.y,
            rot: random(TWO_PI),
            si: (Math.random() * shapes.length) | 0,
            r, g, b
        });
    }
}

// =====================================================
// A) FAST PIXEL SAMPLING
// =====================================================
function samplePixelFast(g, x, y) {
    const w = g.width, h = g.height;
    const xi = x < 0 ? 0 : x >= w ? w - 1 : (x | 0);
    const yi = y < 0 ? 0 : y >= h ? h - 1 : (y | 0);
    const idx = 4 * (yi * w + xi);
    const p = g.pixels;
    return { r: p[idx], g: p[idx + 1], b: p[idx + 2] };
}

// =====================================================
// B) POISSON-DISC SAMPLING (seeded top-left)
// =====================================================
function poissonDiscFromTopLeft(minX, maxX, minY, maxY, r, k) {
    const w = maxX - minX;
    const h = maxY - minY;

    const cell = r / Math.sqrt(2);
    const cols = Math.floor(w / cell) + 1;
    const rows = Math.floor(h / cell) + 1;

    const grid = new Int32Array(cols * rows);
    for (let i = 0; i < grid.length; i++) grid[i] = -1;

    const pts = [];
    const active = [];

    function gridIndex(x, y) {
        const gx = Math.floor((x - minX) / cell);
        const gy = Math.floor((y - minY) / cell);
        return { gx, gy, gi: gx + gy * cols };
    }

    function inBounds(x, y) {
        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }

    function farEnough(x, y) {
        const { gx, gy } = gridIndex(x, y);
        const r2 = r * r;

        for (let yy = gy - 2; yy <= gy + 2; yy++) {
            if (yy < 0 || yy >= rows) continue;
            for (let xx = gx - 2; xx <= gx + 2; xx++) {
                if (xx < 0 || xx >= cols) continue;
                const idx = grid[xx + yy * cols];
                if (idx === -1) continue;
                const p = pts[idx];
                const dx = p.x - x;
                const dy = p.y - y;
                if (dx * dx + dy * dy < r2) return false;
            }
        }
        return true;
    }

    // seed at top-left
    const seed = { x: minX, y: minY };
    pts.push(seed);
    active.push(0);
    grid[gridIndex(seed.x, seed.y).gi] = 0;

    while (active.length > 0) {
        const aIndex = (Math.random() * active.length) | 0;
        const pIndex = active[aIndex];
        const p = pts[pIndex];

        let found = false;

        for (let i = 0; i < k; i++) {
            const ang = Math.random() * Math.PI * 2;
            const mag = r * (1 + Math.random()); // [r, 2r]
            const x = p.x + Math.cos(ang) * mag;
            const y = p.y + Math.sin(ang) * mag;

            if (!inBounds(x, y)) continue;
            if (!farEnough(x, y)) continue;

            const np = { x, y };
            pts.push(np);
            const ni = pts.length - 1;
            active.push(ni);

            grid[gridIndex(x, y).gi] = ni;

            found = true;
            break;
        }

        if (!found) {
            active[aIndex] = active[active.length - 1];
            active.pop();
        }
    }

    // optional: display order top-left -> bottom-right
    pts.sort((a, b) => (a.y === b.y ? a.x - b.x : a.y - b.y));
    return pts;
}

// =====================================================
// C) PAINTING GENERATOR (HSB fixed!)
// =====================================================
function generatePainting(g) {
    // Set BOTH pg and global to HSB while we create colors
    g.push();
    push();

    g.colorMode(HSB, 360, 100, 100, 1);
    colorMode(HSB, 360, 100, 100, 1);

    g.clear();

    const bg = pickBackgroundColor();
    g.background(bg);

    addSoftGrain(g, 0.03);

    const n = floor(random(STROKE_MIN, STROKE_MAX_EXCL));
    const minSide = min(width, height);
    const diag = Math.sqrt(width * width + height * height);

    for (let i = 0; i < n; i++) {
        const targetLen = random(0.5 * minSide, diag);
        const strokeCol = pickHarmoniousStroke(bg);
        drawFlowStroke(g, strokeCol, targetLen);
    }

    g.filter(BLUR, FINAL_BLUR);
    addVignette(g, bg);

    pop();
    g.pop();
}

function pickBackgroundColor() {
    const h = random(0, 360);
    const s = random(12, 35);
    const b = random(70, 92);
    return color(h, s, b, 1); // now interpreted as HSB correctly
}

function pickHarmoniousStroke(bg) {
    const h0 = hue(bg);
    const s0 = saturation(bg);
    const b0 = brightness(bg);

    const scheme = random();
    let dh;
    if (scheme < 0.55) dh = random(-25, 25);
    else if (scheme < 0.85) dh = random([150, -150]) + random(-18, 18);
    else dh = random([60, -60]) + random(-15, 15);

    const h = (h0 + dh + 360) % 360;
    const s = constrain(s0 + random(20, 55), 8, 85);
    const b = constrain(b0 + random(-18, 12), 20, 98);

    return color(h, s, b, random(0.12, 0.22));
}

function drawFlowStroke(g, col, targetLen) {
    g.push();
    g.noFill();

    const baseW = random(25, 50);

    let x = random(width);
    let y = random(height);

    const noiseScale = random(0.0012, 0.0035);
    const step = random(6, 14);
    const maxSteps = Math.floor(targetLen / step);

    const pts = [{ x, y }];

    let len = 0;
    for (let i = 0; i < maxSteps; i++) {
        const a = noise(x * noiseScale, y * noiseScale) * TWO_PI * 3.0;
        x = constrain(x + cos(a) * step, 0, width);
        y = constrain(y + sin(a) * step, 0, height);
        pts.push({ x, y });
        len += step;
        if (len >= targetLen) break;
    }

    // BIG FADE HALO
    const haloLayers = 25;
    const haloSpread = 3.0;
    for (let h = 0; h < haloLayers; h++) {
        const t = h / (haloLayers - 1);
        const w = baseW * (1.0 + t * haloSpread);
        const a = 0.04 * Math.pow(1.0 - t, 1.9);
        g.stroke(color(hue(col), saturation(col), brightness(col), a));
        g.strokeWeight(w);
        g.beginShape();
        for (let k = 0; k < pts.length; k++) g.curveVertex(pts[k].x, pts[k].y);
        g.endShape();
    }

    // Main layered stroke
    const passes = Math.floor(random(55, 90));
    for (let p = 0; p < passes; p++) {
        const jitter = random(0.7, 2.4);
        const w = baseW * (0.55 + 0.75 * noise(p * 0.2)) * random(0.85, 1.1);

        g.strokeWeight(w);
        g.stroke(color(hue(col), saturation(col), brightness(col), random(0.03, 0.085)));

        g.beginShape();
        for (let k = 0; k < pts.length; k++) {
            g.curveVertex(
                pts[k].x + random(-jitter, jitter),
                pts[k].y + random(-jitter, jitter)
            );
        }
        g.endShape();
    }

    // Highlight
    if (random() < 0.65) {
        const hcol = color(
            (hue(col) + random(-10, 10) + 360) % 360,
            constrain(saturation(col) * random(0.6, 0.9), 0, 100),
            constrain(brightness(col) + random(10, 25), 0, 100),
            random(0.06, 0.12)
        );
        g.stroke(hcol);
        g.strokeWeight(baseW * random(0.12, 0.22));
        g.beginShape();
        for (let k = 0; k < pts.length; k++) g.curveVertex(pts[k].x, pts[k].y);
        g.endShape();
    }

    g.pop();
}

function addSoftGrain(g, amount) {
    g.push();
    g.noStroke();
    const n = Math.floor(width * height * 0.0012);
    for (let i = 0; i < n; i++) {
        g.fill(0, 0, random(0, 100), random(0, amount));
        g.circle(random(width), random(height), random(0.5, 1.6));
    }
    g.pop();
}

function addVignette(g, bg) {
    g.push();
    g.noStroke();
    const vcol = color(hue(bg), saturation(bg), max(0, brightness(bg) - 18), 1);
    const layers = 18;

    for (let i = 0; i < layers; i++) {
        const t = i / (layers - 1);
        const a = 0.06 * t;
        g.fill(hue(vcol), saturation(vcol), brightness(vcol), a);
        g.ellipse(width / 2, height / 2, width * (1.05 + t * 0.55), height * (1.05 + t * 0.55));
    }
    g.pop();
}

// =====================================================
// D) TILE SHAPES (cached)
// =====================================================
function generateConvexTile(radius) {
    for (let attempt = 0; attempt < 250; attempt++) {
        const m = Math.floor(random(10, 22));
        const pts = [];

        for (let i = 0; i < m; i++) {
            const a = random(TWO_PI);
            const r = radius * 1.2 * Math.sqrt(random());
            pts.push({ x: r * cos(a), y: r * sin(a) });
        }

        let hull = convexHull(pts);
        if (hull.length < 3 || hull.length > 7) continue;

        hull = normalizeToRadius(hull, radius);
        if (sideRatio(hull) > 5.0) continue;
        if (minSide(hull) < radius * 0.25) continue;

        return hull;
    }

    return regularPoly(Math.floor(random(3, 8)), radius);
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
    for (let i = 0; i < pts.length; i++) { cx += pts[i].x; cy += pts[i].y; }
    return { x: cx / pts.length, y: cy / pts.length };
}

function normalizeToRadius(pts, radius) {
    const c = centroid(pts);
    let maxD = 0;

    const centered = pts.map(p => {
        const x = p.x - c.x;
        const y = p.y - c.y;
        const d = Math.sqrt(x * x + y * y);
        if (d > maxD) maxD = d;
        return { x, y };
    });

    const s = maxD > 0 ? (radius / maxD) : 1;
    return centered.map(p => ({ x: p.x * s, y: p.y * s }));
}

function sideRatio(poly) {
    let mn = Infinity, mx = 0;
    for (let i = 0; i < poly.length; i++) {
        const a = poly[i];
        const b = poly[(i + 1) % poly.length];
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const d = Math.sqrt(dx * dx + dy * dy);
        if (d < mn) mn = d;
        if (d > mx) mx = d;
    }
    return mx / mn;
}

function minSide(poly) {
    let mn = Infinity;
    for (let i = 0; i < poly.length; i++) {
        const a = poly[i];
        const b = poly[(i + 1) % poly.length];
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const d = Math.sqrt(dx * dx + dy * dy);
        if (d < mn) mn = d;
    }
    return mn;
}

// Strict convex hull (no collinear)
function convexHull(points) {
    if (points.length <= 1) return points.slice();

    const pts = points.slice().sort((a, b) => (a.x === b.x ? a.y - b.y : a.x - b.x));
    const cross = (o, a, b) => (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);

    const lower = [];
    for (let i = 0; i < pts.length; i++) {
        const p = pts[i];
        while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) lower.pop();
        lower.push(p);
    }

    const upper = [];
    for (let i = pts.length - 1; i >= 0; i--) {
        const p = pts[i];
        while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) upper.pop();
        upper.push(p);
    }

    lower.pop();
    upper.pop();
    return lower.concat(upper);
}
