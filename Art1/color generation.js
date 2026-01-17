// Random Painting: background + 3–8 harmonious flowing curves with BIG fading edges
// Click or press R to regenerate.

let pg;

function setup() {
    createCanvas(900, 600);
    pixelDensity(2);
    pg = createGraphics(width, height);
    pg.pixelDensity(2);
    noLoop();
    generatePainting();
}

function draw() {
    image(pg, 0, 0);
}

function mousePressed() { regenerate(); }
function keyPressed() { if (key === 'r' || key === 'R') regenerate(); }

function regenerate() {
    generatePainting();
    redraw();
}

function generatePainting() {
    pg.clear();
    pg.colorMode(HSB, 360, 100, 100, 1);
    colorMode(HSB, 360, 100, 100, 1);

    // --- Background (single flat color) ---
    const bg = pickBackgroundColor();
    pg.background(bg);

    // subtle texture so fading looks richer
    addSoftGrain(pg, 0.03);

    // --- Draw 3–8 strokes ---
    const n = floor(random(6, 15));

    const minSide = min(width, height);
    const diag = sqrt(width * width + height * height);

    for (let i = 0; i < n; i++) {
        const targetLen = random(0.5 * minSide, diag);
        const strokeCol = pickHarmoniousStroke(bg);
        drawFlowStroke(pg, strokeCol, targetLen);
    }

    // Stronger blur => bigger fade between background/lines
    pg.filter(BLUR, 15); // try 2.5–6 for more/less

    // subtle vignette
    addVignette(pg, bg);
}

// ---------------- Color helpers ----------------

function pickBackgroundColor() {
    const h = random(0, 360);
    const s = random(12, 35);
    const b = random(70, 92);
    return color(h, s, b, 1);
}

function pickHarmoniousStroke(bg) {
    const h0 = hue(bg);
    const s0 = saturation(bg);
    const b0 = brightness(bg);

    const scheme = random();
    let dh;
    if (scheme < 0.55) dh = random(-25, 25); // analogous
    else if (scheme < 0.85) dh = random([150, -150]) + random(-18, 18); // near-complement
    else dh = random([60, -60]) + random(-15, 15); // split-ish

    const h = (h0 + dh + 360) % 360;
    const s = constrain(s0 + random(20, 55), 8, 85);
    const b = constrain(b0 + random(-18, 12), 20, 98);

    return color(h, s, b, random(0.12, 0.22));
}

// ---------------- Stroke generation ----------------

function drawFlowStroke(g, col, targetLen) {
    g.push();
    g.noFill();

    // Thickness (wide, painterly)
    const baseW = random(10, 32);

    // Random start point
    let x = random(width);
    let y = random(height);

    // Flow field parameters
    const noiseScale = random(0.0012, 0.0035);
    const step = random(6, 14);
    const maxSteps = floor(targetLen / step);

    // Build polyline points following a noise-driven direction field
    const pts = [];
    pts.push({ x, y });

    let len = 0;
    for (let i = 0; i < maxSteps; i++) {
        const a = noise(x * noiseScale, y * noiseScale) * TWO_PI * 3.0;
        const nx = x + cos(a) * step;
        const ny = y + sin(a) * step;

        x = constrain(nx, 0, width);
        y = constrain(ny, 0, height);

        pts.push({ x, y });
        len += step;
        if (len >= targetLen) break;
    }

    // --- BIG FADE HALO (drawn first) ---
    // Bigger fade: increase haloSpread and/or haloLayers
    const haloLayers = 25;     // was 16 in the suggestion
    const haloSpread = 3.0;    // wider fade band
    for (let h = 0; h < haloLayers; h++) {
        const t = h / (haloLayers - 1);
        const w = baseW * (1.0 + t * haloSpread);
        const a = 0.04 * pow(1.0 - t, 1.9); // smooth falloff outward
        g.stroke(color(hue(col), saturation(col), brightness(col), a));
        g.strokeWeight(w);

        g.beginShape();
        for (let k = 0; k < pts.length; k++) g.curveVertex(pts[k].x, pts[k].y);
        g.endShape();
    }

    // Paint the main stroke in many semi-transparent layers for soft edges/fading
    const passes = floor(random(55, 90)); // increased for more buildup
    for (let p = 0; p < passes; p++) {
        const jitter = random(0.7, 2.4); // increased => fluffier boundary
        const w = baseW * (0.55 + 0.75 * noise(p * 0.2)) * random(0.85, 1.1);

        g.strokeWeight(w);
        g.stroke(setAlpha(col, random(0.03, 0.085))); // slightly stronger build-up

        g.beginShape();
        for (let k = 0; k < pts.length; k++) {
            const px = pts[k].x + random(-jitter, jitter);
            const py = pts[k].y + random(-jitter, jitter);
            g.curveVertex(px, py);
        }
        g.endShape();
    }

    // A faint highlight line adds depth
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

function setAlpha(c, a) {
    return color(hue(c), saturation(c), brightness(c), a);
}

// ---------------- Fading/texture helpers ----------------

function addSoftGrain(g, amount) {
    g.push();
    g.noStroke();
    const n = floor(width * height * 0.0012);
    for (let i = 0; i < n; i++) {
        const x = random(width);
        const y = random(height);
        const a = random(0, amount);
        g.fill(0, 0, random(0, 100), a);
        g.circle(x, y, random(0.5, 1.6));
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

        const w = width * (1.05 + t * 0.55);
        const h = height * (1.05 + t * 0.55);
        g.ellipse(width / 2, height / 2, w, h);
    }

    g.pop();
}
