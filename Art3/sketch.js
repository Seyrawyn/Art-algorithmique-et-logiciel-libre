
// p5.js — Trip Walkers (random file per year via manifest)
//
// Art idea:
//   A crowd of tiny "people" continuously walk from START to END points of trips.
//   The points are the beginning and ending of each trajectory (lat/lon).
//
// Requires preprocessed JSON files produced by:
//   preprocess_trip_walkers_auto_manifest.py
//
// Folder layout expected:
//   data/
//     trips/
//       2024/
//         manifest.json
//         ... many *_trips.json
//       2025/
//         manifest.json
//         ... many *_trips.json
//
// Controls:
//   N = pick NEW random file pair (2024 + 2025) and restart
//   R = restart same pair (same files, new seed)
//   S = save PNG
//   I = toggle info box
//   F = toggle slow fade (for motion / ghosting)
//
// Notes:
// - A browser can't list a directory, so we load `manifest.json` to know what files exist.
// - Each trip record is: [startLat, startLon, endLat, endLon, distance_m]

const PREFIX_2024 = "data/trips/2024/";
const PREFIX_2025 = "data/trips/2025/";
const MANIFEST_NAME = "manifest.json";

const CANVAS_SIZE = 1000;
const MARGIN = 45;

// Visual / performance knobs
const MAX_AGENTS_PER_YEAR = 50;  // number of little people per dataset (total ~2x)
const CAP_M = 20000;             // ignore distances above this (should already be filtered by preprocessing)

// Bounds robustness: ignore extreme outliers (e.g., a few bad lon values)
// 0.005 = keep 0.5%..99.5% percentile range
const BOUNDS_Q = 0.005;

// Fade mode makes the animation more "alive" but less "print-like"
let fadeOn = true;         // start with fade ON (still toggle with F)
const FADE_SECONDS = 1.0;  // how long until a line is basically gone
const FADE_TARGET  = 0.05; // how much “ink” remains after FADE_SECONDS (5%)

// Data
let manifest2024 = [];
let manifest2025 = [];
let selectedFile2024 = "";
let selectedFile2025 = "";

let data2024 = null;
let data2025 = null;
let trips2024 = [];
let trips2025 = [];

// We walk through trips in a shuffled order so the drawing uses the whole file,
// not just random sampling.
let nextIdx2024 = 0;
let nextIdx2025 = 0;

// Mapping lat/lon -> canvas
let bounds = null;        // {minLat,maxLat,minLon,maxLon}
let proj = null;          // projection params
let globalMaxM = 10000;   // for log scaling

// Agents (little people)
let agents = [];          // unified list, each agent has a `year` field

// UI state
let infoOn = true;
let ready = false;

function preload() {
  // Load manifests (lists of filenames)
  // NOTE: loadJSON in preload is blocking (p5 waits until loaded).
  manifest2024 = loadJSON(PREFIX_2024 + MANIFEST_NAME);
  manifest2025 = loadJSON(PREFIX_2025 + MANIFEST_NAME);
}

function setup() {
  createCanvas(CANVAS_SIZE, CANVAS_SIZE);
  pixelDensity(2);
  colorMode(HSB, 360, 100, 100, 100);
  background(0, 0, 5);

  // We will start drawing only after the selected trip files are loaded.
  noLoop();

  // Normalize manifests (loadJSON can return object; we want arrays of strings)
  manifest2024 = normalizeManifest(manifest2024);
  manifest2025 = normalizeManifest(manifest2025);

  pickNewPair();
}

function draw() {
  if (!ready) return;

  if (fadeOn) {
    // translucent overlay for motion blur / ghosting
    noStroke();
    const a = 1 - Math.pow(FADE_TARGET, deltaTime / (FADE_SECONDS * 1000));
    fill(0, 0, 5, 100 * a);   // alpha is 0..100 because you set HSB alpha max to 100
    rect(0, 0, width, height);
  }

  blendMode(ADD);

  for (let i = 0; i < agents.length; i++) {
    stepAgent(agents[i]);
  }

  blendMode(BLEND);

  if (infoOn) drawInfo();
}

function keyPressed() {
  if (key === 'n' || key === 'N') pickNewPair();
  if (key === 'r' || key === 'R') restartSamePair();
  if (key === 's' || key === 'S') saveCanvas("trip-walkers", "png");
  if (key === 'i' || key === 'I') infoOn = !infoOn;
  if (key === 'f' || key === 'F') fadeOn = !fadeOn;
}

// -------------------------------
// Loading / restarting
// -------------------------------

function pickNewPair() {
  ready = false;
  noLoop();
  background(0, 0, 5);

  if (!manifest2024.length || !manifest2025.length) {
    console.error("Manifest missing or empty. Check your folder structure.");
    return;
  }

  selectedFile2024 = manifest2024[Math.floor(Math.random() * manifest2024.length)];
  selectedFile2025 = manifest2025[Math.floor(Math.random() * manifest2025.length)];

  loadPair(selectedFile2024, selectedFile2025);
}

function restartSamePair() {
  if (!selectedFile2024 || !selectedFile2025) {
    pickNewPair();
    return;
  }
  ready = false;
  noLoop();
  background(0, 0, 5);

  // Keep same files, but change seeds so the drawing evolves differently
  randomSeed(Math.floor(Math.random() * 1e9));
  noiseSeed(Math.floor(Math.random() * 1e9));

  loadPair(selectedFile2024, selectedFile2025);
}

function loadPair(file2024, file2025) {
  let pending = 2;

  data2024 = null;
  data2025 = null;

  loadJSON(
    PREFIX_2024 + file2024,
    (d) => {
      data2024 = d;
      pending--;
      if (pending === 0) initAfterLoad();
    },
    (err) => console.error("Failed to load 2024 file:", file2024, err)
  );

  loadJSON(
    PREFIX_2025 + file2025,
    (d) => {
      data2025 = d;
      pending--;
      if (pending === 0) initAfterLoad();
    },
    (err) => console.error("Failed to load 2025 file:", file2025, err)
  );
}

function initAfterLoad() {
  trips2024 = extractTrips(data2024);
  trips2025 = extractTrips(data2025);

  if (!trips2024.length || !trips2025.length) {
    console.error("Trips are empty. Check preprocessing output format.");
    return;
  }

  // Shuffle trips so we "use the whole file" but in a random order
  shuffleInPlace(trips2024);
  shuffleInPlace(trips2025);
  nextIdx2024 = 0;
  nextIdx2025 = 0;

  // Robust bounds (quantile-based) to avoid a few outliers ruining the projection
  let b2024 = robustBoundsFromTrips(trips2024, BOUNDS_Q);
  let b2025 = robustBoundsFromTrips(trips2025, BOUNDS_Q);

  bounds = {
    minLat: min(b2024.minLat, b2025.minLat),
    maxLat: max(b2024.maxLat, b2025.maxLat),
    minLon: min(b2024.minLon, b2025.minLon),
    maxLon: max(b2024.maxLon, b2025.maxLon)
  };

  proj = makeProjection(bounds);

  // For scaling
  let md2024 = extractMaxDist(data2024, trips2024);
  let md2025 = extractMaxDist(data2025, trips2025);
  globalMaxM = min(max(md2024, md2025), CAP_M);

  // Fresh seeds for a new run (keeps file selection independent)
  randomSeed(Math.floor(Math.random() * 1e9));
  noiseSeed(Math.floor(Math.random() * 1e9));

  // Create agents
  agents = [];

  let nA = min(MAX_AGENTS_PER_YEAR, Math.max(70, Math.floor(trips2024.length / 70)));
  let nB = min(MAX_AGENTS_PER_YEAR, Math.max(70, Math.floor(trips2025.length / 70)));

  for (let i = 0; i < nA; i++) agents.push(makeAgent(2024));
  for (let i = 0; i < nB; i++) agents.push(makeAgent(2025));

  // Assign initial trips (also draws initial endpoints faintly)
  for (let i = 0; i < agents.length; i++) {
    assignTrip(agents[i]);
  }

  ready = true;
  loop();
}

// -------------------------------
// Agents (little people)
// -------------------------------

function makeAgent(year) {
  return {
    year,

    // current trip (screen-space)
    s: createVector(0, 0),
    e: createVector(0, 0),
    dir: createVector(1, 0),
    perp: createVector(0, 1),
    lenPx: 1,

    // motion
    t: 0,
    speed: 0.01,
    xNorm: 0.5,              // 0..1 based on distance

    // styling
    hueBase: (year === 2024) ? 270 : 70, // Color for each year
    noiseOffset: random(1000),
    walkOffset: random(TWO_PI),
    size: 8
  };
}

function getNextTrip(year) {
  if (year === 2024) {
    if (nextIdx2024 >= trips2024.length) {
      shuffleInPlace(trips2024);
      nextIdx2024 = 0;
    }
    return trips2024[nextIdx2024++];
  } else {
    if (nextIdx2025 >= trips2025.length) {
      shuffleInPlace(trips2025);
      nextIdx2025 = 0;
    }
    return trips2025[nextIdx2025++];
  }
}

function assignTrip(agent) {
  // Try a few times to avoid degenerate / out-of-bounds trips
  for (let attempts = 0; attempts < 60; attempts++) {
    let trip = getNextTrip(agent.year);
    if (!trip || trip.length < 5) continue;

    let slat = trip[0], slon = trip[1], elat = trip[2], elon = trip[3], d = trip[4];
    if (!isFinite(slat) || !isFinite(slon) || !isFinite(elat) || !isFinite(elon) || !isFinite(d)) continue;
    if (d <= 0 || d > CAP_M) continue;

    let s = projectLatLon(slat, slon);
    let e = projectLatLon(elat, elon);

    // If robust bounds clipped outliers, those points end up off-canvas — skip them.
    if (!insideCanvas(s) || !insideCanvas(e)) continue;

    let v = p5.Vector.sub(e, s);
    let lenPx = v.mag();
    if (lenPx < 2) continue;

    agent.s = s;
    agent.e = e;
    agent.dir = v.copy().normalize();
    agent.perp = createVector(-agent.dir.y, agent.dir.x);
    agent.lenPx = lenPx;

    agent.xNorm = normLog(d, globalMaxM);

    // "Walking speed" as a fraction of the segment per frame.
    // Longer trips -> slightly bigger stride, but still take longer overall.
    let stridePx = lerp(0.8, 4.0, agent.xNorm);
    agent.speed = stridePx / agent.lenPx;

    agent.t = 0;

    // Person size relates to distance (subtle)
    agent.size = lerp(6, 11, agent.xNorm);

    // Draw endpoints once (accumulates)
    drawEndpoints(agent);

    return;
  }
}

function stepAgent(agent) {
  let t0 = agent.t;
  let t1 = agent.t + agent.speed;

  // Style by year + distance
  let hue = (agent.hueBase + 70 * agent.xNorm) % 360;
  let sat = 35 + 65 * agent.xNorm;
  let bri = 65 + 35 * (1 - agent.xNorm);

  // Because we redraw the whole partial line every frame, keep alpha modest
  let alpha = 3 + 12 * agent.xNorm;

  stroke(hue, sat, bri, alpha);
  strokeWeight(0.25 + 1.8 * agent.xNorm);

  // If finished: draw the full line once, then stop redrawing it so it fades
  if (t1 >= 1) {
    line(agent.s.x, agent.s.y, agent.e.x, agent.e.y);
    drawPerson(agent, agent.e, hue, sat, bri);
    drawArrival(agent, agent.e);
    assignTrip(agent);
    return;
  }

  // Draw the whole "so far" route (refreshes it, so it doesn't fade mid-walk)
  let onLine = lerpVec(agent.s, agent.e, t1);
  line(agent.s.x, agent.s.y, onLine.x, onLine.y);

  // Person position:
  // Option A (clean): keep person exactly on the line
  // drawPerson(agent, onLine, hue, sat, bri);

  // Option B (if you want the old wiggle back, comment Option A and use this):
  let p1 = onLine.copy();
  let amp = lerp(0.0, 14.0, agent.xNorm);
  let w1 = (noise(agent.noiseOffset, t1 * 3.0) - 0.5) * 2.0;
  p1.add(p5.Vector.mult(agent.perp, w1 * amp));
  drawPerson(agent, p1, hue, sat, bri);

  agent.t = t1;
}

function drawPerson(agent, pos, hue, sat, bri) {
  let heading = atan2(agent.dir.y, agent.dir.x);

  // Walk cycle
  let walk = sin(frameCount * 0.25 + agent.walkOffset);

  push();
  translate(pos.x, pos.y);
  rotate(heading);

  let s = agent.size;

  // Body style: a bright outline + a faint "glow dot"
  noStroke();
  fill(hue, sat, 100, 10);
  circle(0, -s * 0.15, s * 1.2);

  stroke(hue, sat, 100, 42);
  strokeWeight(max(0.7, s * 0.10));
  noFill();

  // Head
  circle(0, -s * 0.55, s * 0.42);

  // Body
  line(0, -s * 0.35, 0, s * 0.25);

  // Arms
  line(0, -s * 0.10,  s * 0.35, s * 0.05 + walk * s * 0.07);
  line(0, -s * 0.10, -s * 0.35, s * 0.05 - walk * s * 0.07);

  // Legs
  line(0, s * 0.25,  s * 0.25, s * 0.62 + walk * s * 0.10);
  line(0, s * 0.25, -s * 0.25, s * 0.62 - walk * s * 0.10);

  pop();
}

function drawEndpoints(agent) {
  let hue = agent.hueBase;
  let sat = 30;
  let bri = 85;

  // Start dot
  noStroke();
  fill(hue, sat, bri, 10);
  circle(agent.s.x, agent.s.y, 2.1);

  // End dot
  fill(hue, sat, bri, 16);
  circle(agent.e.x, agent.e.y, 2.6);

  // Faint baseline segment
  stroke(hue, sat, bri, 4);
  strokeWeight(0.3);
  line(agent.s.x, agent.s.y, agent.e.x, agent.e.y);
}

function drawArrival(agent, pos) {
  let hue = (agent.hueBase + 80 * agent.xNorm) % 360;
  noStroke();
  fill(hue, 80, 100, 28);
  circle(pos.x, pos.y, 1.5 + 8.0 * agent.xNorm);
}

// -------------------------------
// Projection (lat/lon -> canvas)
// -------------------------------

function makeProjection(b) {
  let latPad = (b.maxLat - b.minLat) * 0.05;
  let lonPad = (b.maxLon - b.minLon) * 0.05;

  let minLat = b.minLat - latPad;
  let maxLat = b.maxLat + latPad;
  let minLon = b.minLon - lonPad;
  let maxLon = b.maxLon + lonPad;

  let usableW = width - 2 * MARGIN;
  let usableH = height - 2 * MARGIN;

  let lonRange = maxLon - minLon;
  let latRange = maxLat - minLat;

  let sx = usableW / lonRange;
  let sy = usableH / latRange;
  let scale = min(sx, sy);

  let contentW = lonRange * scale;
  let contentH = latRange * scale;

  let ox = (width - contentW) / 2;
  let oy = (height - contentH) / 2;

  return { minLat, maxLat, minLon, maxLon, scale, ox, oy };
}

function projectLatLon(lat, lon) {
  let x = proj.ox + (lon - proj.minLon) * proj.scale;
  let y = proj.oy + (proj.maxLat - lat) * proj.scale; // invert lat so north is up
  return createVector(x, y);
}

function insideCanvas(v) {
  // A little leeway so points near edges still count
  return (v.x > -10 && v.x < width + 10 && v.y > -10 && v.y < height + 10);
}

// -------------------------------
// Data helpers
// -------------------------------

function normalizeManifest(m) {
  if (Array.isArray(m)) return m;
  if (m && Array.isArray(m.files)) return m.files;
  if (m && Array.isArray(m.manifest)) return m.manifest;

  let out = [];
  if (m && typeof m === "object") {
    for (let k in m) {
      if (typeof m[k] === "string") out.push(m[k]);
    }
  }
  return out;
}

function extractTrips(obj) {
  if (!obj) return [];
  if (Array.isArray(obj)) return obj;
  if (Array.isArray(obj.trips)) return obj.trips;
  if (Array.isArray(obj.data)) return obj.data;
  return [];
}

function extractMaxDist(obj, trips) {
  if (obj && isFinite(obj.maxDistM)) return obj.maxDistM;
  let m = 0;
  for (let i = 0; i < trips.length; i++) {
    let d = trips[i][4];
    if (isFinite(d) && d > m) m = d;
  }
  return m;
}

function normLog(d, maxD) {
  let hi = Math.log(1 + maxD);
  return constrain(Math.log(1 + d) / hi, 0, 1);
}

function lerpVec(a, b, t) {
  return createVector(
    lerp(a.x, b.x, t),
    lerp(a.y, b.y, t)
  );
}

function shuffleInPlace(arr) {
  // Fisher–Yates
  for (let i = arr.length - 1; i > 0; i--) {
    let j = Math.floor(Math.random() * (i + 1));
    let tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
  }
}

function robustBoundsFromTrips(trips, q) {
  // Collect coords
  let lats = [];
  let lons = [];
  for (let i = 0; i < trips.length; i++) {
    let t = trips[i];
    if (!t || t.length < 4) continue;
    let slat = t[0], slon = t[1], elat = t[2], elon = t[3];
    if (!isFinite(slat) || !isFinite(slon) || !isFinite(elat) || !isFinite(elon)) continue;
    lats.push(slat); lats.push(elat);
    lons.push(slon); lons.push(elon);
  }
  lats.sort((a,b) => a-b);
  lons.sort((a,b) => a-b);

  let nLat = lats.length;
  let nLon = lons.length;

  let loLat = lats[Math.floor(q * (nLat - 1))];
  let hiLat = lats[Math.floor((1 - q) * (nLat - 1))];
  let loLon = lons[Math.floor(q * (nLon - 1))];
  let hiLon = lons[Math.floor((1 - q) * (nLon - 1))];

  return { minLat: loLat, maxLat: hiLat, minLon: loLon, maxLon: hiLon };
}

// -------------------------------
// Info overlay
// -------------------------------

function drawInfo() {
  noStroke();
  fill(0, 0, 0, 55);
  rect(18, 18, 720, 132, 12);

  fill(0, 0, 95, 100);
  textSize(13);
  text("Trip Walkers — tiny people walk from trip START to END (N new files, R restart, S save, I info, F fade)", 30, 40);

  text("2024 file: " + selectedFile2024 + "  | trips: " + trips2024.length, 30, 64);
  text("2025 file: " + selectedFile2025 + "  | trips: " + trips2025.length, 30, 84);
  text("Agents: " + agents.length + "  | fade: " + (fadeOn ? "ON" : "OFF") + "  | bounds quantile: " + BOUNDS_Q, 30, 104);

  text(
    "Bounds used: lat[" + nf(bounds.minLat,1,4) + ", " + nf(bounds.maxLat,1,4) + "] lon[" + nf(bounds.minLon,1,4) + ", " + nf(bounds.maxLon,1,4) + "]",
    30, 124
  );
}
