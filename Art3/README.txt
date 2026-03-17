Trip Walkers — p5.js
====================

This sketch draws a crowd of tiny stick-figure "people" walking from the START
point to the END point of each trip (trajectory).

1) Preprocess your CSVs into JSON trip files + manifests
--------------------------------------------------------

Use the included Python script:

  preprocess_trip_walkers_auto_manifest.py

From your project root (where the 'data/' folder lives), run:

  python preprocess_trip_walkers_auto_manifest.py data/output_2024 --recursive --outdir data/trips/2024
  python preprocess_trip_walkers_auto_manifest.py data/output_2025 --recursive --outdir data/trips/2025

This will create:
  data/trips/2024/*.json  + data/trips/2024/manifest.json
  data/trips/2025/*.json  + data/trips/2025/manifest.json

2) Put sketch.js + index.html in your project root
--------------------------------------------------

Expected layout:

  index.html
  sketch.js
  data/
    trips/
      2024/
        manifest.json
        ... *_trips.json
      2025/
        manifest.json
        ... *_trips.json

3) Run from a local web server
------------------------------

Many browsers block JSON loading from file://. From your project root:

  python -m http.server 8000

Then open:
  http://localhost:8000

Controls:
  N = pick NEW random file pair (2024 + 2025)
  R = restart same pair (same files)
  S = save PNG
  I = toggle info box
  F = toggle slow fade
