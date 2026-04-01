# Waterfall Code / Reactive Edition

This version keeps the original smooth procedural waterfall, but adds a second
reactive droplet pass so the water can visibly bounce off the letters of the
scrolling code.

The main differences from the original project are:

- a stronger mask-derived influence field around the text
- a reactive droplet layer with reflection, upward kick, and splash bursts
- slightly sharper default internal rendering for better code legibility
- optional CuPy/CUDA acceleration for the heavy procedural waterfall math
- live background toggle and live foam-opacity control
- 15 projection-ready color palettes that can be cycled during performance
- an independent text-color palette you can change separately from the water palette

## Quick start

```bash
pip install -r requirements.txt
python main.py
```

To start with a different text color from the main waterfall palette:

```bash
python main.py --palette midnight_ice --text-palette ember_lava
```

```bash
python main.py --droplets-only --palette midnight_ice --text-palette ember_lava --fullscreen
```

To start with the procedural waterfall disabled while keeping the scrolling code and droplets visible:

```bash
python main.py --droplets-only --cuda
```

To try CUDA acceleration:

```bash
# First install a CuPy build that matches your NVIDIA/CUDA setup.
python main.py --cuda
```

## Controls

- `Esc` / `Q` — quit
- `Space` — pause / resume
- `R` — restart current file and reset droplets
- `N` / `P` — next / previous code file
- `C` / `Tab` — cycle the visual waterfall / droplet palette
- `T` / `Shift+T` — cycle the text color palette forward / backward independently
- `B` — toggle the procedural waterfall background on/off while keeping the code visible
- `[` / `]` — decrease / increase collision-foam opacity
- `D` — toggle reactive droplets on/off
- `F` — toggle fullscreen
- `L` — reload the samples directory
- `S` — save a screenshot

A small status label appears on screen when you change one of these live
controls, so you can use them even in fullscreen during a presentation.

## What changed

### 1. Text influence is now shared across both render layers
The code mask is converted into:

- `spread` — soft area of influence around letters
- `edge` — contour energy around glyph shapes
- `pressure` — stronger collision field for droplets
- `nx`, `ny` — outward normals used to push droplets away from the code

### 2. A droplet layer was added on top of the procedural waterfall
The original project looked elegant, but the interaction stayed mostly inside
one smooth shader-like field. This edition adds visible water droplets that:

- fall from the top in front of the waterfall
- bend as they enter the influence zone of the text
- bounce away from letter contours using the local normal field
- flare briefly on impact
- spawn small splash bursts on stronger collisions

### 3. The waterfall itself reacts more strongly to the code
The procedural layer now uses stronger text-driven warp, refraction, wake, and
edge sparkle terms. The letters feel more submerged inside the flow.

### 4. Foam can now stay translucent over the code
The collision foam uses an adjustable live opacity control and is blended with
standard alpha instead of only additive bloom, so the code remains readable.

### 5. There are now 15 palettes to cycle live
You can cycle through cool, neutral, warm, and neon variants during the show,
including `midnight_ice`, `aurora_mint`, `violet_noir`, `ember_lava`,
`forest_mist`, and more.

### 6. Text color can now be changed independently
The scrolling code no longer has to inherit its color from the main water
palette. The text uses the same 15 palette families, but on its own control, so
you can keep it matched to the waterfall or deliberately contrast it.

## CUDA / GPU note

This project uses a practical optional path with **CuPy** for the heavy array
math in the procedural waterfall renderer. If CuPy is available and you pass
`--cuda`, the render field runs on the GPU; if not, the project falls back to
plain NumPy automatically.

The droplet pass still draws through pygame, which keeps the project simple and
portable while letting the most expensive full-frame math move to the GPU.

## Good starting presets

For a presentation screen:

- `python main.py --fullscreen --palette midnight_ice`
- `python main.py --fullscreen --palette ink_silver --render-height 360`
- `python main.py --fullscreen --palette midnight_ice --render-height 420 --cuda`
- `python main.py --fullscreen --droplets-only --palette storm_cyan`  # code + droplets on black
- `python main.py --fullscreen --palette ink_silver --text-palette ember_lava`

## Most important tweak points

Open `config.py`.

### To get more bounce on the letters
- `droplet_bounce_restitution`
- `droplet_upward_boost`
- `droplet_repel_strength`
- `droplet_tangent_strength`
- `droplet_splash_probability`

### To control the collision foam visibility
- `foam_opacity`
- `foam_opacity_step`
- `droplet_splash_probability`
- `droplet_splash_count_min`
- `droplet_splash_count_max`

### To make the embedded-code interaction stronger
- `text_distortion`
- `text_reaction`
- `letter_refraction`
- `wake_strength`
- `edge_sparkle`
- `blur_passes`

### To change code readability / projection sharpness
- `render_height`
- `code_font_size`
- `code_scroll_speed`
- `code_wobble_px`

### To trade performance for detail
- `render_height`
- `droplet_max_count`
- `droplet_spawn_rate`
- `mist_particles`

## File structure

```txt
waterfall_code_art_reactive/
├── main.py
├── config.py
├── palettes.py
├── code_scroller.py
├── influence.py
├── waterfall.py
├── droplets.py
├── code_samples/
└── README.md
```
