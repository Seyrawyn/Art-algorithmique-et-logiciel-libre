# Waterfall Code Art

A presentation-ready generative art scene in Python.

The piece shows a dark, high-contrast waterfall made from layered procedural bands and mist. Source code scrolls upward through the scene one file at a time. The water reacts around the letters with local distortion, shimmer, and turbulence so the code feels embedded in the flow instead of simply placed on top.

## Features

- Real-time animated waterfall
- Slow upward source-code scroll, one file at a time
- Local water reaction around the text
- Fullscreen toggle
- Palette cycling
- Reloadable code samples directory
- Screenshot capture
- Clean Python-only project structure

## File structure

```text
waterfall_code_art/
├── README.md
├── requirements.txt
├── config.py
├── palettes.py
├── code_scroller.py
├── waterfall.py
├── main.py
└── code_samples/
    ├── flow_field.py
    ├── palette_lab.py
    ├── signal_chain.py
    └── vector_cache.py
```

## Setup

Create a virtual environment if you want a clean install, then install the dependencies:

```bash
python -m pip install -r requirements.txt
```

## Run

From the project folder:

```bash
python main.py
```

Useful launch options:

```bash
python main.py --fullscreen
python main.py --palette deep_teal
python main.py --render-height 360
python main.py --samples /absolute/path/to/your/code/files
```

## Controls

- `Esc` or `Q` — quit
- `Space` — pause / resume
- `R` — restart current file and reset the animation clock
- `N` / `P` — next / previous file
- `C` or `Tab` — cycle palette
- `F` — toggle fullscreen
- `L` — reload the samples directory
- `S` — save a screenshot into `captures/`

## How the animation works

### 1. Off-screen rendering
The artwork is rendered to a smaller internal surface, then smoothly upscaled to the window or projector resolution. This keeps the frame rate high while giving the final image a soft, cinematic projected feel.

### 2. Waterfall synthesis
The waterfall is not a physics simulation. It is a layered procedural field made from:

- vertical ribbon bands
- finer thread bands
- micro streaks
- a drifting mist layer
- a light moving grain texture

These layers are built from trigonometric wave fields and shaped with exponential falloff. The result reads as flowing water without the overhead of a full fluid solver.

### 3. Code mask generation
Each code file is rendered as anti-aliased text on a transparent surface. The lines scroll upward at a slow, steady speed. A very slight line wobble keeps the text from feeling static.

### 4. Water/text interaction
The alpha mask of the text is converted into a small influence field:

- the mask is blurred to create a local area of influence
- gradients around the letters generate edge energy
- that influence perturbs the waterfall’s horizontal warp
- extra shimmer and turbulence are added around the text edges

This makes the water look agitated where it crosses the code.

### 5. Embedded text compositing
The code is composited inside the procedural image itself. The renderer slightly darkens the water under the text, then adds a cool luminous text color and a subtle halo. Because the distortion and compositing happen in the same render pass, the text feels fused into the waterfall.

## Swapping the displayed code files

Put your own files in `code_samples/`.

By default the project loads these file types:

- `.py`
- `.txt`
- `.md`
- `.json`
- `.toml`
- `.yaml`
- `.yml`
- `.glsl`

If you want different file types, change `supported_extensions` in `config.py`.

## Best settings for projection

For a large white screen, a good starting point is:

- fullscreen mode
- `midnight_ice` or `ink_silver`
- `render_height = 270` for stable performance
- `render_height = 360` if your machine is fast enough and you want a sharper look

## Parameters you can tweak

Open `config.py`.

### To change the overall motion
- `water_speed` — overall downward motion speed
- `vertical_density` — how stretched the waterfall feels vertically
- `side_sway` — broad side-to-side movement

### To change the interaction with the code
- `text_distortion` — how much the letters bend the flow
- `text_reaction` — how turbulent the water becomes around the text
- `blur_passes` — how wide the reaction halo is

### To change the code look
- `code_scroll_speed` — speed of the upward scroll
- `code_font_size` — base font size in the internal render surface
- `code_wobble_px` — how much each line gently wavers
- `code_margin_x`, `code_margin_y` — text placement margins

### To change the visual density
- `ribbon_frequency` — number of broad water bands
- `thread_frequency` — number of finer streaks
- `micro_frequency` — amount of high-frequency detail
- `ribbon_sharpness`, `thread_sharpness`, `micro_sharpness` — edge crispness of those layers

### To change the texture feel
- `grain_amount` — subtle projected grain
- `render_height` — quality/performance trade-off

### To change colors
Use a different preset in `palettes.py`, or create your own palette and set `palette_name` in `config.py`.

## Notes

- The code is intentionally soft and atmospheric rather than syntax-highlighted.
- The whole piece is designed to look elegant and immersive when projected, not like an editor UI.
- If you want a brighter, more monochrome stage look, start with `ink_silver`.
