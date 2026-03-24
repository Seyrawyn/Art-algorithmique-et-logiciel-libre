from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class Settings:
    """Editable project settings.

    Most of the visual character lives here, so you can change the scene without
    digging through the render code.
    """

    # Window / runtime -----------------------------------------------------
    window_size: tuple[int, int] = (1600, 900)
    start_fullscreen: bool = False
    fps: int = 60
    title: str = "Waterfall Code / Generative Projection Scene"

    # Render resolution ----------------------------------------------------
    # The scene is rendered to a smaller off-screen surface and then upscaled.
    # This keeps the animation smooth while giving it a soft projected feel.
    render_height: int = 270

    # Source code overlay --------------------------------------------------
    samples_dir: str = "code_samples"
    supported_extensions: tuple[str, ...] = (
        ".py",
        ".txt",
        ".md",
        ".json",
        ".toml",
        ".yaml",
        ".yml",
        ".glsl",
    )
    code_scroll_speed: float = 18.0   # simulated pixels / second
    code_margin_x: int = 18
    code_margin_y: int = 18
    code_font_size: int = 12
    code_line_spacing: int = 2
    code_wobble_px: float = 1.75
    code_file_gap: int = 24
    code_tab_size: int = 4
    code_edge_fade: int = 18

    mono_font_candidates: tuple[str, ...] = (
        "JetBrains Mono",
        "Fira Code",
        "Cascadia Code",
        "Consolas",
        "DejaVu Sans Mono",
        "Menlo",
        "Monaco",
        "Courier New",
        "Courier",
    )

    # Waterfall motion -----------------------------------------------------
    water_speed: float = 0.95
    vertical_density: float = 4.2
    side_sway: float = 0.065
    text_distortion: float = 0.085
    text_reaction: float = 0.26
    blur_passes: int = 4
    grain_amount: float = 0.018

    ribbon_frequency: float = 18.0
    thread_frequency: float = 36.0
    micro_frequency: float = 72.0
    ribbon_sharpness: float = 8.0
    thread_sharpness: float = 13.0
    micro_sharpness: float = 18.0

    # Look -----------------------------------------------------------------
    palette_name: str = "midnight_ice"
    show_startup_controls: bool = True
    screenshot_dir: str = "captures"

    # Misc -----------------------------------------------------------------
    default_seed: int = 7
    data_root: Path = field(default_factory=lambda: Path(__file__).resolve().parent)

    def render_size(self, output_size: tuple[int, int]) -> tuple[int, int]:
        """Compute the off-screen render size for the current output aspect."""
        output_w, output_h = output_size
        render_h = max(180, min(self.render_height, output_h))
        render_w = max(320, int(round(render_h * (output_w / max(1, output_h)))))
        return render_w, render_h

    def project_path(self, relative: str) -> Path:
        return self.data_root / relative

    def samples_path(self) -> Path:
        return self.project_path(self.samples_dir)

    def screenshot_path(self) -> Path:
        return self.project_path(self.screenshot_dir)
