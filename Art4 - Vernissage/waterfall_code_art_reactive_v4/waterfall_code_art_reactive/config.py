from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class Settings:
    """Editable project settings.

    The defaults below keep the original atmospheric waterfall, but add a
    stronger text interaction layer and a reactive droplet pass inspired by the
    test-claude version.
    """

    # Window / runtime -----------------------------------------------------
    window_size: tuple[int, int] = (1600, 900)
    start_fullscreen: bool = False
    fps: int = 60
    title: str = "Waterfall Code / Reactive Projection Scene"

    # Render resolution ----------------------------------------------------
    # The internal render is still low-ish for a projected, slightly soft look,
    # but raised a bit so the code shapes remain legible when droplets bounce.
    render_height: int = 320

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
    code_scroll_speed: float = 16.0   # simulated pixels / second
    code_margin_x: int = 18
    code_margin_y: int = 18
    code_font_size: int = 13
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
    water_speed: float = 0.98
    vertical_density: float = 4.35
    side_sway: float = 0.072
    text_distortion: float = 0.110
    text_reaction: float = 0.38
    blur_passes: int = 5
    grain_amount: float = 0.018

    ribbon_frequency: float = 18.0
    thread_frequency: float = 36.0
    micro_frequency: float = 80.0
    ribbon_sharpness: float = 8.5
    thread_sharpness: float = 14.0
    micro_sharpness: float = 20.0

    # Extra text-driven shaping -------------------------------------------
    letter_refraction: float = 0.058
    wake_strength: float = 0.17
    edge_sparkle: float = 0.14

    # Reactive droplet layer ----------------------------------------------
    enable_droplets: bool = True
    show_background: bool = True

    droplet_max_count: int = 720
    droplet_spawn_rate: float = 260.0  # droplets / second
    droplet_min_speed: float = 30.0
    droplet_max_speed: float = 82.0
    droplet_min_size: float = 1.1
    droplet_max_size: float = 3.2
    droplet_alpha_min: int = 90
    droplet_alpha_max: int = 220
    droplet_wobble: float = 9.0
    droplet_wind: float = 7.0
    droplet_repel_strength: float = 32.0
    droplet_tangent_strength: float = 14.0
    droplet_drag_in_text: float = 0.34
    droplet_bounce_restitution: float = 0.72
    droplet_upward_boost: float = 28.0
    droplet_collision_cooldown: float = 0.065
    droplet_separation: float = 4.2

    mist_particles: int = 110
    mist_alpha: int = 22
    mist_drift: float = 7.0
    mist_max_size: float = 5.0

    droplet_splash_probability: float = 0.48
    droplet_splash_count_min: int = 2
    droplet_splash_count_max: int = 5
    droplet_splash_speed: float = 25.0

    foam_opacity: float = 0.55
    foam_opacity_step: float = 0.08
    foam_opacity_min: float = 0.0
    foam_opacity_max: float = 1.25

    # GPU / CUDA -----------------------------------------------------------
    prefer_cuda: bool = False

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
