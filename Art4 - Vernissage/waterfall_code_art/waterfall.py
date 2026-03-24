from __future__ import annotations

import numpy as np

from config import Settings
from palettes import Palette


class WaterfallRenderer:
    def __init__(
        self,
        settings: Settings,
        render_size: tuple[int, int],
        palette: Palette,
        seed: int = 7,
    ) -> None:
        self.settings = settings
        self.width, self.height = render_size
        self.palette = palette
        self.seed = seed

        x = np.linspace(0.0, 1.0, self.width, dtype=np.float32)
        y = np.linspace(0.0, 1.0, self.height, dtype=np.float32)
        self.x, self.y = np.meshgrid(x, y)
        self.curtain = 0.82 + 0.18 * np.exp(-((self.x - 0.5) ** 2) / 0.18)

        rng = np.random.default_rng(seed)
        self.grain = rng.random((self.height, self.width), dtype=np.float32)
        self._update_palette_arrays()

    def _color(self, rgb: tuple[int, int, int]) -> np.ndarray:
        return np.array(rgb, dtype=np.float32) / 255.0

    def _update_palette_arrays(self) -> None:
        self.bg_top = self._color(self.palette.bg_top)
        self.bg_bottom = self._color(self.palette.bg_bottom)
        self.water_dark = self._color(self.palette.water_dark)
        self.water_mid = self._color(self.palette.water_mid)
        self.water_bright = self._color(self.palette.water_bright)
        self.code_low = self._color(self.palette.code_low)
        self.code_high = self._color(self.palette.code_high)
        self.mist = self._color(self.palette.mist)
        self.base_background = self.bg_top * (1.0 - self.y[..., None]) + self.bg_bottom * self.y[..., None]

    def set_palette(self, palette: Palette) -> None:
        self.palette = palette
        self._update_palette_arrays()

    def _influence_fields(
        self, mask: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        core = mask.astype(np.float32)
        spread = core.copy()

        for _ in range(self.settings.blur_passes):
            spread = (
                spread * 4.0
                + np.roll(spread, 1, axis=0)
                + np.roll(spread, -1, axis=0)
                + np.roll(spread, 1, axis=1)
                + np.roll(spread, -1, axis=1)
            ) / 8.0

        halo = np.clip(spread - core * 0.55, 0.0, 1.0)
        gx = 0.5 * (np.roll(spread, -1, axis=1) - np.roll(spread, 1, axis=1))
        gy = 0.5 * (np.roll(spread, -1, axis=0) - np.roll(spread, 1, axis=0))
        edge = np.clip(np.sqrt((gx * gx) + (gy * gy)) * 5.0, 0.0, 1.0)
        return core, spread, halo, gx, gy, edge

    def render(self, mask: np.ndarray, t: float) -> np.ndarray:
        core, spread, halo, gx, gy, edge = self._influence_fields(mask)

        time = t * self.settings.water_speed
        x = self.x
        y = self.y

        # Large-scale downward motion and lateral sway.
        warp = (
            self.settings.side_sway * np.sin((6.0 * y) + (1.8 * time))
            + 0.03 * np.sin((18.0 * y) - (0.7 * time) + 2.5 * np.sin((4.2 * x) + (0.2 * time)))
            + 0.018 * np.sin((42.0 * y) + (14.0 * x) - (1.15 * time))
        )
        text_warp = self.settings.text_distortion * (
            (0.75 * gx)
            + 0.45 * spread * np.sin((34.0 * y) - (5.3 * time) + (24.0 * x))
        )

        xw = x + warp + text_warp
        yw = y + (self.settings.text_distortion * 0.35 * gy)
        flow = (yw * self.settings.vertical_density) - time

        ribbon_phase = (xw * self.settings.ribbon_frequency) + 0.72 * np.sin(flow * 1.25)
        thread_phase = (xw * self.settings.thread_frequency) - 0.42 * np.sin((flow * 2.1) + (3.1 * x))
        micro_phase = (xw * self.settings.micro_frequency) + 0.22 * np.sin((flow * 4.1) - (2.0 * x))

        ribbons = np.exp(-self.settings.ribbon_sharpness * np.sin(np.pi * ribbon_phase) ** 2)
        threads = np.exp(-self.settings.thread_sharpness * np.sin(np.pi * thread_phase) ** 2)
        micro = np.exp(-self.settings.micro_sharpness * np.sin(np.pi * micro_phase) ** 2)

        mist_wave = 0.5 + 0.5 * np.sin((flow * 1.55) + 4.8 * np.sin((4.0 * xw) + (1.35 * time)))
        foam_wave = 0.5 + 0.5 * np.sin((35.0 * yw) - (7.0 * time) + 3.0 * np.sin(10.0 * xw))
        shimmer = 0.5 + 0.5 * np.sin((22.0 * xw) + (48.0 * yw) - (6.8 * time))

        agitation = spread * (0.55 + 0.45 * np.sin((18.0 * x) + (48.0 * y) - (7.5 * time)))
        turbulence = edge * (0.5 + 0.5 * np.sin((26.0 * xw) + (30.0 * yw) - (8.5 * time)))

        intensity = (
            0.10
            + 0.54 * ribbons
            + 0.34 * threads * (0.55 + 0.45 * mist_wave)
            + 0.15 * micro * (0.25 + 0.75 * foam_wave)
            + 0.06 * mist_wave
        )
        intensity *= self.curtain
        intensity += self.settings.text_reaction * agitation
        intensity += 0.18 * turbulence
        intensity += 0.06 * shimmer
        intensity = np.clip(intensity, 0.0, 1.35)

        grad = np.clip(intensity, 0.0, 1.0)
        highlight = np.clip((grad - 0.42) * 1.95, 0.0, 1.0)

        water = self.water_dark * (1.0 - grad[..., None]) + self.water_mid * grad[..., None]
        water = water * (1.0 - highlight[..., None]) + self.water_bright * highlight[..., None]

        rgb = self.base_background * (1.0 - 0.88 * grad[..., None])
        rgb += water * (0.28 + 0.94 * grad[..., None])

        halo_glow = halo[..., None] * self.mist * (0.08 + 0.10 * mist_wave[..., None])
        rgb += halo_glow

        # Sink the text slightly into the water, then let it flare at the edges.
        text_mix = 0.24 + 0.38 * grad + 0.20 * (0.5 + 0.5 * np.sin((4.2 * time) + (12.0 * y) + (1.4 * x)))
        text_color = self.code_low * (1.0 - text_mix[..., None]) + self.code_high * text_mix[..., None]
        rgb = rgb * (1.0 - 0.34 * core[..., None]) + text_color * core[..., None] * 0.86 + rgb * core[..., None] * 0.14
        rgb += halo[..., None] * self.code_high * 0.08

        if self.settings.grain_amount > 0.0:
            shift_y = int(time * 28.0) % self.height
            shift_x = int(time * 9.0) % self.width
            grain = np.roll(self.grain, shift=(shift_y, shift_x), axis=(0, 1))
            rgb += (grain[..., None] - 0.5) * self.settings.grain_amount

        rgb = np.clip(rgb, 0.0, 1.0)
        return (rgb * 255.0).astype(np.uint8)
