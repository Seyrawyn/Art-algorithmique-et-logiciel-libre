from __future__ import annotations

import numpy as np

from config import Settings
from influence import InfluenceField
from palettes import Palette

try:  # Optional CUDA backend.
    import cupy as cp  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    cp = None


class ArrayBackend:
    """Tiny adapter that lets the renderer use NumPy or CuPy."""

    def __init__(self, prefer_cuda: bool = False) -> None:
        self.xp = np
        self.use_cuda = False
        self.label = "CPU / NumPy"

        if not prefer_cuda or cp is None:
            return

        try:
            device_count = cp.cuda.runtime.getDeviceCount()
            if device_count > 0:
                self.xp = cp
                self.use_cuda = True
                self.label = f"CUDA / CuPy ({device_count} device{'s' if device_count != 1 else ''})"
        except Exception:
            # Silently fall back to CPU when CuPy is installed but not usable.
            self.xp = np
            self.use_cuda = False
            self.label = "CPU / NumPy"

    def to_device(self, array: np.ndarray):
        if self.use_cuda:
            return cp.asarray(array, dtype=cp.float32)
        return np.asarray(array, dtype=np.float32)

    def to_host(self, array):
        if self.use_cuda:
            return cp.asnumpy(array)
        return np.asarray(array)


class WaterfallRenderer:
    def __init__(
        self,
        settings: Settings,
        render_size: tuple[int, int],
        palette: Palette,
        text_palette: Palette | None = None,
        seed: int = 7,
        prefer_cuda: bool = False,
    ) -> None:
        self.settings = settings
        self.width, self.height = render_size
        self.palette = palette
        self.text_palette = text_palette or palette
        self.seed = seed
        self.backend = ArrayBackend(prefer_cuda or settings.prefer_cuda)
        self.xp = self.backend.xp

        xp = self.xp
        x = xp.linspace(0.0, 1.0, self.width, dtype=xp.float32)
        y = xp.linspace(0.0, 1.0, self.height, dtype=xp.float32)
        self.x, self.y = xp.meshgrid(x, y)
        self.curtain = 0.82 + 0.18 * xp.exp(-((self.x - 0.5) ** 2) / 0.18)

        rng = np.random.default_rng(seed)
        grain = rng.random((self.height, self.width), dtype=np.float32)
        self.grain = self.backend.to_device(grain)
        self._update_palette_arrays()
        self._update_text_palette_arrays()

    def _color(self, rgb: tuple[int, int, int]):
        xp = self.xp
        return xp.array(rgb, dtype=xp.float32) / 255.0

    def _update_palette_arrays(self) -> None:
        self.bg_top = self._color(self.palette.bg_top)
        self.bg_bottom = self._color(self.palette.bg_bottom)
        self.water_dark = self._color(self.palette.water_dark)
        self.water_mid = self._color(self.palette.water_mid)
        self.water_bright = self._color(self.palette.water_bright)
        self.code_low = self._color(self.palette.code_low)
        self.code_high = self._color(self.palette.code_high)
        self.mist = self._color(self.palette.mist)
        self.base_background = (self.bg_top * (1.0 - self.y[..., None])) + (self.bg_bottom * self.y[..., None])

    def _update_text_palette_arrays(self) -> None:
        self.text_low = self._color(self.text_palette.code_low)
        self.text_high = self._color(self.text_palette.code_high)
        self.text_mist = self._color(self.text_palette.mist)
        self.text_edge = (self.text_high * 0.66) + (self.text_mist * 0.34)

    def set_palette(self, palette: Palette) -> None:
        self.palette = palette
        self._update_palette_arrays()

    def set_text_palette(self, palette: Palette) -> None:
        self.text_palette = palette
        self._update_text_palette_arrays()

    def _fields_to_backend(self, fields: InfluenceField):
        return {
            "core": self.backend.to_device(fields.core),
            "spread": self.backend.to_device(fields.spread),
            "halo": self.backend.to_device(fields.halo),
            "gx": self.backend.to_device(fields.gx),
            "gy": self.backend.to_device(fields.gy),
            "edge": self.backend.to_device(fields.edge),
            "pressure": self.backend.to_device(fields.pressure),
        }

    def render_text_layer(self, fields: InfluenceField, t: float) -> np.ndarray:
        """Render just the scrolling code/halo layer on black.

        This is used when the procedural waterfall background is toggled off so
        the source code remains visible underneath the reactive droplets.
        """
        xp = self.xp
        f = self._fields_to_backend(fields)
        core = f["core"]
        halo = f["halo"]
        edge = f["edge"]
        pressure = f["pressure"]

        time = xp.float32(t * self.settings.water_speed)
        x = self.x
        y = self.y

        pulse = 0.5 + (0.5 * xp.sin((4.1 * time) + (12.0 * y) + (1.4 * x)))
        halo_wave = 0.5 + (0.5 * xp.sin((18.0 * x) + (36.0 * y) - (6.2 * time)))
        edge_wave = 0.5 + (0.5 * xp.sin((24.0 * x) + (54.0 * y) - (7.0 * time)))

        text_mix = xp.clip(0.24 + (0.36 * core) + (0.24 * pressure) + (0.16 * pulse), 0.0, 1.0)
        text_color = (self.text_low * (1.0 - text_mix[..., None])) + (self.text_high * text_mix[..., None])

        rgb = xp.zeros_like(self.base_background)
        rgb += halo[..., None] * self.text_low * (0.12 + (0.16 * halo_wave[..., None]))
        rgb += halo[..., None] * self.text_mist * (0.05 + (0.06 * halo_wave[..., None]))
        rgb += edge[..., None] * self.text_edge * (0.03 + (0.05 * edge_wave[..., None]))
        rgb += core[..., None] * text_color * 0.96

        rgb = xp.clip(rgb, 0.0, 1.0)
        return (self.backend.to_host(rgb) * 255.0).astype(np.uint8)

    def render(self, fields: InfluenceField, t: float) -> np.ndarray:
        xp = self.xp
        f = self._fields_to_backend(fields)
        core = f["core"]
        spread = f["spread"]
        halo = f["halo"]
        gx = f["gx"]
        gy = f["gy"]
        edge = f["edge"]
        pressure = f["pressure"]

        time = xp.float32(t * self.settings.water_speed)
        x = self.x
        y = self.y

        # Base downward motion with multiple scales of side-sway.
        warp = (
            self.settings.side_sway * xp.sin((6.0 * y) + (1.8 * time))
            + 0.035 * xp.sin((18.0 * y) - (0.8 * time) + 2.8 * xp.sin((4.4 * x) + (0.24 * time)))
            + 0.020 * xp.sin((44.0 * y) + (15.0 * x) - (1.22 * time))
        )

        # The code mask bends the flow more strongly than before, so the letters
        # feel like submerged obstacles rather than a purely cosmetic overlay.
        text_warp = self.settings.text_distortion * (
            (0.82 * gx)
            + (0.56 * spread * xp.sin((36.0 * y) - (5.7 * time) + (28.0 * x)))
            + (0.22 * edge * xp.sin((54.0 * y) + (18.0 * x) - (7.1 * time)))
        )
        refraction = self.settings.letter_refraction * pressure * (
            (0.70 * gx * xp.cos((28.0 * y) - (6.2 * time)))
            - (0.48 * gy * xp.sin((24.0 * x) + (4.9 * time)))
        )

        xw = x + warp + text_warp + refraction
        yw = y + (self.settings.text_distortion * 0.42 * gy) - (0.012 * pressure * xp.sin((18.0 * x) - (6.0 * time)))
        flow = (yw * self.settings.vertical_density) - time

        ribbon_phase = (xw * self.settings.ribbon_frequency) + (0.72 * xp.sin(flow * 1.25))
        thread_phase = (xw * self.settings.thread_frequency) - (0.44 * xp.sin((flow * 2.1) + (3.2 * x)))
        micro_phase = (xw * self.settings.micro_frequency) + (0.24 * xp.sin((flow * 4.3) - (2.2 * x)))

        ribbons = xp.exp(-self.settings.ribbon_sharpness * xp.sin(np.pi * ribbon_phase) ** 2)
        threads = xp.exp(-self.settings.thread_sharpness * xp.sin(np.pi * thread_phase) ** 2)
        micro = xp.exp(-self.settings.micro_sharpness * xp.sin(np.pi * micro_phase) ** 2)

        mist_wave = 0.5 + (0.5 * xp.sin((flow * 1.55) + 4.8 * xp.sin((4.0 * xw) + (1.35 * time))))
        foam_wave = 0.5 + (0.5 * xp.sin((35.0 * yw) - (7.0 * time) + 3.0 * xp.sin(10.0 * xw)))
        shimmer = 0.5 + (0.5 * xp.sin((22.0 * xw) + (48.0 * yw) - (6.8 * time)))
        collision_wake = pressure * (0.55 + 0.45 * xp.sin((22.0 * x) + (54.0 * y) - (8.6 * time)))
        agitation = spread * (0.55 + 0.45 * xp.sin((18.0 * x) + (48.0 * y) - (7.4 * time)))
        turbulence = edge * (0.5 + 0.5 * xp.sin((26.0 * xw) + (30.0 * yw) - (8.5 * time)))

        intensity = (
            0.10
            + (0.54 * ribbons)
            + (0.34 * threads * (0.55 + (0.45 * mist_wave)))
            + (0.15 * micro * (0.25 + (0.75 * foam_wave)))
            + (0.06 * mist_wave)
        )
        intensity *= self.curtain
        intensity += self.settings.text_reaction * agitation
        intensity += self.settings.wake_strength * collision_wake
        intensity += self.settings.edge_sparkle * turbulence
        intensity += 0.06 * shimmer
        intensity = xp.clip(intensity, 0.0, 1.38)

        grad = xp.clip(intensity, 0.0, 1.0)
        highlight = xp.clip((grad - 0.40) * 2.0, 0.0, 1.0)

        water = (self.water_dark * (1.0 - grad[..., None])) + (self.water_mid * grad[..., None])
        water = (water * (1.0 - highlight[..., None])) + (self.water_bright * highlight[..., None])

        rgb = self.base_background * (1.0 - (0.88 * grad[..., None]))
        rgb += water * (0.28 + (0.95 * grad[..., None]))
        rgb += halo[..., None] * self.mist * (0.08 + (0.10 * mist_wave[..., None]))

        # Sink the text slightly into the water, but flare the edges and local halo.
        text_mix = (
            0.20
            + (0.40 * grad)
            + (0.24 * pressure)
            + (0.16 * (0.5 + 0.5 * xp.sin((4.1 * time) + (12.0 * y) + (1.4 * x))))
        )
        text_color = (self.text_low * (1.0 - text_mix[..., None])) + (self.text_high * text_mix[..., None])
        rgb *= 1.0 - (0.10 * core[..., None])
        rgb = (rgb * (1.0 - (0.40 * core[..., None]))) + (text_color * core[..., None] * 0.92)
        rgb += halo[..., None] * self.text_high * (0.06 + (0.08 * collision_wake[..., None]))
        rgb += edge[..., None] * self.text_edge * (0.03 + (0.06 * shimmer[..., None]))

        if self.settings.grain_amount > 0.0:
            shift_y = int((t * self.settings.water_speed) * 28.0) % self.height
            shift_x = int((t * self.settings.water_speed) * 9.0) % self.width
            grain = xp.roll(self.grain, shift=(shift_y, shift_x), axis=(0, 1))
            rgb += (grain[..., None] - 0.5) * self.settings.grain_amount

        rgb = xp.clip(rgb, 0.0, 1.0)
        return (self.backend.to_host(rgb) * 255.0).astype(np.uint8)
