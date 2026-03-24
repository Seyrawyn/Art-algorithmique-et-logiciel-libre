from __future__ import annotations

from dataclasses import dataclass

import numpy as np


Color = tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class Palette:
    name: str
    bg_top: Color
    bg_bottom: Color
    water_dark: Color
    water_mid: Color
    water_bright: Color
    code_low: Color
    code_high: Color
    mist: Color

    def as_float(self, color: Color) -> np.ndarray:
        return np.array(color, dtype=np.float32) / 255.0


PALETTES: dict[str, Palette] = {
    "midnight_ice": Palette(
        name="midnight_ice",
        bg_top=(3, 6, 14),
        bg_bottom=(0, 0, 0),
        water_dark=(8, 30, 48),
        water_mid=(48, 140, 196),
        water_bright=(228, 247, 255),
        code_low=(60, 128, 162),
        code_high=(236, 250, 255),
        mist=(170, 220, 255),
    ),
    "ink_silver": Palette(
        name="ink_silver",
        bg_top=(10, 10, 12),
        bg_bottom=(0, 0, 0),
        water_dark=(28, 28, 34),
        water_mid=(124, 136, 150),
        water_bright=(242, 246, 250),
        code_low=(108, 116, 126),
        code_high=(250, 252, 255),
        mist=(190, 196, 205),
    ),
    "deep_teal": Palette(
        name="deep_teal",
        bg_top=(1, 10, 12),
        bg_bottom=(0, 0, 0),
        water_dark=(5, 44, 44),
        water_mid=(35, 180, 170),
        water_bright=(230, 255, 250),
        code_low=(65, 158, 150),
        code_high=(238, 255, 250),
        mist=(120, 230, 216),
    ),
}


def palette_names() -> list[str]:
    return list(PALETTES.keys())


def get_palette(name: str) -> Palette:
    return PALETTES.get(name, PALETTES["midnight_ice"])
