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
    "glacier_blue": Palette(
        name="glacier_blue",
        bg_top=(2, 8, 24),
        bg_bottom=(0, 0, 0),
        water_dark=(6, 24, 64),
        water_mid=(68, 156, 255),
        water_bright=(236, 245, 255),
        code_low=(92, 165, 255),
        code_high=(242, 248, 255),
        mist=(152, 204, 255),
    ),
    "aurora_mint": Palette(
        name="aurora_mint",
        bg_top=(1, 12, 14),
        bg_bottom=(0, 0, 0),
        water_dark=(8, 36, 32),
        water_mid=(58, 220, 180),
        water_bright=(236, 255, 248),
        code_low=(112, 212, 184),
        code_high=(244, 255, 250),
        mist=(154, 255, 226),
    ),
    "polar_gold": Palette(
        name="polar_gold",
        bg_top=(10, 8, 4),
        bg_bottom=(0, 0, 0),
        water_dark=(40, 28, 10),
        water_mid=(196, 152, 64),
        water_bright=(255, 243, 218),
        code_low=(210, 178, 98),
        code_high=(255, 248, 230),
        mist=(255, 220, 162),
    ),
    "cobalt_fog": Palette(
        name="cobalt_fog",
        bg_top=(4, 6, 18),
        bg_bottom=(0, 0, 0),
        water_dark=(12, 18, 58),
        water_mid=(92, 112, 220),
        water_bright=(236, 244, 255),
        code_low=(118, 136, 218),
        code_high=(246, 250, 255),
        mist=(176, 192, 255),
    ),
    "violet_noir": Palette(
        name="violet_noir",
        bg_top=(8, 2, 18),
        bg_bottom=(0, 0, 0),
        water_dark=(30, 10, 48),
        water_mid=(150, 78, 220),
        water_bright=(248, 240, 255),
        code_low=(182, 128, 236),
        code_high=(252, 246, 255),
        mist=(216, 178, 255),
    ),
    "ultraviolet_sleet": Palette(
        name="ultraviolet_sleet",
        bg_top=(6, 4, 18),
        bg_bottom=(0, 0, 0),
        water_dark=(18, 16, 54),
        water_mid=(132, 118, 255),
        water_bright=(242, 244, 255),
        code_low=(176, 170, 255),
        code_high=(250, 250, 255),
        mist=(206, 202, 255),
    ),
    "ember_lava": Palette(
        name="ember_lava",
        bg_top=(18, 6, 2),
        bg_bottom=(0, 0, 0),
        water_dark=(52, 18, 6),
        water_mid=(236, 96, 36),
        water_bright=(255, 236, 220),
        code_low=(236, 132, 78),
        code_high=(255, 245, 236),
        mist=(255, 182, 126),
    ),
    "sunset_magenta": Palette(
        name="sunset_magenta",
        bg_top=(14, 2, 12),
        bg_bottom=(0, 0, 0),
        water_dark=(48, 12, 34),
        water_mid=(224, 76, 148),
        water_bright=(255, 236, 248),
        code_low=(234, 120, 176),
        code_high=(255, 244, 250),
        mist=(255, 178, 218),
    ),
    "rose_quartz": Palette(
        name="rose_quartz",
        bg_top=(12, 8, 14),
        bg_bottom=(0, 0, 0),
        water_dark=(42, 24, 48),
        water_mid=(206, 150, 186),
        water_bright=(255, 244, 248),
        code_low=(224, 170, 202),
        code_high=(255, 248, 252),
        mist=(255, 208, 226),
    ),
    "amber_smoke": Palette(
        name="amber_smoke",
        bg_top=(14, 10, 6),
        bg_bottom=(0, 0, 0),
        water_dark=(48, 34, 16),
        water_mid=(186, 132, 68),
        water_bright=(255, 246, 226),
        code_low=(206, 158, 108),
        code_high=(255, 250, 236),
        mist=(246, 208, 152),
    ),
    "forest_mist": Palette(
        name="forest_mist",
        bg_top=(2, 10, 6),
        bg_bottom=(0, 0, 0),
        water_dark=(8, 34, 20),
        water_mid=(74, 180, 110),
        water_bright=(236, 255, 242),
        code_low=(118, 196, 148),
        code_high=(244, 255, 246),
        mist=(176, 240, 196),
    ),
    "storm_cyan": Palette(
        name="storm_cyan",
        bg_top=(2, 10, 16),
        bg_bottom=(0, 0, 0),
        water_dark=(8, 40, 58),
        water_mid=(46, 200, 232),
        water_bright=(232, 252, 255),
        code_low=(86, 208, 224),
        code_high=(242, 252, 255),
        mist=(166, 238, 255),
    ),
}


def palette_names() -> list[str]:
    return list(PALETTES.keys())


def get_palette(name: str) -> Palette:
    return PALETTES.get(name, PALETTES["midnight_ice"])
