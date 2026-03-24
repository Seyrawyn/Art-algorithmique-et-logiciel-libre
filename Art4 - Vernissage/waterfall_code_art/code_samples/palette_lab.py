from __future__ import annotations

from dataclasses import dataclass


Color = tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class Palette:
    background_top: Color
    background_bottom: Color
    water_dark: Color
    water_mid: Color
    water_bright: Color
    mist: Color


def blend(a: Color, b: Color, t: float) -> Color:
    return tuple(int((x * (1.0 - t)) + (y * t)) for x, y in zip(a, b))


MIDNIGHT_ICE = Palette(
    background_top=(3, 5, 12),
    background_bottom=(0, 0, 0),
    water_dark=(7, 28, 44),
    water_mid=(44, 132, 186),
    water_bright=(226, 248, 255),
    mist=(166, 220, 255),
)


INK_SILVER = Palette(
    background_top=(12, 12, 14),
    background_bottom=(0, 0, 0),
    water_dark=(28, 28, 34),
    water_mid=(120, 134, 148),
    water_bright=(244, 247, 250),
    mist=(194, 199, 205),
)


DEEP_TEAL = Palette(
    background_top=(1, 10, 12),
    background_bottom=(0, 0, 0),
    water_dark=(8, 44, 44),
    water_mid=(32, 180, 170),
    water_bright=(232, 255, 250),
    mist=(124, 228, 214),
)


def gradient_row(palette: Palette, steps: int = 8) -> list[Color]:
    row = []
    for index in range(steps):
        t = index / max(1, steps - 1)
        dark_mid = blend(palette.water_dark, palette.water_mid, t)
        row.append(blend(dark_mid, palette.water_bright, t * t))
    return row


if __name__ == "__main__":
    for name, palette in {
        "midnight_ice": MIDNIGHT_ICE,
        "ink_silver": INK_SILVER,
        "deep_teal": DEEP_TEAL,
    }.items():
        print(name, gradient_row(palette))
