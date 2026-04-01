from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(slots=True)
class FlowBand:
    frequency: float
    phase: float
    amplitude: float

    def sample(self, x: float, y: float, t: float) -> float:
        angle = (x * self.frequency) + (y * 0.9) - (t * 0.6) + self.phase
        return math.sin(angle) * self.amplitude


class FlowField:
    def __init__(self) -> None:
        self.bands = [
            FlowBand(frequency=6.0, phase=0.1, amplitude=0.9),
            FlowBand(frequency=11.0, phase=1.8, amplitude=0.4),
            FlowBand(frequency=17.0, phase=2.7, amplitude=0.2),
        ]

    def sway(self, x: float, y: float, t: float) -> float:
        total = 0.0
        for band in self.bands:
            total += band.sample(x, y, t)
        return total

    def velocity(self, x: float, y: float, t: float) -> tuple[float, float]:
        sx = self.sway(x, y, t)
        drift = 0.25 * math.sin((y * 9.0) + (t * 1.4))
        return sx * 0.08, 0.9 + drift


def ribbon(x: float, sharpness: float = 8.0) -> float:
    wave = math.sin(x * math.pi)
    return math.exp(-(wave * wave) * sharpness)


def brightness(x: float, y: float, t: float) -> float:
    field = FlowField()
    vx, vy = field.velocity(x, y, t)
    band_a = ribbon((x + vx) * 14.0 + math.sin((y * 4.0) - t * 0.8))
    band_b = ribbon((x - vx) * 28.0 + math.sin((y * 8.0) - t * 1.2), sharpness=13.0)
    mist = 0.5 + 0.5 * math.sin((y * 13.0) - (t * 2.0) + x * 3.0)
    return 0.1 + (band_a * 0.55) + (band_b * 0.25) + (mist * 0.1) + vy * 0.05


if __name__ == "__main__":
    values = []
    for step in range(12):
        t = step / 12.0
        values.append(brightness(0.37, 0.62, t))
    print([round(v, 3) for v in values])
