from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class Vec2:
    x: float
    y: float

    def __add__(self, other: "Vec2") -> "Vec2":
        return Vec2(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Vec2") -> "Vec2":
        return Vec2(self.x - other.x, self.y - other.y)

    def scale(self, amount: float) -> "Vec2":
        return Vec2(self.x * amount, self.y * amount)

    def tuple(self) -> tuple[float, float]:
        return self.x, self.y


class RibbonCache:
    def __init__(self, width: int) -> None:
        self.width = width
        self.points = [Vec2(index / max(1, width - 1), 0.0) for index in range(width)]

    def disturb(self, force: float) -> None:
        for index, point in enumerate(self.points):
            offset = ((index % 7) - 3) * 0.0015
            self.points[index] = Vec2(point.x + offset * force, point.y)

    def fall(self, amount: float) -> None:
        for index, point in enumerate(self.points):
            self.points[index] = Vec2(point.x, point.y + amount)

    def sample(self, index: int) -> Vec2:
        return self.points[index % len(self.points)]


if __name__ == "__main__":
    cache = RibbonCache(9)
    cache.disturb(0.5)
    cache.fall(0.04)
    print([cache.sample(i).tuple() for i in range(9)])
