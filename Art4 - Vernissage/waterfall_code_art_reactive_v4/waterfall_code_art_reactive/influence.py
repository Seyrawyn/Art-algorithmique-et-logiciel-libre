from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class InfluenceField:
    """Text-derived fields used by both the waterfall and the droplet layer."""

    core: np.ndarray
    spread: np.ndarray
    halo: np.ndarray
    gx: np.ndarray
    gy: np.ndarray
    edge: np.ndarray
    pressure: np.ndarray
    nx: np.ndarray
    ny: np.ndarray


class InfluenceBuilder:
    """Build a small bundle of mask-derived scalar/vector fields.

    The raw code mask is expanded into:
      - `spread`: soft local area of influence
      - `edge`: gradient energy around letter contours
      - `pressure`: stronger collision field used by droplets
      - `nx`, `ny`: outward-facing normals used for bounce/reflection
    """

    def __init__(self, blur_passes: int = 4, edge_gain: float = 6.0) -> None:
        self.blur_passes = max(1, blur_passes)
        self.edge_gain = edge_gain

    def build(self, mask: np.ndarray) -> InfluenceField:
        core = np.clip(mask.astype(np.float32, copy=False), 0.0, 1.0)
        spread = core.copy()

        for _ in range(self.blur_passes):
            spread = (
                (spread * 4.0)
                + np.roll(spread, 1, axis=0)
                + np.roll(spread, -1, axis=0)
                + np.roll(spread, 1, axis=1)
                + np.roll(spread, -1, axis=1)
            ) / 8.0

        halo = np.clip(spread - (core * 0.55), 0.0, 1.0)
        gx = 0.5 * (np.roll(spread, -1, axis=1) - np.roll(spread, 1, axis=1))
        gy = 0.5 * (np.roll(spread, -1, axis=0) - np.roll(spread, 1, axis=0))
        edge = np.clip(np.sqrt((gx * gx) + (gy * gy)) * self.edge_gain, 0.0, 1.0)
        pressure = np.clip((core * 0.85) + (spread * 0.45) + (edge * 0.95), 0.0, 1.0)

        # Outward normals: the negative gradient points away from the text body.
        norm = np.sqrt((gx * gx) + (gy * gy)) + 1.0e-6
        nx = np.where(edge > 0.001, -gx / norm, 0.0).astype(np.float32, copy=False)
        ny = np.where(edge > 0.001, -gy / norm, 0.0).astype(np.float32, copy=False)

        return InfluenceField(
            core=core,
            spread=spread.astype(np.float32, copy=False),
            halo=halo.astype(np.float32, copy=False),
            gx=gx.astype(np.float32, copy=False),
            gy=gy.astype(np.float32, copy=False),
            edge=edge.astype(np.float32, copy=False),
            pressure=pressure.astype(np.float32, copy=False),
            nx=nx,
            ny=ny,
        )
