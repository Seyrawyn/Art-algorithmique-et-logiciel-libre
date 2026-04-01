from __future__ import annotations

import math
import random
from dataclasses import dataclass

import pygame

from config import Settings
from influence import InfluenceField
from palettes import Palette


@dataclass(slots=True)
class Particle:
    x: float
    y: float
    vx: float
    vy: float
    size: float
    color: tuple[int, int, int]
    alpha: int
    base_vy: float
    phase: float
    age: float = 0.0
    cooldown: float = 0.0
    glow: float = 0.0
    prev_x: float = 0.0
    prev_y: float = 0.0


class ReactiveDropletLayer:
    """A droplet/spray layer that bounces off the scrolling code contours.

    The base waterfall remains procedural and smooth, while this layer adds the
    more literal droplet behaviour that was stronger in the test-claude version.
    """

    def __init__(
        self,
        settings: Settings,
        render_size: tuple[int, int],
        palette: Palette,
        seed: int = 7,
    ) -> None:
        self.settings = settings
        self.width, self.height = render_size
        self.scale = self.height / 270.0
        self.palette = palette
        self.rng = random.Random(seed)
        self.time = 0.0
        self.spawn_budget = 0.0
        self.foam_opacity = float(settings.foam_opacity)

        self.drops: list[Particle] = []
        self.splashes: list[Particle] = []
        self.mist: list[Particle] = []
        self._circle_cache: dict[tuple[int, tuple[int, int, int], int], pygame.Surface] = {}

        self._init_mist()

    # ------------------------------------------------------------------
    # Palette / init helpers
    # ------------------------------------------------------------------
    def _mix(self, a: tuple[int, int, int], b: tuple[int, int, int], t: float) -> tuple[int, int, int]:
        return (
            int(a[0] * (1.0 - t) + b[0] * t),
            int(a[1] * (1.0 - t) + b[1] * t),
            int(a[2] * (1.0 - t) + b[2] * t),
        )

    def _water_colors(self) -> list[tuple[int, int, int]]:
        return [
            self.palette.water_mid,
            self.palette.water_bright,
            self._mix(self.palette.water_mid, self.palette.water_bright, 0.5),
            self._mix(self.palette.water_dark, self.palette.water_mid, 0.75),
        ]

    def set_palette(self, palette: Palette) -> None:
        self.palette = palette
        self._circle_cache.clear()
        recolor_bank = self._water_colors()
        for particle in self.drops:
            particle.color = self.rng.choice(recolor_bank)
        for particle in self.splashes:
            particle.color = self.rng.choice(recolor_bank)
        mist_color = self.palette.mist
        for particle in self.mist:
            particle.color = mist_color

    def set_foam_opacity(self, value: float) -> float:
        cfg = self.settings
        self.foam_opacity = max(cfg.foam_opacity_min, min(cfg.foam_opacity_max, float(value)))
        return self.foam_opacity

    def adjust_foam_opacity(self, delta: float) -> float:
        return self.set_foam_opacity(self.foam_opacity + delta)

    def reset(self) -> None:
        self.time = 0.0
        self.spawn_budget = 0.0
        self.drops.clear()
        self.splashes.clear()
        self.mist.clear()
        self._init_mist()

    def _init_mist(self) -> None:
        cfg = self.settings
        mist_count = max(0, int(round(cfg.mist_particles * self.scale)))
        for _ in range(mist_count):
            self.mist.append(
                Particle(
                    x=self.rng.uniform(0.0, self.width),
                    y=self.rng.uniform(0.0, self.height),
                    vx=self.rng.uniform(-cfg.mist_drift, cfg.mist_drift) * 0.12 * self.scale,
                    vy=self.rng.uniform(-cfg.mist_drift, cfg.mist_drift) * 0.06 * self.scale,
                    size=self.rng.uniform(1.8, cfg.mist_max_size) * self.scale,
                    color=self.palette.mist,
                    alpha=self.rng.randint(4, max(5, cfg.mist_alpha)),
                    base_vy=0.0,
                    phase=self.rng.uniform(0.0, math.tau),
                )
            )

    # ------------------------------------------------------------------
    # Spawning
    # ------------------------------------------------------------------
    def _make_drop(self) -> Particle:
        cfg = self.settings
        center = self.width * 0.5
        spread = self.width * 0.32
        x = max(-6.0, min(self.width + 6.0, self.rng.gauss(center, spread)))
        size = self.rng.uniform(cfg.droplet_min_size, cfg.droplet_max_size) * self.scale
        base_vy = self.rng.uniform(cfg.droplet_min_speed, cfg.droplet_max_speed) * self.scale
        color = self.rng.choice(self._water_colors())
        alpha = self.rng.randint(cfg.droplet_alpha_min, cfg.droplet_alpha_max)
        vx = self.rng.uniform(-6.0, 6.0) * self.scale
        vy = base_vy * self.rng.uniform(0.92, 1.08)
        y = self.rng.uniform(-22.0, 2.0) * self.scale
        phase = self.rng.uniform(0.0, math.tau)
        return Particle(
            x=x,
            y=y,
            vx=vx,
            vy=vy,
            size=size,
            color=color,
            alpha=alpha,
            base_vy=base_vy,
            phase=phase,
            prev_x=x,
            prev_y=y,
        )

    def spawn(self, dt: float) -> None:
        if not self.settings.enable_droplets:
            return
        if len(self.drops) >= self.settings.droplet_max_count:
            return

        self.spawn_budget += self.settings.droplet_spawn_rate * dt
        spawn_count = int(self.spawn_budget)
        if spawn_count <= 0:
            return

        self.spawn_budget -= spawn_count
        room = self.settings.droplet_max_count - len(self.drops)
        spawn_count = min(spawn_count, room)
        for _ in range(spawn_count):
            self.drops.append(self._make_drop())

    def _spawn_splash(self, x: float, y: float, nx: float, ny: float, energy: float) -> None:
        cfg = self.settings
        tangent_x, tangent_y = -ny, nx
        count = self.rng.randint(cfg.droplet_splash_count_min, cfg.droplet_splash_count_max)
        count = max(1, count)

        for _ in range(count):
            outward = self.rng.uniform(0.45, 1.05) * cfg.droplet_splash_speed * self.scale * (0.45 + energy)
            tangent = self.rng.uniform(-0.95, 0.95) * cfg.droplet_splash_speed * 0.75 * self.scale
            vx = (nx * outward) + (tangent_x * tangent)
            vy = (ny * outward) + (tangent_y * tangent) - (cfg.droplet_upward_boost * 0.25 * self.scale)
            size = self.rng.uniform(0.45, 1.35) * self.scale
            alpha = self.rng.randint(120, 235)
            color = self.rng.choice(self._water_colors())
            phase = self.rng.uniform(0.0, math.tau)
            self.splashes.append(
                Particle(
                    x=x + self.rng.uniform(-1.8, 1.8) * self.scale,
                    y=y + self.rng.uniform(-1.2, 1.2) * self.scale,
                    vx=vx,
                    vy=vy,
                    size=size,
                    color=color,
                    alpha=alpha,
                    base_vy=0.0,
                    phase=phase,
                    prev_x=x,
                    prev_y=y,
                )
            )

    # ------------------------------------------------------------------
    # Update
    # ------------------------------------------------------------------
    def update(self, fields: InfluenceField, dt: float, t: float) -> None:
        self.time = t
        self.spawn(dt)

        cfg = self.settings
        wind = math.sin((t * 0.68) + 0.4) * cfg.droplet_wind * self.scale
        gravity = 16.0 * self.scale

        alive_drops: list[Particle] = []
        for p in self.drops:
            p.prev_x, p.prev_y = p.x, p.y
            p.age += dt
            p.cooldown = max(0.0, p.cooldown - dt)
            p.glow = max(0.0, p.glow - (dt * 2.6))

            # Relax back toward a slightly wobbly downward motion.
            target_vx = wind + math.sin((p.phase + (t * 5.8) + (p.y * 0.085))) * cfg.droplet_wobble * 0.4 * self.scale
            p.vx += (target_vx - p.vx) * min(1.0, dt * 1.8)
            p.vy += (p.base_vy - p.vy) * min(1.0, dt * 1.1)
            p.vy += gravity * dt

            ix = int(p.x)
            iy = int(p.y)
            if 0 <= ix < self.width and 0 <= iy < self.height:
                core = float(fields.core[iy, ix])
                spread = float(fields.spread[iy, ix])
                edge = float(fields.edge[iy, ix])
                pressure = float(fields.pressure[iy, ix])
                nx = float(fields.nx[iy, ix])
                ny = float(fields.ny[iy, ix])

                if pressure > 0.015:
                    tangent_x, tangent_y = -ny, nx
                    p.vx += nx * cfg.droplet_repel_strength * pressure * dt * self.scale
                    p.vy += ny * cfg.droplet_repel_strength * pressure * 0.72 * dt * self.scale

                    swirl = math.sin((p.y * 0.22) + (t * 8.5) + p.phase)
                    p.vx += tangent_x * cfg.droplet_tangent_strength * spread * swirl * dt * self.scale
                    p.vy += tangent_y * cfg.droplet_tangent_strength * spread * swirl * 0.55 * dt * self.scale

                    drag = 1.0 - (cfg.droplet_drag_in_text * core * min(1.0, dt * 4.0))
                    p.vx *= drag
                    p.vy *= max(0.55, drag)

                    # Bounce against the contour normal when moving into the text body.
                    dot = (p.vx * nx) + (p.vy * ny)
                    hit = pressure > 0.14 and p.cooldown <= 0.0 and dot < -(5.0 * self.scale)
                    if hit:
                        restitution = cfg.droplet_bounce_restitution * (0.6 + (0.4 * edge))
                        p.vx -= (1.0 + restitution) * dot * nx
                        p.vy -= (1.0 + restitution) * dot * ny
                        p.vy -= cfg.droplet_upward_boost * (0.24 + pressure) * self.scale
                        p.x += nx * cfg.droplet_separation * (0.4 + pressure) * self.scale
                        p.y += ny * cfg.droplet_separation * (0.4 + pressure) * self.scale
                        p.glow = max(p.glow, 1.0)
                        p.cooldown = cfg.droplet_collision_cooldown

                        splash_chance = cfg.droplet_splash_probability * min(1.0, 0.35 + pressure)
                        if self.rng.random() < splash_chance:
                            self._spawn_splash(p.x, p.y, nx, ny, pressure)

            p.x += p.vx * dt
            p.y += p.vy * dt

            if -14.0 * self.scale <= p.x <= self.width + (14.0 * self.scale) and p.y <= self.height + (20.0 * self.scale):
                alive_drops.append(p)

        self.drops = alive_drops

        # Spray / splash particles
        splash_alive: list[Particle] = []
        for p in self.splashes:
            p.prev_x, p.prev_y = p.x, p.y
            p.age += dt
            p.vy += 48.0 * self.scale * dt
            p.x += p.vx * dt
            p.y += p.vy * dt
            p.alpha = max(0, p.alpha - int(220 * dt))
            if p.alpha > 0 and p.age < 0.9 and -8 <= p.x <= self.width + 8 and p.y <= self.height + 12:
                splash_alive.append(p)
        self.splashes = splash_alive

        # Mist behind the droplets.
        for p in self.mist:
            p.x += (p.vx + wind * 0.12) * dt
            p.y += p.vy * dt
            p.phase += dt * 0.6
            p.vx += math.sin(p.phase) * cfg.mist_drift * 0.02 * self.scale * dt
            if p.x < -18.0:
                p.x = self.width + 18.0
            elif p.x > self.width + 18.0:
                p.x = -18.0
            if p.y < -18.0:
                p.y = self.height + 18.0
            elif p.y > self.height + 18.0:
                p.y = -18.0

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------
    def _get_circle(self, radius: int, color: tuple[int, int, int], alpha: int) -> pygame.Surface:
        key = (radius, color, alpha)
        cached = self._circle_cache.get(key)
        if cached is not None:
            return cached

        size = radius * 2 + 4
        surf = pygame.Surface((size, size), pygame.SRCALPHA)
        pygame.draw.circle(surf, (*color, alpha), (size // 2, size // 2), radius)
        self._circle_cache[key] = surf
        if len(self._circle_cache) > 1800:
            # Remove a chunk of oldest-insertion keys. In CPython dicts preserve insertion order.
            for old_key in list(self._circle_cache.keys())[:450]:
                self._circle_cache.pop(old_key, None)
        return surf

    def _blit_circle(
        self,
        surface: pygame.Surface,
        x: float,
        y: float,
        radius: float,
        color: tuple[int, int, int],
        alpha: int,
        *,
        additive: bool = True,
    ) -> None:
        alpha = max(0, min(255, alpha))
        if alpha <= 0:
            return
        r = max(1, int(radius))
        circ = self._get_circle(r, color, alpha)
        dest = (int(x) - (circ.get_width() // 2), int(y) - (circ.get_height() // 2))
        if additive:
            surface.blit(circ, dest, special_flags=pygame.BLEND_ADD)
        else:
            surface.blit(circ, dest)

    def draw(self, surface: pygame.Surface) -> None:
        # Mist first.
        for p in self.mist:
            self._blit_circle(surface, p.x, p.y, p.size, p.color, p.alpha)

        foam_scale = max(0.0, self.foam_opacity)

        # Main droplets with short trails.
        for p in self.drops:
            dx = p.x - p.prev_x
            dy = p.y - p.prev_y
            trail_len = math.hypot(dx, dy)
            trail_steps = min(3, max(1, int(trail_len / max(1.3, p.size * 1.7))))

            if p.glow > 0.0 and foam_scale > 0.0:
                glow_alpha = int(72 * p.glow * foam_scale)
                glow_radius = p.size * (1.8 + (0.6 * p.glow))
                glow_color = self._mix(self.palette.water_bright, self.palette.mist, 0.35)
                self._blit_circle(surface, p.x, p.y, glow_radius, glow_color, glow_alpha, additive=False)

            for step in range(trail_steps, 0, -1):
                frac = step / (trail_steps + 1)
                tx = p.x - (dx * frac)
                ty = p.y - (dy * frac)
                trail_alpha = int(p.alpha * (0.12 + 0.34 * (1.0 - frac)))
                trail_radius = max(1.0, p.size * (0.65 + (0.22 * (1.0 - frac))))
                self._blit_circle(surface, tx, ty, trail_radius, p.color, trail_alpha)

            self._blit_circle(surface, p.x, p.y, p.size, p.color, p.alpha)

        # Splash / foam particles on top, drawn with standard alpha so the code stays visible.
        if foam_scale > 0.0:
            for p in self.splashes:
                foam_alpha = int(p.alpha * foam_scale)
                foam_radius = p.size * 1.08
                foam_color = self._mix(p.color, self.palette.mist, 0.42)
                self._blit_circle(surface, p.x, p.y, foam_radius, foam_color, foam_alpha, additive=False)
