from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pygame

from code_scroller import CodeScroller
from config import Settings
from palettes import get_palette, palette_names
from waterfall import WaterfallRenderer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Waterfall-inspired generative art with scrolling code.")
    parser.add_argument("--fullscreen", action="store_true", help="Start in fullscreen mode.")
    parser.add_argument("--width", type=int, default=1600, help="Window width in windowed mode.")
    parser.add_argument("--height", type=int, default=900, help="Window height in windowed mode.")
    parser.add_argument("--palette", default="midnight_ice", choices=palette_names(), help="Palette preset.")
    parser.add_argument("--samples", default="code_samples", help="Directory of code files to scroll.")
    parser.add_argument("--render-height", type=int, default=270, help="Off-screen render height. Higher = sharper, slower.")
    return parser.parse_args()


def create_screen(fullscreen: bool, windowed_size: tuple[int, int]) -> pygame.Surface:
    flags = pygame.FULLSCREEN if fullscreen else pygame.RESIZABLE
    size = pygame.display.get_desktop_sizes()[0] if fullscreen else windowed_size
    try:
        screen = pygame.display.set_mode(size, flags, vsync=1)
    except TypeError:
        screen = pygame.display.set_mode(size, flags)
    return screen


def resolve_samples_path(settings: Settings, raw: str) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        path = settings.project_path(raw)
    return path


def build_scene(
    settings: Settings,
    output_size: tuple[int, int],
    samples_dir: Path,
    palette_name: str,
    seed: int,
) -> tuple[CodeScroller, WaterfallRenderer, pygame.Surface, pygame.Surface]:
    render_size = settings.render_size(output_size)
    scroller = CodeScroller(settings, render_size, samples_dir)
    renderer = WaterfallRenderer(settings, render_size, get_palette(palette_name), seed=seed)
    sim_surface = pygame.Surface(render_size).convert()
    scaled_surface = pygame.Surface(output_size).convert()
    return scroller, renderer, sim_surface, scaled_surface


def save_screenshot(screen: pygame.Surface, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = output_dir / f"waterfall_code_{timestamp}.png"
    pygame.image.save(screen, str(path))
    return path


def print_controls() -> None:
    print("Controls:")
    print("  Esc / Q  quit")
    print("  Space    pause / resume")
    print("  R        restart current file and reset the animation clock")
    print("  N / P    next / previous code file")
    print("  C / Tab  cycle palette")
    print("  F        toggle fullscreen")
    print("  L        reload the samples directory")
    print("  S        save a screenshot")


def main() -> None:
    args = parse_args()

    settings = Settings(
        window_size=(args.width, args.height),
        start_fullscreen=args.fullscreen,
        render_height=args.render_height,
        palette_name=args.palette,
    )

    pygame.init()
    pygame.display.set_caption(settings.title)

    if settings.show_startup_controls:
        print_controls()

    fullscreen = settings.start_fullscreen
    windowed_size = settings.window_size
    screen = create_screen(fullscreen, windowed_size)

    samples_dir = resolve_samples_path(settings, args.samples)
    palette_order = palette_names()
    palette_index = palette_order.index(settings.palette_name)

    scroller, renderer, sim_surface, scaled_surface = build_scene(
        settings,
        screen.get_size(),
        samples_dir,
        palette_order[palette_index],
        settings.default_seed,
    )

    clock = pygame.Clock()
    running = True
    paused = False
    time_accumulator = 0.0

    if fullscreen:
        pygame.mouse.set_visible(False)

    while running:
        dt = clock.tick(settings.fps) / 1000.0

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

            elif event.type == pygame.VIDEORESIZE and not fullscreen:
                windowed_size = (event.w, event.h)
                screen = create_screen(False, windowed_size)
                scroller, renderer, sim_surface, scaled_surface = build_scene(
                    settings,
                    screen.get_size(),
                    samples_dir,
                    palette_order[palette_index],
                    settings.default_seed,
                )

            elif event.type == pygame.KEYDOWN:
                if event.key in (pygame.K_ESCAPE, pygame.K_q):
                    running = False
                elif event.key == pygame.K_SPACE:
                    paused = not paused
                elif event.key == pygame.K_r:
                    time_accumulator = 0.0
                    scroller.restart_current()
                elif event.key == pygame.K_n:
                    scroller.next_file()
                elif event.key == pygame.K_p:
                    scroller.previous_file()
                elif event.key in (pygame.K_c, pygame.K_TAB):
                    palette_index = (palette_index + 1) % len(palette_order)
                    renderer.set_palette(get_palette(palette_order[palette_index]))
                elif event.key == pygame.K_l:
                    scroller.reload_directory()
                elif event.key == pygame.K_s:
                    path = save_screenshot(screen, settings.screenshot_path())
                    print(f"Saved screenshot: {path}")
                elif event.key == pygame.K_f:
                    fullscreen = not fullscreen
                    screen = create_screen(fullscreen, windowed_size)
                    pygame.mouse.set_visible(not fullscreen)
                    scroller, renderer, sim_surface, scaled_surface = build_scene(
                        settings,
                        screen.get_size(),
                        samples_dir,
                        palette_order[palette_index],
                        settings.default_seed,
                    )

        if not paused:
            time_accumulator += dt
            scroller.update(dt)

        mask = scroller.render_mask(time_accumulator)
        frame = renderer.render(mask, time_accumulator)
        pygame.surfarray.blit_array(sim_surface, frame.swapaxes(0, 1))
        pygame.transform.smoothscale(sim_surface, screen.get_size(), scaled_surface)
        screen.blit(scaled_surface, (0, 0))
        pygame.display.flip()

    pygame.quit()


if __name__ == "__main__":
    main()
