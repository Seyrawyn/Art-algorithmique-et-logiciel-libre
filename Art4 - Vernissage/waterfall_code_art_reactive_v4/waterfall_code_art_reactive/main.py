from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pygame

from code_scroller import CodeScroller
from config import Settings
from droplets import ReactiveDropletLayer
from influence import InfluenceBuilder
from palettes import get_palette, palette_names
from waterfall import WaterfallRenderer


STATUS_DURATION = 1.7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Waterfall-inspired generative art with scrolling code.")
    parser.add_argument("--fullscreen", action="store_true", help="Start in fullscreen mode.")
    parser.add_argument("--width", type=int, default=1600, help="Window width in windowed mode.")
    parser.add_argument("--height", type=int, default=900, help="Window height in windowed mode.")
    parser.add_argument("--palette", default="midnight_ice", choices=palette_names(), help="Main waterfall / droplet palette preset.")
    parser.add_argument(
        "--text-palette",
        default=None,
        choices=palette_names(),
        help="Independent text color palette preset. Defaults to the main palette.",
    )
    parser.add_argument("--samples", default="code_samples", help="Directory of code files to scroll.")
    parser.add_argument("--render-height", type=int, default=320, help="Off-screen render height. Higher = sharper, slower.")
    parser.add_argument("--cuda", action="store_true", help="Use CuPy/CUDA for the waterfall field if available.")
    parser.add_argument("--no-droplets", action="store_true", help="Disable the reactive droplet layer.")
    parser.add_argument("--droplets-only", action="store_true", help="Start with the procedural waterfall disabled while keeping the scrolling code visible.")
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
    text_palette_name: str,
    seed: int,
) -> tuple[CodeScroller, InfluenceBuilder, WaterfallRenderer, ReactiveDropletLayer, pygame.Surface, pygame.Surface]:
    render_size = settings.render_size(output_size)
    palette = get_palette(palette_name)
    text_palette = get_palette(text_palette_name)
    scroller = CodeScroller(settings, render_size, samples_dir)
    influence = InfluenceBuilder(settings.blur_passes)
    renderer = WaterfallRenderer(
        settings,
        render_size,
        palette,
        text_palette=text_palette,
        seed=seed,
        prefer_cuda=settings.prefer_cuda,
    )
    droplets = ReactiveDropletLayer(settings, render_size, palette, seed=seed + 101)
    sim_surface = pygame.Surface(render_size).convert()
    scaled_surface = pygame.Surface(output_size).convert()
    return scroller, influence, renderer, droplets, sim_surface, scaled_surface


def save_screenshot(screen: pygame.Surface, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = output_dir / f"waterfall_code_{timestamp}.png"
    pygame.image.save(screen, str(path))
    return path


def print_controls() -> None:
    print("Controls:")
    print("  Esc / Q     quit")
    print("  Space       pause / resume")
    print("  R           restart current file and reset the animation clock")
    print("  N / P       next / previous code file")
    print("  C / Tab     cycle visual palette")
    print("  T / Shift+T cycle text color palette forward / backward")
    print("  B           toggle waterfall background on/off (code stays visible)")
    print("  [ / ]       decrease / increase foam opacity")
    print("  D           toggle droplets")
    print("  F           toggle fullscreen")
    print("  L           reload the samples directory")
    print("  S           save a screenshot")


def draw_status_overlay(screen: pygame.Surface, font: pygame.font.Font, text: str, alpha: float) -> None:
    alpha = max(0.0, min(1.0, alpha))
    if alpha <= 0.0 or not text:
        return

    text_surface = font.render(text, True, (244, 248, 255))
    text_surface.set_alpha(int(255 * alpha))

    pad_x = 12
    pad_y = 8
    box = pygame.Surface((text_surface.get_width() + pad_x * 2, text_surface.get_height() + pad_y * 2), pygame.SRCALPHA)
    box.fill((0, 0, 0, int(150 * alpha)))
    box.blit(text_surface, (pad_x, pad_y))

    dest = (18, screen.get_height() - box.get_height() - 18)
    screen.blit(box, dest)


def main() -> None:
    args = parse_args()

    settings = Settings(
        window_size=(args.width, args.height),
        start_fullscreen=args.fullscreen,
        render_height=args.render_height,
        palette_name=args.palette,
        prefer_cuda=args.cuda,
        enable_droplets=not args.no_droplets,
        show_background=not args.droplets_only,
    )

    pygame.init()
    pygame.display.set_caption(settings.title)

    if settings.show_startup_controls:
        print_controls()

    fullscreen = settings.start_fullscreen
    windowed_size = settings.window_size
    background_enabled = settings.show_background
    screen = create_screen(fullscreen, windowed_size)

    samples_dir = resolve_samples_path(settings, args.samples)
    palette_order = palette_names()
    palette_index = palette_order.index(settings.palette_name)

    initial_text_palette = args.text_palette or settings.palette_name
    text_palette_index = palette_order.index(initial_text_palette)

    scroller, influence_builder, renderer, droplets, sim_surface, scaled_surface = build_scene(
        settings,
        screen.get_size(),
        samples_dir,
        palette_order[palette_index],
        palette_order[text_palette_index],
        settings.default_seed,
    )

    status_font = pygame.font.SysFont("DejaVu Sans", 18)
    status_text = ""
    status_timer = 0.0

    def set_status(message: str) -> None:
        nonlocal status_text, status_timer
        status_text = message
        status_timer = STATUS_DURATION
        print(message)

    def rebuild_scene(output_size: tuple[int, int]) -> None:
        nonlocal scroller, influence_builder, renderer, droplets, sim_surface, scaled_surface
        current_foam = droplets.foam_opacity
        scroller, influence_builder, renderer, droplets, sim_surface, scaled_surface = build_scene(
            settings,
            output_size,
            samples_dir,
            palette_order[palette_index],
            palette_order[text_palette_index],
            settings.default_seed,
        )
        droplets.set_foam_opacity(current_foam)
        print(f"Renderer backend: {renderer.backend.label}")

    print(f"Renderer backend: {renderer.backend.label}")
    print(f"Reactive droplets: {'on' if settings.enable_droplets else 'off'}")
    print(f"Visual palettes available: {len(palette_order)}")
    print(f"Text palettes available: {len(palette_order)}")

    if not background_enabled:
        set_status("Background off — code stays visible")
    else:
        set_status(
            f"Visual: {palette_order[palette_index]} ({palette_index + 1}/{len(palette_order)}) · "
            f"Text: {palette_order[text_palette_index]} ({text_palette_index + 1}/{len(palette_order)})"
        )

    clock = pygame.time.Clock()
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
                rebuild_scene(screen.get_size())
                set_status(f"Resized to {event.w}×{event.h}")

            elif event.type == pygame.KEYDOWN:
                if event.key in (pygame.K_ESCAPE, pygame.K_q):
                    running = False
                elif event.key == pygame.K_SPACE:
                    paused = not paused
                    set_status("Paused" if paused else "Resumed")
                elif event.key == pygame.K_r:
                    time_accumulator = 0.0
                    scroller.restart_current()
                    droplets.reset()
                    set_status("Restarted current file")
                elif event.key == pygame.K_n:
                    scroller.next_file()
                    droplets.reset()
                    set_status("Next code file")
                elif event.key == pygame.K_p:
                    scroller.previous_file()
                    droplets.reset()
                    set_status("Previous code file")
                elif event.key in (pygame.K_c, pygame.K_TAB):
                    palette_index = (palette_index + 1) % len(palette_order)
                    palette_name = palette_order[palette_index]
                    palette = get_palette(palette_name)
                    renderer.set_palette(palette)
                    droplets.set_palette(palette)
                    set_status(f"Visual palette: {palette_name} ({palette_index + 1}/{len(palette_order)})")
                elif event.key == pygame.K_t:
                    direction = -1 if (event.mod & pygame.KMOD_SHIFT) else 1
                    text_palette_index = (text_palette_index + direction) % len(palette_order)
                    text_palette_name = palette_order[text_palette_index]
                    renderer.set_text_palette(get_palette(text_palette_name))
                    set_status(f"Text color: {text_palette_name} ({text_palette_index + 1}/{len(palette_order)})")
                elif event.key == pygame.K_b:
                    background_enabled = not background_enabled
                    label = "Background on" if background_enabled else "Background off — code stays visible"
                    set_status(label)
                elif event.key in (pygame.K_LEFTBRACKET, pygame.K_MINUS):
                    value = droplets.adjust_foam_opacity(-settings.foam_opacity_step)
                    set_status(f"Foam opacity: {int(round(value * 100))}%")
                elif event.key in (pygame.K_RIGHTBRACKET, pygame.K_EQUALS):
                    value = droplets.adjust_foam_opacity(settings.foam_opacity_step)
                    set_status(f"Foam opacity: {int(round(value * 100))}%")
                elif event.key == pygame.K_d:
                    settings.enable_droplets = not settings.enable_droplets
                    if not settings.enable_droplets:
                        droplets.reset()
                    set_status(f"Droplets: {'on' if settings.enable_droplets else 'off'}")
                elif event.key == pygame.K_l:
                    scroller.reload_directory()
                    droplets.reset()
                    set_status("Reloaded samples directory")
                elif event.key == pygame.K_s:
                    path = save_screenshot(screen, settings.screenshot_path())
                    set_status(f"Saved screenshot: {path.name}")
                elif event.key == pygame.K_f:
                    fullscreen = not fullscreen
                    screen = create_screen(fullscreen, windowed_size)
                    pygame.mouse.set_visible(not fullscreen)
                    rebuild_scene(screen.get_size())
                    set_status("Fullscreen on" if fullscreen else "Fullscreen off")

        if not paused:
            time_accumulator += dt
            scroller.update(dt)

        mask = scroller.render_mask(time_accumulator)
        influence = influence_builder.build(mask)
        if settings.enable_droplets and not paused:
            droplets.update(influence, dt, time_accumulator)

        if background_enabled:
            frame = renderer.render(influence, time_accumulator)
        else:
            frame = renderer.render_text_layer(influence, time_accumulator)
        pygame.surfarray.blit_array(sim_surface, frame.swapaxes(0, 1))

        if settings.enable_droplets:
            droplets.draw(sim_surface)

        pygame.transform.smoothscale(sim_surface, screen.get_size(), scaled_surface)
        screen.blit(scaled_surface, (0, 0))

        if status_timer > 0.0 and status_text:
            draw_status_overlay(screen, status_font, status_text, min(1.0, status_timer / STATUS_DURATION))
            status_timer = max(0.0, status_timer - dt)

        pygame.display.flip()

    pygame.quit()


if __name__ == "__main__":
    main()
