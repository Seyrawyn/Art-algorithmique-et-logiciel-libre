#!/usr/bin/env python3
"""
preprocess_trip_walkers_auto_manifest.py

Preprocess trip CSVs into JSON files *for p5.js animation* where each trip is:

  [startLat, startLon, endLat, endLon, distance_m]

It also writes a `manifest.json` by default so a browser sketch can pick a random
file (browsers cannot list directories).

Inputs can be:
  - CSV files
  - directories (it scans for *.csv)
  - glob patterns (expanded by Python, so it works even if your shell doesn't expand '*')

Outputs (per-file mode, default):
  --outdir/
     <name>_trips.json      (JSON object with {"trips":[...], "bounds":..., ...})
     manifest.json          (list of filenames in the folder)

Outputs (merge mode):
  --merge-out somefile.json (one merged JSON object)
  optionally --manifest-out to write a manifest containing only that file

Examples
--------
# Per-file + manifest (recommended):
python preprocess_trip_walkers_auto_manifest.py data/output_2024 --recursive --outdir data/trips/2024
python preprocess_trip_walkers_auto_manifest.py data/output_2025 --recursive --outdir data/trips/2025

# Merge into one JSON per year:
python preprocess_trip_walkers_auto_manifest.py data/output_2024 --recursive --merge-out data/trips/2024_trips.json
python preprocess_trip_walkers_auto_manifest.py data/output_2025 --recursive --merge-out data/trips/2025_trips.json
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


# -----------------------------
# Column detection
# -----------------------------

def find_column(cols: Iterable[str], must_have: List[str]) -> Optional[str]:
    """Find first column whose uppercased name contains ALL substrings in must_have."""
    for c in cols:
        u = str(c).upper()
        if all(token in u for token in must_have):
            return str(c)
    return None


# -----------------------------
# Trip extraction + distance
# -----------------------------

def haversine_m(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    """Vectorized great-circle distance (meters)."""
    R = 6371000.0
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * (np.sin(dl / 2.0) ** 2)
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    return R * c


def compute_trip_rows(df: pd.DataFrame, cap_m: float, min_m: float, round_coords: int, round_dist: int) -> Tuple[List[List[float]], dict, float]:
    """
    Returns:
      trips: list of [slat, slon, elat, elon, dist_m]
      bounds: {minLat,maxLat,minLon,maxLon}
      max_dist: max distance in meters
    """
    # Detect columns (robust-ish)
    start_lat = find_column(df.columns, ["START", "LAT"])
    start_lon = (
        find_column(df.columns, ["START", "LON"])
        or find_column(df.columns, ["START", "LNG"])
        or find_column(df.columns, ["START", "LONG"])
    )
    end_lat = find_column(df.columns, ["END", "LAT"])
    end_lon = (
        find_column(df.columns, ["END", "LON"])
        or find_column(df.columns, ["END", "LNG"])
        or find_column(df.columns, ["END", "LONG"])
    )

    if not all([start_lat, start_lon, end_lat, end_lon]):
        raise ValueError(
            "Could not detect required columns.\n"
            f"Detected: start_lat={start_lat}, start_lon={start_lon}, end_lat={end_lat}, end_lon={end_lon}\n"
            f"Available columns: {list(df.columns)}"
        )

    # Numeric conversion
    for c in [start_lat, start_lon, end_lat, end_lon]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    sub = df[[start_lat, start_lon, end_lat, end_lon]].dropna()
    if sub.empty:
        return [], {"minLat": None, "maxLat": None, "minLon": None, "maxLon": None}, 0.0

    slat = sub[start_lat].to_numpy(dtype=np.float64, copy=False)
    slon = sub[start_lon].to_numpy(dtype=np.float64, copy=False)
    elat = sub[end_lat].to_numpy(dtype=np.float64, copy=False)
    elon = sub[end_lon].to_numpy(dtype=np.float64, copy=False)

    # Sanity masks
    mask = (
        (slat >= -90) & (slat <= 90) &
        (elat >= -90) & (elat <= 90) &
        (slon >= -180) & (slon <= 180) &
        (elon >= -180) & (elon <= 180) &
        ~((slat == 0.0) & (slon == 0.0)) &
        ~((elat == 0.0) & (elon == 0.0))
    )

    slat = slat[mask]; slon = slon[mask]; elat = elat[mask]; elon = elon[mask]
    if slat.size == 0:
        return [], {"minLat": None, "maxLat": None, "minLon": None, "maxLon": None}, 0.0

    d = haversine_m(slat, slon, elat, elon)

    good = np.isfinite(d) & (d >= float(min_m)) & (d <= float(cap_m))
    slat = slat[good]; slon = slon[good]; elat = elat[good]; elon = elon[good]; d = d[good]

    if d.size == 0:
        return [], {"minLat": None, "maxLat": None, "minLon": None, "maxLon": None}, 0.0

    # Rounding (keeps file size down a bit)
    if round_coords is not None and round_coords >= 0:
        slat = np.round(slat, round_coords)
        slon = np.round(slon, round_coords)
        elat = np.round(elat, round_coords)
        elon = np.round(elon, round_coords)

    if round_dist is not None and round_dist >= 0:
        d = np.round(d, round_dist)

    trips = np.column_stack([slat, slon, elat, elon, d]).tolist()

    bounds = {
        "minLat": float(np.min(np.concatenate([slat, elat]))),
        "maxLat": float(np.max(np.concatenate([slat, elat]))),
        "minLon": float(np.min(np.concatenate([slon, elon]))),
        "maxLon": float(np.max(np.concatenate([slon, elon]))),
    }

    max_dist = float(np.max(d))
    return trips, bounds, max_dist


# -----------------------------
# CSV reading
# -----------------------------

def read_csv_auto(path: Path) -> pd.DataFrame:
    """Read CSV with delimiter sniffing."""
    try:
        return pd.read_csv(path, sep=None, engine="python", encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=None, engine="python", encoding="latin-1")


# -----------------------------
# Input expansion (works even if your shell doesn't expand '*')
# -----------------------------

_GLOB_CHARS = set("*?[")
def looks_like_glob(s: str) -> bool:
    return any(ch in s for ch in _GLOB_CHARS) or ("**" in s)


def expand_inputs(inputs: List[str], recursive: bool) -> List[Path]:
    """
    Expand a list of inputs into concrete CSV file Paths.

    Each input can be:
      - a CSV file path
      - a directory (we scan for *.csv)
      - a glob pattern (expanded by Python, so it works even if your shell doesn't expand '*')
    """
    out: List[Path] = []

    for raw in inputs:
        p = Path(raw)

        # 1) Existing directory -> scan
        if p.exists() and p.is_dir():
            pattern = "**/*.csv" if recursive else "*.csv"
            out.extend(sorted(p.glob(pattern)))
            continue

        # 2) Glob expansion (or non-existing path that might be a glob)
        if looks_like_glob(raw) or not p.exists():
            matches = glob.glob(raw, recursive=recursive)
            if matches:
                for m in matches:
                    mp = Path(m)
                    if mp.is_file() and mp.suffix.lower() == ".csv":
                        out.append(mp)
                continue

        # 3) Existing file
        if p.exists() and p.is_file():
            if p.suffix.lower() == ".csv":
                out.append(p)
            else:
                print(f"Warning: skipping non-CSV file: {p}", file=sys.stderr)
            continue

        print(f"Warning: no match for input: {raw}", file=sys.stderr)

    # Dedupe while preserving order
    uniq: List[Path] = []
    seen = set()
    for f in out:
        try:
            key = str(f.resolve())
        except Exception:
            key = str(f)
        if key not in seen:
            seen.add(key)
            uniq.append(f)

    return uniq


# -----------------------------
# Output naming + manifest
# -----------------------------

def safe_slug(s: str) -> str:
    s = s.replace("/", "__").replace("\\", "__")
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def output_name(in_path: Path, name_mode: str) -> str:
    if name_mode == "stem":
        base = in_path.stem
    elif name_mode == "parent_stem":
        base = f"{in_path.parent.name}__{in_path.stem}"
    elif name_mode == "relpath":
        try:
            rel = in_path.resolve().relative_to(Path.cwd().resolve())
            base = safe_slug(str(rel.with_suffix("")))
        except Exception:
            base = safe_slug(str(in_path.with_suffix("")))
    else:
        raise ValueError(f"Unknown name_mode: {name_mode}")

    return f"{base}_trips.json"


def write_manifest(outdir: Path, manifest_name: str, pattern: str = "*_trips*.json") -> Path:
    """Write a manifest (list of filenames) for all trip JSON files in outdir."""
    files = sorted([p.name for p in outdir.glob(pattern) if p.is_file() and p.name != manifest_name])
    manifest_path = outdir / manifest_name
    with open(manifest_path, "w", encoding="utf-8") as fp:
        json.dump(files, fp, indent=2)
    return manifest_path


# -----------------------------
# CLI
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Preprocess trip CSVs into endpoint+distance JSONs for p5.js.")
    ap.add_argument(
        "inputs",
        nargs="+",
        help="CSV files, glob patterns, and/or directories (e.g. data/output_2024)",
    )
    ap.add_argument(
        "--recursive",
        action="store_true",
        help="When an input is a directory or glob, include subfolders",
    )
    ap.add_argument("--cap-m", type=float, default=20000.0, help="Drop trips above this many meters")
    ap.add_argument("--min-m", type=float, default=10.0, help="Drop trips shorter than this many meters")

    ap.add_argument("--round-coords", type=int, default=5, help="Round lat/lon to this many decimals")
    ap.add_argument("--round-dist", type=int, default=1, help="Round distance (meters) to this many decimals")

    # Per-file mode
    ap.add_argument("--outdir", type=str, default=".", help="Output directory (per-file mode)")
    ap.add_argument(
        "--name-mode",
        type=str,
        default="parent_stem",
        choices=["stem", "parent_stem", "relpath"],
        help="How to name outputs in per-file mode (avoid overwrites)",
    )

    # Manifest
    ap.add_argument(
        "--manifest-name",
        type=str,
        default="manifest.json",
        help="Filename of the manifest written into --outdir (per-file mode)",
    )
    ap.add_argument(
        "--no-manifest",
        action="store_true",
        help="Disable writing the manifest in per-file mode",
    )

    # Merge mode
    ap.add_argument(
        "--merge-out",
        type=str,
        default=None,
        help="If set, merge ALL input CSVs into a single JSON output path",
    )
    ap.add_argument(
        "--manifest-out",
        type=str,
        default=None,
        help="Optional: if set in merge mode, write a manifest that contains ONLY the merged file name",
    )

    args = ap.parse_args()

    files = expand_inputs(args.inputs, recursive=bool(args.recursive))
    if not files:
        print("No CSV files found from the provided inputs.", file=sys.stderr)
        sys.exit(2)

    # -----------------
    # Merge mode
    # -----------------
    if args.merge_out:
        merge_path = Path(args.merge_out)
        merge_path.parent.mkdir(parents=True, exist_ok=True)

        merged_trips: List[List[float]] = []
        bounds_all = {"minLat": None, "maxLat": None, "minLon": None, "maxLon": None}
        max_dist_all = 0.0
        ok_files = 0

        for f in files:
            try:
                df = read_csv_auto(f)
                trips, bounds, maxd = compute_trip_rows(
                    df,
                    cap_m=args.cap_m,
                    min_m=args.min_m,
                    round_coords=args.round_coords,
                    round_dist=args.round_dist,
                )
                if trips:
                    merged_trips.extend(trips)
                    max_dist_all = max(max_dist_all, maxd)

                    # Merge bounds
                    if bounds_all["minLat"] is None:
                        bounds_all = bounds
                    else:
                        bounds_all = {
                            "minLat": min(bounds_all["minLat"], bounds["minLat"]),
                            "maxLat": max(bounds_all["maxLat"], bounds["maxLat"]),
                            "minLon": min(bounds_all["minLon"], bounds["minLon"]),
                            "maxLon": max(bounds_all["maxLon"], bounds["maxLon"]),
                        }

                ok_files += 1
                print(f"{f}: +{len(trips)} trips")
            except Exception as e:
                print(f"ERROR processing {f}: {e}", file=sys.stderr)

        payload = {
            "trips": merged_trips,
            "count": len(merged_trips),
            "bounds": bounds_all,
            "maxDistM": float(max_dist_all),
            "source": "MERGED",
        }

        with open(merge_path, "w", encoding="utf-8") as fp:
            json.dump(payload, fp)

        print(f"Merged {ok_files}/{len(files)} files -> {merge_path} ({len(merged_trips)} trips)")

        if args.manifest_out:
            man_path = Path(args.manifest_out)
            man_path.parent.mkdir(parents=True, exist_ok=True)
            with open(man_path, "w", encoding="utf-8") as fp:
                json.dump([merge_path.name], fp, indent=2)
            print(f"Wrote merge manifest -> {man_path}")

        return

    # -----------------
    # Per-file mode
    # -----------------
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for f in files:
        try:
            df = read_csv_auto(f)
            trips, bounds, maxd = compute_trip_rows(
                df,
                cap_m=args.cap_m,
                min_m=args.min_m,
                round_coords=args.round_coords,
                round_dist=args.round_dist,
            )
        except Exception as e:
            print(f"ERROR processing {f}: {e}", file=sys.stderr)
            continue

        out_path = outdir / output_name(f, args.name_mode)

        # Avoid accidental overwrite
        if out_path.exists():
            i = 1
            while True:
                candidate = outdir / out_path.name.replace("_trips.json", f"_trips_{i}.json")
                if not candidate.exists():
                    out_path = candidate
                    break
                i += 1

        payload = {
            "trips": trips,
            "count": len(trips),
            "bounds": bounds,
            "maxDistM": float(maxd),
            "source": f.name,
        }

        with open(out_path, "w", encoding="utf-8") as fp:
            json.dump(payload, fp)

        print(f"{f.name}: wrote {len(trips)} trips -> {out_path}")

    if not args.no_manifest:
        man = write_manifest(outdir, args.manifest_name)
        try:
            n = len(json.loads(man.read_text(encoding="utf-8")))
        except Exception:
            n = 0
        print(f"Wrote manifest -> {man} (files listed: {n})")


if __name__ == "__main__":
    main()
