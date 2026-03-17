#!/usr/bin/env python3
"""
Split a large CSV into ~N MB parts, named 1.csv, 2.csv, 3.csv, ...
Streaming: does NOT load the whole file into memory.

Modes:
  - raw (default): fast, copies lines as-is (assumes each record is one physical line)
  - csv: CSV-aware (handles quoted newlines) but rewrites output formatting

Examples:
  python split_csv.py big.csv out_parts --target-mb 10
  python split_csv.py big.csv out_parts --target-bytes 10000000
  python split_csv.py big.csv out_parts --mode raw --no-repeat-header
  python split_csv.py big.csv out_parts --mode csv --encoding utf-8
"""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path
from typing import List


def split_raw(
        input_path: Path,
        output_dir: Path,
        target_bytes: int,
        has_header: bool = True,
        repeat_header: bool = True,
) -> List[Path]:
    """
    Fast mode: reads the input file as bytes and splits on newline boundaries.
    Preserves original bytes exactly.
    WARNING: If your CSV has embedded newlines inside quoted fields, use split_csv_aware().
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    created: List[Path] = []

    def open_part(part_num: int):
        out_path = output_dir / f"{part_num}.csv"
        out_f = open(out_path, "wb")
        created.append(out_path)
        return out_f

    with open(input_path, "rb") as f:
        first_line = f.readline()
        if first_line == b"":
            return created  # empty input

        if has_header:
            header = first_line
            # file pointer is now after header
        else:
            header = b""
            f.seek(0)  # treat first line as data

        part = 1
        out = open_part(part)
        current_size = 0
        data_lines_in_part = 0

        # Always write header to part 1 if present
        if has_header:
            out.write(header)
            current_size += len(header)

        for line in f:
            # If adding this line would exceed the limit AND we already have some data in this part,
            # rotate to a new part.
            if (current_size + len(line) > target_bytes) and (data_lines_in_part > 0):
                out.close()
                part += 1
                out = open_part(part)
                current_size = 0
                data_lines_in_part = 0

                if has_header and repeat_header:
                    out.write(header)
                    current_size += len(header)

            out.write(line)
            current_size += len(line)
            data_lines_in_part += 1

        out.close()

    return created


def split_csv_aware(
        input_path: Path,
        output_dir: Path,
        target_bytes: int,
        encoding: str = "utf-8",
        has_header: bool = True,
        repeat_header: bool = True,
        delimiter: str = ",",
        quotechar: str = '"',
) -> List[Path]:
    """
    CSV-aware mode: uses csv.reader so it won't break on quoted newlines,
    but output is rewritten (still valid CSV).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    created: List[Path] = []

    def open_part(part_num: int):
        out_path = output_dir / f"{part_num}.csv"
        out_f = open(out_path, "wb")
        created.append(out_path)
        return out_f

    # Reusable buffer to serialize one CSV row at a time
    row_buf = io.StringIO()
    buf_writer = csv.writer(
        row_buf,
        delimiter=delimiter,
        quotechar=quotechar,
        lineterminator="\n",
    )

    def row_to_bytes(row) -> bytes:
        row_buf.seek(0)
        row_buf.truncate(0)
        buf_writer.writerow(row)
        return row_buf.getvalue().encode(encoding)

    with open(input_path, "r", encoding=encoding, newline="") as f:
        reader = csv.reader(f, delimiter=delimiter, quotechar=quotechar)

        header_row = None
        header_bytes = b""
        if has_header:
            try:
                header_row = next(reader)
                header_bytes = row_to_bytes(header_row)
            except StopIteration:
                header_row = None

        part = 1
        out = open_part(part)
        current_size = 0
        rows_in_part = 0

        # Always write header to part 1 if present
        if header_row is not None:
            out.write(header_bytes)
            current_size += len(header_bytes)

        for row in reader:
            row_bytes = row_to_bytes(row)

            if (current_size + len(row_bytes) > target_bytes) and (rows_in_part > 0):
                out.close()
                part += 1
                out = open_part(part)
                current_size = 0
                rows_in_part = 0

                if (header_row is not None) and repeat_header:
                    out.write(header_bytes)
                    current_size += len(header_bytes)

            out.write(row_bytes)
            current_size += len(row_bytes)
            rows_in_part += 1

        out.close()

    return created


def main() -> int:
    p = argparse.ArgumentParser(
        description="Split a large CSV into ~N MB parts named 1.csv, 2.csv, 3.csv, ..."
    )
    p.add_argument("input_csv", type=Path, help="Path to the big CSV file")
    p.add_argument("output_dir", type=Path, help="Directory where parts will be written")

    p.add_argument(
        "--target-mb",
        type=float,
        default=10.0,
        help="Target part size in MiB (1 MiB = 1024*1024 bytes). Default: 10",
    )
    p.add_argument(
        "--target-bytes",
        type=int,
        default=None,
        help="Exact target size in bytes (overrides --target-mb). Example: 10000000",
    )

    p.add_argument(
        "--mode",
        choices=["raw", "csv"],
        default="raw",
        help="raw=fast line-based copy (default). csv=CSV-aware (handles quoted newlines, rewrites output).",
    )

    p.add_argument(
        "--no-header",
        action="store_true",
        help="Treat the first line as data (do not treat as a header row).",
    )
    p.add_argument(
        "--no-repeat-header",
        action="store_true",
        help="If there is a header, write it only in part 1 (not repeated in every part).",
    )

    # Only used for --mode csv
    p.add_argument("--encoding", default="utf-8", help="(csv mode) Input encoding. Default: utf-8")
    p.add_argument("--delimiter", default=",", help="(csv mode) CSV delimiter. Default: ,")
    p.add_argument("--quotechar", default='"', help='(csv mode) CSV quote char. Default: "')

    args = p.parse_args()

    if not args.input_csv.exists():
        raise SystemExit(f"Input file not found: {args.input_csv}")

    target_bytes = args.target_bytes if args.target_bytes is not None else int(args.target_mb * 1024 * 1024)

    has_header = not args.no_header
    repeat_header = has_header and (not args.no_repeat_header)

    if args.mode == "raw":
        parts = split_raw(
            input_path=args.input_csv,
            output_dir=args.output_dir,
            target_bytes=target_bytes,
            has_header=has_header,
            repeat_header=repeat_header,
        )
    else:
        parts = split_csv_aware(
            input_path=args.input_csv,
            output_dir=args.output_dir,
            target_bytes=target_bytes,
            encoding=args.encoding,
            has_header=has_header,
            repeat_header=repeat_header,
            delimiter=args.delimiter,
            quotechar=args.quotechar,
        )

    print(f"Done. Created {len(parts)} file(s) in: {args.output_dir.resolve()}")
    if parts:
        print(f"First: {parts[0].name}   Last: {parts[-1].name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())