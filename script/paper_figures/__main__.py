"""Command-line entry points for data recomputation and saved-data rendering."""
import argparse
from pathlib import Path
from . import PUBLISHED_DATA, recompute_data, render_all


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    compute = commands.add_parser("recompute", help="Recompute the fixed v3 datasets")
    compute.add_argument("--output-dir", type=Path, required=True)
    render = commands.add_parser("render", help="Render Figures 1-9 from saved CSV data")
    render.add_argument("--data-dir", type=Path, default=PUBLISHED_DATA)
    render.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "recompute":
        print(recompute_data(args.output_dir))
    else:
        import matplotlib
        matplotlib.use("Agg")
        for number, paths in render_all(args.output_dir, data_dir=args.data_dir).items():
            print(f"{number}: {paths[0]} | {paths[1]}")


if __name__ == "__main__":
    main()
