"""Reproduce the v3 paper's fixed numerical datasets and Figures 1-9."""
from .data import PUBLISHED_DATA, recompute_data, validate_dataset


def render_figure(number, output_dir, *, data_dir=PUBLISHED_DATA):
    """Render one figure without recomputing datasets."""
    from .plots import render_figure as render
    return render(number, output_dir, data_dir=data_dir)


def render_all(output_dir, *, data_dir=PUBLISHED_DATA):
    """Render all nine figures without recomputing datasets."""
    from .plots import render_all as render
    return render(output_dir, data_dir=data_dir)


__all__ = ["PUBLISHED_DATA", "recompute_data", "validate_dataset", "render_figure", "render_all"]
