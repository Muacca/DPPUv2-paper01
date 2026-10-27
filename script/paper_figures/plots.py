"""Render the paper's nine figures from the supplied v3 CSV datasets.

Figures 1 and 4 are schematics; Figure 4 evaluates the polynomials printed in
its panels. All data plots read saved CSV values. No parameter scan or dataset
recomputation is performed by this module.
"""
from pathlib import Path
import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, FancyBboxPatch
from .data import GEOS, PUBLISHED_DATA, read_rows, validate_dataset
COL = ['#0072B2', '#E69F00', '#D8DDE2']
GEO_COL = ['#0072B2', '#D55E00', '#009E73']
LABEL = {'S3': '$S^3\\times S^1$',
 'T3': '$T^3\\times S^1$',
 'Nil3': '$\\mathrm{Nil}^3\\times S^1$'}
STYLE = {'font.family': 'DejaVu Serif',
 'mathtext.fontset': 'dejavuserif',
 'font.size': 11,
 'axes.labelsize': 11,
 'axes.titlesize': 12,
 'legend.fontsize': 9,
 'xtick.labelsize': 9,
 'ytick.labelsize': 9,
 'axes.spines.top': False,
 'axes.spines.right': False,
 'axes.grid': True,
 'grid.alpha': 0.18,
 'lines.linewidth': 1.8,
 'savefig.facecolor': 'white',
 'svg.fonttype': 'none'}
FIGURE_NAMES = {1: 'fig01_pipeline',
 2: 'fig02_variants',
 3: 'fig03_decomposition',
 4: 'fig04_extrema',
 5: 'fig05_full_phases',
 6: 'fig06_diagnostics_04_05',
 7: 'fig07_diagnostics_08_09',
 8: 'fig08_representatives',
 9: 'fig09_scaling'}

def decorate(ax, xlabel='Scale $s$'):
    ax.set_xlabel(xlabel)
    ax.set_ylabel('$U=-S$')
    ax.axhline(0, color='#808080', lw=0.7)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4), useMathText=True)

def plot_pipeline(data_dir):
    fig, ax = plt.subplots(figsize=(8, 5.8))
    ax.set(xlim=(0, 10), ylim=(0, 10))
    ax.axis('off')

    def box(x, y, w, h, text, color='#F0F4F7', size=10):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=.12', fc=color, ec='#465666', lw=1))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=size)

    def arrow(x, y, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x, y), arrowprops=dict(arrowstyle='->', color='#465666', lw=1.3))
    box(0.2, 8.4, 4.1, 1.2, 'Coframe + geometric torsion T\nCartan connection and curvature')
    box(5.7, 8.4, 4.1, 1.2, 'Global quotient + periods\nVolume, fixed parameters, path')
    arrow(2.25, 8.3, 3.8, 7.65)
    arrow(7.75, 8.3, 6.2, 7.65)
    box(1.3, 6.4, 7.4, 1.1, '$U=-S=P(As^3+Bs^2+Cs),\\quad P>0$' + '\nFULL: B=0   |   TT = positive REE: mixed term')
    arrow(5, 6.3, 5, 5.75)
    box(0.5, 4.4, 9, 1.2, 'Analytic classification on s>0\nA>0: signs of C, B and D=B²-3AC\nA=0: treat the lower-degree polynomial directly', size=10)
    arrow(2, 4.3, 2, 3.7)
    arrow(5, 4.3, 5, 3.7)
    arrow(8, 4.3, 8, 3.7)
    box(0.4, 2.2, 2.9, 1.3, 'Type I\npositive maximum\nthen strict minimum', '#DCEEF8')
    box(3.55, 2.2, 2.9, 1.3, 'Type II\nstrict minimum\nno preceding maximum', '#FFF1D2')
    box(6.7, 2.2, 2.9, 1.3, 'Type III\nno strict minimum\non s>0', '#E9ECEF')
    box(0.5, 0.2, 9, 1.15, 'Independent checks: exterior forms, global integrals, roots/heights\nNumerically unresolved boundaries and input/search failures are reported separately.\nRadius extrema do not establish a full vacuum or dynamical stability.', color='#FAFAFA', size=9)
    return fig

def plot_variants(data_dir):
    data = read_rows(data_dir, 'variant_curves.csv.gz')
    fig, axs = plt.subplots(2, 3, figsize=(8, 5.6), layout='constrained')
    for col, g in enumerate(GEOS):
        dat = [r for r in data if r['geometry'] == g]
        x = np.array([r['s'] for r in dat])
        for row in range(2):
            ax = axs[row, col]
            for variant, color, style in [('FULL', COL[0], '-'), ('TT', COL[1], '-'), ('REE', '#222222', '--')]:
                ax.plot(x, [r['U_' + variant] for r in dat], color=color, ls=style, label=variant if variant != 'REE' else 'positive REE')
            ax.set_xlim(0, 6 if row == 0 else 0.5)
            if row == 1:
                sub = [r for r in dat if r['s'] <= 0.5]
                vals = [r['U_' + v] for r in sub for v in ('FULL', 'TT')]
                ax.set_ylim(min(vals) * 1.2, max(vals) * 1.1)
            decorate(ax, '$r$' if g != 'Nil3' else '$R$')
            if row == 0:
                ax.set_title(LABEL[g], pad=5)
        axs[1, col].set_title('Near-origin detail', fontsize=10)
    axs[0, 0].legend(loc='upper left', fontsize=8)
    fig.suptitle('$V=-4,\\ \\eta=3,\\ \\theta=1,\\ \\kappa=L=1$', fontsize=12)
    return fig

def plot_decomposition(data_dir):
    data = read_rows(data_dir, 'variant_curves.csv.gz')
    fig, axs = plt.subplots(1, 3, figsize=(8, 3.5), layout='constrained')
    for ax, g in zip(axs, GEOS):
        dat = [r for r in data if r['geometry'] == g]
        x = [r['s'] for r in dat]
        for key, label, color, style in [('U_FULL', 'FULL total', COL[0], '-'), ('cubic', 'cubic', '#009E73', '--'), ('linear', 'linear', '#CC79A7', ':'), ('diagnostic_quadratic', 'TT / REE quadratic', COL[1], '-.')]:
            ax.plot(x, [r[key] for r in dat], label=label, color=color, ls=style)
        ax.axhline(0, color='black', ls=':', lw=1, label='FULL NY = 0')
        ax.set_title(LABEL[g])
        decorate(ax, '$r$' if g != 'Nil3' else '$R$')
        ax.set_xlim(0, 6)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.14), fontsize=9)
    fig.suptitle('$V=-4,\\ \\eta=3,\\ \\theta=1,\\ \\kappa=L=1$', fontsize=12)
    return fig

def plot_extrema(data_dir):
    fig, axs = plt.subplots(1, 3, figsize=(8, 3.2), layout='constrained')
    cases = [(1, -3, 2, 'Type I'), (1, 0, -3, 'Type II'), (1, 0, 1, 'Type III')]
    x = np.linspace(0.002, 2.25, 500)
    for i, (ax, (A, B, C, title)) in enumerate(zip(axs, cases)):
        y = A * x ** 3 + B * x * x + C * x
        ax.plot(x, y, color=COL[i] if i < 2 else '#59636D')
        ax.plot(0, 0, 'o', mfc='white', mec='#333333', ms=5, clip_on=False)
        ax.set_title(title)
        ax.text(0.03, 0.94, f'$U=s^3{B:+g}s^2{C:+g}s$' if B else f'$U=s^3{C:+g}s$', transform=ax.transAxes, fontsize=9, va='top')
        if i < 2:
            D = B * B - 3 * A * C
            rp = (-B + math.sqrt(D)) / (3 * A)
            rm = (-B - math.sqrt(D)) / (3 * A)
            low = A * rp ** 3 + B * rp * rp + C * rp
            high = A * rm ** 3 + B * rm * rm + C * rm if i == 0 else 0
            ax.plot(rp, low, 'o', color=COL[i], ms=5)
            if i == 0:
                ax.plot(rm, high, 'o', color=COL[i], ms=5)
            xx = 1.95
            ax.annotate('', xy=(xx, high), xytext=(xx, low), arrowprops=dict(arrowstyle='<->', color='#333333'))
            ax.text(xx - 0.08, (high + low) / 2, '$\\Delta U_{\\mathrm{I}}$' if i == 0 else '$\\Delta U_{\\mathrm{II}}$', ha='right', fontsize=10)
            ax.set_ylim(min(y) * 1.2, max(max(y) * 1.1, high + 0.3))
        else:
            ax.set_ylim(-0.3, max(y) * 1.1)
        ax.set_xlim(0, 2.25)
        decorate(ax, '$s>0$')
    return fig

def phase_panel(ax, dat, g, t, full=False):
    x = np.array(sorted({r['eta'] for r in dat}))
    y = np.array(sorted({r['V'] for r in dat}))
    values = np.array([{'I': 0, 'II': 1, 'III': 2}[r['phase']] for r in dat]).reshape(len(y), len(x))
    ax.pcolormesh(x, y, values, cmap=ListedColormap(COL), norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], 3), shading='nearest', rasterized=False)
    boundaries = []
    if g == 'S3':
        boundaries = [-2, 2]
        if not full and 16 * t * t < 3:
            boundaries.append(math.sqrt(12 / (3 - 16 * t * t)))
    elif g == 'T3' and (not full):
        boundaries = [0]
    elif g == 'Nil3' and (not full) and (64 * t * t > 12):
        boundaries = [1 / math.sqrt(64 * t * t - 12)]
    for b in boundaries:
        if -5 <= b <= 10:
            ax.axvline(b, color='#333333', lw=0.9, ls='--')
    ax.axhline(0, color='#333333', lw=1.2)
    if full:
        for h in [-2, 2] if g == 'S3' else [0] if g == 'T3' else []:
            ax.plot(h, 0, 'o', mfc='white', mec='black', ms=5, clip_on=False)
    ax.set(xlim=(-5, 10), ylim=(-5, 0), xlabel='Geometric $\\eta$', ylabel='Geometric $V$')
    ax.set_xticks([-5, 0, 2, 5, 10])
    ax.set_yticks([-5, -2.5, 0])
    ax.grid(False)

def plot_full_phases(data_dir):
    data = read_rows(data_dir, 'phase_grid.csv.gz')
    fig, axs = plt.subplots(1, 3, figsize=(8, 3.4), layout='constrained')
    for ax, g in zip(axs, GEOS):
        dat = [r for r in data if r['geometry'] == g and r['variant'] == 'FULL']
        phase_panel(ax, dat, g, 1, True)
        ax.set_title(LABEL[g])
    fig.suptitle('FULL: identical at $\\theta=0,1,2$; $\\kappa=L=1$', fontsize=12)
    fig.legend(handles=[Patch(color=COL[1], label='II: strict radius minimum'), Patch(color=COL[2], label='III: no strict minimum')], loc='lower center', bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=9)
    return fig

def plot_diagnostic_phases(data_dir, thetas):
    data = read_rows(data_dir, 'phase_grid.csv.gz')
    fig, axs = plt.subplots(3, 2, figsize=(8, 6.9), layout='constrained')
    for row, g in enumerate(GEOS):
        for col, t in enumerate(thetas):
            ax = axs[row, col]
            dat = [r for r in data if r['geometry'] == g and r['variant'] == 'TT' and (r['theta'] == t)]
            phase_panel(ax, dat, g, t)
            ax.set_title(LABEL[g] + f', $\\theta={t:g}$', fontsize=11)
    fig.suptitle('TT = positive REE diagnostics; $\\kappa=L=1$', fontsize=12)
    fig.legend(handles=[Patch(color=c, label=l) for c, l in zip(COL, ['I: barrier + minimum', 'II: minimum, no barrier', 'III: no strict minimum'])], loc='lower center', bbox_to_anchor=(0.5, -0.035), ncol=3, fontsize=9)
    return fig

def plot_representatives(data_dir):
    data = read_rows(data_dir, 'representative_curves.csv.gz')
    fig, axs = plt.subplots(1, 3, figsize=(8, 3.8), layout='constrained')
    for ax, g in zip(axs, GEOS):
        for h, col in [(3, '#59636D'), (0, COL[1])]:
            dat = [r for r in data if r['geometry'] == g and r['eta'] == h]
            ax.plot([r['s'] for r in dat], [r['U_FULL'] for r in dat], color=col, label='$\\eta=\\pm3$' if h else '$\\eta=0$')
        if g == 'S3':
            e = next((row for row in read_rows(data_dir, 'representative_values.csv.gz') if row['geometry'] == 'S3' and row['variant'] == 'FULL' and (row['eta'] == 0)))
            ax.plot(e['rmin'], e['Umin'], 'o', color=COL[1])
            ax.annotate('$r_0=\\sqrt{3}/2$', xy=(e['rmin'], e['Umin']), xytext=(0.12, 0.4), textcoords='axes fraction', fontsize=9, arrowprops=dict(arrowstyle='->', lw=0.7))
        ax.set_title(LABEL[g])
        decorate(ax, '$r$' if g != 'Nil3' else '$R$')
        ax.set_xlim(0, 2)
        ax.legend(fontsize=9, loc='upper left')
    fig.suptitle('FULL: $V=-4,\\ \\theta=1,\\ \\kappa=L=1$', fontsize=12)
    return fig

def plot_scaling(data_dir):
    data = read_rows(data_dir, 'scaling.csv.gz')
    fig, axs = plt.subplots(2, 2, figsize=(8, 5.8), layout='constrained')
    for col, key in enumerate(('rmin', 'height')):
        dat = [r for r in data if r['variant'] == 'FULL']
        axs[0, col].loglog([-r['V'] for r in dat], [r[key] for r in dat], color=GEO_COL[0])
        axs[0, col].set_title('$S^3$ FULL, $\\eta=0$ (Type II)')
        for g, color in zip(GEOS, GEO_COL):
            dat = [r for r in data if r['variant'] == 'TT' and r['geometry'] == g]
            axs[1, col].loglog([-r['V'] for r in dat], [r[key] for r in dat], color=color, label=LABEL[g])
        axs[1, col].set_title('TT = positive REE, $\\eta=3$ (Type I)')
        for row in range(2):
            axs[row, col].set_xlabel('$|V|$, with $V<0$')
            axs[row, col].set_ylabel('Minimum scale $s_0$' if col == 0 else 'Depth $\\Delta U_{\\mathrm{II}}$' if row == 0 else 'Barrier $\\Delta U_{\\mathrm{I}}$')
            axs[row, col].grid(True, which='both', alpha=0.15)
    axs[1, 1].legend(fontsize=9)
    fig.suptitle('Valid-minimum scaling: $\\theta=\\kappa=L=1$', fontsize=12)
    return fig

def render_figure(number, output_dir, *, data_dir=PUBLISHED_DATA):
    """Save one numbered figure as 360-dpi PNG and SVG, returning both paths.

    Choose a new output directory or unused filenames. Existing figure files
    are not overwritten. The notebook can display the returned PNG directly.
    """
    if number not in FIGURE_NAMES:
        raise ValueError('Figure number must be an integer from 1 to 9')
    validate_dataset(data_dir)
    dest = Path(output_dir).resolve()
    paths = tuple((dest / (FIGURE_NAMES[number] + suffix) for suffix in ('.png', '.svg')))
    if any((path.exists() for path in paths)):
        raise FileExistsError(f'Figure output already exists: {paths[0].stem}')
    dest.mkdir(parents=True, exist_ok=True)
    factories = {1: plot_pipeline, 2: plot_variants, 3: plot_decomposition, 4: plot_extrema, 5: plot_full_phases, 6: lambda source: plot_diagnostic_phases(source, [0.4, 0.5]), 7: lambda source: plot_diagnostic_phases(source, [0.8, 0.9]), 8: plot_representatives, 9: plot_scaling}
    defaults = {key: value for key, value in matplotlib.rcParamsDefault.items() if key != 'backend'}
    with matplotlib.rc_context({**defaults, **STYLE, 'svg.hashsalt': 'paper-v3'}):
        fig = factories[number](Path(data_dir))
        try:
            fig.savefig(paths[0], dpi=360, bbox_inches='tight', pad_inches=0.08)
            fig.savefig(paths[1], bbox_inches='tight', pad_inches=0.08, metadata={'Date': None})
        finally:
            plt.close(fig)
    return paths

def render_all(output_dir, *, data_dir=PUBLISHED_DATA):
    """Render Figures 1-9 from saved data into a new directory."""
    validate_dataset(data_dir)
    dest = Path(output_dir).resolve()
    dest.mkdir(parents=True, exist_ok=False)
    return {number: render_figure(number, dest, data_dir=data_dir) for number in FIGURE_NAMES}
