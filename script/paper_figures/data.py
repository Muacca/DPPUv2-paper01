"""Fixed v3 figure datasets, using geometric torsion and kappa=L=1.

The reduced coefficients are Eqs. M3-03, M3-06 and M3-10; extrema and
phase criteria follow Section 4. This module evaluates those polynomials,
not the full field equations or a general parameter scanner.
"""
import csv
import gzip
import hashlib
import io
import json
import math
import platform
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import numpy as np
PUBLISHED_DATA = Path(__file__).resolve().parents[2] / 'data' / 'paper_v3'
GEOS = ('S3', 'T3', 'Nil3')
PARAMETERS = {'kappa': 1,
 'L': 1,
 'torsion': 'geometric T',
 'code_map': 'V_tau=-V; eta_tau=-eta',
 'scale': 's=r (S3,T3), s=R (Nil3)',
 'T3_path': 'R1=R2=R3=r',
 'S3_physical_radius': 'a=r/2',
 'volumes': {'S3': 'pi^2 L r^3/4', 'T3': '(2 pi)^3 L r^3', 'Nil3': '(2 pi)^4 L R^3'},
 'Nil_lattice': ['2pi', '2pi', '4pi^2']}
GRIDS = {'phase_V': [-5, 0, 51],
 'phase_eta': [-5, 10, 151],
 'theta_FULL_verified': [0, 1, 2],
 'theta_diagnostic': [0.4, 0.5, 0.8, 0.9],
 'curve_s': [0.002, 6, 1201],
 'representative_s': [0.002, 2, 601],
 'scaling_absV': [0.5, 6, 121]}
COUNTS = {'variant_curves': 3603,
 'representative_curves': 5409,
 'phase_grid': 115515,
 'scaling': 484,
 'representative_values': 15,
 'local_invariants': 3}
EQUATIONS = ['M3-03',
 'M3-06',
 'M3-10',
 'M4-01',
 'M4-02',
 'M4-03',
 'M4-04',
 'M4-05',
 'M5-01',
 'M5-02',
 'M5-03',
 'M6-02',
 'A-05']
DATA_FILES = tuple((name + '.csv.gz' for name in COUNTS))

def _coefficients(g, V, h, t, variant='FULL'):
    p = {'S3': math.pi ** 2 / 12, 'T3': 8 * math.pi ** 3 / 3, 'Nil3': 4 * math.pi ** 4 / 3}[g]
    z = 4 if g == 'Nil3' else 1
    return (p, z * V * V, 0 if variant == 'FULL' else z * 12 * V * h * t, {'S3': 9 * h * h - 36, 'T3': 9 * h * h, 'Nil3': 36 * h * h + 3}[g])

@lru_cache(None)
def _phase(g, Vsign, h, t, variant):
    if not Vsign:
        return 'III'
    v = F(Vsign)
    h = F(str(h))
    t = F(str(t))
    _, A, B, C = _coefficients(g, v, h, t, variant)
    D = B * B - 3 * A * C
    if C < 0 or (C == 0 and B < 0):
        return 'II'
    if C > 0 and B < 0 and (D > 0):
        return 'I'
    return 'III'

def _extrema(g, V, h, t, variant='FULL'):
    p, A, B, C = _coefficients(g, V, h, t, variant)
    kind = _phase(g, 0 if V == 0 else 1 if V > 0 else -1, h, t, variant)
    if kind == 'III':
        return dict(phase=kind, rmin=None, rmax=None, Umin=None, height=None)
    d = math.sqrt(B * B - 3 * A * C)
    r = -C / (B + d) if B >= 0 else (-B + d) / (3 * A)
    rmax = (-B - d) / (3 * A) if kind == 'I' else None
    Umin = -p * r * r * (2 * A * r + B)
    height = 4 * p * d ** 3 / (27 * A * A) if kind == 'I' else -Umin
    return dict(phase=kind, rmin=r, rmax=rmax, Umin=Umin, height=height)

def _curve(g, V, h, t, r):
    p, A, B, C = _coefficients(g, V, h, t, 'TT')
    cubic = p * A * r ** 3
    linear = p * C * r
    quad = p * B * r * r
    return dict(geometry=g, V=V, eta=h, V_tau=-V, eta_tau=-h, theta=t, s=r, cubic=cubic, linear=linear, diagnostic_quadratic=quad, U_FULL=cubic + linear, U_TT=cubic + linear + quad, U_REE=cubic + linear + quad)

def _save_csv(dest, name, rows):
    with (dest / name).open('xb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as packed:
            with io.TextIOWrapper(packed, encoding='utf-8', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

def validate_dataset(data_dir=PUBLISHED_DATA):
    """Check the published schema, fixed figure inputs, and CSV file hashes."""
    data_dir = Path(data_dir)
    meta = json.loads((data_dir / 'parameters.json').read_text(encoding='utf-8'))
    for key, expected in (('schema', 'paper-figure-data-v1'), ('parameters', PARAMETERS), ('grids', GRIDS), ('counts', COUNTS)):
        if meta.get(key) != expected:
            raise ValueError(f'{key} differs from the v3 figure specification')
    for name in DATA_FILES:
        digest = hashlib.sha256((data_dir / name).read_bytes()).hexdigest()
        if meta.get('data_sha256', {}).get(name) != digest:
            raise ValueError(f'Dataset hash mismatch: {name}')
    return meta

def read_rows(data_dir, name):
    """Read a numerical figure CSV; empty extrema remain None."""
    with gzip.open(Path(data_dir) / name, 'rt', encoding='utf-8') as stream:
        return [{key: value if key in ('geometry', 'variant', 'phase') else float(value) if value else None for key, value in row.items()} for row in csv.DictReader(stream)]

def recompute_data(output_dir):
    """Recompute all six v3 CSV datasets into a new directory.

    Inputs and grids are fixed to PARAMETERS and GRIDS. Existing directories,
    including the supplied data directory, are never overwritten. Coordinates
    V and eta are geometric; V_tau=-V and eta_tau=-eta are recorded alongside.
    """
    dest = Path(output_dir).resolve()
    dest.mkdir(parents=True, exist_ok=False)
    curves = [_curve(g, -4, 3, 1, float(r)) for g in GEOS for r in np.linspace(0.002, 6, 1201)]
    _save_csv(dest, 'variant_curves.csv.gz', curves)
    representatives = [_curve(g, -4, h, 1, float(r)) for g in GEOS for h in (-3, 0, 3) for r in np.linspace(0.002, 2, 601)]
    _save_csv(dest, 'representative_curves.csv.gz', representatives)
    maps = []
    for g in GEOS:
        for variant, theta in [('FULL', 1), ('TT', 0.4), ('TT', 0.5), ('TT', 0.8), ('TT', 0.9)]:
            for iv in range(-50, 1):
                for ih in range(-50, 101):
                    V = iv / 10
                    h = ih / 10
                    result = _extrema(g, V, h, theta, variant)
                    maps.append(dict(geometry=g, variant=variant, V=V, eta=h, V_tau=-V, eta_tau=-h, theta=theta, **result))
    _save_csv(dest, 'phase_grid.csv.gz', maps)
    scaling = []
    for g, variant, h in [('S3', 'FULL', 0)] + [(g, 'TT', 3) for g in GEOS]:
        for magnitude in np.geomspace(0.5, 6, 121):
            V = -float(magnitude)
            scaling.append(dict(geometry=g, variant=variant, V=V, eta=h, V_tau=-V, eta_tau=-h, theta=1, **_extrema(g, V, h, 1, variant)))
    _save_csv(dest, 'scaling.csv.gz', scaling)
    table = []
    for g in GEOS:
        for variant, h in [('FULL', 3), ('FULL', 2), ('FULL', 0), ('FULL', -3), ('TT', 3)]:
            table.append(dict(geometry=g, variant=variant, V=-4, eta=h, V_tau=4, eta_tau=-h, theta=1, **_extrema(g, -4, h, 1, variant)))
    _save_csv(dest, 'representative_values.csv.gz', table)
    invariants = []
    for g in GEOS:
        Rlc = {'S3': 24, 'T3': 0, 'Nil3': -0.5}[g]
        invariants.append(dict(geometry=g, s=1, V=-4, eta=3, V_tau=4, eta_tau=-3, R_LC=Rlc, R_EC=Rlc - 54 - 32 / 3, torsion_norm=216 + 32 / 3, n_TT=48, n_REE=48, n_FULL=0))
    _save_csv(dest, 'local_invariants.csv.gz', invariants)
    meta = {'schema': 'paper-figure-data-v1', 'parameters': PARAMETERS, 'grids': GRIDS, 'counts': COUNTS, 'equations': EQUATIONS, 'generator': {'module': 'paper_figures.data', 'python': platform.python_version(), 'numpy': np.__version__, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, 'data_sha256': {name: hashlib.sha256((dest / name).read_bytes()).hexdigest() for name in DATA_FILES}, 'notes': ['Geometric torsion T; kappa=L=1; fixed v3 figure grids.', 'TT and positive REE share the diagnostic data.', 'No-minimum fields are blank. Type I is barrier height; Type II is depth.', 'These are the published phase slices, not a full three-dimensional sweep.']}
    (dest / 'parameters.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    return dest
