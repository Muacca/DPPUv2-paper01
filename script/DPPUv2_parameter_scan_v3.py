#!/usr/bin/env python3
"""
DPPUv2 Parameter Scan: Phase Diagram Generation
================================================================

Purpose: Systematically scan (V, η, θ_NY) parameter space for 3 topologies
         and 3 NY variants, generating CSV data for phase diagram visualization.

Topologies: S³xS¹, T³xS¹, Nil³xS¹
NY Variants: FULL, TT, REE

Output: CSV files with stability classification (type-I/type-II/type-III)

Author: Muacca
Date: 2025-12-17
"""

import numpy as np
from scipy.optimize import minimize_scalar
import csv
import os
from datetime import datetime
from multiprocessing import Pool, cpu_count
from functools import partial
import argparse
import time

# ==============================================================================
# Constants
# ==============================================================================

PI = np.pi
KAPPA = 1.0
L = 1.0

# For T³×S¹, we set R2 = R3 = R1 = r (isotropic scaling)
# R2 and R3 constants are no longer needed

# Search bounds for r optimization
R_MIN = 0.01
R_MAX = 1000000.0
R_BOUNDARY_THRESHOLD = 0.02  # If r0 < this or > R_MAX - this, consider unstable

# ==============================================================================
# Potential Functions (from v3 engine logs, 2025-12-14)
# ==============================================================================

# --- S³×S¹ ---

def V_S3_FULL(r, V_param, eta, theta):
    """
    S³xS¹ with FULL Nieh-Yan (N = TT - Ree)
    V(r) = 2*pi²*L*r*(V²r² + 6*V*κ²*θ*(η-4)*r + 9η² + 72η + 108) / (3κ²)
    """
    if r <= 0:
        return 1e50
    term1 = V_param**2 * r**2
    term2 = 6 * V_param * KAPPA**2 * theta * (eta - 4) * r
    term3 = 9 * eta**2 + 72 * eta + 108
    return (2 * PI**2 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_S3_TT(r, V_param, eta, theta):
    """
    S³xS¹ with TT-only Nieh-Yan
    V(r) = 2*pi²*L*r*(V²r² + 12*V*η*κ²*θ*r + 9η² + 72η + 108) / (3κ²)
    """
    if r <= 0:
        return 1e50
    term1 = V_param**2 * r**2
    term2 = 12 * V_param * eta * KAPPA**2 * theta * r
    term3 = 9 * eta**2 + 72 * eta + 108
    return (2 * PI**2 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_S3_REE(r, V_param, eta, theta):
    """
    S³xS¹ with Ree-only Nieh-Yan
    V(r) = 2*pi²*L*r*(V²r² + 6*V*κ²*θ*(η+4)*r + 9η² + 72η + 108) / (3κ²)
    """
    if r <= 0:
        return 1e50
    term1 = V_param**2 * r**2
    term2 = 6 * V_param * KAPPA**2 * theta * (eta + 4) * r
    term3 = 9 * eta**2 + 72 * eta + 108
    return (2 * PI**2 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


# --- T³×S¹ ---
# Note: Isotropic scaling R1 = R2 = R3 = r

def V_T3_FULL(r, V_param, eta, theta):
    """
    T³xS¹ with FULL Nieh-Yan (Isotropic: R1 = R2 = R3 = r)
    
    Original from engine log:
    V(r) = 16*pi⁴*L*R2*R3*(R1²V² + 6*R1*V*η*κ²*θ + 9η²) / (3*R1*κ²)
    
    With R2 = R3 = R1 = r (isotropic):
    V(r) = 16*pi⁴*L*r² * (r²V² + 6*r*V*η*κ²*θ + 9η²) / (3*r*κ²)
         = 16*pi⁴*L / (3κ²) * r * (r²V² + 6*r*V*η*κ²*θ + 9η²)
         = 16*pi⁴*L / (3κ²) * (V²r³ + 6*V*η*κ²*θ*r² + 9η²*r)
    """
    if r <= 0:
        return 1e50
    prefactor = 16 * PI**4 * L / (3 * KAPPA**2)
    term1 = V_param**2 * r**3
    term2 = 6 * V_param * eta * KAPPA**2 * theta * r**2
    term3 = 9 * eta**2 * r
    return prefactor * (term1 + term2 + term3)


def V_T3_TT(r, V_param, eta, theta):
    """
    T³xS¹ with TT-only Nieh-Yan (Isotropic: R1 = R2 = R3 = r)
    
    Original: V(r) = 16*pi⁴*L*R2*R3*(R1²V² + 12*R1*V*η*κ²*θ + 9η²) / (3*R1*κ²)
    
    With R2 = R3 = R1 = r:
    V(r) = 16*pi⁴*L / (3κ²) * (V²r³ + 12*V*η*κ²*θ*r² + 9η²*r)
    """
    if r <= 0:
        return 1e50
    prefactor = 16 * PI**4 * L / (3 * KAPPA**2)
    term1 = V_param**2 * r**3
    term2 = 12 * V_param * eta * KAPPA**2 * theta * r**2
    term3 = 9 * eta**2 * r
    return prefactor * (term1 + term2 + term3)


def V_T3_REE(r, V_param, eta, theta):
    """
    T³xS¹ with Ree-only Nieh-Yan (same as FULL for T³)
    """
    return V_T3_FULL(r, V_param, eta, theta)


# --- Nil³×S¹ ---
# Note: Using R as the variable radius

def V_Nil3_FULL(r, V_param, eta, theta):
    """
    Nil³xS¹ with FULL Nieh-Yan
    V(r) = 4*pi⁴*L*R*(4R²V² + 8*R*V*κ²*θ*(3η+1) + 36η² - 24η - 9) / (3κ²)
    """
    if r <= 0:
        return 1e50
    term1 = 4 * V_param**2 * r**2
    term2 = 8 * V_param * KAPPA**2 * theta * (3 * eta + 1) * r
    term3 = 36 * eta**2 - 24 * eta - 9
    return (4 * PI**4 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_Nil3_TT(r, V_param, eta, theta):
    """
    Nil³xS¹ with TT-only Nieh-Yan
    V(r) = 4*pi⁴*L*R*(4R²V² + 48*R*V*η*κ²*θ + 36η² - 24η - 9) / (3κ²)
    
    Note: From log, the action was:
    -4*pi**4*L*R*(-4*R**2*V**2 - 48*R*V*eta*kappa**2*theta_NY - 36*eta**2 + 24*eta + 9)/(3*kappa**2)
    So V = -S gives:
    4*pi**4*L*R*(4*R**2*V**2 + 48*R*V*eta*kappa**2*theta_NY + 36*eta**2 - 24*eta - 9)/(3*kappa**2)
    """
    if r <= 0:
        return 1e50
    term1 = 4 * V_param**2 * r**2
    term2 = 48 * V_param * eta * KAPPA**2 * theta * r
    term3 = 36 * eta**2 - 24 * eta - 9
    return (4 * PI**4 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


def V_Nil3_REE(r, V_param, eta, theta):
    """
    Nil³xS¹ with Ree-only Nieh-Yan
    V(r) = 4*pi⁴*L*R*(4R²V² + 8*R*V*κ²*θ*(3η-1) + 36η² - 24η - 9) / (3κ²)
    """
    if r <= 0:
        return 1e50
    term1 = 4 * V_param**2 * r**2
    term2 = 8 * V_param * KAPPA**2 * theta * (3 * eta - 1) * r
    term3 = 36 * eta**2 - 24 * eta - 9
    return (4 * PI**4 * L * r * (term1 + term2 + term3)) / (3 * KAPPA**2)


# ==============================================================================
# Potential Function Registry
# ==============================================================================

POTENTIAL_FUNCTIONS = {
    ('S3', 'FULL'): V_S3_FULL,
    ('S3', 'TT'): V_S3_TT,
    ('S3', 'REE'): V_S3_REE,
    ('T3', 'FULL'): V_T3_FULL,
    ('T3', 'TT'): V_T3_TT,
    ('T3', 'REE'): V_T3_REE,
    ('Nil3', 'FULL'): V_Nil3_FULL,
    ('Nil3', 'TT'): V_Nil3_TT,
    ('Nil3', 'REE'): V_Nil3_REE,
}


# ==============================================================================
# Stability Analysis
# ==============================================================================

def analyze_stability(potential_func, V_param, eta, theta):
    """
    Analyze the stability of a potential configuration.
    
    Returns:
        r0: Stable radius (None if unstable)
        delta_V: Barrier height or well depth
        stability_type: 'type-I', 'type-II', or 'type-III'
    
    Stability types:
        - 'type-I': Local minimum with barrier (V increases toward r=0)
        - 'type-II': Local minimum, V decreases from r=0 toward minimum
                    (universe spontaneously rolls to stable radius)
        - 'type-III': No local minimum in physical region
    
    Physical interpretation:
        For potentials of form V(r) ~ r * (A*r² + B*r + C):
        - C > 0: V increases from r=0 (classical barrier)
        - C < 0: V decreases from r=0 (rolling/spontaneous nucleation)
    """
    
    def v(r):
        return potential_func(r, V_param, eta, theta)
    
    # Find minimum
    try:
        res_min = minimize_scalar(v, bounds=(R_MIN, R_MAX), method='bounded')
    except Exception:
        return None, None, 'type-III'
    
    if not res_min.success:
        return None, None, 'type-III'
    
    r0 = res_min.x
    v_min = res_min.fun
    
    # Check if hitting bounds
    if r0 < R_BOUNDARY_THRESHOLD or r0 > R_MAX - R_BOUNDARY_THRESHOLD:
        return None, None, 'type-III'
    
    # Check curvature (second derivative) - must be positive for minimum
    h = 1e-5
    try:
        d2v = (v(r0 + h) - 2 * v(r0) + v(r0 - h)) / h**2
    except Exception:
        return None, None, 'type-III'
    
    if d2v <= 0:
        return None, None, 'type-III'
    
    # === Determine barrier type based on slope at r=0 ===
    # Check the derivative dV/dr at small r
    # If dV/dr > 0 near r=0: potential increases -> classical barrier
    # If dV/dr < 0 near r=0: potential decreases -> rolling
    
    r_test = R_MIN * 2  # Small but not too small
    dr = R_MIN * 0.1
    
    try:
        v1 = v(r_test)
        v2 = v(r_test + dr)
        slope_near_zero = (v2 - v1) / dr
    except Exception:
        slope_near_zero = 0
    
    # Also sample to find barrier height
    r_samples = np.linspace(R_MIN, r0, 30)
    v_samples = [v(r) for r in r_samples]
    v_max_before_min = max(v_samples)
    
    if slope_near_zero > 0:
        # Potential increases from r=0: classical barrier exists
        stability_type = 'type-I'
        # Barrier height is max value before minimum minus minimum
        delta_V = v_max_before_min - v_min
    else:
        # Potential decreases from r=0: rolling behavior
        # No barrier at origin, universe spontaneously nucleates
        stability_type = 'type-II'
        # "Depth" is from V(r~0) down to V(r0)
        # Since potential goes V(0)~0 -> decreases -> v_min -> increases
        delta_V = v_samples[0] - v_min  # Height from near-origin to bottom
        if delta_V < 0:
            delta_V = abs(v_min)  # Fallback: just report well depth
    
    return r0, delta_V, stability_type


# ==============================================================================
# Single Point Calculation (for multiprocessing)
# ==============================================================================

def calculate_single_point(args):
    """
    Calculate stability for a single parameter point.
    Used for multiprocessing.
    """
    topology, ny_variant, V_param, eta, theta = args
    V_param = round(V_param, 2)
    eta = round(eta, 2)
    theta = round(theta, 2)

    potential_func = POTENTIAL_FUNCTIONS[(topology, ny_variant)]
    r0, delta_V, stability_type = analyze_stability(potential_func, V_param, eta, theta)
    
    return {
        'topology': topology,
        'ny_variant': ny_variant,
        'V': V_param,
        'eta': eta,
        'theta': theta,
        'r0': r0,
        'delta_V': delta_V,
        'stability_type': stability_type
    }


# ==============================================================================
# Parameter Grid Generation
# ==============================================================================

def generate_parameter_grid(V_range, eta_range, theta_range, topologies, ny_variants):
    """
    Generate all parameter combinations for scanning.
    """
    grid = []
    for topology in topologies:
        for ny_variant in ny_variants:
            for V_param in V_range:
                for eta in eta_range:
                    for theta in theta_range:
                        grid.append((topology, ny_variant, V_param, eta, theta))
    return grid


# ==============================================================================
# Main Scan Function
# ==============================================================================

def run_scan(V_points=6, eta_points=11, theta_points=6,
             V_min=0.0, V_max=5.0,
             eta_min=-5.0, eta_max=5.0,
             theta_min=0.0, theta_max=5.0,
             topologies=None, ny_variants=None,
             output_dir='output', n_workers=None):
    """
    Run the full parameter scan.
    """
    
    if topologies is None:
        topologies = ['S3', 'T3', 'Nil3']
    if ny_variants is None:
        ny_variants = ['FULL', 'TT', 'REE']
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)
    
    # Generate parameter ranges
    # Note: Avoid V=0 exactly to prevent division issues
    V_range = np.linspace(V_min, V_max, V_points)
    eta_range = np.linspace(eta_min, eta_max, eta_points)
    theta_range = np.linspace(theta_min, theta_max, theta_points)
    
    print("=" * 70)
    print("DPPUv2 Parameter Scan")
    print("=" * 70)
    print(f"V range: [{V_min}, {V_max}] with {V_points} points")
    print(f"η range: [{eta_min}, {eta_max}] with {eta_points} points")
    print(f"θ range: [{theta_min}, {theta_max}] with {theta_points} points")
    print(f"Topologies: {topologies}")
    print(f"NY Variants: {ny_variants}")
    print(f"Workers: {n_workers}")
    
    # Generate grid
    print("\nGenerating parameter grid...")
    grid = generate_parameter_grid(V_range, eta_range, theta_range, topologies, ny_variants)
    total_points = len(grid)
    print(f"Total points: {total_points:,}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Run scan with multiprocessing
    print(f"\nStarting scan with {n_workers} workers...")
    start_time = time.time()
    
    results = []
    
    with Pool(n_workers) as pool:
        # Process in chunks for progress reporting
        chunk_size = 10000
        for i in range(0, total_points, chunk_size):
            chunk = grid[i:i + chunk_size]
            chunk_results = pool.map(calculate_single_point, chunk)
            results.extend(chunk_results)
            
            elapsed = time.time() - start_time
            progress = (i + len(chunk)) / total_points * 100
            eta_remaining = elapsed / (i + len(chunk)) * (total_points - i - len(chunk))
            print(f"  Progress: {progress:.1f}% ({i + len(chunk):,}/{total_points:,}) "
                  f"- Elapsed: {elapsed:.1f}s - ETA: {eta_remaining:.1f}s")
    
    total_time = time.time() - start_time
    print(f"\nScan complete. Total time: {total_time:.1f}s")
    
    # Save results
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Save combined CSV
    combined_file = os.path.join(output_dir, f'dppu_scan_all_{timestamp}.csv')
    print(f"\nSaving combined results to: {combined_file}")
    
    with open(combined_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['topology', 'ny_variant', 'V', 'eta', 'theta',
                                                'r0', 'delta_V', 'stability_type'])
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    
    # Save separate files per topology/variant
    for topology in topologies:
        for ny_variant in ny_variants:
            subset = [r for r in results if r['topology'] == topology and r['ny_variant'] == ny_variant]
            filename = os.path.join(output_dir, f'dppu_scan_{topology}_{ny_variant}_{timestamp}.csv')
            
            with open(filename, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['topology', 'ny_variant', 'V', 'eta', 'theta', 'r0', 'delta_V', 'stability_type'])
                writer.writeheader()
                for row in subset:
                    writer.writerow(row)
            
            # Count statistics
            n_stable = sum(1 for r in subset if r['stability_type'] in ['type-I', 'type-II'])
            n_type_I = sum(1 for r in subset if r['stability_type'] == 'type-I')
            n_type_II = sum(1 for r in subset if r['stability_type'] == 'type-II')
            n_type_III = sum(1 for r in subset if r['stability_type'] == 'type-III')
            
            print(f"  {topology}-{ny_variant}: "
                  f"stable={n_stable:,} (type-I={n_type_I:,}, type-II={n_type_II:,}), "
                  f"type-III={n_type_III:,}")
    
    print(f"\nOutput files saved to: {output_dir}/")
    print("Done.")
    
    return results


# ==============================================================================
# CLI Entry Point
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='DPPUv2 High-Resolution Parameter Scan',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full scan (default parameters from specification)
  python DPPUv2_parameter_scan_v3.py
  
  # Quick test with reduced resolution
  python DPPUv2_parameter_scan_v3.py --V-points 20 --eta-points 30 --theta-points 10
  
  # Single topology
  python DPPUv2_parameter_scan_v3.py --topologies S3
  
  # Custom output directory
  python DPPUv2_parameter_scan_v3.py --output-dir ./my_results
        """
    )
    
    # Parameter range arguments
    parser.add_argument('--V-points', type=int, default=51, help='Number of V points (default: 51)')
    parser.add_argument('--eta-points', type=int, default=151, help='Number of η points (default: 151)')
    parser.add_argument('--theta-points', type=int, default=51, help='Number of θ points (default: 51)')
    
    parser.add_argument('--V-min', type=float, default=0.0, help='Minimum V (default: 0.0)')
    parser.add_argument('--V-max', type=float, default=5.0, help='Maximum V (default: 5.0)')
    parser.add_argument('--eta-min', type=float, default=-10.0, help='Minimum η (default: -10.0)')
    parser.add_argument('--eta-max', type=float, default=5.0, help='Maximum η (default: 5.0)')
    parser.add_argument('--theta-min', type=float, default=0.0, help='Minimum θ (default: 0.0)')
    parser.add_argument('--theta-max', type=float, default=5.0, help='Maximum θ (default: 5.0)')
    
    # Topology/variant selection
    parser.add_argument('--topologies', nargs='+', default=['S3', 'T3', 'Nil3'],
                        choices=['S3', 'T3', 'Nil3'], help='Topologies to scan')
    parser.add_argument('--ny-variants', nargs='+', default=['FULL', 'TT', 'REE'],
                        choices=['FULL', 'TT', 'REE'], help='NY variants to scan')
    
    # Output and parallelization
    parser.add_argument('--output-dir', type=str, default='output', help='Output directory')
    parser.add_argument('--workers', type=int, default=None, help='Number of parallel workers')
    
    args = parser.parse_args()
    
    run_scan(
        V_points=args.V_points,
        eta_points=args.eta_points,
        theta_points=args.theta_points,
        V_min=args.V_min,
        V_max=args.V_max,
        eta_min=args.eta_min,
        eta_max=args.eta_max,
        theta_min=args.theta_min,
        theta_max=args.theta_max,
        topologies=args.topologies,
        ny_variants=args.ny_variants,
        output_dir=args.output_dir,
        n_workers=args.workers
    )


if __name__ == '__main__':
    main()
