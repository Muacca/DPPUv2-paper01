#!/usr/bin/env python3
"""
DPPUv2 Runner for S3xS1 Topology (v3.0)
========================================

Purpose: Implements topology-specific setup for S3xS1 manifold with v3 engine.

Topology: S3 (3-sphere) x S1 (circle)
Structure Constants: C^i_jk = (4/r)eps_ijk  (SU(2) Lie algebra)
Volume: V = 2pi^2 L r^3
Curvature: R_LC = 24/r^2 (for eta=0)

Changes from v2.0:
- Import from DPPUv2_engine_core_v3 (strict antisymmetry check)
- CLI: --mode AX/VT/MX (NO backward compatibility with 1C/1D)
- CLI: --ny-variant TT/REE/FULL (BOTH removed)
- Optional checkpoint via --checkpoint-dir

Author: Muacca
Version: 3.0
Date: 2025-12-14
"""

import sys
from pathlib import Path
import argparse

from sympy import symbols, Matrix, S, pi
from sympy.tensor.array import MutableDenseNDimArray

# Import v3 engine
from DPPUv2_engine_core_v3 import (
    BaseFrameEngine, Mode, NyVariant,
    ComputationLogger, CheckpointManager, epsilon_symbol
)

# ============================================================
# S³×S¹ Engine Implementation
# ============================================================

class S3S1Engine(BaseFrameEngine):
    """
    S3xS1 Topology Implementation (v3.0)
    Fully compatible with SPEC v3 requirements
    """
    
    def step_E4_1_setup(self):
        """Setup parameters and symbolic variables for S3xS1"""
        self._log_info(f"Mode: {self.mode.value}")
        self._log_info(f"NY Variant: {self.ny_variant.value}")
        self._log_info("Topology: S3xS1")
        self._log_info("Basis: Strictly Orthonormal Frame (Anholonomic)")
        
        r = symbols('r', positive=True, real=True)
        L = symbols('L', positive=True, real=True)
        kappa = symbols('kappa', positive=True, real=True)
        theta_NY = symbols('theta_NY', real=True)
        
        # Torsion parameters based on mode
        if self.mode == Mode.AX:
            # Axial-only: S^μ ≠ 0, T_μ = 0
            eta = symbols('eta', real=True)
            V = S.Zero
            q = 2 * eta
            self._log_info("Mode AX: Axial-only (T1)")
        elif self.mode == Mode.VT:
            # Vector-trace-only: T_μ ≠ 0, S^μ = 0
            eta = S.Zero
            V = symbols('V', real=True, positive=True)
            q = S.Zero
            self._log_info("Mode VT: Vector-trace-only (T2)")
        elif self.mode == Mode.MX:
            # Mixed: both non-zero
            eta = symbols('eta', real=True)
            V = symbols('V', real=True, positive=True)
            q = 2 * eta
            self._log_info("Mode MX: Mixed (T1 + T2)")
        else:
            raise ValueError(f"Unknown mode: {self.mode}")
        
        self.data['params'] = {
            'r': r, 'L': L, 'kappa': kappa, 'theta_NY': theta_NY,
            'eta': eta, 'V': V, 'q': q
        }
        self.data['dim'] = 4
        self._log_success("Setup complete")
    
    def step_E4_2_metric_and_frame(self):
        """Define Frame Metric, Volume, and Structure Constants for S3xS1"""
        self._log_info("Setting up Frame Metric and Structure Constants...")
        
        dim = self.data['dim']
        r = self.data['params']['r']
        L = self.data['params']['L']
        
        # Frame metric (Identity for orthonormal frame)
        metric = Matrix.eye(dim)
        self.data['metric_frame'] = metric
        
        # Total volume: V = 2π²Lr³
        total_volume = 2 * pi**2 * L * r**3
        self.data['total_volume'] = total_volume
        
        # Structure constants: C^i_jk = (4/r)ε_ijk for i,j,k ∈ {0,1,2}
        C = MutableDenseNDimArray.zeros(dim, dim, dim)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    eps = epsilon_symbol(i, j, k)
                    if eps != 0:
                        C[i, j, k] = 4 * eps / r
        
        self.data['structure_constants'] = C
        self._log_success("Metric and structure constants defined")


# ============================================================
# Main CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="DPPUv2 S3xS1 Runner v3.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python DPPUv2_runner_S3S1_v3.py --mode AX --ny-variant TT
  python DPPUv2_runner_S3S1_v3.py --mode MX --ny-variant FULL --checkpoint-dir checkpoints
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--mode',
        type=str,
        required=True,
        choices=['AX', 'VT', 'MX'],
        help='Torsion mode: AX (Axial), VT (Vector-trace), MX (Mixed)'
    )
    
    parser.add_argument(
        '--ny-variant',
        type=str,
        required=True,
        choices=['TT', 'REE', 'FULL'],
        help='Nieh-Yan variant: TT (TT-only), REE (Ree-only), FULL (TT-Ree)'
    )
    
    # Optional arguments
    parser.add_argument(
        '--checkpoint-dir',
        type=str,
        default=None,
        help='Enable checkpoints and save to this directory (default: OFF)'
    )
    
    parser.add_argument(
        '--log-file',
        type=str,
        default='dppu_s3s1_v3.log',
        help='Log file path (default: dppu_s3s1_v3.log)'
    )
    
    parser.add_argument(
        '--start-step',
        type=str,
        default='E4.1',
        help='Starting step (default: E4.1)'
    )
    
    args = parser.parse_args()
    
    # Convert string arguments to Enums
    mode = Mode[args.mode]
    ny_variant = NyVariant[args.ny_variant]
    
    # Setup logger
    logger = ComputationLogger(log_file=args.log_file)
    
    # Setup checkpoint manager
    checkpoint_enabled = (args.checkpoint_dir is not None)
    ckpt_mgr = CheckpointManager(
        output_dir=args.checkpoint_dir if checkpoint_enabled else "checkpoints_s3s1_v3",
        enabled=checkpoint_enabled
    )
    
    # Display configuration
    print("="*70)
    print("DPPUv2 S3xS1 Runner v3.0")
    print("="*70)
    print(f"Mode:         {mode.value}")
    print(f"NY Variant:   {ny_variant.value}")
    print(f"Checkpoints:  {'ENABLED' if checkpoint_enabled else 'DISABLED'}")
    if checkpoint_enabled:
        print(f"  Directory:  {args.checkpoint_dir}")
    print(f"Log File:     {args.log_file}")
    print("="*70)
    print()
    
    # Create and run engine
    try:
        engine = S3S1Engine(
            mode=mode,
            ny_variant=ny_variant,
            logger=logger,
            checkpoint_mgr=ckpt_mgr
        )
        
        engine.run(start_step=args.start_step)
        
        print("\n" + "="*70)
        print("SUCCESS: Computation completed without errors")
        print("="*70)
        print(f"Results saved to: {args.log_file}")
        if checkpoint_enabled:
            print(f"Checkpoints saved to: {args.checkpoint_dir}/")
        
    except Exception as e:
        print("\n" + "="*70)
        print("FAILED: Computation terminated with error")
        print("="*70)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
