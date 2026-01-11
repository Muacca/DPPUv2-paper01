#!/usr/bin/env python3
"""
DPPUv2 Runner for T3xS1 Topology (v3.0)
========================================

Purpose: Implements topology-specific setup for T3xS1 (flat) manifold with v3 engine.

Topology: T3 (3-torus) x S1 (circle)
Structure Constants: All zero (Abelian)
Volume: V = (2π)⁴ L R₁R₂R₃
Curvature: R_LC = 0 (flat)

Changes from v2.0:
- Import from DPPUv2_engine_core_v3 (strict antisymmetry check)
- CLI: --mode AX/VT/MX (NO backward compatibility)
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
    ComputationLogger, CheckpointManager
)

# ============================================================
# T³×S¹ Engine Implementation
# ============================================================

class T3S1Engine(BaseFrameEngine):
    """
    T3xS1 Topology Implementation (v3.0)
    Flat manifold with zero structure constants
    """
    
    def step_E4_1_setup(self):
        """Setup parameters and symbolic variables for T3xS1"""
        self._log_info(f"Mode: {self.mode.value}")
        self._log_info(f"NY Variant: {self.ny_variant.value}")
        self._log_info("Topology: T3xS1")
        self._log_info("Basis: Orthonormal Frame (Abelian)")
        
        R1 = symbols('R1', positive=True, real=True)
        R2 = symbols('R2', positive=True, real=True)
        R3 = symbols('R3', positive=True, real=True)
        L = symbols('L', positive=True, real=True)
        kappa = symbols('kappa', positive=True, real=True)
        theta_NY = symbols('theta_NY', real=True)
        
        # Torsion parameters based on mode
        if self.mode == Mode.AX:
            eta = symbols('eta', real=True)
            V = S.Zero
            q = 2 * eta
            self._log_info("Mode AX: Axial-only (T1)")
        elif self.mode == Mode.VT:
            eta = S.Zero
            V = symbols('V', real=True, positive=True)
            q = S.Zero
            self._log_info("Mode VT: Vector-trace-only (T2)")
        elif self.mode == Mode.MX:
            eta = symbols('eta', real=True)
            V = symbols('V', real=True, positive=True)
            q = 2 * eta
            self._log_info("Mode MX: Mixed (T1 + T2)")
        else:
            raise ValueError(f"Unknown mode: {self.mode}")
        
        self.data['params'] = {
            'R1': R1, 'R2': R2, 'R3': R3, 'L': L,
            'kappa': kappa, 'theta_NY': theta_NY,
            'eta': eta, 'V': V, 'q': q,
            'r': R1  # For compatibility with potential function
        }
        self.data['dim'] = 4
        self._log_success("Setup complete")
    
    def step_E4_2_metric_and_frame(self):
        """Define Frame Metric, Volume, and Structure Constants for T3xS1"""
        self._log_info("Setting up Frame Metric and Structure Constants...")
        
        dim = self.data['dim']
        R1 = self.data['params']['R1']
        R2 = self.data['params']['R2']
        R3 = self.data['params']['R3']
        L = self.data['params']['L']
        
        # Frame metric (Identity for orthonormal frame)
        metric = Matrix.eye(dim)
        self.data['metric_frame'] = metric
        
        # Total volume: V = (2π)⁴ L R₁R₂R₃
        total_volume = (2*pi)**4 * L * R1 * R2 * R3
        self.data['total_volume'] = total_volume
        
        # Structure constants: All zero (Abelian/flat)
        C = MutableDenseNDimArray.zeros(dim, dim, dim)
        self.data['structure_constants'] = C
        
        self._log_success("Metric and structure constants defined (flat manifold)")


# ============================================================
# Main CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="DPPUv2 T³ x S¹ Runner v3.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python DPPUv2_runner_T3S1_v3.py --mode AX --ny-variant TT
  python DPPUv2_runner_T3S1_v3.py --mode MX --ny-variant FULL --checkpoint-dir checkpoints
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
        default='dppu_t3s1_v3.log',
        help='Log file path (default: dppu_t3s1_v3.log)'
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
        output_dir=args.checkpoint_dir if checkpoint_enabled else "checkpoints_t3s1_v3",
        enabled=checkpoint_enabled
    )
    
    # Display configuration
    print("="*70)
    print("DPPUv2 T3xS1 Runner v3.0")
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
        engine = T3S1Engine(
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
