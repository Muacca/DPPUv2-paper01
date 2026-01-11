#!/usr/bin/env python3
"""
DPPUv2 Runner for Nil³ x S¹ Topology (v3.0) - CONVENTIONS Compliant
==================================================================

Purpose: Implements topology-specific setup for Nil³ x S¹ (Heisenberg) manifold with v3 engine.

Topology: Nil³ (Heisenberg nilmanifold / Bianchi II) x S¹ (circle)

CONVENTIONS Compliance:
-----------------------
This runner implements CONVENTIONS-compliant definitions:

1. Structure Constants (CONVENTIONS 3.2):
   [E_b, E_c] = -C^a_{bc} E_a
   
   For Heisenberg: [E_0, E_1] = (1/R) E_2
   => C^2_{01} = -1/R, C^2_{10} = +1/R

2. Levi-Civita Connection (CONVENTIONS 5):
   Left-invariant frame allows the general Koszul formula to correctly
   compute the Levi-Civita connection:
   
   Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} - C^b_{ac})
   
   This formula works for Nil³ (non bi-invariant) because it's derived
   from the general Koszul identity without assuming bi-invariance.

3. Metric Compatibility (CONVENTIONS 4):
   The computed connection automatically satisfies Γ_{abc} + Γ_{bac} = 0
   (verified by engine's internal check).

4. Riemann Antisymmetry (CONVENTIONS 8.3):
   With the general Koszul connection, R_{abcd} = -R_{bacd} holds.

Mathematical Background:
------------------------
Nil³ is the 3-dimensional Heisenberg group with left-invariant metric.
Unlike S³ (SU(2)), Nil³ does NOT admit a bi-invariant metric because
the Heisenberg group is nilpotent (not semisimple).

The general Koszul formula (CONVENTIONS 5) correctly handles this case
because it works for any left-invariant frame, not just bi-invariant ones.

Curvature:
---------
For torsion-free case (eta=0, V=0), the Levi-Civita Ricci scalar is negative.
With torsion, the total Ricci scalar R depends on torsion parameters.
See engine output (step E4.8) for the explicit formula.

Author: Muacca
Version: 3.0 (CONVENTIONS compliant, general Koszul)
Date: 2025-12-14
"""

import sys
from pathlib import Path
import argparse

from sympy import symbols, Matrix, S, pi, Rational, simplify
from sympy.tensor.array import MutableDenseNDimArray

# Import v3 engine
from DPPUv2_engine_core_v3 import (
    BaseFrameEngine, Mode, NyVariant,
    ComputationLogger, CheckpointManager
)

# ============================================================
# Nil³×S¹ Engine Implementation (CONVENTIONS Compliant)
# ============================================================

class Nil3S1Engine(BaseFrameEngine):
    """
    Nil³×S¹ Topology Implementation (v3.0) - CONVENTIONS Compliant
    
    Heisenberg nilmanifold using general Koszul formula from engine.
    
    Key Features:
    - Nil³ is NOT bi-invariant (nilpotent group)
    - Left-invariant frame allows general Koszul formula (CONVENTIONS 5)
    - No step_E4_3 override needed - engine handles it correctly
    """
    
    def step_E4_1_setup(self):
        """Setup parameters and symbolic variables for Nil3xS1"""
        self._log_info(f"Mode: {self.mode.value}")
        self._log_info(f"NY Variant: {self.ny_variant.value}")
        self._log_info("Topology: Nil3xS1 (Heisenberg / Bianchi II)")
        self._log_info("Basis: Orthonormal Frame (Left-Invariant)")
        self._log_info("")
        self._log_info("CONVENTIONS Compliance:")
        self._log_info("  - Structure constants: CONVENTIONS 3.2")
        self._log_info("  - Connection: General Koszul (CONVENTIONS 5)")
        self._log_info("  - Metric compatibility: Auto-verified by engine")
        
        # Scale parameter R for Nil³
        R = symbols('R', positive=True, real=True)
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
            'R': R, 'L': L, 'kappa': kappa, 'theta_NY': theta_NY,
            'eta': eta, 'V': V, 'q': q,
            'r': R  # Alias for compatibility
        }
        self.data['dim'] = 4
        self._log_success("Setup complete")
    
    def step_E4_2_metric_and_frame(self):
        """Define Frame Metric, Volume, and Structure Constants for Nil3xS1"""
        self._log_info("Setting up Frame Metric and Structure Constants...")
        
        dim = self.data['dim']
        R = self.data['params']['R']
        L = self.data['params']['L']
        
        # Frame metric (Identity for orthonormal frame)
        metric = Matrix.eye(dim)
        self.data['metric_frame'] = metric
        self._log_info("Frame metric: eta_ab = diag(1,1,1,1)")
        
        # Total volume for Nil3xS1
        # Standard compact quotient: (2pi)^4 L R^3
        total_volume = (2*pi)**4 * L * R**3
        self.data['total_volume'] = total_volume
        self._log_info(f"Total volume: (2pi)^4 L R^3")
        
        # Structure constants (CONVENTIONS 3.2)
        # [E_b, E_c] = -C^a_{bc} E_a
        # Heisenberg: [E_0, E_1] = (1/R) E_2
        # => -C^2_{01} = 1/R => C^2_{01} = -1/R
        
        C = MutableDenseNDimArray.zeros(dim, dim, dim)
        lam = 1 / R
        
        # CONVENTIONS 3.2 compliant
        C[2, 0, 1] = -lam   # C^2_{01} = -1/R
        C[2, 1, 0] = +lam   # C^2_{10} = +1/R (antisymmetry)
        
        self.data['structure_constants'] = C
        
        self._log_info("Structure constants (CONVENTIONS 3.2):")
        self._log_info(f"  [E_0, E_1] = (1/R) E_2")
        self._log_info(f"  => C^2_01 = -1/R, C^2_10 = +1/R")
        
        # Verify antisymmetry
        self._verify_structure_constant_antisymmetry(C, dim)
        
        self._log_success("Metric and structure constants defined")
    
    def _verify_structure_constant_antisymmetry(self, C, dim):
        """Verify C^a_{bc} = -C^a_{cb} (CONVENTIONS 8.1)"""
        self._log_info("Verifying structure constant antisymmetry...")
        violations = 0
        for a in range(dim):
            for b in range(dim):
                for c in range(b+1, dim):
                    check = simplify(C[a,b,c] + C[a,c,b])
                    if check != 0:
                        self._log_warning(f"C^{a}_{b}{c} + C^{a}_{c}{b} = {check}")
                        violations += 1
        if violations == 0:
            self._log_success("Structure constant antisymmetry: PASSED")
        else:
            raise ValueError(f"Structure constant antisymmetry violated: {violations} cases")


# ============================================================
# Main CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="DPPUv2 Nil³ x S¹ (Heisenberg) Runner v3.0 - CONVENTIONS Compliant",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
CONVENTIONS Compliance:
  This runner implements CONVENTIONS-compliant definitions for Nil³.
  
  Left-invariant frame allows general Koszul formula (CONVENTIONS 5):
  - Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} - C^b_{ac})
  - Works for Nil³ despite being non bi-invariant
  - No manual override needed - engine handles it correctly
  - Metric compatibility verified automatically

Topology Properties:
  - Heisenberg group (Bianchi Type II)
  - Structure: [E_0, E_1] = (1/R) E_2
  - Background curvature: R_LC = -1/(2R²) < 0 (negative)
  - Parallelizable (global frame exists)

Examples:
  python DPPUv2_runner_Nil3S1_v3.py --mode AX --ny-variant TT
  python DPPUv2_runner_Nil3S1_v3.py --mode VT --ny-variant FULL
  python DPPUv2_runner_Nil3S1_v3.py --mode MX --ny-variant FULL --checkpoint-dir checkpoints
        """
    )
    
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
    
    parser.add_argument(
        '--checkpoint-dir',
        type=str,
        default=None,
        help='Enable checkpoints and save to this directory (default: OFF)'
    )
    
    parser.add_argument(
        '--log-file',
        type=str,
        default=None,
        help='Log file path (default: auto-generated)'
    )
    
    parser.add_argument(
        '--start-step',
        type=str,
        default='E4.1',
        help='Starting step (default: E4.1)'
    )
    
    args = parser.parse_args()
    
    mode = Mode[args.mode]
    ny_variant = NyVariant[args.ny_variant]
    
    if args.log_file is None:
        from datetime import datetime
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_filename = f"dppu_nil3s1_{args.mode}_{args.ny_variant}_{timestamp}.log"
    else:
        log_filename = args.log_file
    
    logger = ComputationLogger(log_file=log_filename)
    
    checkpoint_enabled = (args.checkpoint_dir is not None)
    ckpt_mgr = CheckpointManager(
        output_dir=args.checkpoint_dir if checkpoint_enabled else "checkpoints_nil3s1_v3",
        enabled=checkpoint_enabled
    )
    
    print("="*70)
    print("DPPUv2 Nil3xS1 Runner v3.0 (CONVENTIONS Compliant)")
    print("="*70)
    print(f"Mode:         {mode.value}")
    print(f"NY Variant:   {ny_variant.value}")
    print(f"Checkpoints:  {'ENABLED' if checkpoint_enabled else 'DISABLED'}")
    print(f"Log File:     {log_filename}")
    print("="*70)
    print()
    print("CONVENTIONS Compliance:")
    print("  - Structure constants: sec.3.2 (antisymmetric)")
    print("  - Connection: General Koszul (sec.5)")
    print("  - Metric compatibility: sec.4 (auto-verified)")
    print("  - Riemann antisymmetry: sec.8.3 (verified)")
    print()
    print("Topology Properties:")
    print("  - Heisenberg group (Bianchi II, nilpotent)")
    print("  - NOT bi-invariant (unlike S3)")
    print("  - R_LC = -1/(2R^2) < 0 (negative curvature)")
    print("="*70)
    print()
    
    try:
        engine = Nil3S1Engine(
            mode=mode,
            ny_variant=ny_variant,
            logger=logger,
            checkpoint_mgr=ckpt_mgr
        )
        
        engine.run(start_step=args.start_step)
        
        print("\n" + "="*70)
        print("SUCCESS: Computation completed without errors")
        print("="*70)
        print(f"Results saved to: {log_filename}")
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
