#!/usr/bin/env python3
"""
DPPUv2 Engine Core v3.0: Refactored for Paper Publication
==========================================================

MAJOR CHANGES from v2.0:
1. Mode names: AX (Axial), VT (Vector-trace), MX (Mixed) 
2. NY variant: TT, REE, FULL
3. Riemann antisymmetry check: STRICT 3-level verification (PROVED_ZERO/WITNESS_NONZERO/UNPROVED)
4. Infrastructure separation: logger/checkpoint are optional
5. Enum-based mode management

Purpose: Provides computation logic for Einstein-Cartan + Nieh-Yan
         with strict self-checking for paper publication.

References:
- Chandia & Zanelli (1997): Phys. Rev. D 55, 7580 [hep-th/9702025]
- Hehl et al. (1976): Rev. Mod. Phys. 48, 393
- Nieh & Yan (1982): J. Math. Phys. 23, 373

Author: Muacca
Version: 3.0
Date: 2025-12-14
"""

import sys
import os
from pathlib import Path
import pickle
from datetime import datetime
import traceback
from typing import Dict, Any, Tuple, Optional, List
from enum import Enum

import numpy as np

from sympy import (symbols, Matrix, simplify, sqrt, pi, S, Rational, expand,
                   cancel, together, trigsimp, factor, diff, lambdify, Symbol)
from sympy.tensor.array import MutableDenseNDimArray
from sympy.simplify.fu import fu

# ============================================================
# Enums for Type-Safe Mode Management
# ============================================================

class Mode(Enum):
    """Torsion ansatz mode"""
    AX = "AX"    # Axial-only (S^μ ≠ 0, T_μ = 0)
    VT = "VT"    # Vector-trace-only (T_μ ≠ 0, S^μ = 0)
    MX = "MX"    # Mixed (both non-zero)


class NyVariant(Enum):
    """Nieh-Yan variant selection"""
    TT = "TT"      # TT-only
    REE = "REE"    # Ree-only
    FULL = "FULL"  # TT - Ree


# ============================================================
# Custom Exceptions
# ============================================================

class RiemannAntisymmetryError(Exception):
    """Raised when Riemann tensor violates antisymmetry."""
    
    def __init__(self, violation_type: str, violations: List[Dict]):
        self.violation_type = violation_type
        self.violations = violations
        
        # Build detailed error message
        msg_lines = ["Riemann tensor antisymmetry violated", ""]
        msg_lines.append(f"Violation Type: {violation_type}")
        
        if violations:
            msg_lines.append(f"Violated components (first 3):")
            for v in violations[:3]:
                indices = v['indices']
                antisym = v['antisymmetry']
                residual = v['residual']
                
                msg_lines.append(f"  [{indices[0]},{indices[1]},{indices[2]},{indices[3]}]: "
                               f"residual = {residual}")
                
                if 'witness_point' in v:
                    witness = v['witness_point']
                    value = v['witness_value']
                    point_str = ', '.join([f"{k}={v:.3f}" for k, v in witness.items()])
                    msg_lines.append(f"    (witness at {point_str} → value={value:.3e})")
            
            msg_lines.append("")
            msg_lines.append(f"Total violations: {len(violations)}")
        
        super().__init__('\n'.join(msg_lines))


# ============================================================
# Helper Functions
# ============================================================

def epsilon_symbol(i: int, j: int, k: int) -> int:
    """Levi-Civita symbol for indices 0,1,2"""
    if (i, j, k) in [(0,1,2), (1,2,0), (2,0,1)]:
        return 1
    elif (i, j, k) in [(2,1,0), (0,2,1), (1,0,2)]:
        return -1
    else:
        return 0


def levi_civita_4d(mu: int, nu: int, rho: int, sigma: int) -> int:
    """Levi-Civita symbol for 4D indices (0,1,2,3)"""
    indices = [mu, nu, rho, sigma]
    if len(set(indices)) != 4:
        return 0
    inversions = 0
    for i in range(4):
        for j in range(i+1, 4):
            if indices[i] > indices[j]:
                inversions += 1
    return 1 if inversions % 2 == 0 else -1


# ============================================================
# Riemann Antisymmetry Verification (3-Level)
# ============================================================

def prove_zero(expr, assumptions_dict=None, timeout_seconds=10) -> Tuple[bool, str]:
    """
    Try to prove expr == 0 symbolically using normalization pipeline.
    
    Returns:
        (True, method_name) if PROVED_ZERO
        (False, "UNPROVED") otherwise
    """
    if expr == 0:
        return True, "PROVED_ZERO (trivial)"
    
    if assumptions_dict:
        expr = expr.subs(assumptions_dict)
    
    # Stage 1: Lightweight simplifications
    for name, func in [("simplify", simplify), ("expand", expand)]:
        try:
            result = func(expr)
            if result == 0:
                return True, f"PROVED_ZERO (via {name})"
        except Exception:
            continue
    
    # Stage 2: Heavy simplifications
    for name, func in [("factor", factor), ("cancel", cancel), 
                       ("together", together), ("trigsimp", trigsimp)]:
        try:
            result = func(expr)
            if result == 0:
                return True, f"PROVED_ZERO (via {name})"
        except Exception:
            continue
    
    # Stage 3: Combined pipeline
    try:
        result = trigsimp(simplify(factor(together(expr))))
        if result == 0:
            return True, "PROVED_ZERO (via pipeline)"
    except Exception:
        pass
    
    return False, "UNPROVED"


def find_nonzero_witness(expr, symbols_list, n_points=10, precision=50) -> Tuple[bool, Optional[Dict], Optional[float]]:
    """
    Try to find numerical witness that expr != 0.
    
    Returns:
        (True, witness_point, witness_value) if WITNESS_NONZERO found
        (False, None, None) if all test points give ~0
    """
    try:
        import mpmath
        f = lambdify(symbols_list, expr, modules=['mpmath'])
    except Exception:
        return False, None, None
    
    test_points = generate_test_points(symbols_list, n_points)
    mpmath.mp.dps = precision
    threshold = mpmath.mpf(10) ** (-precision + 10)
    
    for point in test_points:
        try:
            val = f(*point.values())
            if abs(val) > threshold:
                return True, point, float(val)  # Witness found
        except (ValueError, ZeroDivisionError, TypeError):
            continue  # Skip singular points
    
    return False, None, None  # All points ~0


def generate_test_points(symbols_list, n_points):
    """Generate diverse test points avoiding special values."""
    points = []
    for _ in range(n_points):
        point = {}
        for sym in symbols_list:
            if sym.is_positive:
                point[sym] = np.random.uniform(0.1, 10.0)
            else:
                point[sym] = np.random.uniform(-10.0, 10.0)
        points.append(point)
    return points


def verify_antisymmetry_strict(R_abcd, dim, symbols_list=None, logger=None) -> None:
    """
    Strictly verify Riemann antisymmetry with 3-level judgment.
    
    Args:
        R_abcd: Riemann tensor array
        dim: Dimension
        symbols_list: (Deprecated, unused) Symbol list for witness search
        logger: Optional logger
    
    Raises RiemannAntisymmetryError if any violation found.
    
    Note: symbols_list parameter is deprecated and ignored.
          Each residual's free_symbols are used instead for safer witness search.
    """
    if logger:
        logger.info("Verifying Riemann antisymmetry (STRICT 3-level mode)...")
    
    violations = []
    debug_mode = os.getenv('DPPU_DEBUG_RIEMANN', '0') == '1'
    
    # Check ab antisymmetry: R_abcd = -R_bacd
    for a in range(dim):
        for b in range(a+1, dim):
            for c in range(dim):
                for d in range(dim):
                    residual = R_abcd[a,b,c,d] + R_abcd[b,a,c,d]
                    
                    # Stage 1: Try symbolic proof
                    proved, method = prove_zero(residual)
                    if proved:
                        continue  # PASSED
                    
                    # Stage 2: Try numerical witness
                    # Use residual's free_symbols instead of params (safer)
                    residual_symbols = sorted(residual.free_symbols, key=str)
                    has_witness, witness, value = find_nonzero_witness(residual, residual_symbols)
                    
                    if has_witness:
                        violations.append({
                            'type': 'WITNESS_NONZERO',
                            'indices': (a,b,c,d),
                            'antisymmetry': 'ab',
                            'residual': residual,
                            'witness_point': witness,
                            'witness_value': value
                        })
                    else:
                        violations.append({
                            'type': 'UNPROVED',
                            'indices': (a,b,c,d),
                            'antisymmetry': 'ab',
                            'residual': residual
                        })
    
    # Check cd antisymmetry: R_abcd = -R_abdc
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(c+1, dim):
                    residual = R_abcd[a,b,c,d] + R_abcd[a,b,d,c]
                    
                    proved, method = prove_zero(residual)
                    if proved:
                        continue
                    
                    # Use residual's free_symbols instead of params (safer)
                    residual_symbols = sorted(residual.free_symbols, key=str)
                    has_witness, witness, value = find_nonzero_witness(residual, residual_symbols)
                    
                    if has_witness:
                        violations.append({
                            'type': 'WITNESS_NONZERO',
                            'indices': (a,b,c,d),
                            'antisymmetry': 'cd',
                            'residual': residual,
                            'witness_point': witness,
                            'witness_value': value
                        })
                    else:
                        violations.append({
                            'type': 'UNPROVED',
                            'indices': (a,b,c,d),
                            'antisymmetry': 'cd',
                            'residual': residual
                        })
    
    # Report results
    if violations:
        violation_type = violations[0]['type']  # All should be same type
        
        if debug_mode:
            # Development mode: warning only
            if logger:
                logger.warning(f"Riemann violations: {len(violations)} (DEBUG_MODE: continuing)")
                for v in violations[:3]:
                    logger.warning(f"  {v}")
        else:
            # Production mode: raise exception
            raise RiemannAntisymmetryError(violation_type, violations)
    else:
        if logger:
            logger.success("Riemann antisymmetry: PASSED (all PROVED_ZERO)")


# ============================================================
# Infrastructure Classes (Optional)
# ============================================================

class ComputationLogger:
    """Manages step-by-step logging with timestamps"""
    
    def __init__(self, log_file: str = "dppu_engine_v3.log"):
        self.log_file = Path(log_file)
        self.current_step = None
        self.step_start_time = None
        self.total_start_time = datetime.now()
        
        with open(self.log_file, 'w', encoding='utf-8') as f:
            f.write("=" * 70 + "\n")
            f.write("DPPUv2 Engine Core v3.0\n")
            f.write("=" * 70 + "\n")
            f.write(f"Started: {self.total_start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    def step(self, step_id: str, description: str) -> str:
        if self.step_start_time:
            duration = (datetime.now() - self.step_start_time).total_seconds()
            self._log(f"  Duration: {duration:.2f}s\n")
        
        self.current_step = step_id
        self.step_start_time = datetime.now()
        timestamp = self.step_start_time.strftime("%H:%M:%S")
        
        desc_str = description.strip() if description else "No description"
        
        msg = f"[{timestamp}] STEP {step_id}: {desc_str}"
        print("=" * 70)
        print(msg)
        print("=" * 70)
        self._log(msg + "\n")
        
        return step_id
    
    def info(self, message: str, indent: int = 2):
        prefix = " " * indent
        msg = f"{prefix}{message}"
        print(msg)
        self._log(msg + "\n")
    
    def warning(self, message: str):
        msg = f"  [WARNING] {message}"
        print(msg)
        self._log(msg + "\n")
    
    def error(self, message: str):
        msg = f"  [ERROR] {message}"
        print(msg)
        self._log(msg + "\n")
    
    def success(self, message: str):
        msg = f"  [SUCCESS] {message}"
        print(msg)
        self._log(msg + "\n")
    
    def _log(self, message: str):
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(message)
    
    def finalize(self):
        if self.step_start_time:
            duration = (datetime.now() - self.step_start_time).total_seconds()
            self._log(f"  Duration: {duration:.2f}s\n")
        
        total_duration = (datetime.now() - self.total_start_time).total_seconds()
        msg = f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        msg += f"Total runtime: {total_duration:.1f}s ({total_duration/60:.1f} min)\n"
        print(msg)
        self._log(msg)


class CheckpointManager:
    """Manages saving and loading of computation checkpoints"""
    
    def __init__(self, output_dir: str = "checkpoints_v3", enabled: bool = False):
        self.output_dir = Path(output_dir)
        self.enabled = enabled
        
        if self.enabled:
            self.output_dir.mkdir(exist_ok=True)
    
    def save(self, step_id: str, data: Dict, metadata: Optional[Dict] = None) -> Optional[Path]:
        if not self.enabled:
            return None
        
        checkpoint = {
            'step_id': step_id,
            'timestamp': datetime.now().isoformat(),
            'data': data,
            'metadata': metadata or {}
        }
        
        filepath = self.output_dir / f"checkpoint_{step_id.replace('.', '_')}.pkl"
        
        try:
            with open(filepath, 'wb') as f:
                pickle.dump(checkpoint, f, protocol=pickle.HIGHEST_PROTOCOL)
            return filepath
        except Exception as e:
            print(f"  [ERROR] Failed to save checkpoint: {e}")
            return None
    
    def load(self, step_id: str) -> Tuple[Dict, Dict]:
        filepath = self.output_dir / f"checkpoint_{step_id.replace('.', '_')}.pkl"
        
        if not filepath.exists():
            raise FileNotFoundError(f"Checkpoint {step_id} not found at {filepath}")
        
        with open(filepath, 'rb') as f:
            checkpoint = pickle.load(f)
        
        print(f"  [OK] Loaded checkpoint: {step_id}")
        return checkpoint['data'], checkpoint['metadata']


# ============================================================
# Abstract Base Engine v3
# ============================================================

class BaseFrameEngine:
    """
    Abstract Base Class for Einstein-Cartan + Nieh-Yan Calculations (v3.0)
    
    NEW in v3.0:
    - Mode enum: Mode.AX/VT/MX (no backward compatibility)
    - NY variant enum: NyVariant.TT/REE/FULL
    - Strict Riemann antisymmetry check (3-level verification)
    - Optional logger/checkpoint (can be None)
    
    Topology-specific runners must implement:
    - step_E4_1_setup
    - step_E4_2_metric_and_frame
    """
    
    def __init__(self, 
                 mode: Mode,
                 ny_variant: NyVariant,
                 logger: Optional[ComputationLogger] = None,
                 checkpoint_mgr: Optional[CheckpointManager] = None):
        
        # Validate types
        if not isinstance(mode, Mode):
            raise TypeError(f"mode must be Mode enum, got {type(mode)}")
        if not isinstance(ny_variant, NyVariant):
            raise TypeError(f"ny_variant must be NyVariant enum, got {type(ny_variant)}")
        
        self.mode = mode
        self.ny_variant = ny_variant
        self.logger = logger
        self.ckpt = checkpoint_mgr
        self.data = {}
        
        self.steps = [
            ("E4.1", self.step_E4_1_setup),
            ("E4.2", self.step_E4_2_metric_and_frame),
            ("E4.3", self.step_E4_3_christoffel_frame),
            ("E4.4", self.step_E4_4_torsion_ansatz_frame),
            ("E4.5", self.step_E4_5_contortion_frame),
            ("E4.6", self.step_E4_6_ec_connection_frame),
            ("E4.7", self.step_E4_7_riemann_tensor_frame),
            ("E4.8", self.step_E4_8_ricci_scalar_frame),
            ("E4.9", self.step_E4_9_torsion_scalar_frame),
            ("E4.10", self.step_E4_10_nieh_yan_frame),
            ("E4.11", self.step_E4_11_lagrangian),
            ("E4.12", self.step_E4_12_angular_integration),
            ("E4.13", self.step_E4_13_effective_potential),
            ("E4.14", self.step_E4_14_stability_analysis),
            ("E4.15", self.step_E4_15_summary),
        ]
    
    def _log_step(self, step_id: str, description: str):
        """Log step start (optional)"""
        if self.logger:
            self.logger.step(step_id, description)
    
    def _log_info(self, message: str, indent: int = 2):
        """Log info message (optional)"""
        if self.logger:
            self.logger.info(message, indent)
    
    def _log_success(self, message: str):
        """Log success message (optional)"""
        if self.logger:
            self.logger.success(message)
    
    def _log_warning(self, message: str):
        """Log warning message (optional)"""
        if self.logger:
            self.logger.warning(message)
    
    def _log_error(self, message: str):
        """Log error message (optional)"""
        if self.logger:
            self.logger.error(message)
    
    def _save_checkpoint(self, step_id: str, metadata: Optional[Dict] = None):
        """Save checkpoint (optional)"""
        if self.ckpt:
            full_metadata = {
                'step_name': step_id,
                'mode': self.mode.value,
                'ny_variant': self.ny_variant.value
            }
            if metadata:
                full_metadata.update(metadata)
            self.ckpt.save(step_id, self.data, full_metadata)
    
    def run(self, start_step: str = "E4.1"):
        """Execute computation pipeline"""
        # Display configuration
        self._log_info("="*70)
        self._log_info("EXECUTION CONFIGURATION")
        self._log_info("="*70)
        self._log_info(f"Mode: {self.mode.value}")
        self._log_info(f"NY Variant: {self.ny_variant.value}")
        self._log_info("="*70)
        
        try:
            start_idx = next(i for i, (sid, _) in enumerate(self.steps) if sid == start_step)
        except StopIteration:
            raise ValueError(f"Invalid start step: {start_step}")
        
        for step_id, step_func in self.steps[start_idx:]:
            try:
                doc = step_func.__doc__
                desc = doc.strip() if doc else "Executing step..."
                self._log_step(step_id, desc)
                
                step_func()
                self._save_checkpoint(step_id)
            except Exception as e:
                self._log_error(f"Step {step_id} failed: {e}")
                self._log_error(traceback.format_exc())
                raise
        
        if self.logger:
            self.logger.finalize()
    
    # ========================================
    # Abstract Methods (Must be Implemented)
    # ========================================
    
    def step_E4_1_setup(self):
        """Setup parameters and symbolic variables"""
        raise NotImplementedError("This method must be implemented by the topology runner.")
    
    def step_E4_2_metric_and_frame(self):
        """Define Frame Metric, Volume, and Structure Constants"""
        raise NotImplementedError("This method must be implemented by the topology runner.")
    
    # ========================================
    # Step E4.3: Christoffel (CORRECTED Koszul Formula)
    # ========================================
    
    def step_E4_3_christoffel_frame(self):
        """Compute Frame Connection (Spin Connection) from Structure Constants via GENERAL Koszul Formula"""
        self._log_info("Computing Frame Connection Γ^a_bc via general Koszul formula...")
        
        dim = self.data['dim']
        C = self.data['structure_constants']
        
        Gamma_LC = MutableDenseNDimArray.zeros(dim, dim, dim)
        
        # GENERAL Koszul formula for orthonormal frame (CONVENTIONS 3.2 consistent):
        # Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} - C^b_{ac})
        # This works for general left-invariant frames, including Nil³×S¹
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    term1 = C[a, b, c]
                    term2 = C[c, b, a]  # Changed from C[c, a, b]
                    term3 = C[b, a, c]  # Changed from C[b, c, a]
                    
                    val = (term1 + term2 - term3) / S(2)
                    if val != S.Zero:
                        Gamma_LC[a, b, c] = simplify(val)
        
        self.data['connection_LC'] = Gamma_LC
        
        # Log non-zero components
        nonzero_count = 0
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    if Gamma_LC[a, b, c] != S.Zero:
                        nonzero_count += 1
                        if nonzero_count <= 10:  # Show first 10 components
                            self._log_info(f"  Gamma^{a}_{{{b}{c}}} = {Gamma_LC[a,b,c]}")
        if nonzero_count > 10:
            self._log_info(f"  ... and {nonzero_count - 10} more non-zero components")
        
        # Verify metric compatibility as a safety check
        self._verify_metric_compatibility()
        
        self._log_success("Frame Connection computed")
    
    def _verify_metric_compatibility(self):
        """Verify Γ_{abc} + Γ_{bac} = 0 (metric compatibility for orthonormal frame)"""
        dim = self.data['dim']
        Gamma = self.data['connection_LC']
        
        violations = []
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    # Γ_{abc} + Γ_{bac} should be 0
                    check = simplify(Gamma[a, b, c] + Gamma[b, a, c])
                    if check != S.Zero:
                        violations.append((a, b, c, check))
        
        if violations:
            self._log_error("Metric compatibility VIOLATED:")
            for a, b, c, val in violations[:5]:  # Show first 5 violations
                self._log_error(f"  Γ_{{{a}{b}{c}}} + Γ_{{{b}{a}{c}}} = {val}")
            raise ValueError(f"Metric compatibility check failed: {len(violations)} violations")
        else:
            self._log_info("  Metric compatibility (Γ_abc + Γ_bac = 0): VERIFIED")
    
    # ========================================
    # Step E4.4: Torsion Ansatz (Frame Basis)
    # ========================================
    
    def step_E4_4_torsion_ansatz_frame(self):
        """Construct Torsion Tensor in Frame Basis"""
        self._log_info(f"Constructing Torsion T_abc (Mode: {self.mode.value})...")
        
        dim = self.data['dim']
        r = self.data['params']['r']
        eta = self.data['params']['eta']
        V = self.data['params']['V']
        
        T_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)
        
        # T1: Totally Antisymmetric (Axial)
        if self.mode in [Mode.AX, Mode.MX]:
            if eta != S.Zero:
                self._log_info("  Adding T1 (Axial/Antisymmetric)...")
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            eps = epsilon_symbol(i, j, k)
                            if eps != 0:
                                T_tensor[i, j, k] = 2 * eta * eps / r
        
        # T2: Vector Trace
        if self.mode in [Mode.VT, Mode.MX]:
            if V != S.Zero:
                self._log_info("  Adding T2 (Vector Trace)...")
                T_vec = [S.Zero, S.Zero, S.Zero, V]
                metric = self.data['metric_frame']
                
                for a in range(dim):
                    for b in range(dim):
                        for c in range(dim):
                            val = (metric[a,c] * T_vec[b] - metric[a,b] * T_vec[c]) * Rational(1, 3)
                            if val != S.Zero:
                                T_tensor[a, b, c] = T_tensor[a, b, c] + val
        
        self.data['torsion_tensor'] = T_tensor
        self._log_success("Torsion constructed")
    
    # ========================================
    # Step E4.5: Contortion (Frame Basis)
    # ========================================
    
    def step_E4_5_contortion_frame(self):
        """Compute Contortion K_abc in Frame Basis"""
        self._log_info("Computing Contortion...")
        dim = self.data['dim']
        T = self.data['torsion_tensor']
        
        K = MutableDenseNDimArray.zeros(dim, dim, dim)
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    val = (T[a,b,c] + T[b,c,a] - T[c,a,b]) * Rational(1, 2)
                    if val != 0:
                        K[a,b,c] = simplify(val)
        
        self.data['contortion'] = K
        self._log_success("Contortion computed")
    
    # ========================================
    # Step E4.6: EC Connection (Frame Basis)
    # ========================================
    
    def step_E4_6_ec_connection_frame(self):
        """Compute EC Connection"""
        self._log_info("Computing EC Connection...")
        dim = self.data['dim']
        Gamma_LC = self.data['connection_LC']
        K = self.data['contortion']
        
        Gamma_EC = MutableDenseNDimArray.zeros(dim, dim, dim)
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    Gamma_EC[a,b,c] = Gamma_LC[a,b,c] + K[a,b,c]
        
        self.data['connection_EC'] = Gamma_EC
        self._log_success("EC Connection computed")
    
    # ========================================
    # Step E4.7: Riemann Tensor (Frame Basis)
    # ========================================
    
    def step_E4_7_riemann_tensor_frame(self):
        """Compute Riemann Tensor R^a_bcd in Frame Basis"""
        self._log_info("Computing Riemann Tensor...")
        dim = self.data['dim']
        Gamma = self.data['connection_EC']
        C = self.data['structure_constants']
        
        Riemann = MutableDenseNDimArray.zeros(dim, dim, dim, dim)
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        term = S.Zero
                        for e in range(dim):
                            term += Gamma[a,e,c] * Gamma[e,b,d]
                            term -= Gamma[a,e,d] * Gamma[e,b,c]
                        for e in range(dim):
                            term -= Gamma[a,b,e] * C[e,c,d]
                        
                        if term != 0:
                            Riemann[a,b,c,d] = simplify(term)
        
        self.data['riemann'] = Riemann
        
        # Mandatory: Strict Riemann antisymmetry check
        self._verify_riemann_antisymmetry_strict()
        
        self._log_success("Riemann Tensor computed and verified")
    
    def _verify_riemann_antisymmetry_strict(self):
        """Verify R_{abcd} antisymmetries with STRICT 3-level check"""
        Riemann = self.data['riemann']
        eta = self.data['metric_frame']
        dim = self.data['dim']
        
        # Lower first index: R_abcd = η_ae R^e_bcd
        R_abcd = MutableDenseNDimArray.zeros(dim, dim, dim, dim)
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        val = S.Zero
                        for e in range(dim):
                            val += eta[a, e] * Riemann[e, b, c, d]
                        R_abcd[a, b, c, d] = val
        
        # STRICT verification (raises exception if violated)
        # Note: symbols_list is now extracted per-residual inside verify_antisymmetry_strict
        verify_antisymmetry_strict(R_abcd, dim, None, self.logger)
        
        # Store for NY calculation
        self.data['riemann_abcd'] = R_abcd
    
    # ========================================
    # Remaining steps (E4.8-E4.15) - Similar refactoring
    # ========================================
    
    def step_E4_8_ricci_scalar_frame(self):
        """Compute Ricci Scalar R = R^a_baa"""
        self._log_info("Computing Ricci Scalar...")
        dim = self.data['dim']
        Riemann = self.data['riemann']
        
        Ricci = Matrix.zeros(dim, dim)
        for b in range(dim):
            for d in range(dim):
                val = S.Zero
                for a in range(dim):
                    val += Riemann[a,b,a,d]
                Ricci[b,d] = simplify(val)
        
        R_scalar = S.Zero
        for a in range(dim):
            R_scalar += Ricci[a,a]
        
        R_scalar = simplify(R_scalar)
        self.data['ricci_scalar'] = R_scalar
        
        self._log_info(f"Ricci Scalar R = {R_scalar}")
    
    def step_E4_9_torsion_scalar_frame(self):
        """Compute Torsion Scalar T = T_{abc} T^{abc}"""
        self._log_info("Computing Torsion Scalar...")
        dim = self.data['dim']
        T = self.data['torsion_tensor']
        
        val = S.Zero
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    val += T[a,b,c] * T[a,b,c]
        
        self.data['torsion_scalar'] = simplify(val)
        self._log_info(f"Torsion Scalar T = {self.data['torsion_scalar']}")
    
    def step_E4_10_nieh_yan_frame(self):
        """Compute NY density (all 3 variants, then select)"""
        self._log_info(f"Computing Nieh-Yan Density (variant: {self.ny_variant.value})...")
        dim = self.data['dim']
        T = self.data['torsion_tensor']
        R_abcd = self.data['riemann_abcd']
        
        # 1. TT-only term
        N_TT = S.Zero
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        eps = levi_civita_4d(a,b,c,d)
                        if eps != 0:
                            sum_term = S.Zero
                            for e in range(dim):
                                sum_term += T[e,a,b] * T[e,c,d]
                            N_TT += Rational(1,4) * eps * sum_term
        
        self.data['ny_density_TT'] = simplify(N_TT)
        self._log_info(f"  N_TT = {self.data['ny_density_TT']}")
        
        # 2. Ree term
        N_Ree = S.Zero
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        eps = levi_civita_4d(a,b,c,d)
                        if eps != 0:
                            N_Ree += Rational(1,4) * eps * R_abcd[a,b,c,d]
        
        self.data['ny_density_Ree'] = simplify(N_Ree)
        self._log_info(f"  N_Ree = {self.data['ny_density_Ree']}")
        
        # 3. Full NY
        self.data['ny_density_full'] = simplify(self.data['ny_density_TT'] - self.data['ny_density_Ree'])
        self._log_info(f"  N_full = N_TT - N_Ree = {self.data['ny_density_full']}")
        
        # 4. Select adopted NY
        if self.ny_variant == NyVariant.TT:
            self.data['nieh_yan_density'] = self.data['ny_density_TT']
            self._log_info(f">>> ADOPTED_NY: N_TT")
        elif self.ny_variant == NyVariant.REE:
            self.data['nieh_yan_density'] = self.data['ny_density_Ree']
            self._log_info(f">>> ADOPTED_NY: N_Ree")
        elif self.ny_variant == NyVariant.FULL:
            self.data['nieh_yan_density'] = self.data['ny_density_full']
            self._log_info(f">>> ADOPTED_NY: N_full")
        
        self._log_success(f"NY computation complete (variant: {self.ny_variant.value})")
    
    def step_E4_11_lagrangian(self):
        """Construct Lagrangian"""
        self._log_info("Constructing Lagrangian...")
        
        R = self.data['ricci_scalar']
        N = self.data['nieh_yan_density']
        k = self.data['params']['kappa']
        t = self.data['params']['theta_NY']
        
        L = R / (2*k**2) + t * N
        self.data['lagrangian'] = simplify(L)
        
        self._log_info(f"L = R/(2*kappa^2) + theta_NY*N")
    
    def step_E4_12_angular_integration(self):
        """Integrate Action"""
        self._log_info("Integrating Action...")
        Vol = self.data['total_volume']
        Lag = self.data['lagrangian']
        Action = simplify(Lag * Vol)
        self.data['action'] = Action
        self._log_info(f"Total Action = {Action}")
    
    def step_E4_13_effective_potential(self):
        """Extract Potential V(r)"""
        self._log_info("Extracting Potential...")
        Action = self.data['action']
        V = -Action
        self.data['potential'] = V
        self._log_info(f"Effective Potential V(r) = {V}")
    
    def step_E4_14_stability_analysis(self):
        """Analyze Stability"""
        self._log_info("Stability Analysis...")
        V = self.data['potential']
        r = self.data['params']['r']
        dV = simplify(diff(V, r))
        d2V = simplify(diff(dV, r))
        
        self._log_info(f"dV/dr = {dV}")
        self.data['equilibria'] = []
    
    def step_E4_15_summary(self):
        """Summary"""
        self._log_info("Computation complete. Check log for details.")
    
    def get_effective_potential_function(self):
        """Returns a fast, numpy-compatible function for V(r)"""
        if 'potential' not in self.data:
            raise RuntimeError("Potential not computed yet. Run engine.run() first.")
        
        r = self.data['params']['r']
        V = self.data['params']['V']
        eta = self.data['params']['eta']
        theta = self.data['params']['theta_NY']
        L = self.data['params']['L']
        kappa = self.data['params']['kappa']
        
        sym_expr = self.data['potential']
        
        fast_func = lambdify([r, V, eta, theta, L, kappa], sym_expr, modules='numpy')
        
        return fast_func

    def get_potential_decomposition(self):
        """
        Decompose effective potential by powers of r.
        
        Returns:
            dict: {
                'r^3': SymPy expression (coefficient of r³),
                'r^2': SymPy expression (coefficient of r²),
                'r^1': SymPy expression (coefficient of r),
                'r^0': SymPy expression (constant term),
                'r^-1': SymPy expression (coefficient of 1/r),
                'total': SymPy expression (full potential),
                'params': dict of parameter symbols
            }
        
        Note:
            - r^3 term: typically from Vector torsion (V²)
            - r^2 term: typically from Nieh-Yan coupling (V×η×θ)
            - r^1 term: typically from Curvature (η²)
            - Other terms may appear depending on topology
        """
        from sympy import Poly, symbols, expand, collect, S
        
        if 'potential' not in self.data:
            raise RuntimeError("Potential not computed yet. Run engine.run() first.")
        
        r = self.data['params']['r']
        potential = self.data['potential']
        
        # Expand and collect by powers of r
        expanded = expand(potential)
        collected = collect(expanded, r, evaluate=False)
        
        # Extract coefficients for each power of r
        decomposition = {
            'r^3': S.Zero,
            'r^2': S.Zero,
            'r^1': S.Zero,
            'r^0': S.Zero,
            'r^-1': S.Zero,
            'total': potential,
            'params': self.data['params'].copy()
        }
        
        for term, coeff in collected.items():
            # Determine the power of r in this term
            if term == r**3:
                decomposition['r^3'] = coeff
            elif term == r**2:
                decomposition['r^2'] = coeff
            elif term == r:
                decomposition['r^1'] = coeff
            elif term == 1 or term == S.One:
                decomposition['r^0'] = coeff
            elif term == 1/r or term == r**(-1):
                decomposition['r^-1'] = coeff
        
        return decomposition

    def get_potential_term_functions(self):
        """
        Get lambdified functions for each potential term.
        
        Returns:
            dict: {
                'r^3': callable(r, V, eta, theta, L, kappa) -> numpy array,
                'r^2': callable(...),
                'r^1': callable(...),
                'total': callable(...),
                'labels': dict mapping term keys to display labels
            }
        """
        decomp = self.get_potential_decomposition()
        
        r = self.data['params']['r']
        V = self.data['params']['V']
        eta = self.data['params']['eta']
        theta = self.data['params']['theta_NY']
        L = self.data['params']['L']
        kappa = self.data['params']['kappa']
        
        args = [r, V, eta, theta, L, kappa]
        
        result = {
            'labels': {
                'r^3': r'$V^2$ term ($r^3$)',
                'r^2': r'$V \cdot \eta \cdot \theta$ coupling ($r^2$)',
                'r^1': r'$\eta^2$ + const ($r$)',
                'r^0': r'Constant',
                'r^-1': r'$1/r$ term',
                'total': r'Total $V(r)$'
            }
        }
        
        # Create lambdified functions for each term
        for key in ['r^3', 'r^2', 'r^1', 'r^0', 'r^-1', 'total']:
            coeff = decomp[key]
            if key == 'total':
                # Total is already the full expression
                expr = coeff
            else:
                # Multiply coefficient by appropriate power of r
                if key == 'r^3':
                    expr = coeff * r**3
                elif key == 'r^2':
                    expr = coeff * r**2
                elif key == 'r^1':
                    expr = coeff * r
                elif key == 'r^0':
                    expr = coeff
                elif key == 'r^-1':
                    expr = coeff / r
            
            if expr == S.Zero:
                # Return a function that returns zeros
                result[key] = lambda r_arr, *args, e=expr: np.zeros_like(r_arr)
            else:
                result[key] = lambdify(args, expr, modules='numpy')
        
        return result
