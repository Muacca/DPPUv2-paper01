# =============================================================================
# Cell 3: Interactive Phase Diagram & Potential Viewer (v3 with T³ Anisotropy)
# =============================================================================
# Prerequisites:
#   - Cell 0: %matplotlib widget
#   - Cell 1: Basic imports (numpy, matplotlib, sys, os)
#   - Engine files (DPPUv2_engine_core_v3.py, DPPUv2_runner_*.py) in path
#
# Features:
#   - Topology selection (S³×S¹, T³×S¹, Nil³×S¹)
#   - Torsion mode selection (MX/VT/AX) with on-demand engine generation
#   - NY variant selection (FULL/TT/REE)
#   - T³×S¹ anisotropy parameters (α=R2/R1, β=R3/R1)
#   - Cached potential functions (computed once per configuration)
#   - Interactive phase diagram with click-to-select points
#   - Multi-point potential comparison (up to 3 points)
#   - Axis scale switching (linear/log/symlog)
#
# v3 Changes:
#   - Added α and β sliders for T³ anisotropy (conditionally visible)
#   - Extended cache key to include (α, β) for T³
#   - Modified get_potential_function to accept and apply anisotropy
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import ipywidgets as widgets
from ipywidgets import Layout, HBox, VBox, Output
from IPython.display import display, clear_output
from scipy.optimize import minimize_scalar
import sys
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Engine Imports
# =============================================================================
# Ensure engine files are in path
sys.path.append(os.getcwd())

try:
    from DPPUv2_engine_core_v3 import Mode, NyVariant
    from DPPUv2_runner_S3S1_v3 import S3S1Engine
    from DPPUv2_runner_Nil3S1_v3 import Nil3S1Engine
    from DPPUv2_runner_T3S1_v3 import T3S1Engine
    print("✓ Engine modules imported successfully.")
except ImportError as e:
    print(f"✗ Error importing engines: {e}")
    print("  Make sure DPPUv2_engine_core_v3.py and runner files are in the current directory.")
    raise

# =============================================================================
# Constants
# =============================================================================
PI = np.pi
KAPPA = 1.0
L_SCALE = 1.0

# Phase classification constants
R_MIN = 0.01
R_MAX = 1000000.0
R_BOUNDARY_THRESHOLD = 0.02

# =============================================================================
# Engine Registry
# =============================================================================
ENGINE_CLASSES = {
    'S3': S3S1Engine,
    'T3': T3S1Engine,
    'Nil3': Nil3S1Engine,
}

MODE_ENUMS = {
    'MX': Mode.MX,
    'VT': Mode.VT,
    'AX': Mode.AX,
}

NY_ENUMS = {
    'FULL': NyVariant.FULL,
    'TT': NyVariant.TT,
    'REE': NyVariant.REE,
}

# =============================================================================
# Potential Function Cache
# =============================================================================
_potential_cache = {}  # key: (topology, mode, ny_variant[, alpha, beta]) → value: callable

def get_potential_function(topology, mode, ny_variant, alpha_T3=1.0, beta_T3=1.0, status_callback=None):
    """
    Get potential function with caching.
    
    Args:
        topology: 'S3', 'T3', or 'Nil3'
        mode: 'MX', 'VT', or 'AX'
        ny_variant: 'FULL', 'TT', or 'REE'
        alpha_T3: R2/R1 ratio (T3 only, default=1.0 for isotropy)
        beta_T3: R3/R1 ratio (T3 only, default=1.0 for isotropy)
        status_callback: Optional function to display status messages
    
    Returns:
        Callable function f(r, V, eta, theta_NY) -> V_eff
    """
    from sympy import symbols, lambdify, S
    
    # Cache key includes anisotropy parameters for T3
    if topology == 'T3':
        key = (topology, mode, ny_variant, alpha_T3, beta_T3)
    else:
        key = (topology, mode, ny_variant)
    
    if key in _potential_cache:
        if status_callback:
            if topology == 'T3':
                status_callback(f"Using cached: {topology}×S¹ / {mode} / {ny_variant} / α={alpha_T3:.2f}, β={beta_T3:.2f}")
            else:
                status_callback(f"Using cached: {topology}×S¹ / {mode} / {ny_variant}")
        return _potential_cache[key]
    
    # Generate new function
    if status_callback:
        if topology == 'T3':
            status_callback(f"Generating: {topology}×S¹ / {mode} / {ny_variant} / α={alpha_T3:.2f}, β={beta_T3:.2f}... (this may take a few seconds)")
        else:
            status_callback(f"Generating: {topology}×S¹ / {mode} / {ny_variant}... (this may take a few seconds)")
    
    engine_class = ENGINE_CLASSES[topology]
    mode_enum = MODE_ENUMS[mode]
    ny_enum = NY_ENUMS[ny_variant]
    
    # Suppress engine output during generation
    original_stdout = sys.stdout
    try:
        sys.stdout = open(os.devnull, 'w')
        
        engine = engine_class(mode=mode_enum, ny_variant=ny_enum)
        engine.run()
        
        # === T3 Anisotropy: R2 = α*r, R3 = β*r ===
        if topology == 'T3':
            params = engine.data['params']
            if 'R2' in params and 'R3' in params:
                r_sym = params['r']  # This is R1
                R2_sym = params['R2']
                R3_sym = params['R3']
                engine.data['potential'] = engine.data['potential'].subs({
                    R2_sym: alpha_T3 * r_sym,
                    R3_sym: beta_T3 * r_sym
                })
        
        # === Custom lambdify to handle AX/VT modes ===
        # In AX mode: V = 0 (not a symbol)
        # In VT mode: eta = 0 (not a symbol)
        # We need to create consistent interface f(r, V, eta, theta_NY)
        
        params = engine.data['params']
        sym_expr = engine.data['potential']
        
        r_sym = params['r']
        L_sym = params['L']
        kappa_sym = params['kappa']
        theta_sym = params['theta_NY']
        V_sym = params['V']
        eta_sym = params['eta']
        
        # Check which parameters are actual symbols vs fixed values
        V_is_symbol = hasattr(V_sym, 'is_Symbol') and V_sym.is_Symbol
        eta_is_symbol = hasattr(eta_sym, 'is_Symbol') and eta_sym.is_Symbol
        
        # Build lambdify argument list with only actual symbols
        lambda_args = [r_sym]
        if V_is_symbol:
            lambda_args.append(V_sym)
        if eta_is_symbol:
            lambda_args.append(eta_sym)
        lambda_args.extend([theta_sym, L_sym, kappa_sym])
        
        # Create lambdified function
        raw_func = lambdify(lambda_args, sym_expr, modules='numpy')
        
    finally:
        sys.stdout.close()
        sys.stdout = original_stdout
    
    # Create wrapper with consistent interface: f(r, V_param, eta_param, theta_NY)
    if mode == 'MX':
        # Both V and eta are symbols
        def wrapped_func(r, V_param, eta_param, theta_NY):
            return raw_func(r, V_param, eta_param, theta_NY, L_SCALE, KAPPA)
    elif mode == 'AX':
        # V = 0 (fixed), only eta is symbol
        def wrapped_func(r, V_param, eta_param, theta_NY):
            # V_param is ignored (always 0 in this mode)
            return raw_func(r, eta_param, theta_NY, L_SCALE, KAPPA)
    elif mode == 'VT':
        # eta = 0 (fixed), only V is symbol
        def wrapped_func(r, V_param, eta_param, theta_NY):
            # eta_param is ignored (always 0 in this mode)
            return raw_func(r, V_param, theta_NY, L_SCALE, KAPPA)
    
    # Cache it
    _potential_cache[key] = wrapped_func
    
    if status_callback:
        if topology == 'T3':
            status_callback(f"✓ Generated and cached: {topology}×S¹ / {mode} / {ny_variant} / α={alpha_T3:.2f}, β={beta_T3:.2f}")
        else:
            status_callback(f"✓ Generated and cached: {topology}×S¹ / {mode} / {ny_variant}")
    
    return wrapped_func


def get_cache_status():
    """Return current cache contents as a formatted string."""
    if not _potential_cache:
        return "Cache is empty."
    
    lines = ["Cached configurations:"]
    for key in sorted(_potential_cache.keys()):
        if len(key) == 5:  # T3 with anisotropy
            lines.append(f"  • {key[0]}×S¹ / {key[1]} / {key[2]} / α={key[3]:.2f}, β={key[4]:.2f}")
        else:  # S3 or Nil3 without anisotropy
            lines.append(f"  • {key[0]}×S¹ / {key[1]} / {key[2]}")
    return "\n".join(lines)


def clear_cache():
    """Clear the potential function cache."""
    global _potential_cache
    _potential_cache = {}
    print("Cache cleared.")

# =============================================================================
# Phase Classification (from parameter_scan_v3)
# =============================================================================

def classify_phase(potential_func, V_param, eta, theta):
    """
    Classify the stability type of the potential.
    
    Returns:
        1: Type I  - Local minimum with barrier (V increases toward r=0)
        2: Type II - Rolling (V decreases from r=0 toward minimum)
        3: Type III - No local minimum in physical region
    """
    def v(r):
        return potential_func(r, V_param, eta, theta)
    
    # Find minimum
    try:
        res_min = minimize_scalar(v, bounds=(R_MIN, R_MAX), method='bounded')
    except Exception:
        return 3
    
    if not res_min.success:
        return 3
    
    r0 = res_min.x
    
    # Check if hitting bounds
    if r0 < R_BOUNDARY_THRESHOLD or r0 > R_MAX - R_BOUNDARY_THRESHOLD:
        return 3
    
    # Check curvature (second derivative) - must be positive for minimum
    h = 1e-5
    try:
        d2v = (v(r0 + h) - 2 * v(r0) + v(r0 - h)) / h**2
    except Exception:
        return 3
    
    if d2v <= 0:
        return 3
    
    # Check slope near r=0
    r_test = R_MIN * 2
    dr = R_MIN * 0.1
    
    try:
        v1 = v(r_test)
        v2 = v(r_test + dr)
        slope_near_zero = (v2 - v1) / dr
    except Exception:
        slope_near_zero = 0
    
    if slope_near_zero > 0:
        return 1  # Type I: barrier at origin
    else:
        return 2  # Type II: rolling

# =============================================================================
# Interactive Viewer Class
# =============================================================================

class DPPUv2InteractiveViewer:
    """
    Interactive phase diagram and potential viewer for DPPUv2.
    Engine-direct version with on-demand caching.
    v3: Added T³ anisotropy parameters (α, β).
    
    Usage in notebook:
        %run DPPUv2_interactive_viewer_v3.py
        
        # Optional: Customize slider ranges before instantiation
        DPPUv2InteractiveViewer.SLIDER_V_MAX_LIMIT = 30.0
        DPPUv2InteractiveViewer.SLIDER_ETA_MIN_LIMIT = -30.0
        
        viewer = DPPUv2InteractiveViewer()
        viewer.display()
    """
    
    # =========================================================================
    # Class-level Configuration (can be modified before instantiation)
    # =========================================================================
    
    # Phase Diagram Range - Slider Limits
    SLIDER_V_MAX_MIN = 1.0
    SLIDER_V_MAX_MAX = 50.0
    SLIDER_V_MAX_DEFAULT = 10.0
    
    SLIDER_ETA_MIN_MIN = -50.0
    SLIDER_ETA_MIN_MAX = 0.0
    SLIDER_ETA_MIN_DEFAULT = -15.0
    
    SLIDER_ETA_MAX_MIN = 0.0
    SLIDER_ETA_MAX_MAX = 50.0
    SLIDER_ETA_MAX_DEFAULT = 5.0
    
    # Potential Plot Range - Slider Limits
    SLIDER_R_MAX_MIN = 1  # log scale (10^1)
    SLIDER_R_MAX_MAX = 6  # log scale (10^6)
    SLIDER_R_MAX_DEFAULT = 100.0
    
    SLIDER_VEFF_MIN_MIN = 0  # log scale (10^0)
    SLIDER_VEFF_MIN_MAX = 8  # log scale (10^8)
    SLIDER_VEFF_MIN_DEFAULT = 1e4
    
    SLIDER_VEFF_MAX_MIN = 2  # log scale (10^2)
    SLIDER_VEFF_MAX_MAX = 6  # log scale (10^6)
    SLIDER_VEFF_MAX_DEFAULT = 1e3
    
    def __init__(self):
        # State
        self.topology = 'S3'
        self.torsion_mode = 'MX'
        self.ny_variant = 'FULL'
        self.theta_NY = 1.0
        
        # T3 anisotropy parameters (R2 = α*r, R3 = β*r)
        self.alpha_T3 = 1.0  # Default: isotropic
        self.beta_T3 = 1.0   # Default: isotropic
        
        # Phase diagram range
        self.V_min = 0.0
        self.V_max = 10.0
        self.eta_min = -10.0
        self.eta_max = 10.0
        
        # Potential plot range
        self.r_max = 20.0
        self.Veff_min = -1e6
        self.Veff_max = 1e6
        self.r_log_scale = False
        self.Veff_log_scale = False
        
        # Selected points (up to 3)
        self.selected_points = [
            {'V': 1.0, 'eta': -4.0, 'active': True},
            {'V': 2.0, 'eta': -4.0, 'active': False},
            {'V': 2.0, 'eta': -0.0, 'active': False},
        ]
        self.current_point_idx = 0
        
        # Cached phase data
        self.phase_data = None
        self.V_grid = None
        self.eta_grid = None
        
        # Current potential function
        self._current_func = None
        self._current_config = None
        
        # Figure and axes
        self.fig = None
        self.ax_phase = None
        self.ax_potential = None
        self.click_cid = None
        
        # Create widgets
        self._create_widgets()
    
    def _create_widgets(self):
        """Create all ipywidgets controls."""
        style = {'description_width': '120px'}
        layout_slider = Layout(width='350px')
        layout_dropdown = Layout(width='200px')
        
        # --- Configuration Section ---
        self.w_topology = widgets.Dropdown(
            options=[('S³×S¹', 'S3'), ('T³×S¹', 'T3'), ('Nil³×S¹', 'Nil3')],
            value='S3',
            description='Topology:',
            style=style,
            layout=layout_dropdown
        )
        
        self.w_torsion_mode = widgets.Dropdown(
            options=[('Mixed (MX)', 'MX'), ('Vector-Trace (VT)', 'VT'), ('Axial (AX)', 'AX')],
            value='MX',
            description='Torsion Mode:',
            style=style,
            layout=layout_dropdown
        )
        
        self.w_ny_variant = widgets.Dropdown(
            options=[('FULL (TT−Ree)', 'FULL'), ('TT only', 'TT'), ('Ree only', 'REE')],
            value='FULL',
            description='NY Variant:',
            style=style,
            layout=layout_dropdown
        )
        
        self.w_theta_NY = widgets.FloatSlider(
            value=1.0, min=0.0, max=10.0, step=0.1,
            description='θ_NY:',
            style=style,
            layout=layout_slider,
            readout_format='.1f'
        )
        
        # --- T3 Anisotropy Parameters (conditionally visible) ---
        self.w_alpha_T3 = widgets.FloatSlider(
            value=1.0, min=0.001, max=1.0, step=0.001,
            description='α (R2/R1):',
            style=style,
            layout=layout_slider,
            readout_format='.2f'
        )
        
        self.w_beta_T3 = widgets.FloatSlider(
            value=1.0, min=0.001, max=1.0, step=0.001,
            description='β (R3/R1):',
            style=style,
            layout=layout_slider,
            readout_format='.2f'
        )
        
        # Group anisotropy controls into a VBox for easier show/hide
        self.w_anisotropy_box = VBox([
            widgets.HTML('<h5 style="margin-top:10px;">T³ Anisotropy (R2=α·R1, R3=β·R1)</h5>'),
            self.w_alpha_T3,
            self.w_beta_T3,
        ])
        
        # Initially hidden (only for T3)
        self.w_anisotropy_box.layout.display = 'none'
        
        # --- Phase Diagram Range ---
        self.w_V_max = widgets.FloatSlider(
            value=self.SLIDER_V_MAX_DEFAULT,
            min=self.SLIDER_V_MAX_MIN,
            max=self.SLIDER_V_MAX_MAX,
            step=1.0,
            description='V max:',
            style=style,
            layout=layout_slider
        )
        
        self.w_eta_min = widgets.FloatSlider(
            value=self.SLIDER_ETA_MIN_DEFAULT,
            min=self.SLIDER_ETA_MIN_MIN,
            max=self.SLIDER_ETA_MIN_MAX,
            step=1.0,
            description='η min:',
            style=style,
            layout=layout_slider
        )
        
        self.w_eta_max = widgets.FloatSlider(
            value=self.SLIDER_ETA_MAX_DEFAULT,
            min=self.SLIDER_ETA_MAX_MIN,
            max=self.SLIDER_ETA_MAX_MAX,
            step=1.0,
            description='η max:',
            style=style,
            layout=layout_slider
        )
        
        # --- Potential Plot Range ---
        self.w_r_max = widgets.FloatLogSlider(
            value=self.SLIDER_R_MAX_DEFAULT,
            base=10,
            min=self.SLIDER_R_MAX_MIN,
            max=self.SLIDER_R_MAX_MAX,
            step=0.1,
            description='r max:',
            style=style,
            layout=layout_slider,
            readout_format='.0e'
        )
        
        self.w_Veff_min = widgets.FloatLogSlider(
            value=self.SLIDER_VEFF_MIN_DEFAULT,
            base=10,
            min=self.SLIDER_VEFF_MIN_MIN,
            max=self.SLIDER_VEFF_MIN_MAX,
            step=0.1,
            description='|V_eff| min:',
            style=style,
            layout=layout_slider,
            readout_format='.0e'
        )
        
        self.w_Veff_max = widgets.FloatLogSlider(
            value=self.SLIDER_VEFF_MAX_DEFAULT,
            base=10,
            min=self.SLIDER_VEFF_MAX_MIN,
            max=self.SLIDER_VEFF_MAX_MAX,
            step=0.1,
            description='V_eff max:',
            style=style,
            layout=layout_slider,
            readout_format='.0e'
        )
        
        self.w_r_scale = widgets.ToggleButtons(
            options=['Linear', 'Log'],
            value='Linear',
            description='r axis:',
            style=style
        )
        
        self.w_Veff_scale = widgets.ToggleButtons(
            options=['Linear', 'Symlog'],
            value='Linear',
            description='V_eff axis:',
            style=style
        )
        
        # --- Point Selection ---
        self.w_point_select = widgets.Dropdown(
            options=[('Point 1 (Black)', 0), ('Point 2 (Green)', 1), ('Point 3 (Red)', 2)],
            value=0,
            description='Active Point:',
            style=style,
            layout=layout_dropdown
        )
        
        self.w_point_V = widgets.FloatSlider(
            value=3.0, min=0.0, max=20.0, step=0.1,
            description='V:',
            style=style,
            layout=layout_slider
        )
        
        self.w_point_eta = widgets.FloatSlider(
            value=-5.0, min=-20.0, max=20.0, step=0.1,
            description='η:',
            style=style,
            layout=layout_slider
        )
        
        self.w_point_active = widgets.Checkbox(
            value=True,
            description='Show this point',
            style=style
        )
        
        # --- Status Display ---
        self.w_status = widgets.HTML(
            value='<i style="color:gray;">Ready. Click "Draw Phase Diagram" to start.</i>',
            layout=Layout(width='100%')
        )
        
        self.w_cache_status = widgets.HTML(
            value=f'<pre style="font-size:11px; color:gray;">{get_cache_status()}</pre>',
            layout=Layout(width='100%')
        )
        
        # --- Buttons ---
        self.w_draw_button = widgets.Button(
            description='Draw Phase Diagram',
            button_style='primary',
            layout=Layout(width='180px', height='40px')
        )
        self.w_draw_button.on_click(self._on_draw_clicked)
        
        self.w_update_potential_button = widgets.Button(
            description='Update Potential',
            button_style='success',
            layout=Layout(width='180px', height='40px')
        )
        self.w_update_potential_button.on_click(self._on_update_potential_clicked)
        
        self.w_clear_cache_button = widgets.Button(
            description='Clear Cache',
            button_style='warning',
            layout=Layout(width='120px', height='30px')
        )
        self.w_clear_cache_button.on_click(self._on_clear_cache_clicked)
        
        # --- Output Area ---
        self.output = Output()
        
        # --- Event Handlers ---
        self.w_topology.observe(self._on_topology_changed, names='value')
        self.w_point_select.observe(self._on_point_select_changed, names='value')
        self.w_point_V.observe(self._on_point_slider_changed, names='value')
        self.w_point_eta.observe(self._on_point_slider_changed, names='value')
        self.w_point_active.observe(self._on_point_active_changed, names='value')
    
    def _set_status(self, msg, color='black'):
        """Update status display."""
        self.w_status.value = f'<span style="color:{color};">{msg}</span>'
    
    def _update_cache_display(self):
        """Update cache status display."""
        self.w_cache_status.value = f'<pre style="font-size:11px; color:gray;">{get_cache_status()}</pre>'
    
    def _on_topology_changed(self, change):
        """Show/hide T3-specific anisotropy sliders."""
        if change['new'] == 'T3':
            self.w_anisotropy_box.layout.display = 'flex'
        else:
            self.w_anisotropy_box.layout.display = 'none'
    
    def _on_point_select_changed(self, change):
        """Handle point selection dropdown change."""
        idx = change['new']
        self.current_point_idx = idx
        pt = self.selected_points[idx]
        
        # Update sliders without triggering callbacks
        self.w_point_V.unobserve(self._on_point_slider_changed, names='value')
        self.w_point_eta.unobserve(self._on_point_slider_changed, names='value')
        self.w_point_active.unobserve(self._on_point_active_changed, names='value')
        
        self.w_point_V.value = pt['V']
        self.w_point_eta.value = pt['eta']
        self.w_point_active.value = pt['active']
        
        self.w_point_V.observe(self._on_point_slider_changed, names='value')
        self.w_point_eta.observe(self._on_point_slider_changed, names='value')
        self.w_point_active.observe(self._on_point_active_changed, names='value')
    
    def _on_point_slider_changed(self, change):
        """Handle point V/eta slider change."""
        idx = self.current_point_idx
        self.selected_points[idx]['V'] = self.w_point_V.value
        self.selected_points[idx]['eta'] = self.w_point_eta.value
        self._update_point_markers()
    
    def _on_point_active_changed(self, change):
        """Handle point active checkbox change."""
        idx = self.current_point_idx
        self.selected_points[idx]['active'] = change['new']
        self._update_point_markers()
        self._update_potential_plot()
    
    def _on_draw_clicked(self, b):
        """Handle Draw Phase Diagram button click."""
        self._read_widget_values()
        self._compute_phase_diagram()
        self._draw_all()
        self._update_cache_display()
    
    def _on_update_potential_clicked(self, b):
        """Handle Update Potential button click."""
        self._read_widget_values()
        self._update_potential_plot()
    
    def _on_clear_cache_clicked(self, b):
        """Handle Clear Cache button click."""
        clear_cache()
        self._current_func = None
        self._current_config = None
        self._update_cache_display()
        self._set_status("Cache cleared.", color='orange')
    
    def _read_widget_values(self):
        """Read current widget values into state."""
        self.topology = self.w_topology.value
        self.torsion_mode = self.w_torsion_mode.value
        self.ny_variant = self.w_ny_variant.value
        self.theta_NY = self.w_theta_NY.value
        
        # T3 anisotropy
        if self.topology == 'T3':
            self.alpha_T3 = self.w_alpha_T3.value
            self.beta_T3 = self.w_beta_T3.value
        
        self.V_min = 0.0
        self.V_max = self.w_V_max.value
        self.eta_min = self.w_eta_min.value
        self.eta_max = self.w_eta_max.value
        
        self.r_max = self.w_r_max.value
        self.Veff_min = -self.w_Veff_min.value
        self.Veff_max = self.w_Veff_max.value
        self.r_log_scale = (self.w_r_scale.value == 'Log')
        self.Veff_log_scale = (self.w_Veff_scale.value == 'Symlog')
    
    def _get_potential_func(self):
        """Get the appropriate potential function (with caching)."""
        # Config includes anisotropy for T3
        if self.topology == 'T3':
            config = (self.topology, self.torsion_mode, self.ny_variant, self.alpha_T3, self.beta_T3)
        else:
            config = (self.topology, self.torsion_mode, self.ny_variant)
        
        if self._current_config != config:
            if self.topology == 'T3':
                self._current_func = get_potential_function(
                    self.topology, 
                    self.torsion_mode, 
                    self.ny_variant,
                    alpha_T3=self.alpha_T3,
                    beta_T3=self.beta_T3,
                    status_callback=self._set_status
                )
            else:
                self._current_func = get_potential_function(
                    self.topology, 
                    self.torsion_mode, 
                    self.ny_variant,
                    status_callback=self._set_status
                )
            self._current_config = config
        
        return self._current_func
    
    def _compute_phase_diagram(self):
        """Compute phase classification grid."""
        # Resolution: 100 x 100
        n_points = 100
        
        V_range = np.linspace(self.V_min, self.V_max, n_points)
        eta_range = np.linspace(self.eta_min, self.eta_max, n_points)
        
        self.V_grid, self.eta_grid = np.meshgrid(V_range, eta_range)
        self.phase_data = np.zeros_like(self.V_grid)
        
        func = self._get_potential_func()
        
        # Status message includes anisotropy for T3
        if self.topology == 'T3':
            status_str = f"{self.topology}×S¹ / {self.torsion_mode} / {self.ny_variant} / α={self.alpha_T3:.2f}, β={self.beta_T3:.2f}"
        else:
            status_str = f"{self.topology}×S¹ / {self.torsion_mode} / {self.ny_variant}"
        
        self._set_status(f"Computing phase diagram for {status_str}...", color='blue')
        
        total = n_points * n_points
        for i in range(n_points):
            for j in range(n_points):
                V_val = self.V_grid[i, j]
                eta_val = self.eta_grid[i, j]
                self.phase_data[i, j] = classify_phase(func, V_val, eta_val, self.theta_NY)
            
            # Progress indicator (update every 10%)
            if (i + 1) % 10 == 0:
                pct = 100 * (i + 1) / n_points
                self._set_status(f"Computing phase diagram... {pct:.0f}%", color='blue')
        
        self._set_status(f"✓ Phase diagram computed: {status_str} / θ_NY={self.theta_NY:.1f}", color='green')
    
    def _draw_all(self):
        """Draw both phase diagram and potential plot."""
        with self.output:
            clear_output(wait=True)
            
            # Create figure
            if self.fig is not None:
                plt.close(self.fig)
            
            self.fig, (self.ax_phase, self.ax_potential) = plt.subplots(
                1, 2, figsize=(14, 5), constrained_layout=True
            )
            
            self._draw_phase_diagram()
            self._draw_potential_plot()
            
            # Connect click event
            if self.click_cid is not None:
                self.fig.canvas.mpl_disconnect(self.click_cid)
            self.click_cid = self.fig.canvas.mpl_connect('button_press_event', self._on_click)
            
            plt.show()
    
    def _draw_phase_diagram(self):
        """Draw the phase diagram."""
        ax = self.ax_phase
        ax.clear()
        
        if self.phase_data is None:
            ax.text(0.5, 0.5, 'Click "Draw Phase Diagram" to compute', 
                    ha='center', va='center', transform=ax.transAxes)
            return
        
        # Color definitions
        color_type_I = '#90EE90'    # Light green (Stable Well)
        color_type_II = '#FFFF99'   # Light yellow (Rolling)
        color_type_III = '#E8E8E8'  # Light gray (Unstable) - distinguishable from white
        color_error = '#FFFFFF'     # White (error/undefined)
        
        # Colormap indexed by classify_phase return values:
        # 0 = error, 1 = Type I, 2 = Type II, 3 = Type III
        colors = [color_error, color_type_I, color_type_II, color_type_III]
        cmap = mcolors.ListedColormap(colors)
        
        # Bounds: [-0.5, 0.5, 1.5, 2.5, 3.5]
        bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
        norm = mcolors.BoundaryNorm(bounds, cmap.N)
        
        # Plot phase map
        ax.pcolormesh(self.V_grid, self.eta_grid, self.phase_data, 
                      cmap=cmap, norm=norm, shading='auto')
        
        # Contour lines at phase boundaries
        ax.contour(self.V_grid, self.eta_grid, self.phase_data, 
                   levels=[1.5, 2.5], colors=['darkgreen', 'darkorange'], linewidths=1.5)
        
        # Plot selected points
        self._draw_point_markers(ax)
        
        # Labels and title
        ax.set_xlabel(r'$V$ (Vector Torsion)', fontsize=11)
        ax.set_ylabel(r'$\eta$ (Axial Torsion)', fontsize=11)
        
        # Title includes anisotropy for T3
        if self.topology == 'T3':
            title_str = f'{self.topology}×S¹ Phase Diagram (α={self.alpha_T3:.2f}, β={self.beta_T3:.2f})\n'
        else:
            title_str = f'{self.topology}×S¹ Phase Diagram\n'
        title_str += f'({self.torsion_mode} / {self.ny_variant}, θ_NY={self.theta_NY:.1f})'
        ax.set_title(title_str, fontsize=12)
        
        ax.set_xlim(self.V_min, self.V_max)
        ax.set_ylim(self.eta_min, self.eta_max)
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Legend
        legend_elements = [
            Patch(facecolor=color_type_I, edgecolor='darkgreen', label='Type I (Stable Well)'),
            Patch(facecolor=color_type_II, edgecolor='darkorange', label='Type II (Rolling)'),
            Patch(facecolor=color_type_III, edgecolor='gray', label='Type III (Unstable)')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    def _draw_point_markers(self, ax):
        """Draw markers for selected points on phase diagram."""
        markers = ['o', 's', '^']  # Circle, Square, Triangle
        colors = ['black', 'green', 'red']
        marker_sizes = [100, 80, 80]
        
        for i, (pt, marker, color, ms) in enumerate(zip(self.selected_points, markers, colors, marker_sizes)):
            if pt['active']:
                ax.scatter(pt['V'], pt['eta'], marker=marker, s=ms, 
                           facecolors='none', edgecolors=color, linewidths=2,
                           zorder=10, label=f'Point {i+1}')
    
    def _update_point_markers(self):
        """Update point markers without recomputing phase diagram."""
        if self.ax_phase is not None and self.phase_data is not None:
            self._draw_phase_diagram()
            self.fig.canvas.draw_idle()
    
    def _draw_potential_plot(self):
        """Draw the potential plot for selected points."""
        ax = self.ax_potential
        ax.clear()
        
        func = self._get_potential_func()
        
        # r range
        if self.r_log_scale:
            r_arr = np.logspace(np.log10(0.01), np.log10(self.r_max), 1000)
        else:
            r_arr = np.linspace(0.01, self.r_max, 1000)
        
        colors = ['black', 'green', 'red']
        linestyles = ['-', '--', ':']
        labels = ['Point 1', 'Point 2', 'Point 3']
        
        has_plot = False
        for i, pt in enumerate(self.selected_points):
            if pt['active']:
                V_arr = np.array([func(r, pt['V'], pt['eta'], self.theta_NY) for r in r_arr])
                phase_type = classify_phase(func, pt['V'], pt['eta'], self.theta_NY)
                type_str = ['Error', 'Type I', 'Type II', 'Type III'][phase_type]
                
                ax.plot(r_arr, V_arr, color=colors[i], linestyle=linestyles[i], 
                        linewidth=2, label=f'{labels[i]}: V={pt["V"]:.1f}, η={pt["eta"]:.1f} ({type_str})')
                has_plot = True
        
        if not has_plot:
            ax.text(0.5, 0.5, 'No points selected', ha='center', va='center', 
                    transform=ax.transAxes)
            return
        
        # Zero line
        ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        
        # Axis scales
        if self.r_log_scale:
            ax.set_xscale('log')
        
        if self.Veff_log_scale:
            ax.set_yscale('symlog', linthresh=100)
        
        # Axis limits
        ax.set_xlim(0.01 if self.r_log_scale else 0, self.r_max)
        ax.set_ylim(self.Veff_min, self.Veff_max)
        
        # Labels and title
        ax.set_xlabel(r'$r$ (Radius)', fontsize=11)
        ax.set_ylabel(r'$V_{\mathrm{eff}}(r)$', fontsize=11)
        
        # Title includes anisotropy for T3
        if self.topology == 'T3':
            title_str = f'Effective Potential (θ_NY={self.theta_NY:.1f})\n(α={self.alpha_T3:.2f}, β={self.beta_T3:.2f})'
        else:
            title_str = f'Effective Potential\n(θ_NY={self.theta_NY:.1f})'
        ax.set_title(title_str, fontsize=12)
        
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3, linestyle='--')
    
    def _update_potential_plot(self):
        """Update only the potential plot."""
        if self.ax_potential is not None:
            self._read_widget_values()
            self._draw_potential_plot()
            self.fig.canvas.draw_idle()
    
    def _on_click(self, event):
        """Handle click on phase diagram."""
        if event.inaxes != self.ax_phase:
            return
        
        if event.button == 1:  # Left click
            idx = self.current_point_idx
            V_clicked = event.xdata
            eta_clicked = event.ydata
            
            if V_clicked is not None and eta_clicked is not None:
                # Clamp to range
                V_clicked = np.clip(V_clicked, self.V_min, self.V_max)
                eta_clicked = np.clip(eta_clicked, self.eta_min, self.eta_max)
                
                self.selected_points[idx]['V'] = V_clicked
                self.selected_points[idx]['eta'] = eta_clicked
                self.selected_points[idx]['active'] = True
                
                # Update sliders
                self.w_point_V.unobserve(self._on_point_slider_changed, names='value')
                self.w_point_eta.unobserve(self._on_point_slider_changed, names='value')
                self.w_point_active.unobserve(self._on_point_active_changed, names='value')
                
                self.w_point_V.value = V_clicked
                self.w_point_eta.value = eta_clicked
                self.w_point_active.value = True
                
                self.w_point_V.observe(self._on_point_slider_changed, names='value')
                self.w_point_eta.observe(self._on_point_slider_changed, names='value')
                self.w_point_active.observe(self._on_point_active_changed, names='value')
                
                # Redraw
                self._draw_phase_diagram()
                self._draw_potential_plot()
                self.fig.canvas.draw_idle()
    
    def display(self):
        """Display the interactive viewer."""
        # Clear any existing output (prevents multiple figures when re-running)
        from IPython.display import clear_output
        clear_output(wait=False)
        
        # =====================================================================
        # Top Section: Configuration & Status
        # =====================================================================
        config_box = VBox([
            widgets.HTML('<h4>Configuration</h4>'),
            HBox([self.w_topology, self.w_torsion_mode, self.w_ny_variant]),
            self.w_theta_NY,
            self.w_anisotropy_box,
        ])
        
        status_box = VBox([
            widgets.HTML('<h4>Status</h4>'),
            self.w_status,
            HBox([self.w_clear_cache_button]),
            self.w_cache_status,
        ])
        
        top_section = HBox([config_box, status_box])
        
        # =====================================================================
        # Bottom Section: Phase Diagram Controls (Left) & Potential Controls (Right)
        # =====================================================================
        
        # --- Left: Phase Diagram Controls ---
        phase_diagram_controls = VBox([
            widgets.HTML('<h4 style="border-bottom: 2px solid #4CAF50; padding-bottom: 5px;">Phase Diagram Controls</h4>'),
            widgets.HTML('<h5>Range</h5>'),
            self.w_V_max,
            HBox([self.w_eta_min, self.w_eta_max]),
            widgets.HTML('<h5 style="margin-top:10px;">Point Selection</h5>'),
            widgets.HTML('<p style="font-size:11px; color:gray;">Click on phase diagram or use sliders</p>'),
            self.w_point_select,
            HBox([self.w_point_V, self.w_point_eta]),
            self.w_point_active,
            widgets.HTML('<div style="margin-top:15px;"></div>'),
            self.w_draw_button,
        ])
        
        # --- Right: Potential Plot Controls ---
        potential_plot_controls = VBox([
            widgets.HTML('<h4 style="border-bottom: 2px solid #2196F3; padding-bottom: 5px;">Potential Plot Controls</h4>'),
            widgets.HTML('<h5>Range</h5>'),
            self.w_r_max,
            HBox([self.w_Veff_min, self.w_Veff_max]),
            widgets.HTML('<h5 style="margin-top:10px;">Axis Scale</h5>'),
            HBox([self.w_r_scale, self.w_Veff_scale]),
            widgets.HTML('<div style="margin-top:15px;"></div>'),
            self.w_update_potential_button,
        ])
        
        bottom_section = HBox([phase_diagram_controls, potential_plot_controls])
        
        # =====================================================================
        # Combine all sections
        # =====================================================================
        controls = VBox([top_section, bottom_section])
        
        display(controls)
        display(self.output)
        
        # Initial empty plot
        with self.output:
            self.fig, (self.ax_phase, self.ax_potential) = plt.subplots(
                1, 2, figsize=(14, 5), constrained_layout=True
            )
            self.ax_phase.text(0.5, 0.5, 'Click "Draw Phase Diagram" to start', 
                               ha='center', va='center', transform=self.ax_phase.transAxes,
                               fontsize=12, color='gray')
            self.ax_phase.set_title('Phase Diagram')
            self.ax_potential.text(0.5, 0.5, 'Potential will appear here', 
                                   ha='center', va='center', transform=self.ax_potential.transAxes,
                                   fontsize=12, color='gray')
            self.ax_potential.set_title('Effective Potential')
            
            self.click_cid = self.fig.canvas.mpl_connect('button_press_event', self._on_click)
            plt.show()

# =============================================================================
# Usage Example (for standalone execution)
# =============================================================================
# To use in a Jupyter notebook:
#   %run DPPUv2_interactive_viewer_v3.py
#   viewer = DPPUv2InteractiveViewer()
#   viewer.display()
#
# To customize slider ranges:
#   DPPUv2InteractiveViewer.SLIDER_V_MAX_MAX = 30.0
#   DPPUv2InteractiveViewer.SLIDER_ETA_MIN_MIN = -30.0
#   viewer = DPPUv2InteractiveViewer()
#   viewer.display()
