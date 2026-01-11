import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
import sys
import re
import os
from matplotlib.lines import Line2D

def extract_topology_variant(df):
    """
    Extract topology and NY variant from CSV data
    Assumes all rows have the same topology and variant
    """
    if 'topology' in df.columns and 'ny_variant' in df.columns:
        # Get from CSV columns (assumes uniform data)
        topology = df['topology'].iloc[0] if len(df) > 0 else None
        variant = df['ny_variant'].iloc[0] if len(df) > 0 else None
    else:
        # Fallback: try to extract from file name if columns don't exist
        topology = None
        variant = None
    
    return topology, variant

def create_phase_diagrams(csv_file):
    # 1. Load data
    df = pd.read_csv(csv_file)
    
    # Extract topology and variant from CSV data
    topology, variant = extract_topology_variant(df)
    
    # Get unique theta values ​​and sort
    thetas = np.sort(df['theta'].unique())
    
    for idx, theta_val in enumerate(thetas):
        output_filename = f"phase_diagram_{idx:04d}.png"
        plot_refined_phase_diagram(df, theta_val, output_filename, topology, variant)
        print(f"Saved: {output_filename} (theta={theta_val})")

def plot_refined_phase_diagram(df, theta_val, output_filename, topology=None, variant=None):
    # Extract data with a specific theta
    df_t = df[df['theta'] == theta_val].copy()
    
    # Create a grid
    v_unique = np.sort(df_t['V'].unique())
    eta_unique = np.sort(df_t['eta'].unique())
    
    if len(v_unique) < 2 or len(eta_unique) < 2:
        return # If data is missing
    
    V_grid, Eta_grid = np.meshgrid(v_unique, eta_unique)
    
    # Pivot data
    r0_grid = df_t.pivot(index='eta', columns='V', values='r0').values
    dv_grid = df_t.pivot(index='eta', columns='V', values='delta_V').values
    stab_grid = df_t.pivot(index='eta', columns='V', values='stability_type').values

    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 1. Set background to white
    ax.set_facecolor('white')

    # Common color range (uses log10(r0))
    valid_mask = (stab_grid == 'type-I') | (stab_grid == 'type-II')
    
    # Fixed range: -2 to 5
    vmin, vmax = -2, 5
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap('viridis')
    
    # Adjust this to slightly reduce the color map brightness
    cmap_colors = cmap(np.linspace(0, 1, 256))
    cmap_colors[:, :3] *= 1.0  # Reduce the RGB values ​​to reduce brightness
    cmap_adjusted = ListedColormap(cmap_colors)
    
    if not np.any(valid_mask):
        # If all data is type-III (unstable), display a message
        ax.text(0.5, 0.5, 'All type-III (Unstable)', ha='center', va='center',
               transform=ax.transAxes, fontsize=14, fontweight='bold')
        
        # Set axis limits to match the data range
        ax.set_xlim(V_grid.min(), V_grid.max())
        ax.set_ylim(Eta_grid.min(), Eta_grid.max())
        
        # Add topology and variant information to the title
        title_parts = ["Phase Diagram"]
        if topology and variant:
            title_parts.append(f"{topology} - {variant}")
        title_parts.append(f"$\\theta = {theta_val}$")
        ax.set_title(" | ".join(title_parts), fontsize=14, fontweight='bold')
        ax.set_xlabel('$V$ (vector trace torsion)')
        ax.set_ylabel('$\\eta$ (axial torsion)')
        
        # Add a color bar even when all type-III
        sm = ScalarMappable(norm=norm, cmap=cmap_adjusted)
        cbar = plt.colorbar(sm, ax=ax, ticks=np.arange(-2, 5.5, 0.5))
        cbar.set_label('$\\log_{10}(r_0)$', fontsize=12)
        
        # Create a legend
        legend_elements = [
            Patch(facecolor=cmap_adjusted(0.5), label='type-I (stable)'),
            Patch(facecolor=cmap_adjusted(0.5), hatch='////', label='type-II (rolling)'),
            Patch(facecolor='white', edgecolor='black', label='type-III (unstable)'),
            Line2D([0], [0], color='gray', linewidth=1, label='Contours (white): $\\log_{10}(\\Delta V)$')
        ]
        ax.legend(handles=legend_elements, loc='lower center', ncol=4, 
                 frameon=True, fontsize=9, bbox_to_anchor=(0.5, -0.15))
        
        plt.tight_layout()
        plt.savefig(output_filename, bbox_inches='tight')
        plt.close()
        return
    
    # Calculate log10 of r0
    r0_log = np.log10(np.where(valid_mask & (r0_grid > 0), r0_grid, np.nan))

    # 2. Draw both type-I and type-II heatmaps
    r0_valid = np.where(valid_mask, r0_log, np.nan)
    ax.pcolormesh(V_grid, Eta_grid, r0_valid, cmap=cmap_adjusted, norm=norm, shading='auto', alpha=1.0)

    # 3. Add hatching to the type-II region
    rolling_mask = np.where(stab_grid == 'type-II', 1, 0)
    if np.any(rolling_mask > 0):
        # Temporarily change RcParams to make the hatching lines thinner
        import matplotlib as mpl
        original_hatch_linewidth = mpl.rcParams['hatch.linewidth']
        mpl.rcParams['hatch.linewidth'] = 0.3  # default=1.0
        ax.contourf(V_grid, Eta_grid, rolling_mask, levels=[0.5, 1.5], 
                   colors='none', hatches=['////'], alpha=0)
        mpl.rcParams['hatch.linewidth'] = original_hatch_linewidth  # Revert

    # 4. Add a color bar (0.5 ticks)
    sm = ScalarMappable(norm=norm, cmap=cmap_adjusted)
    cbar = plt.colorbar(sm, ax=ax, ticks=np.arange(-2, 5.5, 0.5))
    cbar.set_label('$\\log_{10}(r_0)$', fontsize=12)

    # 5. Log10(delta_V) contours (1.0 intervals)
    dv_log = np.log10(np.where(valid_mask & (dv_grid > 0), dv_grid, np.nan))
    if not np.all(np.isnan(dv_log)):
        # Fix contour levels to integer values
        dv_min = np.floor(np.nanmin(dv_log))
        dv_max = np.ceil(np.nanmax(dv_log))
        levels = np.arange(dv_min, dv_max + 1, 1.0)
        cs = ax.contour(V_grid, Eta_grid, dv_log, levels=levels,
                      colors='white', alpha=0.5, linewidths=0.6)
        ax.clabel(cs, inline=True, fontsize=8, fmt='%.0f')

    # Add topology and variant information to the title
    title_parts = ["Phase Diagram"]
    if topology and variant:
        title_parts.append(f"{topology} - {variant}")
    title_parts.append(f"$\\theta = {theta_val}$")
    ax.set_title(" | ".join(title_parts), fontsize=14, fontweight='bold')
    ax.set_xlabel('$V$ (vector trace torsion)')
    ax.set_ylabel('$\\eta$ (axial torsion)')
    
    # Create a legend
    legend_elements = [
        Patch(facecolor=cmap_adjusted(0.5), label='type-I (stable)'),
        Patch(facecolor=cmap_adjusted(0.5), hatch='////', label='type-II (rolling)'),
        Patch(facecolor='white', edgecolor='black', label='type-III (unstable)'),
        Line2D([0], [0], color='gray', linewidth=1, label='Contours (white): $\\log_{10}(\\Delta V)$')
    ]
    ax.legend(handles=legend_elements, loc='lower center', ncol=4, 
             frameon=True, fontsize=9, bbox_to_anchor=(0.5, -0.15))
    
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: CSV file path is required.")
        print("Usage: python DPPUv2_visualize_phasemap_v3.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    print(f"Processing: {csv_file}")
    create_phase_diagrams(csv_file) 

