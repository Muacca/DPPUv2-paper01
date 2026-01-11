import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def create_phase_matrix_diagrams(csv_file):
    """
    Generate a 3x3 matrix phase diagram (topology x variant)
    Generate one image for each theta value
    """
    # Load data
    df = pd.read_csv(csv_file)
    
    # Get and sort unique theta values
    thetas = np.sort(df['theta'].unique())
    
    # Define topology and variant order
    topologies = ['S3', 'T3', 'Nil3']
    variants = ['FULL', 'TT', 'REE']
    
    # Generate a 3x3 matrix image for each theta value
    for idx, theta_val in enumerate(thetas):
        output_filename = f"phase_matrix_{idx:04d}.png"
        plot_phase_matrix(df, theta_val, topologies, variants, output_filename)
        print(f"Saved: {output_filename} (theta={theta_val})")

def plot_phase_matrix(df, theta_val, topologies, variants, output_filename):
    """
    Draw a 3x3 matrix for a specific theta value
    """
    # Create a 3x3 subplot
    fig, axes = plt.subplots(3, 3, figsize=(16, 13), sharex=True, sharey=True)
    
    # Colormap settings (same for all panels)
    vmin, vmax = -2, 5
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap('viridis')
    
    # Adjust this to slightly reduce the brightness of the colormap
    cmap_colors = cmap(np.linspace(0, 1, 256))
    cmap_colors[:, :3] *= 1.0  # Lower the RGB values ​​to reduce the brightness
    cmap_adjusted = ListedColormap(cmap_colors)
    
    # Draw each panel
    for i, topology in enumerate(topologies):
        for j, variant in enumerate(variants):
            ax = axes[i, j]
            ax.set_facecolor('white')  # Set background to white
            
            # Filter data
            df_filtered = df[
                (df['topology'] == topology) & 
                (df['ny_variant'] == variant) & 
                (df['theta'] == theta_val)
            ].copy()
            
            if len(df_filtered) < 4:
                # If data is missing
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', 
                       transform=ax.transAxes, fontsize=12)
                ax.set_title(f"{topology} - {variant}", fontsize=10, fontweight='bold')
                continue
            
            # Create a grid
            v_unique = np.sort(df_filtered['V'].unique())
            eta_unique = np.sort(df_filtered['eta'].unique())
            
            if len(v_unique) < 2 or len(eta_unique) < 2:
                ax.text(0.5, 0.5, 'Insufficient Data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=10)
                ax.set_title(f"{topology} - {variant}", fontsize=10, fontweight='bold')
                continue
            
            V_grid, Eta_grid = np.meshgrid(v_unique, eta_unique)
            
            # Pivot data
            r0_grid = df_filtered.pivot(index='eta', columns='V', values='r0').values
            dv_grid = df_filtered.pivot(index='eta', columns='V', values='delta_V').values
            stab_grid = df_filtered.pivot(index='eta', columns='V', values='stability_type').values
            
            # Mask valid data (type-I or type-II)
            valid_mask = (stab_grid == 'type-I') | (stab_grid == 'type-II')
            
            if not np.any(valid_mask):
                ax.text(0.5, 0.5, 'All type-III (Unstable)', ha='center', va='center',
                       transform=ax.transAxes, fontsize=10)
                ax.set_title(f"{topology} - {variant}", fontsize=10, fontweight='bold')
                continue
            
            # Calculate log10 of r0
            r0_log = np.log10(np.where(valid_mask & (r0_grid > 0), r0_grid, np.nan))
            
            # Draw both type-I and type-II heatmaps
            r0_valid = np.where(valid_mask, r0_log, np.nan)
            mesh = ax.pcolormesh(V_grid, Eta_grid, r0_valid, cmap=cmap_adjusted, norm=norm, 
                                shading='auto', alpha=1.0)
            
            # Add hatching to the type-II region
            rolling_mask = np.where(stab_grid == 'type-II', 1, 0)
            if np.any(rolling_mask > 0):
                # Temporarily change RcParams to make the hatching lines thinner
                import matplotlib as mpl
                original_hatch_linewidth = mpl.rcParams['hatch.linewidth']
                mpl.rcParams['hatch.linewidth'] = 0.3  # default=1.0
                ax.contourf(V_grid, Eta_grid, rolling_mask, levels=[0.5, 1.5], 
                           colors='none', hatches=['////'], alpha=0)
                mpl.rcParams['hatch.linewidth'] = original_hatch_linewidth  # Revert
            
            # Log10(delta_V) contours (1.0 intervals)
            dv_log = np.log10(np.where(valid_mask & (dv_grid > 0), dv_grid, np.nan))
            if not np.all(np.isnan(dv_log)):
                # Fix contour levels to integer values
                dv_min = np.floor(np.nanmin(dv_log))
                dv_max = np.ceil(np.nanmax(dv_log))
                levels = np.arange(dv_min, dv_max + 1, 1.0)
                cs = ax.contour(V_grid, Eta_grid, dv_log, levels=levels,
                              colors='white', alpha=0.5, linewidths=0.6)
                ax.clabel(cs, inline=True, fontsize=8, fmt='%.0f')
            
            # set title
            ax.set_title(f"{topology} - {variant}", fontsize=10, fontweight='bold')
            
            # Axis labels (outer panels only)
            if i == 2:  # Bottom row
                ax.set_xlabel('$V$ (vector trace torsion)', fontsize=9)
            if j == 0:  # Left column
                ax.set_ylabel('$\\eta$ (axial torsion)', fontsize=9)
    
    # Add a color bar to the right
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
    sm = ScalarMappable(norm=norm, cmap=cmap_adjusted)
    cbar = fig.colorbar(sm, cax=cbar_ax, ticks=np.arange(-2, 5.5, 0.5))
    cbar.set_label('$\\log_{10}(r_0)$', fontsize=12)
    
    # Create a legend
    legend_elements = [
        Patch(facecolor=cmap_adjusted(0.5), label='type-I (stable)'),
        Patch(facecolor=cmap_adjusted(0.5), hatch='////', label='type-II (rolling)'),
        Patch(facecolor='white', edgecolor='black', label='type-III (unstable)'),
        Line2D([0], [0], color='gray', linewidth=1, label='Contours (white): $\\log_{10}(\\Delta V)$')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, 
              frameon=True, fontsize=11, bbox_to_anchor=(0.5, 0.02))
    
    # Overall title
    fig.suptitle(f'Phase Diagram Matrix: Topology x NY Variant ($\\theta = {theta_val}$)',
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Error: CSV file path is required.")
        print("Usage: python DPPUv2_visualize_phasematrix_v3.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    print(f"Processing: {csv_file}")
    create_phase_matrix_diagrams(csv_file)
