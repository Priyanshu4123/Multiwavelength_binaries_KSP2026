# =============================================================================
# plot_lightcurves.py — Generate light curve plots for all classified sources
# =============================================================================
# Produces:
#   - For Gold sources: 4-panel plot (MAXI 2-4keV, Swift 15-50keV or MAXI 4-10keV,
#                                     ZTF g/r, ATLAS c/o)
#   - For all sources:  single combined light curve with all available bands
#
# Usage:
#   python plot_lightcurves.py
#
# Reads from:
#   results/crossmatch_results.csv
#   data/lightcurves/{maxi,swift,ztf,atlas}/<source_name>.csv
# =============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

from config import LIGHTCURVE_DIR, RESULTS_DIR, PLOT_DIR

# ── Plot style ────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.size'       : 11,
    'axes.titlesize'  : 11,
    'axes.labelsize'  : 11,
    'legend.fontsize' : 9,
    'figure.dpi'      : 120,
})

BAND_COLORS = {
    'maxi_2_4'  : '#e63946',   # red
    'maxi_4_10' : '#f4a261',   # orange
    'swift'     : '#457b9d',   # blue
    'ztf_g'     : '#2a9d8f',   # teal
    'ztf_r'     : '#e76f51',   # salmon
    'ztf_i'     : '#8338ec',   # purple
    'atlas_c'   : '#06d6a0',   # cyan
    'atlas_o'   : '#fb8500',   # amber
}

BAND_LABELS = {
    'maxi_2_4'  : 'MAXI 2–4 keV',
    'maxi_4_10' : 'MAXI 4–10 keV',
    'swift'     : 'Swift/BAT 15–50 keV',
    'ztf_g'     : 'ZTF g-band',
    'ztf_r'     : 'ZTF r-band',
    'ztf_i'     : 'ZTF i-band',
    'atlas_c'   : 'ATLAS cyan',
    'atlas_o'   : 'ATLAS orange',
}


# =============================================================================
# DATA LOADERS
# =============================================================================

def load_maxi(source_name):
    path = os.path.join(LIGHTCURVE_DIR, "maxi", f"{source_name.replace(' ', '_')}.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df = df.dropna(subset=['MJD'])
    df['MJD'] = pd.to_numeric(df['MJD'], errors='coerce')
    return df


def load_swift(source_name):
    path = os.path.join(LIGHTCURVE_DIR, "swift", f"{source_name.replace(' ', '_')}.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df = df.dropna(subset=['MJD'])
    df['MJD'] = pd.to_numeric(df['MJD'], errors='coerce')
    return df


def load_ztf(source_name):
    path = os.path.join(LIGHTCURVE_DIR, "ztf", f"{source_name.replace(' ', '_')}.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    # ALeRCE uses 'mjd', fid: 1=g, 2=r, 3=i
    if 'mjd' in df.columns:
        df = df.rename(columns={'mjd': 'MJD'})
    df['MJD'] = pd.to_numeric(df['MJD'], errors='coerce')
    return df


def load_atlas(source_name):
    path = os.path.join(LIGHTCURVE_DIR, "atlas", f"{source_name.replace(' ', '_')}.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    # ATLAS uses 'MJD' and filter column 'F' with values 'c' or 'o'
    df['MJD'] = pd.to_numeric(df['MJD'], errors='coerce')
    return df


# =============================================================================
# PANEL HELPERS
# =============================================================================

def _errorbar_panel(ax, mjd, rate, err, color, label, ylabel, title=None):
    """Draw a single light curve panel with error bars."""
    mask = np.isfinite(mjd) & np.isfinite(rate) & np.isfinite(err)
    if mask.sum() == 0:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                ha='center', va='center', color='grey')
    else:
        ax.errorbar(mjd[mask], rate[mask], yerr=err[mask],
                    fmt='.', ms=3, lw=0.6, color=color,
                    alpha=0.7, label=label)
        ax.axhline(0, color='grey', lw=0.5, ls='--')
        ax.legend(loc='upper right', framealpha=0.5)

    ax.set_ylabel(ylabel)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if title:
        ax.set_title(title, fontsize=10)


# =============================================================================
# PLOT 1 — 4-PANEL PLOT (Gold sources)
# =============================================================================

def plot_4panel(source_name, has_outburst, outdir=None):
    """
    4-panel figure:
      Panel 1: MAXI 2–4 keV
      Panel 2: MAXI 4–10 keV  OR  Swift 15–50 keV (whichever is available)
      Panel 3: ZTF g + r bands (or ATLAS c if no ZTF)
      Panel 4: ZTF i band      (or ATLAS o if no ZTF)
    """
    maxi_lc  = load_maxi(source_name)
    swift_lc = load_swift(source_name)
    ztf_lc   = load_ztf(source_name)
    atlas_lc = load_atlas(source_name)

    outdir = outdir or os.path.join(PLOT_DIR, "gold_4panel")
    os.makedirs(outdir, exist_ok=True)

    fig = plt.figure(figsize=(12, 10))
    fig.suptitle(
        f"{source_name}  —  Gold sample"
        f"  ({'with outburst' if has_outburst else 'quiescent'})",
        fontsize=13, fontweight='bold', y=0.98
    )
    gs = gridspec.GridSpec(4, 1, hspace=0.08, figure=fig)
    axes = [fig.add_subplot(gs[i]) for i in range(4)]

    # ── Panel 1: MAXI 2–4 keV ────────────────────────────────────────────────
    if maxi_lc is not None and 'rate_2_4keV' in maxi_lc.columns:
        _errorbar_panel(axes[0],
                        maxi_lc['MJD'].values,
                        maxi_lc['rate_2_4keV'].values,
                        maxi_lc['err_2_4keV'].values,
                        BAND_COLORS['maxi_2_4'],
                        BAND_LABELS['maxi_2_4'],
                        ylabel='Count rate\n(ct/s)')
    else:
        axes[0].text(0.5, 0.5, 'MAXI 2–4 keV — no data',
                     transform=axes[0].transAxes, ha='center', va='center', color='grey')
        axes[0].set_ylabel('Count rate\n(ct/s)')

    # ── Panel 2: Swift 15–50 keV (preferred) or MAXI 4–10 keV ───────────────
    if swift_lc is not None and 'rate_15_50keV' in swift_lc.columns:
        _errorbar_panel(axes[1],
                        swift_lc['MJD'].values,
                        swift_lc['rate_15_50keV'].values,
                        swift_lc['err_15_50keV'].values,
                        BAND_COLORS['swift'],
                        BAND_LABELS['swift'],
                        ylabel='Count rate\n(ct/s)')
    elif maxi_lc is not None and 'rate_4_10keV' in maxi_lc.columns:
        _errorbar_panel(axes[1],
                        maxi_lc['MJD'].values,
                        maxi_lc['rate_4_10keV'].values,
                        maxi_lc['err_4_10keV'].values,
                        BAND_COLORS['maxi_4_10'],
                        BAND_LABELS['maxi_4_10'],
                        ylabel='Count rate\n(ct/s)')
    else:
        axes[1].text(0.5, 0.5, 'No second X-ray band',
                     transform=axes[1].transAxes, ha='center', va='center', color='grey')
        axes[1].set_ylabel('Count rate\n(ct/s)')

    # ── Panel 3: ZTF g + r  or  ATLAS cyan ───────────────────────────────────
    axes[2].axhline(0, color='grey', lw=0.5, ls='--')
    if ztf_lc is not None and 'fid' in ztf_lc.columns:
        for fid, key, label in [(1, 'ztf_g', BAND_LABELS['ztf_g']),
                                  (2, 'ztf_r', BAND_LABELS['ztf_r'])]:
            sub = ztf_lc[ztf_lc['fid'] == fid]
            if len(sub) > 0 and 'magpsf' in sub.columns:
                axes[2].errorbar(sub['MJD'].values, sub['magpsf'].values,
                                 yerr=sub['sigmapsf'].values if 'sigmapsf' in sub.columns else None,
                                 fmt='.', ms=3, lw=0.6,
                                 color=BAND_COLORS[key], alpha=0.7, label=label)
        axes[2].invert_yaxis()
        axes[2].set_ylabel('Magnitude')
        axes[2].legend(loc='upper right', framealpha=0.5)
    elif atlas_lc is not None and 'F' in atlas_lc.columns:
        sub = atlas_lc[atlas_lc['F'] == 'c']
        if len(sub) > 0:
            axes[2].errorbar(sub['MJD'].values, sub['uJy'].values,
                             yerr=sub['duJy'].values if 'duJy' in sub.columns else None,
                             fmt='.', ms=3, lw=0.6,
                             color=BAND_COLORS['atlas_c'], alpha=0.7,
                             label=BAND_LABELS['atlas_c'])
        axes[2].set_ylabel('Flux (µJy)')
        axes[2].legend(loc='upper right', framealpha=0.5)
    else:
        axes[2].text(0.5, 0.5, 'No optical band 1 data',
                     transform=axes[2].transAxes, ha='center', va='center', color='grey')
        axes[2].set_ylabel('Flux / Mag')

    # ── Panel 4: ZTF i  or  ATLAS orange ─────────────────────────────────────
    axes[3].axhline(0, color='grey', lw=0.5, ls='--')
    if ztf_lc is not None and 'fid' in ztf_lc.columns:
        sub = ztf_lc[ztf_lc['fid'] == 3]
        if len(sub) > 0 and 'magpsf' in sub.columns:
            axes[3].errorbar(sub['MJD'].values, sub['magpsf'].values,
                             yerr=sub['sigmapsf'].values if 'sigmapsf' in sub.columns else None,
                             fmt='.', ms=3, lw=0.6,
                             color=BAND_COLORS['ztf_i'], alpha=0.7,
                             label=BAND_LABELS['ztf_i'])
            axes[3].invert_yaxis()
        axes[3].set_ylabel('Magnitude')
        axes[3].legend(loc='upper right', framealpha=0.5)
    elif atlas_lc is not None and 'F' in atlas_lc.columns:
        sub = atlas_lc[atlas_lc['F'] == 'o']
        if len(sub) > 0:
            axes[3].errorbar(sub['MJD'].values, sub['uJy'].values,
                             yerr=sub['duJy'].values if 'duJy' in sub.columns else None,
                             fmt='.', ms=3, lw=0.6,
                             color=BAND_COLORS['atlas_o'], alpha=0.7,
                             label=BAND_LABELS['atlas_o'])
        axes[3].set_ylabel('Flux (µJy)')
        axes[3].legend(loc='upper right', framealpha=0.5)
    else:
        axes[3].text(0.5, 0.5, 'No optical band 2 data',
                     transform=axes[3].transAxes, ha='center', va='center', color='grey')
        axes[3].set_ylabel('Flux / Mag')

    # Shared x-axis
    for ax in axes[:-1]:
        ax.set_xticklabels([])
    axes[-1].set_xlabel('MJD')

    fname = os.path.join(outdir, f"{source_name.replace(' ', '_')}_4panel.pdf")
    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {fname}")


# =============================================================================
# PLOT 2 — ALL-BAND LIGHT CURVE (all sources)
# =============================================================================

def plot_all_bands(source_name, has_outburst, classification, outdir=None):
    """
    Single figure with all available bands overplotted on a shared time axis,
    using twin y-axes for X-ray (count rate) and optical (magnitude/flux).
    """
    maxi_lc  = load_maxi(source_name)
    swift_lc = load_swift(source_name)
    ztf_lc   = load_ztf(source_name)
    atlas_lc = load_atlas(source_name)

    outdir = outdir or os.path.join(PLOT_DIR, "all_bands")
    os.makedirs(outdir, exist_ok=True)

    fig, ax1 = plt.subplots(figsize=(12, 4))
    ax2 = ax1.twinx()   # optical axis

    plotted = False

    # ── X-ray panels ──────────────────────────────────────────────────────────
    if maxi_lc is not None:
        for col, err_col, key in [
            ('rate_2_4keV',  'err_2_4keV',  'maxi_2_4'),
            ('rate_4_10keV', 'err_4_10keV', 'maxi_4_10'),
        ]:
            if col in maxi_lc.columns:
                m = np.isfinite(maxi_lc['MJD']) & np.isfinite(maxi_lc[col])
                ax1.errorbar(maxi_lc['MJD'][m], maxi_lc[col][m],
                             yerr=maxi_lc[err_col][m],
                             fmt='.', ms=2, lw=0.5,
                             color=BAND_COLORS[key], alpha=0.6,
                             label=BAND_LABELS[key])
                plotted = True

    if swift_lc is not None and 'rate_15_50keV' in swift_lc.columns:
        m = np.isfinite(swift_lc['MJD']) & np.isfinite(swift_lc['rate_15_50keV'])
        ax1.errorbar(swift_lc['MJD'][m], swift_lc['rate_15_50keV'][m],
                     yerr=swift_lc['err_15_50keV'][m],
                     fmt='.', ms=2, lw=0.5,
                     color=BAND_COLORS['swift'], alpha=0.6,
                     label=BAND_LABELS['swift'])
        plotted = True

    # ── Optical ───────────────────────────────────────────────────────────────
    if ztf_lc is not None and 'fid' in ztf_lc.columns:
        for fid, key in [(1, 'ztf_g'), (2, 'ztf_r'), (3, 'ztf_i')]:
            sub = ztf_lc[ztf_lc['fid'] == fid]
            if len(sub) > 0 and 'magpsf' in sub.columns:
                m = np.isfinite(sub['MJD']) & np.isfinite(sub['magpsf'])
                ax2.errorbar(sub['MJD'][m], sub['magpsf'][m],
                             yerr=sub['sigmapsf'][m] if 'sigmapsf' in sub.columns else None,
                             fmt='.', ms=2, lw=0.5,
                             color=BAND_COLORS[key], alpha=0.6,
                             label=BAND_LABELS[key])
                plotted = True
        ax2.invert_yaxis()
        ax2.set_ylabel('Magnitude (ZTF)')

    elif atlas_lc is not None and 'F' in atlas_lc.columns:
        for filt, key in [('c', 'atlas_c'), ('o', 'atlas_o')]:
            sub = atlas_lc[atlas_lc['F'] == filt]
            if len(sub) > 0 and 'uJy' in sub.columns:
                m = np.isfinite(sub['MJD']) & np.isfinite(sub['uJy'])
                ax2.errorbar(sub['MJD'][m], sub['uJy'][m],
                             yerr=sub['duJy'][m] if 'duJy' in sub.columns else None,
                             fmt='.', ms=2, lw=0.5,
                             color=BAND_COLORS[key], alpha=0.6,
                             label=BAND_LABELS[key])
                plotted = True
        ax2.set_ylabel('Flux µJy (ATLAS)')

    ax1.set_xlabel('MJD')
    ax1.set_ylabel('X-ray count rate (ct/s)')
    ax1.axhline(0, color='grey', lw=0.5, ls='--')

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc='upper left', fontsize=8, framealpha=0.5)

    outburst_str = 'outburst' if has_outburst else 'quiescent'
    ax1.set_title(f"{source_name}  —  {classification}  ({outburst_str})",
                  fontsize=11, fontweight='bold')

    if not plotted:
        ax1.text(0.5, 0.5, 'No light curve data available',
                 transform=ax1.transAxes, ha='center', va='center', color='grey')

    fname = os.path.join(outdir, f"{source_name.replace(' ', '_')}_allbands.pdf")
    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {fname}")


# =============================================================================
# MAIN
# =============================================================================

def run_plots():
    results_path = os.path.join(RESULTS_DIR, "crossmatch_results.csv")
    if not os.path.exists(results_path):
        raise FileNotFoundError(
            f"{results_path} not found. Run pipeline.py first."
        )

    df = pd.read_csv(results_path)
    print(f"Generating plots for {len(df)} sources...")

    for _, row in df.iterrows():
        name           = row['source_name']
        classification = row['classification']
        has_outburst   = bool(row['has_outburst'])

        print(f"\n{name}  [{classification}]")

        # 4-panel plot — Gold sources only
        if classification == 'Gold':
            plot_4panel(name, has_outburst)

        # All-band light curve — every source
        plot_all_bands(name, has_outburst, classification)

    print("\nAll plots done.")
    print(f"  4-panel plots → {os.path.join(PLOT_DIR, 'gold_4panel')}/")
    print(f"  All-band plots → {os.path.join(PLOT_DIR, 'all_bands')}/")


if __name__ == "__main__":
    run_plots()
