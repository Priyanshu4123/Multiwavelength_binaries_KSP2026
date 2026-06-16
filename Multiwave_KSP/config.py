# =============================================================================
# config.py — Configuration and shared utilities for the crossmatch pipeline
# =============================================================================
# Edit the credentials and paths below before running the pipeline.
# All other scripts import from this file.

import os

# All paths are resolved relative to this config.py file, so the pipeline
# works correctly regardless of which directory you run Python from.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Credentials ───────────────────────────────────────────────────────────────
ATLAS_USER = "#######"   # Register at fallingstar-data.com/forcedphot
ATLAS_PASS = "#######"

# ── Survey start dates (MJD) ──────────────────────────────────────────────────
# Restricting the ATLAS search window dramatically reduces query time.
# ZTF began: 2018-03-20  → MJD 58196
# ATLAS began: 2015-07-01 → MJD 57204 (but good all-sky from ~2017)
ZTF_START_MJD   = 58196   # 2018-03-20
ATLAS_START_MJD = 57940   # 2017-06-01  (conservative start for reliable coverage)

# ── Search radius ─────────────────────────────────────────────────────────────
CONE_RADIUS_ARCSEC = 5     # arcseconds — tune if needed for your source types

# ── Outburst detection ────────────────────────────────────────────────────────
# A source is flagged as "has outburst" if its peak flux exceeds this multiple
# of its median flux in the MAXI or Swift light curve.
OUTBURST_SIGMA_THRESHOLD = 5   # sigma above median

# ── File paths ────────────────────────────────────────────────────────────────
# CSV catalogs — expected in the same folder as this config.py
SWIFT_CSV  = os.path.join(BASE_DIR, "Swift_BAT_Transient_Sources.csv")
MAXI_CSV   = os.path.join(BASE_DIR, "MAXI_Sources.csv")

DATA_DIR        = os.path.join(BASE_DIR, "data")
RESULTS_DIR     = os.path.join(BASE_DIR, "results")
LIGHTCURVE_DIR  = os.path.join(DATA_DIR, "lightcurves")
PLOT_DIR        = os.path.join(RESULTS_DIR, "plots")

for d in [DATA_DIR, RESULTS_DIR, LIGHTCURVE_DIR, PLOT_DIR,
          os.path.join(LIGHTCURVE_DIR, "maxi"),
          os.path.join(LIGHTCURVE_DIR, "swift"),
          os.path.join(LIGHTCURVE_DIR, "ztf"),
          os.path.join(LIGHTCURVE_DIR, "atlas")]:
    os.makedirs(d, exist_ok=True)

# ── ATLAS API base URL ────────────────────────────────────────────────────────
ATLAS_BASEURL = "https://fallingstar-data.com/forcedphot"
