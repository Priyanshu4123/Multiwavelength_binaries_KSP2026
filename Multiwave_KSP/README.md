# Multiwavelength Crossmatch Pipeline
## MAXI × Swift × ZTF × ATLAS

Reproducible pipeline for crossmatching X-ray transient sources across
X-ray (MAXI, Swift/BAT) and optical (ZTF, ATLAS) catalogs, generating
classified source lists and multi-band light curve plots.

---

## Setup

### 1. Install dependencies
```bash
pip install alerce requests pandas numpy astropy matplotlib
```

### 2. Register for ATLAS
Create a free account at https://fallingstar-data.com/forcedphot/register/
Then edit `config.py` and fill in your credentials:
```python
ATLAS_USER = "your_username"
ATLAS_PASS = "your_password"
```

### 3. Place catalog files in the working directory
- `Swift_BAT_Transient_Sources.csv`
- `MAXI_Sources.csv`

---

## Usage

### Step 1 — Run the crossmatch pipeline
```bash
python pipeline.py
```
This will:
- Crossmatch every MAXI source against Swift (cone search, 5 arcsec)
- Fetch and save MAXI and Swift/BAT light curves locally
- Query ZTF via ALeRCE for each source
- Query ATLAS only where ZTF has no data (run overnight for full catalog)
- Detect outbursts from X-ray light curves
- Classify each source as Gold / Silver / Bronze
- Save incremental results to `results/crossmatch_results.csv`

### Step 2 — Generate plots
```bash
python plot_lightcurves.py
```
This will:
- Generate a 4-panel plot for every Gold source
- Generate an all-band light curve for every source

---

## Output structure

```
results/
  crossmatch_results.csv      ← main classification table
  plots/
    gold_4panel/              ← 4-panel PDFs for Gold sources
    all_bands/                ← all-band light curve PDFs for all sources

data/
  lightcurves/
    maxi/                     ← MAXI light curves (.csv per source)
    swift/                    ← Swift/BAT light curves (.csv per source)
    ztf/                      ← ZTF light curves (.csv per source)
    atlas/                    ← ATLAS forced photometry (.csv per source)
```

All fetched data is saved locally — re-running the pipeline skips
already-downloaded sources, so queries are never repeated.

---

## Classification scheme

### Sources with outbursts
| Class  | Criteria |
|--------|----------|
| Gold   | Outburst captured by MAXI ± Swift AND optical coverage (ZTF or ATLAS) with ≥2 X-ray bands AND ≥2 optical bands |
| Silver | Both X-ray and optical coverage but incomplete band count |
| Bronze | Data in only one wavelength regime (X-ray only or optical only) |

### Sources without outbursts (quiescent)
| Class  | Criteria |
|--------|----------|
| Gold   | MAXI light curve present ± Swift AND optical coverage with ≥2 X-ray bands AND ≥2 optical bands |
| Silver | Both X-ray and optical data but incomplete band count |
| Bronze | Data in only one wavelength regime |

---

## Key settings (config.py)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CONE_RADIUS_ARCSEC` | 5 | Search radius for all cone searches |
| `OUTBURST_SIGMA_THRESHOLD` | 5 | Peak / median MAD threshold for outburst detection |
| `ZTF_START_MJD` | 58196 | ZTF survey start (2018-03-20) |
| `ATLAS_START_MJD` | 57940 | ATLAS search window start (2017-06-01) |

---

## Survey details

| Survey | Band | Coverage |
|--------|------|----------|
| MAXI | 2–4 keV, 4–10 keV | Full sky, continuous since 2009 |
| Swift/BAT | 15–50 keV | Full sky transient monitor |
| ZTF | g, r, i (optical) | Dec > −30°, since 2018 |
| ATLAS | cyan, orange (optical) | Dec > −50°, since ~2017 |

---

## Notes
- ATLAS queries are asynchronous and take ~5–30 seconds per source.
  For 400+ sources, run overnight: `nohup python pipeline.py > pipeline.log &`
- Progress is saved after every source, so the pipeline can be safely
  interrupted and restarted.
- If neither ZTF nor ATLAS finds a match, the source is flagged
  `optical_vizier_needed = True` for manual follow-up via VizieR.
