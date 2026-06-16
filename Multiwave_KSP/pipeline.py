# =============================================================================
# pipeline.py — Main crossmatch + data retrieval pipeline
# =============================================================================
# Workflow:
#   1. Load MAXI + Swift catalogs
#   2. For each MAXI source, cone-search Swift (X-ray crossmatch)
#   3. Fetch MAXI and Swift light curves and save locally
#   4. Detect outbursts from light curves
#   5. Optical crossmatch: ZTF first, ATLAS only if ZTF misses
#   6. Save optical light curves locally
#   7. Classify each source as Gold / Silver / Bronze
#      (separately for sources with and without outbursts)
#   8. Save final crossmatch table
#
# Usage:
#   python pipeline.py
#
# Requirements:
#   pip install alerce requests pandas numpy astropy matplotlib
# =============================================================================

import os
import ssl
import time
import numpy as np
import pandas as pd
from io import StringIO
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from astropy.coordinates import SkyCoord
import astropy.units as u
from alerce.core import Alerce

from config import (
    ATLAS_USER, ATLAS_PASS, ATLAS_BASEURL,
    ZTF_START_MJD, ATLAS_START_MJD,
    CONE_RADIUS_ARCSEC, OUTBURST_SIGMA_THRESHOLD,
    SWIFT_CSV, MAXI_CSV,
    LIGHTCURVE_DIR, RESULTS_DIR
)

alerce_client = Alerce()


def http_request(method, url, headers=None, data=None, timeout=30):
    import ssl, certifi
    ssl_ctx = ssl.create_default_context(cafile=certifi.where())
    req_headers = headers or {}
    req_data = None
    if data is not None:
        req_data = urlencode(data).encode("utf-8")
        req_headers = {**req_headers, "Content-Type": "application/x-www-form-urlencoded"}
    req = Request(url, data=req_data, headers=req_headers, method=method)
    with urlopen(req, context=ssl_ctx, timeout=timeout) as resp:
        body = resp.read().decode("utf-8")
        return resp.status, body, dict(resp.headers)


# =============================================================================
# 1. LOAD CATALOGS
# =============================================================================

def load_catalogs():
    swift = pd.read_csv(SWIFT_CSV)
    maxi  = pd.read_csv(MAXI_CSV)

    for df in [swift, maxi]:
        df['RA J2000 Degs']  = pd.to_numeric(df['RA J2000 Degs'],  errors='coerce')
        df['Dec J2000 Degs'] = pd.to_numeric(df['Dec J2000 Degs'], errors='coerce')

    swift = swift.dropna(subset=['RA J2000 Degs', 'Dec J2000 Degs']).reset_index(drop=True)
    maxi  = maxi.dropna(subset=['RA J2000 Degs', 'Dec J2000 Degs']).reset_index(drop=True)

    print(f"Loaded {len(maxi)} MAXI sources and {len(swift)} Swift sources.")
    return maxi, swift


# =============================================================================
# 2. SWIFT CONE SEARCH (X-ray crossmatch)
# =============================================================================

def crossmatch_swift(maxi_ra, maxi_dec, swift_df, radius_arcsec=CONE_RADIUS_ARCSEC):
    """
    Find the closest Swift source within radius_arcsec of (maxi_ra, maxi_dec).
    Returns the matching Swift row (Series) or None.
    """
    maxi_coord   = SkyCoord(ra=maxi_ra*u.deg, dec=maxi_dec*u.deg)
    swift_coords = SkyCoord(
        ra=swift_df['RA J2000 Degs'].values*u.deg,
        dec=swift_df['Dec J2000 Degs'].values*u.deg
    )
    sep = maxi_coord.separation(swift_coords).arcsec
    idx = sep.argmin()
    if sep[idx] <= radius_arcsec:
        return swift_df.iloc[idx], sep[idx]
    return None, None


# =============================================================================
# 3. LIGHT CURVE FETCHING — MAXI
# =============================================================================

def radec_to_maxi_id(ra, dec):
    coord   = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    ra_hms  = coord.ra.hms
    dec_dms = coord.dec.dms

    hh   = int(ra_hms.h)
    mm   = int(ra_hms.m)
    dd   = abs(int(dec_dms.d))    # 2 digits, not 3
    sign = '+' if dec >= 0 else '-'

    return f"J{hh:02d}{mm:02d}{sign}{dd:02d}"  # ← :02d not :03d


def fetch_maxi_lightcurve(source_name, ra, dec):
    save_path = os.path.join(LIGHTCURVE_DIR, "maxi", f"{source_name.replace(' ', '_')}.csv")
    if os.path.exists(save_path):
        return pd.read_csv(save_path)

    maxi_id  = radec_to_maxi_id(ra, dec)

    for suffix in ['_lc_1day_all.dat', '_lc_all.dat']:
        url = f"http://maxi.riken.jp/star_data/{maxi_id}/{maxi_id}{suffix}"
        try:
            req  = Request(url, method="GET")
            with urlopen(req, timeout=30) as resp:   # no ssl_ctx — plain HTTP
                text = resp.read().decode('utf-8')
            if len(text) > 100:
                df = pd.read_csv(
                    StringIO(text), sep=r'\s+', comment='!',
                    names=['MJD', 'rate_2_20keV', 'err_2_20keV',
                           'rate_2_4keV',   'err_2_4keV',
                           'rate_4_10keV',  'err_4_10keV',
                           'rate_10_20keV', 'err_10_20keV']
                )
                df.to_csv(save_path, index=False)
                return df
        except Exception as e:
            print(f"  MAXI fetch error for {source_name} ({url}): {e}")

    return None


# =============================================================================
# 4. LIGHT CURVE FETCHING — SWIFT/BAT
# =============================================================================

import certifi
import urllib.request

def fetch_swift_lightcurve(source_name):
    import ssl, certifi
    save_path = os.path.join(LIGHTCURVE_DIR, "swift", f"{source_name.replace(' ', '_')}.csv")
    if os.path.exists(save_path):
        return pd.read_csv(save_path)

    ssl_ctx  = ssl.create_default_context(cafile=certifi.where())
    name_fmt = source_name.replace(' ', '')

    for subdir in ['', 'weak/']:
        for suffix in ['.lc.txt', '.orbit.lc.txt']:
            url = f"https://swift.gsfc.nasa.gov/results/transients/{subdir}{name_fmt}{suffix}"
            try:
                req = Request(url, method="GET")
                with urlopen(req, context=ssl_ctx, timeout=30) as resp:
                    text = resp.read().decode('utf-8')
                if len(text) > 100:
                    df = pd.read_csv(
                        StringIO(text), sep=r'\s+', comment='!',
                        header=None        # let pandas infer, no fixed names
                    )
                    # Column 0=MJD, 1=rate, 2=error regardless of extra columns
                    df = df.iloc[:, [0, 1, 2]]
                    df.columns = ['MJD', 'rate_15_50keV', 'err_15_50keV']
                    df = df.apply(pd.to_numeric, errors='coerce').dropna()
                    df.to_csv(save_path, index=False)
                    return df
            except Exception as e:
                print(f"  Swift fetch error for {source_name} ({url}): {e}")

    return None

# =============================================================================
# 5. OUTBURST DETECTION
# =============================================================================

def detect_outburst(lc_df, rate_col, err_col, sigma=OUTBURST_SIGMA_THRESHOLD):
    """
    Simple outburst detector: flag if peak flux > median + sigma * MAD.
    Returns True/False.
    """
    if lc_df is None or len(lc_df) < 10:
        return False

    rates = lc_df[rate_col].dropna()
    rates = rates[rates > 0]
    if len(rates) < 5:
        return False

    median = rates.median()
    mad    = (rates - median).abs().median()
    return float(rates.max()) > (median + sigma * mad)


# =============================================================================
# 6. OPTICAL CROSSMATCH — ZTF via ALeRCE
# =============================================================================

def query_ztf(ra, dec, source_name, radius_arcsec=CONE_RADIUS_ARCSEC):
    """
    Query ZTF via ALeRCE. If a match is found, fetch and save the light curve.
    Returns (found: bool|None, ztf_oid: str|None, n_bands: int)
    """
    save_path = os.path.join(LIGHTCURVE_DIR, "ztf", f"{source_name.replace(' ', '_')}.csv")

    # If already fetched, skip API call
    if os.path.exists(save_path):
        df = pd.read_csv(save_path)
        n_bands = df['fid'].nunique() if 'fid' in df.columns else 1
        return True, df['oid'].iloc[0] if 'oid' in df.columns else None, n_bands

    try:
        result = alerce_client.query_objects(
            survey="ztf", ra=ra, dec=dec, radius=radius_arcsec
        )
        time.sleep(0.2)

        if result is None or len(result) == 0:
            return False, None, 0

        oid = result.iloc[0]['oid']

        # Fetch full light curve for the best match
        lc = alerce_client.query_lightcurve(oid, format='pandas')
        if lc is not None and len(lc) > 0:
            lc['oid'] = oid
            lc.to_csv(save_path, index=False)
            # ZTF fid: 1=g-band, 2=r-band, 3=i-band
            n_bands = lc['fid'].nunique() if 'fid' in lc.columns else 1
            return True, oid, n_bands

        return True, oid, 0  # found object but no light curve data

    except Exception as e:
        print(f"  ZTF error for {source_name}: {e}")
        return None, None, 0


# =============================================================================
# 7. OPTICAL CROSSMATCH — ATLAS
# =============================================================================

import ssl
import certifi

def get_atlas_token():
    ssl_ctx = ssl.create_default_context(cafile=certifi.where())
    # Patch requests to use certifi certificates
    import requests.adapters
    import urllib3
    session = requests.Session()
    session.verify = certifi.where()
    
    resp = session.post(
        f"{ATLAS_BASEURL}/api-token-auth/",
        data={"username": ATLAS_USER, "password": ATLAS_PASS},
        timeout=30
    )
    resp.raise_for_status()
    return resp.json()["token"]


def query_atlas(ra, dec, source_name, token, radius_arcsec=CONE_RADIUS_ARCSEC):
    """
    Submit an ATLAS forced photometry job, poll until complete, save results.
    Restricts MJD to ATLAS_START_MJD onwards to reduce query time.
    Returns (found: bool|None, n_bands: int)
    None = query failed or timed out; False = queried, no detections.
    """
    save_path = os.path.join(LIGHTCURVE_DIR, "atlas", f"{source_name.replace(' ', '_')}.csv")

    if os.path.exists(save_path):
        df = pd.read_csv(save_path)
        n_bands = df['F'].nunique() if 'F' in df.columns else 1
        return True, n_bands

    # ATLAS doesn't cover dec < -50
    if dec < -50:
        return None, 0

    headers = {"Authorization": f"Token {token}", "Accept": "application/json"}

    try:
        status, body, _ = http_request(
            "POST",
            f"{ATLAS_BASEURL}/queue/",
            headers=headers,
            data={
                "ra": ra, "dec": dec,
                "mjd_min": ATLAS_START_MJD,
                "send_email": False
            },
            timeout=30
        )
        if status != 201:
            print(f"  ATLAS submit failed for {source_name}: {status}")
            return None, 0

        task_url = pd.read_json(StringIO(body), typ="series")["url"]

        # Poll until done (max ~5 min)
        for _ in range(120):   # wait up to 10 minutes instead of 5
            time.sleep(5)
            _, result_body, _ = http_request("GET", task_url, headers=headers, timeout=30)
            result = pd.read_json(StringIO(result_body), typ="series").to_dict()
            if result.get("finishtimestamp"):
                phot_url = result.get("result_url")
                if not phot_url:
                    return False, 0

                _, phot_text, _ = http_request("GET", phot_url, headers=headers, timeout=60)
                df = pd.read_csv(StringIO(phot_text), sep=r'\s+', comment='#')

                if df.empty:
                    return False, 0

                # Filter to S/N > 3 detections
                if 'uJy' in df.columns and 'duJy' in df.columns:
                    df = df[df['uJy'] / df['duJy'] > 3]

                if len(df) == 0:
                    return False, 0

                df.to_csv(save_path, index=False)
                # ATLAS filters: 'o' (orange) and 'c' (cyan) stored in column 'F'
                n_bands = df['F'].nunique() if 'F' in df.columns else 1
                return True, n_bands

        print(f"  ATLAS timed out for {source_name}")
        return None, 0

    except Exception as e:
        print(f"  ATLAS error for {source_name}: {e}")
        return None, 0


# =============================================================================
# 8. CLASSIFICATION
# =============================================================================

def classify(row):
    """
    Classify a source as Gold / Silver / Bronze based on band coverage.

    X-ray bands:   MAXI (2–4 keV, 4–10 keV) + Swift (15–50 keV)
    Optical bands: ZTF (g, r, i) + ATLAS (c, o)

    Gold:   ≥2 X-ray bands AND ≥2 optical bands
    Silver: ≥1 X-ray band  AND ≥1 optical band  (but not Gold)
    Bronze: data in only one wavelength regime
    """
    xray_bands   = row['n_xray_bands']
    optical_bands = row['n_optical_bands']

    in_xray    = xray_bands > 0
    in_optical = optical_bands > 0

    if xray_bands >= 2 and optical_bands >= 2:
        return 'Gold'
    elif in_xray and in_optical:
        return 'Silver'
    elif in_xray or in_optical:
        return 'Bronze'
    else:
        return 'No match'


# =============================================================================
# 9. MAIN PIPELINE
# =============================================================================

def run_pipeline():
    maxi_df, swift_df = load_catalogs()

    # Get ATLAS token once
    try:
        atlas_token = get_atlas_token()
        print("ATLAS token obtained.")
    except Exception as e:
        atlas_token = None
        print(f"Warning: Could not get ATLAS token ({e}). ATLAS queries will be skipped.")

    results = []

    for i, row in maxi_df.iterrows():
        name = row['Source Name']
        ra   = row['RA J2000 Degs']
        dec  = row['Dec J2000 Degs']

        print(f"\n[{i+1}/{len(maxi_df)}] {name}  (RA={ra:.3f}, Dec={dec:.3f})")

        # ── X-ray: MAXI ──────────────────────────────────────────────────────
        maxi_lc = fetch_maxi_lightcurve(name, ra, dec)
        # MAXI has 2 science bands: 2–4 keV and 4–10 keV
        n_maxi_bands = 2 if maxi_lc is not None and len(maxi_lc) > 0 else 0

        # Outburst detection (use 4–10 keV band as primary)
        has_outburst = False
        if maxi_lc is not None and 'rate_4_10keV' in maxi_lc.columns:
            has_outburst = detect_outburst(maxi_lc, 'rate_4_10keV', 'err_4_10keV')

        # ── X-ray: Swift crossmatch ───────────────────────────────────────────
        swift_match, swift_sep = crossmatch_swift(ra, dec, swift_df)
        swift_name = swift_match['Source Name'] if swift_match is not None else None
        n_swift_bands = 0

        if swift_match is not None:
            swift_lc = fetch_swift_lightcurve(swift_name)
            n_swift_bands = 1 if swift_lc is not None and len(swift_lc) > 0 else 0
            print(f"  Swift match: {swift_name} ({swift_sep:.1f} arcsec)")
        else:
            print(f"  Swift: no match")

        n_xray_bands = n_maxi_bands + n_swift_bands

        # ── Optical: ZTF ─────────────────────────────────────────────────────
        ztf_found, ztf_oid, n_ztf_bands = query_ztf(ra, dec, name)
        print(f"  ZTF: {'found' if ztf_found else 'not found'}"
              + (f" ({n_ztf_bands} bands)" if ztf_found else ""))

        # ── Optical: ATLAS (only if ZTF missed) ──────────────────────────────
        atlas_found, n_atlas_bands = False, 0
        if not ztf_found and atlas_token:
            atlas_found, n_atlas_bands = query_atlas(ra, dec, name, atlas_token)
            print(f"  ATLAS: {'found' if atlas_found else 'not found'}"
                  + (f" ({n_atlas_bands} bands)" if atlas_found else ""))
        elif ztf_found:
            print(f"  ATLAS: skipped (ZTF data available)")

        n_optical_bands = n_ztf_bands + n_atlas_bands

        # ── Classify ─────────────────────────────────────────────────────────
        result = {
            'source_name'   : name,
            'ra'            : ra,
            'dec'           : dec,
            # X-ray
            'maxi_found'    : n_maxi_bands > 0,
            'n_maxi_bands'  : n_maxi_bands,
            'swift_match'   : swift_name,
            'swift_sep_arcsec': swift_sep,
            'n_swift_bands' : n_swift_bands,
            'n_xray_bands'  : n_xray_bands,
            # Outburst
            'has_outburst'  : has_outburst,
            # Optical
            'ztf_found'     : ztf_found,
            'ztf_oid'       : ztf_oid,
            'n_ztf_bands'   : n_ztf_bands,
            'atlas_found'   : atlas_found,
            'n_atlas_bands' : n_atlas_bands,
            'n_optical_bands': n_optical_bands,
        }

        result['classification'] = classify(result)
        print(f"  → {result['classification']} | outburst={has_outburst} "
              f"| X-ray bands={n_xray_bands} | optical bands={n_optical_bands}")

        results.append(result)

        # Save incrementally after every source so progress isn't lost
        pd.DataFrame(results).to_csv(
            os.path.join(RESULTS_DIR, "crossmatch_results.csv"), index=False
        )

    # ── Final summary ──────────────────────────────────────────────────────────
    df = pd.DataFrame(results)
    out_path = os.path.join(RESULTS_DIR, "crossmatch_results.csv")
    df.to_csv(out_path, index=False)

    print("\n" + "="*60)
    print("PIPELINE COMPLETE")
    print("="*60)
    print(f"\nTotal sources: {len(df)}")
    print("\nClassification breakdown:")
    print(df['classification'].value_counts().to_string())
    print("\nWith outbursts:")
    print(df[df['has_outburst']]['classification'].value_counts().to_string())
    print("\nWithout outbursts:")
    print(df[~df['has_outburst']]['classification'].value_counts().to_string())
    print(f"\nResults saved to: {out_path}")

    return df


if __name__ == "__main__":
    run_pipeline()
