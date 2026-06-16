# =============================================================================
# poll_atlas.py — Standalone ATLAS result collector
# =============================================================================
# Run this AFTER pipeline_3.py Phase 1 completes.
# It reads the submitted job URLs from results/atlas_jobs.json,
# polls each one until complete, saves the light curves, and updates
# results/crossmatch_results.csv with the new optical band counts
# and reclassified Gold/Silver/Bronze labels.
#
# Can be safely re-run multiple times — already-collected sources are skipped.
#
# Usage:
#   python poll_atlas.py
# =============================================================================

import os
import sys
import json
import time
import socket
import pandas as pd
from io import StringIO

# Import shared utilities from pipeline_2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_3 import http_request, get_atlas_token

from config import (
    ATLAS_BASEURL,
    LIGHTCURVE_DIR,
    RESULTS_DIR
)

ATLAS_JOBS_PATH    = os.path.join(RESULTS_DIR, "atlas_jobs.json")
RESULTS_PATH       = os.path.join(RESULTS_DIR, "crossmatch_results.csv")
ATLAS_LIGHTCURVE_DIR = os.path.join(LIGHTCURVE_DIR, "atlas")

os.makedirs(ATLAS_LIGHTCURVE_DIR, exist_ok=True)


# =============================================================================
# CLASSIFICATION (mirrors pipeline_3.py)
# =============================================================================

def classify(n_xray_bands, n_optical_bands):
    in_xray    = n_xray_bands > 0
    in_optical = n_optical_bands > 0
    if n_xray_bands >= 2 and n_optical_bands >= 2:
        return 'Gold'
    elif in_xray and in_optical:
        return 'Silver'
    elif in_xray or in_optical:
        return 'Bronze'
    else:
        return 'No match'


# =============================================================================
# ATLAS RESULT COLLECTOR
# =============================================================================

def collect_atlas_result(task_url, source_name, token, save_path):
    """
    Poll a single ATLAS job until it finishes, then save results.
    Returns (found: bool|None, n_bands: int).
    None = timed out or network error; False = queried, no detections.
    """
    print(f"    Polling: {task_url}")
    headers = {"Authorization": f"Token {token}", "Accept": "application/json"}

    # check status only once
    for attempt in range(1):  

        # Force fresh DNS resolution each attempt
        try:
            socket.getaddrinfo('fallingstar-data.com', 443)
        except Exception:
            pass

        try:
            status, result_body, _ = http_request("GET", task_url, headers=headers, timeout=30)
        except Exception as e:
            print(f"    Poll error (attempt {attempt+1}): {e}")
            continue

        if status != 200:
            continue

        result = json.loads(result_body)
        queuepos = result.get("queuepos")
        finished = result.get("finished")

        print(
            f"    queue={queuepos} "
            f"finished={finished}"
        )

        if not result.get("finishtimestamp"):
            return None, 0

        # Job finished
        phot_url = result.get("result_url")
        if not phot_url:
            print(f"  Job finished but no result_url for {source_name}")
            return False, 0

        try:
            _, phot_text, _ = http_request("GET", phot_url, headers=headers, timeout=60)
        except Exception as e:
            print(f"  Failed to fetch results for {source_name}: {e}")
            return False, 0

        df = pd.read_csv(StringIO(phot_text), sep=r'\s+', comment='#')

        if df.empty:
            print(f"  Empty result for {source_name}")
            return False, 0

        if 'uJy' in df.columns and 'duJy' in df.columns:
            df = df[df['uJy'] / df['duJy'] > 3]

        if len(df) == 0:
            print(f"  No detections above S/N>3 for {source_name}")
            return False, 0

        df.to_csv(save_path, index=False)
        n_bands = int(df['F'].nunique()) if 'F' in df.columns else 1
        return True, n_bands

    print(f"  Timed out for {source_name}")
    return False, 0


# =============================================================================
# MAIN
# =============================================================================

def run():
    # ── Load job list ─────────────────────────────────────────────────────────
    if not os.path.exists(ATLAS_JOBS_PATH):
        print(f"ERROR: {ATLAS_JOBS_PATH} not found.")
        print("Run pipeline_3.py first to submit ATLAS jobs and generate this file.")
        sys.exit(1)

    with open(ATLAS_JOBS_PATH) as f:
        atlas_jobs = json.load(f)
    print(f"Loaded {len(atlas_jobs)} ATLAS jobs from {ATLAS_JOBS_PATH}")

    # ── Load existing results ─────────────────────────────────────────────────
    if not os.path.exists(RESULTS_PATH):
        print(f"ERROR: {RESULTS_PATH} not found.")
        print("Run pipeline_3.py first to generate the crossmatch results table.")
        sys.exit(1)

    df = pd.read_csv(RESULTS_PATH)

    # allow True / False / None
    df['atlas_found'] = df['atlas_found'].astype('object')

    df = df.set_index('source_name')

    # ── Get fresh ATLAS token ─────────────────────────────────────────────────
    try:
        token = get_atlas_token()
        print("ATLAS token obtained.\n")
    except Exception as e:
        print(f"ERROR: Could not get ATLAS token: {e}")
        sys.exit(1)

    # ── Poll each job ─────────────────────────────────────────────────────────
    total     = len(atlas_jobs)
    completed = 0
    skipped   = 0

    for i, (name, task_url) in enumerate(atlas_jobs.items()):
        save_path = os.path.join(ATLAS_LIGHTCURVE_DIR, f"{name.replace(' ', '_')}.csv")

        print(f"\n[{i+1}/{total}] {name}")

        # Skip if already collected
        if os.path.exists(save_path):
            existing  = pd.read_csv(save_path)
            n_bands   = int(existing['F'].nunique()) if 'F' in existing.columns else 1
            print(f"  Already cached ({n_bands} bands) — skipping")
            skipped += 1

            if name in df.index:
                df.at[name, 'atlas_found']     = True
                df.at[name, 'n_atlas_bands']   = n_bands
                n_optical = int(df.at[name, 'n_ztf_bands'] or 0) + n_bands
                df.at[name, 'n_optical_bands']  = n_optical
                df.at[name, 'classification']   = classify(
                    int(df.at[name, 'n_xray_bands'] or 0), n_optical
                )
            continue

        # Poll the job
        found, n_bands = collect_atlas_result(task_url, name, token, save_path)
        if found is None:
            print("  → still queued")
            if name in df.index:
                df.at[name, 'atlas_found'] = pd.NA
            continue

        n_bands = n_bands or 0
        completed += 1

        print(f"  → {'found' if found else 'not found'} ({n_bands} bands)")

        # Update results table
        if name in df.index:
            df.at[name, 'atlas_found']    = found
            df.at[name, 'n_atlas_bands']  = n_bands
            n_optical = int(df.at[name, 'n_ztf_bands'] or 0) + n_bands
            df.at[name, 'n_optical_bands'] = n_optical
            df.at[name, 'classification']  = classify(
                int(df.at[name, 'n_xray_bands'] or 0), n_optical
            )
        else:
            print(f"  WARNING: {name} not found in crossmatch_results.csv")

        # Save progress after every source so a crash doesn't lose everything
        df.reset_index().to_csv(RESULTS_PATH, index=False)

    # ── Final summary ─────────────────────────────────────────────────────────
    df_final = df.reset_index()
    df_final.to_csv(RESULTS_PATH, index=False)

    print("\n" + "="*60)
    print("ATLAS POLLING COMPLETE")
    print("="*60)
    print(f"  Newly polled : {completed}")
    print(f"  From cache   : {skipped}")
    print(f"  Total jobs   : {total}")
    print("\nUpdated classification breakdown:")
    print(df_final['classification'].value_counts().to_string())
    print(f"\nResults saved to: {RESULTS_PATH}")


if __name__ == "__main__":
    run()