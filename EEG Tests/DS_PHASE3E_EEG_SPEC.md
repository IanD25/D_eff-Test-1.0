# DS Phase 3E: EEG Seizure Classification Transfer — Simplified Spec

**Version:** 1.0 | **Date:** 2026-03-06 | **Platform:** Python (local)
**Runtime:** ≤15 min on M4 MacBook Air | **Branch:** `phase3e-eeg`

## Objective

Test whether the three information-geometric states manifest during epileptic seizures — a genuine synchronization phase transition in neural tissue. Identical FIM pipeline to Phase 2/3B, zero neuroscience-specific tuning.

**Mapping:** Pre-ictal (normal brain) ≈ State 1 (short-range, modular). Ictal (seizure) ≈ State 2 (system-spanning synchronization). Post-ictal (suppression) ≈ State 3 (flat correlations, noise-dominated).

## Data Source

**CHB-MIT Scalp EEG Database** (PhysioNet, free, no registration)
- 23 pediatric patients, 198 annotated seizures
- 23 EEG channels, 256 Hz sampling
- Annotations: seizure start/end times per file
- URL: https://physionet.org/content/chbmit/1.0.0/
- Python access: `mne` + `wfdb` or direct download

**We use patients chb01–chb05** (5 patients, ~30 seizures) for the initial test. Enough for signal detection; extend to all 23 if positive.

## Kill Condition

**P3E-1:** SV₂/SV₁ must increase by ≥0.05 (absolute) during ictal vs. pre-ictal baseline in at least 3 of 5 patients. If fewer than 3 → KILL.

## Design Summary

1. For each seizure file: load 23-channel EEG, segment into rolling windows (2-second windows, 1-second step)
2. Per window: compute 23×23 channel cross-correlation matrix → k-NN graph (k=5) → BFS distances → C(r) → FIM → SVD → diagnostics
3. Epoch each seizure into: pre-ictal (120s before onset), ictal (seizure duration), post-ictal (120s after end)
4. Output: one CSV row per window with all diagnostics + epoch label + patient/seizure ID
5. Analysis done in chat

## Output Format

Single file: `phase3e_results.csv`

```
patient,seizure_id,file,time_sec,epoch,sv2sv1,rank,eta,pr,n_valid
```

Where `epoch` ∈ {pre_ictal, ictal, post_ictal, interictal}.

Plus `phase3e_summary.txt` — per-seizure means for each epoch.

## Pre-Registered Predictions

| ID | Prediction | PASS | KILL |
|---|---|---|---|
| **P3E-1** | SV₂/SV₁ ictal > pre-ictal by ≥0.05 in ≥3/5 patients | Elevation detected | **<3 patients → KILL** |
| P3E-2 | SV₂/SV₁ post-ictal < ictal (State 3 suppression) | Drop after seizure | — |
| P3E-3 | Rank decreases toward 1 during ictal (single-mode dominance) | Rank drops | — |
| P3E-4 | η peaks during or immediately after ictal | η extremum at seizure | — |
| P3E-5 | SV₂/SV₁ begins rising ≥10s before annotated seizure onset | Pre-ictal lead | — |

## Implementation Notes

- `mne` handles EDF file reading and channel selection
- Use only standard 10-20 EEG channels (exclude ECG/EMG if present)
- Cross-correlation computed via `np.corrcoef` on windowed segments
- k=5 for k-NN (23 channels is small — k=5 gives diameter ~3-4)
- Window = 512 samples (2s at 256Hz), step = 256 samples (1s)
- Bandpass filter 1–50 Hz before correlation (remove DC drift and line noise)

## Dependencies

```
pip install mne wfdb numpy scipy matplotlib
```

---

## Complete Code

```python
"""
DS Phase 3E: EEG Seizure Classification Transfer
Fisher Information Matrix applied to EEG channel correlations during epileptic seizures.

Usage:
    python eeg_fisher_phase3e.py

Downloads CHB-MIT data automatically via wfdb. Results saved to phase3e_results/.
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.signal import butter, filtfilt
import os, csv, time, warnings
warnings.filterwarnings("ignore")

# ============================================================
# Configuration
# ============================================================
PATIENTS = ["chb01", "chb02", "chb03", "chb04", "chb05"]
DATA_DIR = "chb-mit-data"
OUTPUT_DIR = "phase3e_results"
SFREQ = 256          # CHB-MIT sampling rate
WIN_SAMPLES = 512    # 2 seconds
STEP_SAMPLES = 256   # 1 second step
PRE_ICTAL_SEC = 120  # 2 min before seizure
POST_ICTAL_SEC = 120 # 2 min after seizure
K_NN = 5
N_FISHER_SAMPLES = 10
N_CR_SAMPLES = 15
BANDPASS = (1.0, 50.0)  # Hz

# Standard 10-20 channel names to keep (CHB-MIT uses these with dashes)
STANDARD_CHANNELS = [
    "FP1-F7", "F7-T7", "T7-P7", "P7-O1",
    "FP1-F3", "F3-C3", "C3-P3", "P3-O1",
    "FP2-F4", "F4-C4", "C4-P4", "P4-O2",
    "FP2-F8", "F8-T8", "T8-P8", "P8-O2",
    "FZ-CZ", "CZ-PZ",
    "P7-T7", "T7-FT9", "FT9-FT10", "FT10-T8", "T8-P8",
]

# ============================================================
# Data loading
# ============================================================
def download_patient_data(patient):
    """Download seizure files for one patient using wfdb."""
    import wfdb
    patient_dir = os.path.join(DATA_DIR, patient)
    os.makedirs(patient_dir, exist_ok=True)

    # Parse summary file for seizure annotations
    summary_url = f"chbmit/{patient}/{patient}-summary.txt"
    try:
        # Download summary
        import urllib.request
        summary_path = os.path.join(patient_dir, f"{patient}-summary.txt")
        if not os.path.exists(summary_path):
            url = f"https://physionet.org/files/chbmit/1.0.0/{patient}/{patient}-summary.txt"
            urllib.request.urlretrieve(url, summary_path)

        seizures = parse_summary(summary_path)
        print(f"  {patient}: {len(seizures)} seizure files found")

        # Download each seizure file
        for sz in seizures:
            fname = sz["file"].replace(".edf", "")
            rec_path = os.path.join(patient_dir, fname)
            if not os.path.exists(rec_path + ".hea"):
                try:
                    wfdb.dl_database(f"chbmit/{patient}", patient_dir, records=[fname])
                except Exception as e:
                    print(f"    Download failed for {fname}: {e}")
        return seizures
    except Exception as e:
        print(f"  Error loading {patient}: {e}")
        return []


def parse_summary(summary_path):
    """Parse CHB-MIT summary file for seizure start/end times."""
    seizures = []
    with open(summary_path) as f:
        lines = f.readlines()

    current_file = None
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("File Name:"):
            current_file = line.split(":")[1].strip()
        elif line.startswith("Seizure") and "Start" in line and "Time" in line:
            # Handle both "Seizure Start Time:" and "Seizure 1 Start Time:"
            start_sec = int(line.split(":")[-1].strip().split()[0])
            i += 1
            end_line = lines[i].strip()
            end_sec = int(end_line.split(":")[-1].strip().split()[0])
            if current_file:
                seizures.append({
                    "file": current_file,
                    "start": start_sec,
                    "end": end_sec,
                })
        i += 1
    return seizures


def load_eeg_segment(patient_dir, filename, start_sec, duration_sec):
    """Load EEG segment, return (n_channels, n_samples) array and channel names."""
    import wfdb
    fname = filename.replace(".edf", "")
    rec_path = os.path.join(patient_dir, fname)

    try:
        record = wfdb.rdrecord(rec_path,
                               sampfrom=int(start_sec * SFREQ),
                               sampto=int((start_sec + duration_sec) * SFREQ))
        signals = record.p_signal.T  # (n_channels, n_samples)
        ch_names = [ch.upper().replace(" ", "").replace(".", "-")
                    for ch in record.sig_name]
        return signals, ch_names
    except Exception as e:
        return None, None


# ============================================================
# Signal processing
# ============================================================
def bandpass_filter(data, low, high, fs, order=4):
    """Apply bandpass filter to each channel."""
    b, a = butter(order, [low / (fs / 2), high / (fs / 2)], btype="band")
    filtered = np.zeros_like(data)
    for i in range(data.shape[0]):
        try:
            filtered[i] = filtfilt(b, a, data[i])
        except:
            filtered[i] = data[i]
    return filtered


def select_eeg_channels(signals, ch_names):
    """Keep only standard 10-20 EEG channels."""
    indices = []
    kept_names = []
    upper_names = [ch.upper().replace(" ", "") for ch in ch_names]
    for i, name in enumerate(upper_names):
        # Match against standard channels (flexible matching)
        for std in STANDARD_CHANNELS:
            if std.replace("-", "") in name.replace("-", "") or name.replace("-", "") in std.replace("-", ""):
                indices.append(i)
                kept_names.append(name)
                break
    if len(indices) < 8:
        # Fallback: use all non-ECG/EMG channels
        for i, name in enumerate(upper_names):
            if not any(x in name for x in ["ECG", "EMG", "VNS", "STI", "DC", "--"]):
                if i not in indices:
                    indices.append(i)
                    kept_names.append(name)
    return signals[indices], kept_names


# ============================================================
# Fisher pipeline (identical to Phase 2/3B)
# ============================================================
def build_knn(corr, k):
    n = corr.shape[0]
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        row = corr[i].copy()
        row[i] = -np.inf
        top_k = np.argsort(row)[-k:]
        for j in top_k:
            adj[i, j] = True
            adj[j, i] = True
    return adj


def compute_cr(corr, adj, n):
    adj_sp = csr_matrix(adj.astype(float))
    samples = np.random.choice(n, min(N_CR_SAMPLES, n), replace=False)
    dist = shortest_path(adj_sp, method="D", indices=samples, directed=False)
    max_r = int(np.nanmax(dist[np.isfinite(dist)]))
    if max_r < 2:
        return None
    max_r = min(max_r, 10)
    C_r = np.zeros(max_r + 1)
    counts = np.zeros(max_r + 1)
    for idx_i, src in enumerate(samples):
        for dst in range(n):
            if src == dst:
                continue
            d = dist[idx_i, dst]
            if np.isfinite(d) and int(d) <= max_r:
                C_r[int(d)] += corr[src, dst]
                counts[int(d)] += 1
    mask = counts > 0
    C_r[mask] /= counts[mask]
    return C_r


def kernel(C_r, distances, n):
    weights = np.zeros(n)
    for u in range(n):
        d = distances[u]
        if np.isfinite(d) and int(d) < len(C_r):
            weights[u] = abs(C_r[int(d)])
    total = np.sum(weights)
    if total < 1e-15:
        return None
    return weights / total


def fisher_diagnostics(C_r, adj, n):
    adj_sp = csr_matrix(adj.astype(float))
    neighbors = {i: list(np.where(adj[i])[0]) for i in range(n)}
    samples = np.random.choice(n, min(N_FISHER_SAMPLES, n), replace=False)
    ranks, etas, prs, sv2sv1s = [], [], [], []

    for v0 in samples:
        nbrs = neighbors.get(v0, [])
        if len(nbrs) < 2:
            continue
        sources = [v0] + list(nbrs)
        try:
            dists = shortest_path(adj_sp, method="D", indices=sources, directed=False)
        except:
            continue
        p_v0 = kernel(C_r, dists[0], n)
        if p_v0 is None:
            continue
        k = len(nbrs)
        scores = np.zeros((k, n))
        valid = True
        for j in range(k):
            p_wj = kernel(C_r, dists[j + 1], n)
            if p_wj is None:
                valid = False
                break
            with np.errstate(divide="ignore", invalid="ignore"):
                log_ratio = np.log(p_wj) - np.log(p_v0)
                log_ratio = np.nan_to_num(log_ratio, nan=0.0, posinf=0.0, neginf=0.0)
            scores[j] = log_ratio
        if not valid:
            continue
        weighted = scores * np.sqrt(p_v0)[np.newaxis, :]
        F = weighted @ weighted.T
        sv = np.linalg.svd(F, compute_uv=False)
        if sv[0] < 1e-15:
            continue
        sv = sv / sv[0]
        gaps = sv[:-1] / (sv[1:] + 1e-15)
        rank = int(np.argmax(gaps) + 1)
        eta = float(sv[min(rank, len(sv) - 1)] / (sv[rank - 1] + 1e-15))
        pr = float(np.sum(sv) ** 2 / (np.sum(sv ** 2) + 1e-15))
        s2s1 = float(sv[1] / (sv[0] + 1e-15)) if len(sv) > 1 else 0.0
        ranks.append(rank)
        etas.append(eta)
        prs.append(pr)
        sv2sv1s.append(s2s1)

    if len(ranks) == 0:
        return None
    return {
        "sv2sv1": float(np.mean(sv2sv1s)),
        "rank": float(np.mean(ranks)),
        "eta": float(np.mean(etas)),
        "pr": float(np.mean(prs)),
        "n_valid": len(ranks),
    }


# ============================================================
# Main analysis loop
# ============================================================
def process_seizure(patient, seizure_info, patient_dir):
    """Process one seizure file: compute Fisher diagnostics in rolling windows."""
    fname = seizure_info["file"]
    sz_start = seizure_info["start"]
    sz_end = seizure_info["end"]
    sz_duration = sz_end - sz_start

    # Load segment: pre-ictal (120s) + ictal + post-ictal (120s)
    load_start = max(0, sz_start - PRE_ICTAL_SEC)
    load_end = sz_end + POST_ICTAL_SEC
    load_duration = load_end - load_start

    signals, ch_names = load_eeg_segment(patient_dir, fname, load_start, load_duration)
    if signals is None:
        return []

    signals, ch_names = select_eeg_channels(signals, ch_names)
    n_channels = signals.shape[0]
    if n_channels < 8:
        print(f"    Skipping {fname}: only {n_channels} channels")
        return []

    # Bandpass filter
    signals = bandpass_filter(signals, BANDPASS[0], BANDPASS[1], SFREQ)

    n_samples = signals.shape[1]
    results = []

    # Sliding window
    win_start = 0
    while win_start + WIN_SAMPLES <= n_samples:
        window_data = signals[:, win_start:win_start + WIN_SAMPLES]

        # Time relative to seizure onset
        abs_time = load_start + win_start / SFREQ
        rel_time = abs_time - sz_start  # negative = pre-ictal

        # Epoch label
        if abs_time < sz_start:
            if abs_time >= sz_start - PRE_ICTAL_SEC:
                epoch = "pre_ictal"
            else:
                epoch = "interictal"
        elif abs_time <= sz_end:
            epoch = "ictal"
        else:
            epoch = "post_ictal"

        # Correlation matrix
        corr = np.corrcoef(window_data)
        if np.any(np.isnan(corr)):
            corr = np.nan_to_num(corr, nan=0.0)
        np.fill_diagonal(corr, 1.0)
        n = corr.shape[0]

        # Fisher pipeline
        adj = build_knn(corr, min(K_NN, n - 1))
        C_r = compute_cr(corr, adj, n)

        if C_r is not None and len(C_r) >= 2:
            diag = fisher_diagnostics(C_r, adj, n)
        else:
            diag = None

        if diag:
            results.append({
                "patient": patient,
                "seizure_file": fname,
                "time_sec": round(abs_time, 1),
                "rel_onset_sec": round(rel_time, 1),
                "epoch": epoch,
                "sv2sv1": round(diag["sv2sv1"], 4),
                "rank": round(diag["rank"], 2),
                "eta": round(diag["eta"], 4),
                "pr": round(diag["pr"], 3),
                "n_valid": diag["n_valid"],
                "n_channels": n_channels,
            })

        win_start += STEP_SAMPLES

    return results


def write_summary(all_results, summary_path):
    """Write per-patient, per-epoch summary statistics."""
    from collections import defaultdict

    # Group by patient and epoch
    groups = defaultdict(lambda: defaultdict(list))
    for r in all_results:
        groups[r["patient"]][r["epoch"]].append(r)

    with open(summary_path, "w") as f:
        f.write("DS Phase 3E: EEG Seizure Results Summary\n")
        f.write("=" * 60 + "\n\n")

        # Per-patient summary
        for patient in sorted(groups.keys()):
            f.write(f"\n--- {patient} ---\n")
            f.write(f"{'Epoch':<15} {'N':>5} {'SV2/SV1':>10} {'Rank':>8} {'Eta':>8} {'PR':>8}\n")
            for epoch in ["pre_ictal", "ictal", "post_ictal"]:
                rows = groups[patient].get(epoch, [])
                if rows:
                    sv = [r["sv2sv1"] for r in rows]
                    rk = [r["rank"] for r in rows]
                    et = [r["eta"] for r in rows]
                    pr = [r["pr"] for r in rows]
                    f.write(f"{epoch:<15} {len(sv):>5} "
                            f"{np.mean(sv):>10.4f} {np.mean(rk):>8.2f} "
                            f"{np.mean(et):>8.4f} {np.mean(pr):>8.3f}\n")
            f.write("\n")

        # P3E-1 kill test
        f.write("\n" + "=" * 60 + "\n")
        f.write("P3E-1 KILL TEST: SV2/SV1 ictal > pre-ictal by >= 0.05\n")
        f.write("-" * 60 + "\n")
        pass_count = 0
        for patient in sorted(groups.keys()):
            pre = groups[patient].get("pre_ictal", [])
            ict = groups[patient].get("ictal", [])
            if pre and ict:
                pre_mean = np.mean([r["sv2sv1"] for r in pre])
                ict_mean = np.mean([r["sv2sv1"] for r in ict])
                delta = ict_mean - pre_mean
                verdict = "PASS" if delta >= 0.05 else "FAIL"
                if verdict == "PASS":
                    pass_count += 1
                f.write(f"  {patient}: pre={pre_mean:.4f} ictal={ict_mean:.4f} "
                        f"delta={delta:+.4f} → {verdict}\n")
        f.write(f"\nP3E-1 OVERALL: {pass_count}/5 patients pass. "
                f"{'KILL TEST PASSES' if pass_count >= 3 else 'KILL TEST FAILS'}\n")

        # P3E-2 through P3E-4
        f.write("\n" + "=" * 60 + "\n")
        f.write("P3E-2: SV2/SV1 post-ictal < ictal\n")
        for patient in sorted(groups.keys()):
            ict = groups[patient].get("ictal", [])
            post = groups[patient].get("post_ictal", [])
            if ict and post:
                ict_m = np.mean([r["sv2sv1"] for r in ict])
                post_m = np.mean([r["sv2sv1"] for r in post])
                f.write(f"  {patient}: ictal={ict_m:.4f} post={post_m:.4f} "
                        f"→ {'PASS' if post_m < ict_m else 'FAIL'}\n")

        f.write("\nP3E-3: Rank decreases during ictal\n")
        for patient in sorted(groups.keys()):
            pre = groups[patient].get("pre_ictal", [])
            ict = groups[patient].get("ictal", [])
            if pre and ict:
                pre_r = np.mean([r["rank"] for r in pre])
                ict_r = np.mean([r["rank"] for r in ict])
                f.write(f"  {patient}: pre_rank={pre_r:.2f} ictal_rank={ict_r:.2f} "
                        f"→ {'PASS' if ict_r < pre_r else 'FAIL'}\n")

        f.write("\nP3E-4: Eta peaks during/after ictal\n")
        for patient in sorted(groups.keys()):
            pre = groups[patient].get("pre_ictal", [])
            ict = groups[patient].get("ictal", [])
            post = groups[patient].get("post_ictal", [])
            if pre and ict:
                pre_e = np.mean([r["eta"] for r in pre])
                ict_e = np.mean([r["eta"] for r in ict])
                post_e = np.mean([r["eta"] for r in post]) if post else 0
                peak = max(ict_e, post_e)
                f.write(f"  {patient}: pre_eta={pre_e:.4f} ictal_eta={ict_e:.4f} "
                        f"post_eta={post_e:.4f} → {'PASS' if peak > pre_e else 'FAIL'}\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)

    all_results = []
    t0 = time.time()

    for patient in PATIENTS:
        print(f"\nProcessing {patient}...")
        patient_dir = os.path.join(DATA_DIR, patient)
        seizures = download_patient_data(patient)

        if not seizures:
            print(f"  No seizures found for {patient}")
            continue

        for si, sz in enumerate(seizures):
            print(f"  Seizure {si+1}/{len(seizures)}: {sz['file']} "
                  f"[{sz['start']}s - {sz['end']}s]")
            results = process_seizure(patient, sz, patient_dir)
            all_results.extend(results)
            print(f"    → {len(results)} windows processed")

    # Write CSV
    csv_path = os.path.join(OUTPUT_DIR, "phase3e_results.csv")
    if all_results:
        fields = ["patient", "seizure_file", "time_sec", "rel_onset_sec",
                   "epoch", "sv2sv1", "rank", "eta", "pr", "n_valid", "n_channels"]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for row in all_results:
                writer.writerow(row)
        print(f"\nWrote {len(all_results)} rows to {csv_path}")

    # Write summary
    summary_path = os.path.join(OUTPUT_DIR, "phase3e_summary.txt")
    write_summary(all_results, summary_path)
    print(f"Wrote summary to {summary_path}")

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.1f}s")
    print(f"Upload {csv_path} and {summary_path} to chat for analysis.")


if __name__ == "__main__":
    main()
```

---

## Validation Gates (checked in code)

- **Gate A:** ≥8 EEG channels after filtering non-EEG → else skip seizure
- **Gate B:** C(r) has ≥2 distance values → else skip window
- **Gate C:** ≥1 valid FIM from n_fisher_samples → else skip window

## Handback

Upload `phase3e_results.csv` and `phase3e_summary.txt` to chat. I will evaluate P3E-1 through P3E-5, generate time-resolved plots around seizure onsets, and assess three-state classification transfer.

## Key Differences from Phase 3B

| Aspect | Phase 3B (Financial) | Phase 3E (EEG) |
|--------|---------------------|----------------|
| N nodes | ~90 stocks | ~18-23 channels |
| k-NN k | 10 | 5 |
| Window | 30/60/90 trading days | 2 seconds (512 samples) |
| Step | 5 days | 1 second |
| Temporal resolution | Weekly | Sub-second |
| Ground truth | Crisis dates (approximate) | Seizure annotations (exact) |
| Expected State 2 | SV₂/SV₁ ≈ 0.65-0.72 | Unknown — may approach 1.0 |
| Pipeline | Identical | Identical |
