# DS Phase 3E-B: EEG Seizure — Frequency-Band Coherence Multiplex FIM
## Resolution Fix via Multiplex Graph Construction

**Version:** 1.0
**Date:** 2026-03-25
**Platform:** Python (local)
**Runtime:** ≤30 min on M4 MacBook Air, 24GB RAM
**Branch:** `phase3e-b-multiplex`
**Same data:** CHB-MIT chb01–chb05 (no new downloads if already cached)
**Depends on:** phase3e results (INCONCLUSIVE — resolution-limited)

---

## Why This Spec Exists

Phase 3E applied the FIM to raw 23-channel correlation matrices and
returned INCONCLUSIVE. Root cause: graph diameter ~2 gives C(r) only
2-3 distance bins — below the minimum viable resolution for the FIM
pipeline to extract spectral structure.

This spec fixes the resolution problem without new data by constructing
a **frequency-band coherence multiplex graph** from the same 23
channels. Four frequency bands × 23 channels = 92 effective nodes with
a graph diameter of ~5-8, sufficient for C(r) to have 5+ bins.

The fix is scientifically motivated, not post-hoc tuning: seizure
dynamics are inherently band-specific (theta, alpha, beta, gamma each
carry different information about synchronization state). The multiplex
graph captures functional connectivity in frequency space — the
biologically appropriate object for seizure analysis.

---

## The Key Change from Phase 3E

| Aspect | Phase 3E (v1) | Phase 3E-B (this spec) |
|--------|--------------|----------------------|
| Input | 23×23 raw correlation | 92×92 multiplex coherence |
| Graph nodes | 23 channels | 92 (23 channels × 4 bands) |
| k-NN k | 5 | 6 |
| Graph diameter | ~2 hops | ~5-8 hops |
| C(r) bins | 2-3 | 5-8 |
| FIM matrix | 5×5 | 6×6 |
| Kernel | |ρ_raw| | |coh_band| |
| Resolution | BELOW minimum | ABOVE minimum |

Everything else is identical to Phase 3E: same data, same patients,
same seizure annotations, same FIM pipeline (build_knn → compute_cr →
fisher_diagnostics), same output format, same kill condition.

---

## The Multiplex Graph Construction

### Four frequency bands

| Band | Range (Hz) | Seizure relevance |
|------|-----------|-------------------|
| Theta | 4–8 Hz | Pre-ictal build-up, hippocampal involvement |
| Alpha | 8–13 Hz | Suppressed during ictal, recovers post-ictal |
| Beta | 13–30 Hz | Focal onset marker, desynchronizes ictal |
| Gamma | 30–50 Hz | High-frequency oscillations, ictal onset |

### Node labeling

Node index = band × 23 + channel_index

    Nodes 0–22:   theta band, channels 0–22
    Nodes 23–45:  alpha band, channels 0–22
    Nodes 46–68:  beta band, channels 0–22
    Nodes 69–91:  gamma band, channels 0–22

### Edge types

**Within-band edges** (spatial): coherence between channels i and j
in the same band. Computed as magnitude-squared coherence (MSC):

    MSC_band(i, j) = |Cxy_band|² / (Cxx_band × Cyy_band)

where Cxy is the cross-spectral density in the band. Range [0, 1].

**Cross-band edges** (spectral): coherence between the SAME channel
across adjacent bands (theta-alpha, alpha-beta, beta-gamma). This
captures the cross-frequency coupling known to change at seizure onset.

    Cross(i_band1, i_band2) = |Cxy_cross| — same channel, adjacent bands

Cross-band edges connect node (band_a, ch_i) to node (band_b, ch_i)
where |band_a - band_b| = 1.

### Why this gives diameter ~5-8

Within-band graph: 23 channels, k=6 NN → diameter ~3-4.
Cross-band edges add "vertical" connections: a path can go
theta→alpha→beta→gamma, traversing 3 extra hops.
Total effective diameter: within-band diameter + cross-band hops
≈ 3 + 3 = ~6. Well above the minimum of 5.

---

## Data Source

Same as Phase 3E. If already downloaded:
    chb-mit-data/chb01/ ... chb05/

If not downloaded, the script reuses the same download logic.

URL: https://physionet.org/content/chbmit/1.0.0/

---

## Kill Condition

**P3E-B-1 (primary kill):** SV₂/SV₁ ictal > pre-ictal by ≥0.05
(absolute) in at least 3 of 5 patients. Same as Phase 3E.

**P3E-B-0 (resolution gate, checked first):** C(r) must have ≥5
distance bins in ≥80% of windows for ≥4/5 patients. If this gate
fails, the resolution fix didn't work — report and stop.

If P3E-B-0 passes and P3E-B-1 fails: genuine KILL. The FIM signal
is not present in EEG seizure data even at adequate resolution.

If P3E-B-0 passes and P3E-B-1 passes: PASS. Resolution was the
issue in Phase 3E; the physics hypothesis is supported.

---

## Pre-Registered Predictions

| ID | Prediction | PASS | KILL/FAIL |
|----|-----------|------|-----------|
| P3E-B-0 | C(r) bins ≥5 in ≥80% of windows | Graph resolved | Stop — design failed |
| **P3E-B-1** | SV₂/SV₁ ictal > pre-ictal ≥0.05 in ≥3/5 patients | Signal present | **KILL** |
| P3E-B-2 | SV₂/SV₁ post-ictal < ictal | State 3 drop | — |
| P3E-B-3 | Rank decreases during ictal | Single-mode dominance | — |
| P3E-B-4 | η peaks during/after ictal | Disorder extremum | — |
| P3E-B-5 | SV₂/SV₁ rises ≥10s before onset | Pre-ictal lead | — |
| P3E-B-6 | Gamma band coherence drives SV swap | Gamma-specific | — |

P3E-B-6 is new: the gamma band is the first to show high-frequency
oscillations at ictal onset. If the SV swap is driven by the gamma
subgraph specifically (gamma within-band edges dominating the FIM
score vectors), that provides mechanistic specificity beyond a scalar
classification.

---

## Output Format

```
phase3e_b_results.csv
```

Same schema as Phase 3E plus two new columns:

```
patient, seizure_file, time_sec, rel_onset_sec, epoch,
sv2sv1, rank, eta, pr, n_valid, n_channels,
n_cr_bins,        ← NEW: number of C(r) distance bins (resolution gate)
dominant_band     ← NEW: which band contributes most to FIM score variance
```

Plus `phase3e_b_summary.txt` — per-patient epoch means with
resolution gate check.

---

## Complete Code

```python
"""
DS Phase 3E-B: EEG Seizure — Frequency-Band Coherence Multiplex FIM

Fixes the resolution problem from Phase 3E by constructing a 92-node
multiplex graph (23 channels × 4 frequency bands) instead of the
23-node raw correlation graph.

Usage:
    python eeg_fisher_phase3e_b.py

Reuses CHB-MIT data from chb-mit-data/ if already downloaded.
Results saved to phase3e_b_results/.
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.signal import butter, filtfilt, welch, csd
import os, csv, time, warnings
warnings.filterwarnings("ignore")

# ============================================================
# Configuration
# ============================================================
PATIENTS = ["chb01", "chb02", "chb03", "chb04", "chb05"]
DATA_DIR = "chb-mit-data"
OUTPUT_DIR = "phase3e_b_results"
SFREQ = 256
WIN_SAMPLES = 512    # 2 seconds
STEP_SAMPLES = 256   # 1 second
PRE_ICTAL_SEC = 120
POST_ICTAL_SEC = 120
K_NN = 6             # Increased from 5 — larger graph needs more edges
N_FISHER_SAMPLES = 10
N_CR_SAMPLES = 20    # More samples for larger graph
BANDPASS_GLOBAL = (1.0, 50.0)

# Frequency bands (Hz)
BANDS = {
    "theta": (4.0, 8.0),
    "alpha": (8.0, 13.0),
    "beta":  (13.0, 30.0),
    "gamma": (30.0, 50.0),
}
BAND_NAMES = ["theta", "alpha", "beta", "gamma"]
N_BANDS = len(BAND_NAMES)

# Minimum C(r) bins for resolution gate P3E-B-0
MIN_CR_BINS = 5

# Standard 10-20 channels
STANDARD_CHANNELS = [
    "FP1-F7", "F7-T7", "T7-P7", "P7-O1",
    "FP1-F3", "F3-C3", "C3-P3", "P3-O1",
    "FP2-F4", "F4-C4", "C4-P4", "P4-O2",
    "FP2-F8", "F8-T8", "T8-P8", "P8-O2",
    "FZ-CZ", "CZ-PZ",
    "P7-T7", "T7-FT9", "FT9-FT10", "FT10-T8", "T8-P8",
]


# ============================================================
# Data loading (reused from Phase 3E)
# ============================================================
def download_patient_data(patient):
    """Download seizure files and parse summary."""
    import urllib.request
    patient_dir = os.path.join(DATA_DIR, patient)
    os.makedirs(patient_dir, exist_ok=True)

    summary_path = os.path.join(patient_dir, f"{patient}-summary.txt")
    if not os.path.exists(summary_path):
        url = (f"https://physionet.org/files/chbmit/1.0.0/"
               f"{patient}/{patient}-summary.txt")
        try:
            urllib.request.urlretrieve(url, summary_path)
        except Exception as e:
            print(f"  Summary download failed for {patient}: {e}")
            return []

    seizures = parse_summary(summary_path)
    print(f"  {patient}: {len(seizures)} seizure files found")

    import wfdb
    for sz in seizures:
        fname = sz["file"].replace(".edf", "")
        rec_path = os.path.join(patient_dir, fname)
        if not os.path.exists(rec_path + ".hea"):
            try:
                wfdb.dl_database(f"chbmit/{patient}",
                                 patient_dir, records=[fname])
            except Exception as e:
                print(f"    Download failed for {fname}: {e}")
    return seizures


def parse_summary(summary_path):
    """Parse CHB-MIT summary file."""
    seizures = []
    with open(summary_path) as f:
        lines = f.readlines()
    current_file = None
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("File Name:"):
            current_file = line.split(":")[1].strip()
        elif "Seizure" in line and "Start" in line and "Time" in line:
            try:
                start_sec = int(line.split(":")[-1].strip().split()[0])
                i += 1
                end_sec = int(lines[i].strip().split(":")[-1].strip().split()[0])
                if current_file:
                    seizures.append({
                        "file": current_file,
                        "start": start_sec,
                        "end": end_sec,
                    })
            except (ValueError, IndexError):
                pass
        i += 1
    return seizures


def load_eeg_segment(patient_dir, filename, start_sec, duration_sec):
    """Load EEG segment."""
    import wfdb
    fname = filename.replace(".edf", "")
    rec_path = os.path.join(patient_dir, fname)
    try:
        record = wfdb.rdrecord(
            rec_path,
            sampfrom=int(start_sec * SFREQ),
            sampto=int((start_sec + duration_sec) * SFREQ)
        )
        signals = record.p_signal.T
        ch_names = [ch.upper().replace(" ", "").replace(".", "-")
                    for ch in record.sig_name]
        return signals, ch_names
    except Exception:
        return None, None


def select_eeg_channels(signals, ch_names):
    """Keep standard 10-20 EEG channels."""
    indices = []
    upper = [c.replace("-", "").upper() for c in ch_names]
    for i, name in enumerate(upper):
        for std in STANDARD_CHANNELS:
            s = std.replace("-", "").upper()
            if s in name or name in s:
                indices.append(i)
                break
    if len(indices) < 8:
        for i, name in enumerate(upper):
            if not any(x in name for x in
                       ["ECG", "EMG", "VNS", "STI", "DC", "--"]):
                if i not in indices:
                    indices.append(i)
    return signals[indices], [ch_names[i] for i in indices]


def bandpass_filter(data, low, high, fs, order=4):
    """Bandpass filter each channel."""
    b, a = butter(order, [low / (fs / 2), high / (fs / 2)], btype="band")
    out = np.zeros_like(data)
    for i in range(data.shape[0]):
        try:
            out[i] = filtfilt(b, a, data[i])
        except Exception:
            out[i] = data[i]
    return out


# ============================================================
# Multiplex coherence graph construction
# ============================================================
def compute_band_coherence(signals, fs, band_low, band_high):
    """
    Compute magnitude-squared coherence (MSC) between all channel pairs
    in a given frequency band.

    Returns: (n_ch, n_ch) coherence matrix, values in [0, 1].
    """
    n_ch, n_samp = signals.shape
    nperseg = min(256, n_samp // 2)
    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    band_mask = (freqs >= band_low) & (freqs <= band_high)

    if band_mask.sum() == 0:
        return np.eye(n_ch) * 0.5

    coh = np.zeros((n_ch, n_ch))
    for i in range(n_ch):
        for j in range(i, n_ch):
            if i == j:
                coh[i, j] = 1.0
                continue
            try:
                f, Pxy = csd(signals[i], signals[j],
                             fs=fs, nperseg=nperseg)
                f, Pxx = welch(signals[i], fs=fs, nperseg=nperseg)
                f, Pyy = welch(signals[j], fs=fs, nperseg=nperseg)
                mask = (f >= band_low) & (f <= band_high)
                if mask.sum() == 0:
                    coh[i, j] = coh[j, i] = 0.0
                    continue
                num = np.abs(Pxy[mask]) ** 2
                den = Pxx[mask] * Pyy[mask]
                msc = np.mean(num / (den + 1e-15))
                coh[i, j] = coh[j, i] = float(np.clip(msc, 0.0, 1.0))
            except Exception:
                coh[i, j] = coh[j, i] = 0.0

    return coh


def build_multiplex_coherence(signals, fs):
    """
    Build 92×92 multiplex coherence matrix from 23 channels × 4 bands.

    Node index: band_idx * n_ch + ch_idx
        0-22:  theta
        23-45: alpha
        46-68: beta
        69-91: gamma

    Within-band edges: MSC between channels in same band.
    Cross-band edges: same channel across adjacent bands
                      (theta↔alpha, alpha↔beta, beta↔gamma).

    Returns: (N_total × N_total) coherence matrix where N_total = n_ch × n_bands
    """
    n_ch = signals.shape[0]
    n_total = n_ch * N_BANDS
    C = np.zeros((n_total, n_total))

    # Within-band coherence
    band_cohs = {}
    for b_idx, band_name in enumerate(BAND_NAMES):
        low, high = BANDS[band_name]
        coh = compute_band_coherence(signals, fs, low, high)
        band_cohs[b_idx] = coh
        # Fill diagonal block
        i_start = b_idx * n_ch
        C[i_start:i_start + n_ch, i_start:i_start + n_ch] = coh

    # Cross-band edges: same channel, adjacent bands
    # Weight = geometric mean of within-band coherence at that channel
    for b_idx in range(N_BANDS - 1):
        b_next = b_idx + 1
        coh_a = band_cohs[b_idx]
        coh_b = band_cohs[b_next]
        for ch in range(n_ch):
            # Cross-band coupling weight for channel ch:
            # mean coherence of channel ch with all others in each band
            w_a = np.mean(np.delete(coh_a[ch], ch))
            w_b = np.mean(np.delete(coh_b[ch], ch))
            cross_weight = float(np.sqrt(w_a * w_b))
            # Node indices
            node_a = b_idx * n_ch + ch
            node_b = b_next * n_ch + ch
            C[node_a, node_b] = cross_weight
            C[node_b, node_a] = cross_weight

    # Zero diagonal
    np.fill_diagonal(C, 0.0)
    return C, n_total


# ============================================================
# FIM pipeline (identical logic to Phase 3E, operates on
# the multiplex coherence matrix instead of raw correlation)
# ============================================================
def build_knn(coh_matrix, k, n):
    """Build k-NN adjacency from coherence matrix."""
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        row = coh_matrix[i].copy()
        row[i] = -np.inf
        top_k = np.argsort(row)[-k:]
        for j in top_k:
            adj[i, j] = True
            adj[j, i] = True
    return adj


def compute_cr(coh_matrix, adj, n):
    """Compute C(r): mean coherence at each graph distance."""
    adj_sp = csr_matrix(adj.astype(float))
    samples = np.random.choice(n, min(N_CR_SAMPLES, n), replace=False)
    dist = shortest_path(adj_sp, method="D", indices=samples,
                         directed=False)
    finite_dists = dist[np.isfinite(dist) & (dist > 0)]
    if len(finite_dists) == 0:
        return None, 0
    max_r = int(np.max(finite_dists))
    max_r = min(max_r, 15)
    if max_r < 2:
        return None, max_r
    C_r = np.zeros(max_r + 1)
    counts = np.zeros(max_r + 1)
    for idx_i, src in enumerate(samples):
        for dst in range(n):
            if src == dst:
                continue
            d = dist[idx_i, dst]
            if np.isfinite(d) and 1 <= int(d) <= max_r:
                C_r[int(d)] += coh_matrix[src, dst]
                counts[int(d)] += 1
    mask = counts > 0
    C_r[mask] /= counts[mask]
    return C_r, max_r


def kernel(C_r, distances, n):
    """Probability kernel from C(r)."""
    weights = np.zeros(n)
    for u in range(n):
        d = distances[u]
        if np.isfinite(d) and 1 <= int(d) < len(C_r):
            weights[u] = abs(C_r[int(d)])
    total = np.sum(weights)
    if total < 1e-15:
        return None
    return weights / total


def fisher_diagnostics(C_r, adj, n, coh_matrix):
    """
    Compute FIM diagnostics from multiplex coherence graph.

    Additional output: dominant_band — which band's nodes contribute
    most variance to the FIM score vectors.
    """
    adj_sp = csr_matrix(adj.astype(float))
    neighbors = {i: list(np.where(adj[i])[0]) for i in range(n)}
    samples = np.random.choice(n, min(N_FISHER_SAMPLES, n), replace=False)
    ranks, etas, prs, sv2sv1s, n_cr_bins_list = [], [], [], [], []
    band_variances = {b: [] for b in BAND_NAMES}

    n_ch_per_band = n // N_BANDS

    for v0 in samples:
        nbrs = neighbors.get(v0, [])
        if len(nbrs) < 2:
            continue
        sources = [v0] + list(nbrs)
        try:
            dists = shortest_path(adj_sp, method="D",
                                  indices=sources, directed=False)
        except Exception:
            continue
        p_v0 = kernel(C_r, dists[0], n)
        if p_v0 is None:
            continue
        k_nbrs = len(nbrs)
        scores = np.zeros((k_nbrs, n))
        valid = True
        for j in range(k_nbrs):
            p_wj = kernel(C_r, dists[j + 1], n)
            if p_wj is None:
                valid = False
                break
            with np.errstate(divide="ignore", invalid="ignore"):
                lr = np.log(p_wj + 1e-15) - np.log(p_v0 + 1e-15)
                lr = np.nan_to_num(lr, nan=0.0, posinf=0.0, neginf=0.0)
            scores[j] = lr
        if not valid:
            continue
        weighted = scores * np.sqrt(p_v0)[np.newaxis, :]
        F = weighted @ weighted.T
        sv = np.linalg.svd(F, compute_uv=False)
        if sv[0] < 1e-15:
            continue
        sv_norm = sv / sv[0]
        gaps = sv_norm[:-1] / (sv_norm[1:] + 1e-15)
        rank = int(np.argmax(gaps) + 1)
        eta = float(sv_norm[min(rank, len(sv_norm)-1)] /
                    (sv_norm[rank-1] + 1e-15))
        pr = float(np.sum(sv_norm)**2 / (np.sum(sv_norm**2) + 1e-15))
        s2s1 = float(sv_norm[1]) if len(sv_norm) > 1 else 0.0
        ranks.append(rank)
        etas.append(eta)
        prs.append(pr)
        sv2sv1s.append(s2s1)
        n_cr_bins_list.append(np.sum(C_r > 0))

        # Per-band score variance for P3E-B-6
        for b_idx, band_name in enumerate(BAND_NAMES):
            b_start = b_idx * n_ch_per_band
            b_end = b_start + n_ch_per_band
            band_var = float(np.var(weighted[:, b_start:b_end]))
            band_variances[band_name].append(band_var)

    if len(ranks) == 0:
        return None

    # Dominant band
    mean_vars = {b: np.mean(v) if v else 0.0
                 for b, v in band_variances.items()}
    dominant = max(mean_vars, key=mean_vars.get)

    return {
        "sv2sv1":       float(np.mean(sv2sv1s)),
        "rank":         float(np.mean(ranks)),
        "eta":          float(np.mean(etas)),
        "pr":           float(np.mean(prs)),
        "n_valid":      len(ranks),
        "n_cr_bins":    float(np.mean(n_cr_bins_list)),
        "dominant_band": dominant,
    }


# ============================================================
# Main seizure processing loop
# ============================================================
def process_seizure(patient, seizure_info, patient_dir):
    """Process one seizure: sliding window multiplex FIM."""
    fname = seizure_info["file"]
    sz_start = seizure_info["start"]
    sz_end = seizure_info["end"]

    load_start = max(0, sz_start - PRE_ICTAL_SEC)
    load_end = sz_end + POST_ICTAL_SEC
    load_dur = load_end - load_start

    signals, ch_names = load_eeg_segment(
        patient_dir, fname, load_start, load_dur)
    if signals is None:
        return []

    signals, ch_names = select_eeg_channels(signals, ch_names)
    n_ch = signals.shape[0]
    if n_ch < 8:
        print(f"    Skipping {fname}: only {n_ch} channels")
        return []

    # Global bandpass
    signals = bandpass_filter(
        signals, BANDPASS_GLOBAL[0], BANDPASS_GLOBAL[1], SFREQ)

    results = []
    win_start = 0
    n_samp = signals.shape[1]

    while win_start + WIN_SAMPLES <= n_samp:
        window = signals[:, win_start:win_start + WIN_SAMPLES]

        abs_time = load_start + win_start / SFREQ
        rel_time = abs_time - sz_start

        # Epoch
        if abs_time < sz_start:
            epoch = "pre_ictal" if abs_time >= sz_start - PRE_ICTAL_SEC \
                    else "interictal"
        elif abs_time <= sz_end:
            epoch = "ictal"
        else:
            epoch = "post_ictal"

        # Build multiplex coherence matrix
        try:
            C_multi, n_total = build_multiplex_coherence(window, SFREQ)
        except Exception as e:
            win_start += STEP_SAMPLES
            continue

        # Build k-NN graph on multiplex
        adj = build_knn(C_multi, min(K_NN, n_total - 1), n_total)

        # C(r) from multiplex
        C_r, n_bins = compute_cr(C_multi, adj, n_total)

        if C_r is None or n_bins < 2:
            win_start += STEP_SAMPLES
            continue

        # FIM diagnostics
        diag = fisher_diagnostics(C_r, adj, n_total, C_multi)

        if diag:
            results.append({
                "patient":       patient,
                "seizure_file":  fname,
                "time_sec":      round(abs_time, 1),
                "rel_onset_sec": round(rel_time, 1),
                "epoch":         epoch,
                "sv2sv1":        round(diag["sv2sv1"], 4),
                "rank":          round(diag["rank"], 2),
                "eta":           round(diag["eta"], 4),
                "pr":            round(diag["pr"], 3),
                "n_valid":       diag["n_valid"],
                "n_channels":    n_ch,
                "n_cr_bins":     round(diag["n_cr_bins"], 1),
                "dominant_band": diag["dominant_band"],
            })

        win_start += STEP_SAMPLES

    return results


# ============================================================
# Summary and output
# ============================================================
def write_summary(all_results, summary_path):
    """Per-patient epoch summary with resolution gate check."""
    from collections import defaultdict
    groups = defaultdict(lambda: defaultdict(list))
    for r in all_results:
        groups[r["patient"]][r["epoch"]].append(r)

    with open(summary_path, "w") as f:
        f.write("DS Phase 3E-B: EEG Multiplex FIM Results Summary\n")
        f.write("=" * 65 + "\n\n")

        # Resolution gate P3E-B-0
        f.write("P3E-B-0 RESOLUTION GATE: C(r) bins >= 5 in >= 80% windows\n")
        f.write("-" * 65 + "\n")
        gate_pass_count = 0
        for patient in sorted(groups.keys()):
            all_bins = []
            for epoch_rows in groups[patient].values():
                all_bins.extend([r["n_cr_bins"] for r in epoch_rows])
            if all_bins:
                pct_above = np.mean(np.array(all_bins) >= MIN_CR_BINS)
                gate_ok = pct_above >= 0.80
                if gate_ok:
                    gate_pass_count += 1
                f.write(f"  {patient}: mean_bins={np.mean(all_bins):.1f}, "
                        f"pct_above_5={pct_above:.1%} → "
                        f"{'PASS' if gate_ok else 'FAIL'}\n")
        f.write(f"\nResolution gate: {gate_pass_count}/5 patients PASS\n")
        if gate_pass_count < 4:
            f.write("WARNING: Resolution fix insufficient — "
                    "diameter still too small\n")
        f.write("\n")

        # Per-patient epoch summary
        for patient in sorted(groups.keys()):
            f.write(f"\n--- {patient} ---\n")
            f.write(f"{'Epoch':<15} {'N':>5} {'SV2/SV1':>10} "
                    f"{'Rank':>8} {'Eta':>8} {'PR':>8} {'Bins':>6}\n")
            for epoch in ["pre_ictal", "ictal", "post_ictal"]:
                rows = groups[patient].get(epoch, [])
                if rows:
                    sv  = [r["sv2sv1"] for r in rows]
                    rk  = [r["rank"]   for r in rows]
                    et  = [r["eta"]    for r in rows]
                    pr  = [r["pr"]     for r in rows]
                    bn  = [r["n_cr_bins"] for r in rows]
                    f.write(f"{epoch:<15} {len(sv):>5} "
                            f"{np.mean(sv):>10.4f} {np.mean(rk):>8.2f} "
                            f"{np.mean(et):>8.4f} {np.mean(pr):>8.3f} "
                            f"{np.mean(bn):>6.1f}\n")
            f.write("\n")

        # P3E-B-1 kill test
        f.write("=" * 65 + "\n")
        f.write("P3E-B-1 KILL TEST: SV2/SV1 ictal > pre-ictal by >= 0.05\n")
        f.write("-" * 65 + "\n")
        pass_count = 0
        for patient in sorted(groups.keys()):
            pre = groups[patient].get("pre_ictal", [])
            ict = groups[patient].get("ictal", [])
            if pre and ict:
                pre_m = np.mean([r["sv2sv1"] for r in pre])
                ict_m = np.mean([r["sv2sv1"] for r in ict])
                delta = ict_m - pre_m
                ok = delta >= 0.05
                if ok:
                    pass_count += 1
                f.write(f"  {patient}: pre={pre_m:.4f} ictal={ict_m:.4f} "
                        f"delta={delta:+.4f} → {'PASS' if ok else 'FAIL'}\n")
        f.write(f"\nP3E-B-1: {pass_count}/5 patients pass. "
                f"{'KILL TEST PASSES' if pass_count >= 3 else 'KILL TEST FAILS'}\n")

        # P3E-B-2 through P3E-B-4
        f.write("\nP3E-B-2: SV2/SV1 post-ictal < ictal\n")
        for patient in sorted(groups.keys()):
            ict  = groups[patient].get("ictal", [])
            post = groups[patient].get("post_ictal", [])
            if ict and post:
                im = np.mean([r["sv2sv1"] for r in ict])
                pm = np.mean([r["sv2sv1"] for r in post])
                f.write(f"  {patient}: ictal={im:.4f} post={pm:.4f} "
                        f"→ {'PASS' if pm < im else 'FAIL'}\n")

        f.write("\nP3E-B-3: Rank decreases during ictal\n")
        for patient in sorted(groups.keys()):
            pre = groups[patient].get("pre_ictal", [])
            ict = groups[patient].get("ictal", [])
            if pre and ict:
                pr_r = np.mean([r["rank"] for r in pre])
                ic_r = np.mean([r["rank"] for r in ict])
                f.write(f"  {patient}: pre={pr_r:.2f} ictal={ic_r:.2f} "
                        f"→ {'PASS' if ic_r < pr_r else 'FAIL'}\n")

        f.write("\nP3E-B-4: Eta peaks during/after ictal\n")
        for patient in sorted(groups.keys()):
            pre  = groups[patient].get("pre_ictal", [])
            ict  = groups[patient].get("ictal", [])
            post = groups[patient].get("post_ictal", [])
            if pre and ict:
                pe  = np.mean([r["eta"] for r in pre])
                ie  = np.mean([r["eta"] for r in ict])
                poe = np.mean([r["eta"] for r in post]) if post else 0.0
                peak = max(ie, poe)
                f.write(f"  {patient}: pre={pe:.4f} ictal={ie:.4f} "
                        f"post={poe:.4f} → {'PASS' if peak > pe else 'FAIL'}\n")

        f.write("\nP3E-B-6: Dominant band during ictal vs pre-ictal\n")
        for patient in sorted(groups.keys()):
            pre = groups[patient].get("pre_ictal", [])
            ict = groups[patient].get("ictal", [])
            if pre and ict:
                from collections import Counter
                pre_dom = Counter([r["dominant_band"] for r in pre]).most_common(1)
                ict_dom = Counter([r["dominant_band"] for r in ict]).most_common(1)
                f.write(f"  {patient}: pre_dominant={pre_dom[0][0] if pre_dom else 'n/a'} "
                        f"ictal_dominant={ict_dom[0][0] if ict_dom else 'n/a'}\n")


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
            continue
        for si, sz in enumerate(seizures):
            print(f"  Seizure {si+1}/{len(seizures)}: {sz['file']} "
                  f"[{sz['start']}s–{sz['end']}s]")
            results = process_seizure(patient, sz, patient_dir)
            all_results.extend(results)
            print(f"    → {len(results)} windows")

    # Write CSV
    csv_path = os.path.join(OUTPUT_DIR, "phase3e_b_results.csv")
    if all_results:
        fields = ["patient", "seizure_file", "time_sec", "rel_onset_sec",
                  "epoch", "sv2sv1", "rank", "eta", "pr", "n_valid",
                  "n_channels", "n_cr_bins", "dominant_band"]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for row in all_results:
                writer.writerow(row)
        print(f"\nWrote {len(all_results)} rows → {csv_path}")

    # Write summary
    summary_path = os.path.join(OUTPUT_DIR, "phase3e_b_summary.txt")
    write_summary(all_results, summary_path)
    print(f"Wrote summary → {summary_path}")

    elapsed = time.time() - t0
    print(f"\nRuntime: {elapsed:.1f}s")
    print("Upload phase3e_b_results.csv and phase3e_b_summary.txt to chat.")


if __name__ == "__main__":
    main()
```

---

## Validation Gates

- **Gate A:** ≥8 EEG channels after filtering → else skip file
- **Gate B:** C(r) bins ≥2 → else skip window (logged in n_cr_bins)
- **Gate C:** ≥1 valid FIM → else skip window
- **Gate P3E-B-0:** C(r) bins ≥5 in ≥80% windows for ≥4/5 patients
  → if this fails, report and STOP (resolution fix failed)

## Runtime Note

Coherence computation (scipy.signal.csd × 4 bands × 23² pairs × N
windows) is slower than raw correlation. Expected runtime:
- 15-25 min for chb01–chb05 on M4 depending on seizure file count
- If too slow: reduce N_CR_SAMPLES from 20 to 10

## Handback

Upload `phase3e_b_results.csv` and `phase3e_b_summary.txt` to chat.
Analysis will:
1. Verify P3E-B-0 resolution gate passed
2. Evaluate P3E-B-1 kill condition
3. Check P3E-B-2 through P3E-B-6
4. If pass: plot time-resolved SV₂/SV₁ around seizure onsets
5. Report verdict: PASS / KILL / INCONCLUSIVE
