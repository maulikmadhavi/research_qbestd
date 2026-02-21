# QbE-STD: Query by Example Spoken Term Detection

## Overview

**Query-by-Example Spoken Term Detection (QbE-STD)** is a speech processing task where you search for a query phrase (by voice) within a large audio database to find similar spoken occurrences, without requiring a manual transcript. This repository contains research implementations spanning a decade of work, backed by a PhD thesis and clean Python reimplementations for experimentation.

The repository is split into two language implementations:
- **MATLAB**: High-performance DTW kernels compiled as MEX (C++) for efficient distance computation
- **Python**: Clean reimplementation and algorithmic experimentation with visualization tools

---

## Repository Layout

```
qbe-std/
├── README.md                          # This file – project overview and navigation
├── matlab/                            # MATLAB + C MEX implementations of DTW kernels
│   ├── Fx_do_SDTW.m                   # Wrapper for segmental DTW
│   ├── DTW_c_skel_nobt.cpp            # Core DTW DP recurrence (C++ MEX)
│   └── [other NSDTW/GTTS variants]    # Additional DTW implementations
├── python/                            # Python experiments & algorithms
│   ├── readme.md                      # Python-specific setup guide
│   ├── WARP.md                        # Development notes for WARP IDE
│   ├── subsequence_dtw.py             # End-to-end DTW on real audio (MFCCs)
│   ├── test_subsequence_dtw.py        # Synthetic test harness for DTW validation
│   ├── k-means.py                     # Vector quantization + clustering demo
│   ├── pixi.toml                      # Pixi environment & dependencies
│   ├── utils/                         # Helper modules
│   │   ├── gen_data.py                # Synthetic data generation
│   │   └── ploting.py                 # Plotting utilities
│   └── data/                          # Audio input files (query + reference WAVs)
└── .pixi/                             # Pixi environment cache
```

---

## MATLAB Implementation

The MATLAB module contains high-performance DTW kernels compiled as C++ MEX functions, developed as part of the PhD research work. These are optimized for production-scale distance matrix computation.

### Purpose
- Fast, production-grade dynamic time warping (DTW) for spoken term detection
- Supports multiple DTW variants: basic, nonsegmental, global time-scale stretching
- Massively optimized C++ implementations with MEX interface

### Compilation
To build MEX files:
```matlab
% In the matlab/ directory
mex DTW_c_skel_nobt.cpp
mex NSDTW_c_skel.cpp
% ... and other .cpp files as needed
```

### DTW Kernel Variants

| Kernel | Description |
|--------|-------------|
| `DTW_c_skel_nobt` | Standard DTW with step constraints (no backtracking) |
| `DTW_c_basic_skel_nobt` | Simplified DTW for full-sequence alignment |
| `NSDTW_c_skel` | Nonsegmental DTW (variants: `_2`, `_4`, `_5`) |
| `NSDTW_c_skel_online` | Online variant with incremental computation |
| `newNSDTW_c_skel` | Improved nonsegmental DTW |
| `newNSDTW_c_skel_online` | Online variant of improved NSDTW |
| `GTTS_DTW_c_skel` | Global Time-Scale Stretching DTW |
| `GTTS_DTW_c_skel_online` | Online variant of GTTS |
| `sub_DTW_c_skel_online` | Subsequence DTW (online) |

### Entry Point
**`Fx_do_SDTW.m`** is the main callable wrapper for **Segmental DTW**:
```matlab
[dist, startpos, endpos, DistMtrx] = Fx_do_SDTW(refcoef, qrycoef, Type_localdist, Warping_adjustment)
```
- **Input**: Reference features (ND × N_ref), query features (ND × N_query)
- **Output**: DTW distance, start/end positions in reference, full distance matrix
- **Supports**: Euclidean ('s'), inner product ('i'), KL divergence ('k'), Bhattacharyya distance ('b')

For details, see [`matlab/README.md`](matlab/README.md).

---

## Python Implementation

A clean reimplementation of core QbE-STD algorithms in Python for experimentation, visualization, and educational purposes.

### Purpose
- Algorithmic prototyping and testing before optimization
- Interactive visualization of DTW alignment paths
- Integration with modern Python libraries (NumPy, Librosa, Matplotlib)

### Environment Setup
Install [Pixi](https://pixi.sh/latest/installation/), then:
```bash
cd python
pixi shell
```
All commands below assume you are inside a `pixi shell`.

**Dependencies** (from `pixi.toml`):
- Python 3.12
- NumPy ≥ 2.2
- Matplotlib ≥ 3.10
- Librosa ≥ 0.11
- Einops ≥ 0.8

### Main Scripts

#### `subsequence_dtw.py` – End-to-End DTW on Real Audio
Runs subsequence DTW on MFCC features extracted from WAV files.

**Pipeline**:
1. Extract 40-dimensional MFCC features from query and reference audio (normalized: zero mean, unit variance)
2. Compute pairwise Euclidean distance matrix between query and reference frames
3. Build accumulated cost matrix `S`, auxiliary length matrix `T`, and 3D backpointer tensor `BP`
4. Backtrack from the best-matching end position to recover the optimal alignment path
5. Visualize and save result to `S_matrix.png`

**Usage**:
```bash
python subsequence_dtw.py
```

**Data**: Expects audio files under `data/`:
- Query: `data/d124cd78-78ce-4603-8ffa-673f35887e61.wav` (or modify hard-coded filename)
- Reference: `data/reference.wav`

**Output**: `S_matrix.png` showing the accumulated cost matrix with the recovered alignment path overlaid.

#### `test_subsequence_dtw.py` – Synthetic DTW Test Harness
Self-contained validation script for the subsequence DTW logic.

**Purpose**: Verify that DP recurrence and backtracking produce valid monotonic alignment paths before deploying to real audio.

**Usage**:
```bash
python test_subsequence_dtw.py
```

**Output**: `test_S_matrix.png` with detailed step-by-step path validation printed to console.

#### `k-means.py` – Vector Quantization + Clustering
Demonstration of k-means clustering on synthetic 2D data.

**Pipeline**:
1. Generate simulated clustered data (4 clusters from 2D Gaussian mixtures)
2. Run iterative k-means with custom convergence and cluster-splitting logic
3. Compute final distances to cluster centers and assign each point
4. Compare assignments to ground-truth labels

**Usage**:
```bash
python k-means.py
```

### Utility Modules

- **`utils/gen_data.py`**: Signal generation helpers
  - `get_20hz_sine_wave`, `get_10hz_sine_wave` – noisy sine waves
  - `generated_simulated_data` – synthetic clustered 2D data for k-means
  - Run directly: `python utils/gen_data.py` for quick visualization

- **`utils/ploting.py`**: Generic plotting helper
  - `plot_wave(t, wave, freq)` – time-domain signal visualization
  - Run directly: `python utils/ploting.py` for sample waveforms

### Linting
```bash
ruff .
```
(Configured via `ruff.toml` to exclude `data/` and `.pixi/` directories)

### Output Artifacts
- `S_matrix.png` – DTW cost matrix with alignment path from `subsequence_dtw.py`
- `test_S_matrix.png` – Synthetic test cost matrix from `test_subsequence_dtw.py`
- Cluster visualizations from `k-means.py`

For detailed development notes, see [`python/readme.md`](python/readme.md) and [`python/WARP.md`](python/WARP.md).

---

## Algorithm Overview: QbE-STD Pipeline

Query-by-Example Spoken Term Detection follows this general pipeline:

1. **Feature Extraction**
   - Extract acoustic features (MFCCs, spectrograms, etc.) from query and reference audio
   - Normalize features (zero mean, unit variance)

2. **Distance Computation**
   - Compute pairwise local distance matrix between query and reference frames (Euclidean, inner product, or KL divergence)

3. **Dynamic Programming (DP)**
   - Build accumulated cost matrix via standard DTW recurrence with step constraints
   - Track start/end positions and backpointer information

4. **Backtracking**
   - Identify the minimum-cost end position in the last row
   - Walk backwards through backpointers to recover the optimal alignment path

5. **Detection**
   - The aligned path spans a contiguous segment of the reference audio
   - Cost score indicates confidence of match

**Variants** (implemented in MATLAB):
- **Nonsegmental DTW (NSDTW)**: Constrain warping window to enforce time-scale limits
- **Segmental DTW**: Search all possible segments of reference for best match
- **Online variants**: Incremental computation as reference grows
- **Global Time-Scale Stretching (GTTS)**: Allow controlled global tempo variations

---

## Reference

For the full research context and algorithmic details, see the PhD thesis:

**Thesis**: *Query-by-Example Spoken Term Detection* (PhD work)
**Link**: [View PDF](http://drsr.daiict.ac.in/jspui/bitstream/123456789/649/1/201121003.pdf)

---

## Getting Started

### Quick Start – Python
```bash
cd python
pixi shell
python subsequence_dtw.py  # Requires audio in data/
python test_subsequence_dtw.py  # Self-contained test
```

### Quick Start – MATLAB
```matlab
cd matlab
mex DTW_c_skel_nobt.cpp  % Compile kernels
[dist, sp, ep, ~] = Fx_do_SDTW(ref_features, query_features, 's', 0);
```

---

## License & Attribution

This project contains research code from a PhD thesis. Both the MATLAB kernels and Python reimplementation are part of the same scientific work.

---

**Repository**: `/mnt/d/speech_research/qbe-std`
**Languages**: MATLAB, Python, C++
**Status**: Research & Experimentation
