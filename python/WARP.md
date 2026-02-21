# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Environment and toolchain
- Python project managed via [`pixi`](https://pixi.sh/latest/installation/).
- Core dependencies are defined in `pixi.toml` (Python 3.12, `numpy`, `matplotlib`, `librosa`, `einops`).
- Lint configuration is in `ruff.toml` (excludes `data/*` and `.pixi/*`).

To enter the project environment:
- From the repo root (`qbe-std`):
  - `pixi shell`

All Python commands below assume you are inside a `pixi shell`.

## Common commands

### Run experiments / scripts
From the repo root:
- Vector quantization + k-means demo:
  - `python k-means.py`
- Subsequence DTW on real audio (expects audio files under `data/`):
  - `python subsequence_dtw.py`
- Synthetic DTW validation script (acts as a test harness for the DP + backtracking logic):
  - `python test_subsequence_dtw.py`
- Utility data-generation or plotting modules (when run directly):
  - `python utils/gen_data.py`
  - `python utils/ploting.py`

### Linting
If `ruff` is installed in your environment, you can lint the codebase with:
- `ruff .`

This uses the settings in `ruff.toml`, which already exclude large/generated content under `data/` and environment files under `.pixi/`.

## High-level architecture

### Top-level scripts
- `k-means.py`
  - Demonstration script for vector quantization followed by k-means clustering on synthetic 2D data.
  - Pulls simulated data from `utils/gen_data.generated_simulated_data` and performs iterative k-means with custom convergence logic and cluster splitting until a requested number of clusters is reached.
  - At the end, computes final distances to cluster centers, assigns each point, and compares assignments to expected cluster labels from the generator.

- `subsequence_dtw.py`
  - End-to-end implementation of subsequence DTW on real audio.
  - Pipeline:
    1. **Feature extraction**: `extract_mfcc` uses `librosa` to load audio, compute 40-dim MFCCs, and normalize them (zero mean, unit variance across time).
    2. **Distance computation**: `compute_distances_einsum` computes a full pairwise Euclidean distance matrix between query and reference MFCC frames using vectorized `numpy.einsum`.
    3. **DP matrices**: builds accumulated cost matrix `S`, auxiliary matrices `T` (length) and `P` (start positions), plus a 3D backpointer tensor `BP` that tracks predecessors for each DP cell.
    4. **Backtracking**: finds the best-matching end column in the last row of `S`, then walks backwards through `BP` to recover the optimal subsequence alignment path.
    5. **Visualization**: plots `S` with the recovered path overlaid and saves `S_matrix.png`.
  - Assumes the presence of audio files under `data/` (e.g. `data/d124cd78-78ce-4603-8ffa-673f35887e61.wav` and `data/reference.wav`).

- `test_subsequence_dtw.py`
  - Self-contained synthetic testbed for the same subsequence DTW core logic used in `subsequence_dtw.py`.
  - Generates random query/reference feature sequences, computes the distance matrix with the same `compute_distances_einsum` routine, and reimplements the DP and backtracking loop.
  - Prints detailed path information (step-by-step validation of allowed moves) and visualizes the accumulated cost matrix with the path overlay, saving `test_S_matrix.png`.
  - Useful as a quick way to verify that algorithmic changes to the DTW logic produce valid monotonic paths before running on real audio.

### Utility modules
- `utils/gen_data.py`
  - Provides signal-generation helpers:
    - `get_20hz_sine_wave` and `get_10hz_sine_wave` for noisy sine waves at fixed frequencies (used for simple signal experiments/plots).
    - `generated_simulated_data` to construct synthetic clustered 2D data used by `k-means.py`:
      - Stacks several 1D Gaussian sequences into a 2×N feature matrix and produces corresponding integer labels for four clusters.
  - When run as a script, plots the generated 2D data for quick visual sanity checks.

- `utils/ploting.py`
  - Contains generic plotting helper `plot_wave(t, wave, freq)` for visualizing time-domain signals.
  - When run as a script, generates two sample waveforms using `get_20hz_sine_wave` and `get_10hz_sine_wave` and displays them.

### Data and outputs
- `data/`
  - Expected to contain input audio for `subsequence_dtw.py` (query and reference WAV files). The exact filenames currently referenced are hard-coded in `subsequence_dtw.py`.
- Generated artifacts:
  - `S_matrix.png` from `subsequence_dtw.py`.
  - `test_S_matrix.png` from `test_subsequence_dtw.py`.

These visualizations are part of the intended workflow for inspecting DTW behavior and alignment paths.
