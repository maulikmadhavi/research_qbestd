# QbE-STD experiments

This repo contains small experiments around Query-by-Example Spoken Term Detection (QbE-STD), reimplementing parts of the PhD work referenced in the thesis: [link](http://drsr.daiict.ac.in/jspui/bitstream/123456789/649/1/201121003.pdf).

## Environment
- Install [pixi](https://pixi.sh/latest/installation/).
- From the repo root, start the environment with:
  - `pixi shell`

All commands below assume you are inside a `pixi shell`.

## Main scripts
- `k-means.py` – vector quantization + k-means clustering on synthetic 2D data generated in `utils/gen_data.py`.
- `subsequence_dtw.py` – subsequence DTW on MFCCs from real audio in `data/`; saves `S_matrix.png` with the alignment path.
- `test_subsequence_dtw.py` – synthetic DTW testbed validating the DP + backtracking path; saves `test_S_matrix.png`.

## Data
- Place query and reference WAV files under `data/`.
- Filenames currently expected by `subsequence_dtw.py` are hard-coded (e.g. `data/d124cd78-78ce-4603-8ffa-673f35887e61.wav`, `data/reference.wav`).
