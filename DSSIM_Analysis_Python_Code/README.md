# DSSIM Analysis — Python Code

Python implementation of structural dissimilarity (DSSIM) analysis, as described in
[Mulvey et al., *Ultramicroscopy* **257** (2024) 113894](https://doi.org/10.1016/j.ultramic.2023.113894).

DSSIM turns a video into a per-pixel map of where structure changed between frames:

```
DSSIM = (1 - SSIM) / 2
```

DSSIM highlights movement as well as structural change. Drift correcting the data is advantageous towards isolating structural change with DSSIM.

## Files

| File | What it is |
|---|---|
| `DSSIM_Analysis_notebook.ipynb` | **Start here.** The analysis and its documentation in one place. Set your parameters in a single cell, run it, and collect the videos and CSVs it writes |
| `dssim.py` | Every function the notebook calls. Import it to run DSSIM from your own script or pipeline, with no notebook involved |
| `test_dssim.py` | Confirms the install works and the maths is right — cross-checks the SSIM calculation against scikit-image's reference and exercises the edge cases. Worth running once after installing |
| `requirements.txt` | The packages needed and their minimum versions |

Example data lives in [`../example_data/`](../example_data/), shared with the other
implementations.

## Getting started

```bash
pip install -r requirements.txt
jupyter notebook DSSIM_Analysis_notebook.ipynb
```

It runs on the bundled example out of the box. Run Jupyter from this folder so the relative
paths resolve, then point `DATA_PATH` at your own data.

**The notebook is the documentation.** It explains every parameter, how to choose it for your
dataset, what each output file contains, ect.

To verify the install: `python test_dssim.py`

## Scripting it instead

```python
import dssim

data = dssim.load_data("my_movie.avi")          # video, TIFF stack, TIFF folder, or array
data = dssim.gaussian_denoise(data, sigma=1)

result = dssim.dssim_analysis(data, times="times.csv", frame_offset=1, radius=3)

result.dssim     # (frames, h, w) float32 - the quantitative output
result.means     # mean DSSIM per frame - the change-over-time curve
result.stats     # DataFrame of frame numbers, times, and mean values

dssim.save_results(result, run_name="run1", output_dir="out", data_frames=data)
```

Every function has a docstring — `help(dssim.dssim_analysis)`.

Frame numbers are 0-based everywhere, including the exported CSV, so any number in the output
indexes straight back into your data.

## Speed

The frame loop is single-threaded. The bundled traffic example — 173 frames of 1280 x 720,
638 MB once loaded as float32 — takes roughly half a minute at `RADIUS = 3`, and longer at
larger radii.


## Contact

Justin T. Mulvey — jtmulvey1@gmail.com

Additional DSSIM examples can be viewed at [justintmulvey.com](https://justintmulvey.com).
