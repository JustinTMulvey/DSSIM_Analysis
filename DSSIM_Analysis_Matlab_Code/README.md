# DSSIM Analysis — MATLAB Code

MATLAB implementation of structural dissimilarity (DSSIM) analysis, as described in
[Mulvey et al., *Ultramicroscopy* **257** (2024) 113894](https://doi.org/10.1016/j.ultramic.2023.113894).

DSSIM turns a video into a per-pixel map of where structure changed between frames:

```
DSSIM = (1 - SSIM) / 2
```

Each frame is compared with the one `frame_offset` positions later. In the DSSIM image,
dark purple regions indicate nothing changed there; bright yellow regions indicate the
local structure changed. No segmentation, tracking or thresholding required.

DSSIM highlights movement as well as structural change. Drift correcting the data is advantageous towards isolating structural change with DSSIM.

## Files

| File | What it is |
|---|---|
| `dssim_analysis_example_script.m` | **Start here.** Edit the Inputs block at the top, run it, and collect the three files it writes |
| `dssim_analysis.m` | The stand-alone analysis function. Its header documents every parameter in detail — `help dssim_analysis` |

Example data lives in [`../example_data/`](../example_data/), shared with the other
implementations.

## Running it

Open `dssim_analysis_example_script.m`, set your Current Folder to this directory so the
relative paths resolve, and run. It writes three files into `output_dir` (`dssim_output/`
by default, created if it doesn't exist):

- `<run_name>_output_DSSIM_video.avi` — original data beside the DSSIM map
- `<run_name>_DSSIM_values.csv` — per-frame times and mean DSSIM
- `<run_name>_DSSIM_run_info.csv` — every parameter used, so the run is reproducible

**The script is the documentation.** The Inputs block explains each parameter and how to
choose it for your dataset; `dssim_analysis.m`'s header covers the algorithm, the SSIM
components and the outputs. Nothing here repeats them.

I strongly recommend starting with a highly binned dataset (<1 GB) and working backward to
full resolution or uncropped data.

Written and tested in MATLAB 2020b. Requires the Image Processing Toolbox for `ssim`,
`imgaussfilt` and `mat2gray`.

## The four parameters that matter

Set in the Inputs block; each is documented in place.

| Parameter | What it does |
|---|---|
| `gauss_filt_std` | Denoising blur applied before analysis, in pixels. DSSIM reacts to *any* frame-to-frame difference and noise contributes heavily, so for low-dose data this is what separates real structural change from counting noise. `0` disables it |
| `frame_offset` | How far apart the compared frames are. `1` = consecutive, which catches the fastest changes. Larger values are more sensitive to slow events, but high values severely limit the temporal resolution (Nyquist limit) |
| `radius` | Standard deviation of the Gaussian weighting, in pixels — **not** the window size. The neighborhood is `2*ceil(3*radius) + 1` px per side, so `radius = 3` gives 19 × 19 px. Tune it to the size of the feature of interest |
| `exponents` | The DSSIM coefficients `[alpha beta gamma]`, weighting the mean, variance and normalized cross-correlation components. Typically left at `[1 1 1]` |

## Memory and speed

The whole dataset is held in RAM. This implementation is memory intensive and was
developed on a workstation with a large amount of RAM; 32 GB is a sensible minimum for
large datasets. Bin or crop first if you run out.

The DSSIM loop itself is single-threaded — the only `parfor` is in the per-frame display
contrast step, which the example script does not use by default.

## Relationship to the Python version

The two produce comparable results. Both denoise with the same Gaussian filter, normalise
the same way, and compute SSIM from the same three components.

Two differences worth knowing:

1. **Frame numbering.** This version writes 1-based frame numbers, matching MATLAB
   indexing. The Python version writes 0-based, matching Python indexing. The two CSVs are
   therefore offset by one.
2. **Neighborhood size at unusual radii.** MATLAB's `ssim` uses `ceil(3*radius)` while
   SciPy uses `round(3*radius)`. These agree at integer and half-integer radii — including
   the default of 3 — and differ by 2 px at values like 2.4 or 3.4.

## Contact

Justin T. Mulvey — jtmulvey1@gmail.com

Additional DSSIM examples can be viewed at [justintmulvey.com](https://justintmulvey.com).
