# Example data

Two datasets, shared by all three implementations. Every example script and both Python
notebooks point here by relative path, so nothing needs configuring to run the demos.

## Traffic

| File | What it is |
|---|---|
| `Traffic_example_for_DSSIM.avi` | 173 frames, 1280 x 720, greyscale traffic-camera footage |
| `time_data.csv` | One column of frame times: 173 rows, 0.0 to 68.8 s at 0.4 s spacing |

Traffic footage is used rather than microscopy because the structural change is obvious
by eye — you can check the DSSIM map against what you can see moving. The moving
vehicles and pedestrians light up; the road, markings and buildings stay dark. Run this
one first.

Used by `DSSIM_Analysis_notebook.ipynb` and `dssim_analysis_example_script.m`.

## Liquid-cell TEM

| File | What it is |
|---|---|
| `LC_TEM_Tifs/` | 30 frames, 512 x 512, 32-bit float `.tif` files (`0001.tif` … `0030.tif`) |

The microscopy case the method was built for: noisy, low-dose data where the change is
*not* obvious by eye. It also exercises the other input paths — a folder of TIFFs rather
than a video file, and no time file, so frames are spaced 1 second apart.

Used by `DSSIM_Analysis_notebook_LC_TEM_example.ipynb`.

The two need different parameters — σ 1 / offset 1 / radius 3 / `(1,1,1)` for traffic
versus 3 / 3 / 1.5 / `(1,0,0)` here — which is the clearest illustration in the repo that
these values follow from the dataset rather than being universal defaults.

## A note on `time_data.csv`

This file originally contained 346 rows at 0.2 s spacing while the AVI decodes to 173
frames — exactly twice as many. Both the MATLAB and Python implementations detect the
mismatch and fall back to 1 second per frame, so the bundled example silently ran with
the wrong time axis.

The AVI appears to be a 2x decimation of whatever produced the original CSV, so the file
has been decimated to match: every other row kept, giving 173 rows at 0.4 s spacing.

If you regenerate this file, the row count must equal the frame count. Both
implementations warn when it doesn't; the Python one names the two numbers.

## Using your own data

Don't overwrite these files — point the scripts at your own data instead.

- **Python**: set `DATA_PATH` and `TIMES` in the notebook's Settings cell
- **MATLAB**: set `paras_dssim.data` and `paras_dssim.times` at the top of
  `dssim_analysis_example_script.m`

Both accept a video file or a folder of TIFFs; the Python version also accepts a
multi-page TIFF stack or a NumPy array.
