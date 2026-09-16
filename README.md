# DSSIM Analysis

ImageJ, MATLAB and Python implementations of structural dissimilarity (DSSIM) analysis,
as described in Mulvey et al., *Ultramicroscopy* **257** (2024) 113894:
https://doi.org/10.1016/j.ultramic.2023.113894

DSSIM turns a video into a per-pixel map of **where structure changed** between frames:

```
DSSIM = (1 - SSIM) / 2
```

DSSIM highlights movement as well as structural change. Drift correcting the data is advantageous towards isolating structural change with DSSIM.

Highlights:

- Structural dissimilarity analysis provides structural dynamics maps from videos.
- It is a segmentation-free method, but gives results comparable to segmentation analysis.
- It is simple and computationally efficient.

https://github.com/JustinTMulvey/DSSIM_Analysis/assets/141819477/ab4b17b4-8638-4f88-bbce-10bd322c4a05

## Repository structure

```
DSSIM_Analysis/
├── example_data/                    Two example datasets, shared by all three
│   ├── Traffic_example_for_DSSIM.avi + time_data.csv
│   └── LC_TEM_Tifs/
├── DSSIM_Analysis_Python_Code/      One notebook per example + dssim.py module
├── DSSIM_Analysis_Matlab_Code/      Function + example script
└── DSSIM_Analysis_ImageJ_Plugin/    Java plugin + PDF user guide
```

## The two examples

Both are in [`example_data/`](example_data/) and both run out of the box.

| | Data | Notebook |
|---|---|---|
| **Traffic** | 173 frames, 1280 × 720 AVI, with a matching time file | `DSSIM_Analysis_notebook.ipynb` |
| **Liquid-cell TEM** | 30 frames, 512 × 512, a folder of `.tif` files | `DSSIM_Analysis_notebook_LC_TEM_example.ipynb` |

**Start with the traffic example.** The structural change is obvious by eye, so you can
check the DSSIM map against what you can see moving — vehicles and pedestrians light up,
the road and buildings stay dark.

**Then look at the LC-TEM example.** It is the microscopy case the method was built for,
and should be easier to adapt to real dataset.

The point of having both is that the four main parameters are dataset-dependent, not
universal. The traffic example runs at `GAUSS_BLUR_SIGMA = 1`, `FRAME_OFFSET = 1`,
`RADIUS = 3`, `EXPONENTS = (1, 1, 1)`; the LC-TEM example needs `3`, `3`, `1.5` and
`(1, 0, 0)`. Comparing the two settings cells is the fastest way to get a feel for how to
choose them for your own data.

The MATLAB script ships pointed at the traffic example. `paras_dssim.data` also accepts a
folder of TIFFs, so it can be pointed at the other dataset.

## Which implementation should I use?

| | Best for | Requires |
|---|---|---|
| **ImageJ plugin** | Small datasets (<200 MB), point-and-click use, no coding | ImageJ / Fiji |
| **Python** | Any size, scripting, integration with other analysis | Python 3.9+ |
| **MATLAB** | Large datasets if MATLAB is already your environment | MATLAB 2020b+ |

Start with the ImageJ plugin to see what DSSIM does to your data. Use Python or MATLAB when
you need to script it or tune parameters systematically.

## Quick start

**Python** — the notebooks document every parameter and output:

```bash
cd DSSIM_Analysis_Python_Code
pip install -r requirements.txt
jupyter notebook DSSIM_Analysis_notebook.ipynb
```

**MATLAB** — open `DSSIM_Analysis_Matlab_Code/dssim_analysis_example_script.m` and run it.
Modify the parameters at the top for your own data.


## Citation

If you use this software in published work, please cite:

> Mulvey, J. T.; Iyer, K. P.; Ortega, T.; Merham, J. G.; Pivak, Y.; Sun, H.;
> Hochbaum, A. I.; Patterson, J. P. "Correlating electrochemical stimulus to structural
> change in liquid electron microscopy videos using the structural dissimilarity metric."
> *Ultramicroscopy* **257** (2024) 113894. https://doi.org/10.1016/j.ultramic.2023.113894

## License

MIT for the MATLAB and Python code and the example data — see [LICENSE](LICENSE).

The ImageJ plugin is derived from third-party code restricted to educational and
research use; see [DSSIM_Analysis_ImageJ_Plugin/LICENSE](DSSIM_Analysis_ImageJ_Plugin/LICENSE).

## Found a bug?

If you found an issue or would like to submit an improvement, please email
jtmulvey1@gmail.com or open an issue on GitHub.

Additional DSSIM examples can be viewed at [justintmulvey.com](https://justintmulvey.com).
