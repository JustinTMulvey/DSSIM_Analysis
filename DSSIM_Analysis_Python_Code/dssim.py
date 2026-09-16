"""
DSSIM Analysis
==============

Python implementation of structural dissimilarity (DSSIM) analysis, as described in:

    Mulvey, J. T. et al., Ultramicroscopy 257 (2024) 113894
    https://doi.org/10.1016/j.ultramic.2023.113894

DSSIM turns a video into a per-pixel map of structural *change* between frames:

    DSSIM = (1 - SSIM) / 2

Frame ``i`` is compared against frame ``i + frame_offset``, so a video of N frames
produces N - frame_offset DSSIM maps. Values near 0 mean "nothing changed here";
larger values mean the local structure changed between frames.

DSSIM highlights movement as well as structural change.
Drift correcting the data is advantageous towards isolating structural change with DSSIM.

The frame loop is single-threaded.

Typical use::

    import dssim
    stack = dssim.load_data("movie.avi")
    result = dssim.dssim_analysis(stack, frame_offset=1, radius=3)
    dssim.write_video(result.rgb_frames, "out.avi", fps=10)
    result.stats.to_csv("out_values.csv", index=False)

Author: Justin T. Mulvey (jtmulvey1@gmail.com), Patterson Group, UC Irvine
Repository: https://github.com/JustinTMulvey/DSSIM_Analysis
"""

from __future__ import annotations

import math
import os
import re
import sys
import time
from dataclasses import dataclass, field
from glob import glob
from typing import Sequence

import cv2
import matplotlib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

__all__ = [
    "load_data",
    "gaussian_denoise",
    "ProgressBar",
    "ssim_map",
    "dssim_analysis",
    "DssimResult",
    "apply_colormap",
    "to_display_rgb",
    "write_video",
    "side_by_side",
    "load_times",
    "save_results",
]

# SSIM stabilising constants from Wang et al. 2004. These should not normally
# be changed.
K1 = 0.01
K2 = 0.03

# Width of the progress bar, in characters.
PROGRESS_WIDTH = 32

# Gaussian window is truncated at 3 standard deviations. Changing this changes the
# neighbourhood size.
_TRUNCATE = 3.0


def _filter_radius(radius: float) -> int:
    """How far the Gaussian window actually reaches, in pixels.

    Mirrors scipy.ndimage.gaussian_filter, which uses int(truncate * sigma + 0.5).
    Deriving it any other way makes the reported window size, and the width of the
    blanked border, disagree with the filter that was really applied.
    """
    return int(_TRUNCATE * radius + 0.5)


# --------------------------------------------------------------------------------
# Progress reporting
# --------------------------------------------------------------------------------

def _format_time(seconds: float) -> str:
    """Seconds as m:ss, or h:mm:ss once it runs past an hour."""
    seconds = int(round(seconds))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes}:{secs:02d}"


class ProgressBar:
    """A single-line progress bar that redraws itself in place.

    Prints something like::

        DSSIM  (##############------------------)  78/172   45%  0:16 elapsed, ~0:19 left

    Uses a carriage return to overwrite one line rather than printing a new line per
    step, so a long run leaves one line of output instead of hundreds. Works in Jupyter,
    VS Code notebooks and a plain terminal.

    ``total`` may be 0 or unknown, in which case only a count is shown.
    """

    def __init__(self, total: int | None, label: str = "", width: int = PROGRESS_WIDTH,
                 enabled: bool = True, stream=None):
        self.total = int(total) if total else None
        self.label = label
        self.width = width
        self.enabled = enabled
        self.stream = stream if stream is not None else sys.stdout
        self.start = time.perf_counter()
        self._last_filled = -1
        self._last_drawn = 0.0
        self._last_text = None
        self._closed = False
        if self.enabled:
            self._draw(0)

    def update(self, done: int) -> None:
        """Redraw for ``done`` completed items, throttled to avoid flooding the output."""
        if not self.enabled or self._closed:
            return
        now = time.perf_counter()
        filled = self._filled(done)
        # Redraw only when the bar visibly changes, or once a second so the timer moves.
        # Every redraw is stored in the saved notebook, so keep the rate modest.
        if filled != self._last_filled or now - self._last_drawn > 1.0:
            self._draw(done)

    def close(self, done: int | None = None) -> None:
        """Draw the final state and move to a new line."""
        if not self.enabled or self._closed:
            return
        self._draw(self.total if done is None else done, final=True)
        self.stream.write("\n")
        self.stream.flush()
        self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def _filled(self, done: int) -> int:
        if not self.total:
            return 0
        return int(self.width * min(done / self.total, 1.0))

    def _draw(self, done: int, final: bool = False) -> None:
        elapsed = time.perf_counter() - self.start
        if self.total:
            frac = min(done / self.total, 1.0)
            filled = self._filled(done)
            bar = "#" * filled + "-" * (self.width - filled)
            text = (f"\r{self.label}({bar}) {done}/{self.total} {frac * 100:4.0f}%  "
                    f"{_format_time(elapsed)} elapsed")
            if 0 < frac < 1 and not final:
                text += f", ~{_format_time(elapsed / frac - elapsed)} left"
            self._last_filled = filled
        else:
            text = f"\r{self.label}{done} frames  {_format_time(elapsed)} elapsed"
        # Skip a redraw that would produce exactly what is already on screen.
        if text == self._last_text:
            return
        # Pad so a shorter line never leaves fragments of the previous one behind.
        self.stream.write(text.ljust(96))
        self.stream.flush()
        self._last_text = text
        self._last_drawn = time.perf_counter()


# --------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------

def load_data(data, max_frames: int | None = None, verbose: bool = True) -> np.ndarray:
    """Load a dataset into a single-channel float32 volume of shape (frames, h, w).

    Accepts any of:

    - a path to a video file (``.avi``, ``.mp4``, ``.mov``, ``.mkv``)
    - a path to a multi-page TIFF stack (``.tif`` / ``.tiff``)
    - a path to a directory of single-page TIFFs (sorted by filename)
    - a NumPy array already shaped (frames, height, width) or (frames, h, w, channels)
    - a list/tuple of 2D arrays

    Colour frames are collapsed to one channel by averaging (to a grayscale image).
    Pixel values are *not* rescaled here; ``dssim_analysis`` normalises internally.

    ``max_frames`` loads only the first N frames. Use it to try parameters quickly on
    a large dataset before committing to a full run.
    """
    if isinstance(data, np.ndarray):
        vol = _to_single_channel(data)
    elif isinstance(data, (list, tuple)):
        vol = _to_single_channel(np.stack([np.asarray(f) for f in data]))
    elif isinstance(data, (str, os.PathLike)):
        path = os.path.normpath(os.fspath(data))
        if not os.path.exists(path):
            raise FileNotFoundError(f"No such file or directory: {path}")
        if os.path.isdir(path):
            vol = _load_tiff_folder(path, max_frames, verbose)
        else:
            ext = os.path.splitext(path)[1].lower()
            if ext in (".avi", ".mp4", ".mov", ".mkv", ".m4v"):
                vol = _load_video(path, max_frames, verbose)
            elif ext in (".tif", ".tiff"):
                vol = _load_tiff_stack(path)
            else:
                raise ValueError(
                    f"Unsupported file type '{ext}'. Expected a video (.avi/.mp4/.mov/.mkv), "
                    "a TIFF stack (.tif/.tiff), or a directory of TIFFs."
                )
    else:
        raise TypeError(
            "data must be a path, a NumPy array, or a list of 2D arrays; "
            f"got {type(data).__name__}"
        )

    if vol.ndim != 3:
        raise ValueError(f"Expected a 3D volume (frames, height, width), got shape {vol.shape}")
    if max_frames is not None:
        vol = vol[:max_frames]
    if vol.shape[0] < 2:
        raise ValueError(f"Need at least 2 frames to compute DSSIM, got {vol.shape[0]}")

    if verbose:
        n, h, w = vol.shape
        print(f"Loaded {n} frames of {w} x {h} px  ({vol.nbytes / 1e6:.1f} MB as float32)")
    return vol


def _natural_key(path: str):
    """Sort key that orders frame_2 before frame_10.

    A plain lexicographic sort puts frame_10 between frame_1 and frame_2, which
    silently reorders an unpadded image sequence and makes DSSIM compare frames
    that are not adjacent in time.
    """
    name = os.path.basename(path)
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", name)]


def _to_single_channel(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 4:
        arr = arr.mean(axis=3)
    return np.ascontiguousarray(arr, dtype=np.float32)


def _load_video(path: str, max_frames: int | None = None, verbose: bool = True) -> np.ndarray:
    cap = cv2.VideoCapture(path)
    if not cap.isOpened():
        raise IOError(f"Could not open video file: {path}")

    # Container metadata is only a hint - some files report 0 or a wrong count - so the
    # bar falls back to a plain frame counter when it looks unreliable.
    reported = int(cap.get(cv2.CAP_PROP_FRAME_COUNT) or 0)
    expected = min(reported, max_frames) if max_frames else reported
    bar = ProgressBar(expected if expected > 0 else None, label="Loading    ", enabled=verbose)

    frames = []
    try:
        while max_frames is None or len(frames) < max_frames:
            ok, frame = cap.read()
            if not ok:
                break
            if frame.ndim == 3:
                # Unweighted channel mean, to match _to_single_channel and the
                # documented behaviour. cv2.COLOR_BGR2GRAY applies luma weights
                # instead, which would make video input disagree with array input.
                frame = frame.mean(axis=2)
            frames.append(frame.astype(np.float32))
            bar.update(len(frames))
    finally:
        cap.release()
        bar.close(len(frames))
    if not frames:
        raise IOError(f"No frames could be read from {path}. The codec may be unsupported.")
    return np.stack(frames)


def _load_tiff_stack(path: str) -> np.ndarray:
    import tifffile
    # Ask the file how many pages it has rather than guessing from the array shape:
    # a single RGB image is (h, w, 3), which is indistinguishable by shape alone
    # from a 3-frame stack and would otherwise be read as a stack of rows.
    with tifffile.TiffFile(path) as handle:
        n_pages = len(handle.pages)
        arr = np.asarray(handle.asarray())
    if n_pages < 2:
        raise ValueError(
            f"{os.path.basename(path)} is a single image, not a stack. "
            "Point load_data() at the containing folder to load a sequence of TIFFs."
        )
    return _to_single_channel(arr)


def _load_tiff_folder(folder: str, max_frames: int | None = None, verbose: bool = True) -> np.ndarray:
    import tifffile
    paths = sorted(set(glob(os.path.join(folder, "*.tif")) + glob(os.path.join(folder, "*.tiff"))),
                   key=_natural_key)
    if not paths:
        raise FileNotFoundError(f"No .tif or .tiff files found in {folder}")
    if max_frames is not None:
        paths = paths[:max_frames]
    frames = []
    bar = ProgressBar(len(paths), label="Loading    ", enabled=verbose)
    try:
        for p in paths:
            frames.append(_to_single_channel(np.asarray(tifffile.imread(p))[None])[0])
            bar.update(len(frames))
    finally:
        bar.close()
    shapes = {f.shape for f in frames}
    if len(shapes) > 1:
        raise ValueError(f"TIFFs in {folder} have inconsistent shapes: {sorted(shapes)}")
    return np.stack(frames)


def load_times(times, n_frames: int, allow_truncate: bool = False,
               verbose: bool = True) -> np.ndarray:
    """Resolve a time vector for ``n_frames`` frames.

    ``times`` may be a path to a one-column CSV, a sequence of numbers, or None.
    If it is None, unreadable, or the wrong length, falls back to 1 second per
    frame and says so.

    ``allow_truncate`` handles the case where you deliberately loaded a subset of the
    data with ``max_frames``: a longer time vector is then trimmed to its first
    ``n_frames`` entries instead of being rejected. Leave it False otherwise, so that
    a mismatch between your times and your data is reported.
    """
    fallback = np.arange(n_frames, dtype=float)

    if times is None:
        if verbose:
            print("No time vector given; assuming 1 second per frame.")
        return fallback

    if isinstance(times, (str, os.PathLike)):
        path = os.path.normpath(os.fspath(times))
        try:
            # A one-column CSV may or may not have a header row. If the first line
            # parses as a number it is data, not a header - reading it as a header
            # would silently drop the first time point.
            with open(path, "r") as handle:
                first_field = handle.readline().split(",")[0].strip().strip('"')
            try:
                float(first_field)
                header = None
            except ValueError:
                header = "infer"
            table = pd.read_csv(path, header=header)
            values = np.asarray(table.iloc[:, 0], dtype=float)
        except Exception as exc:
            print(f"Warning: could not read times from '{path}' ({exc}). Assuming 1 second per frame.")
            return fallback
    else:
        values = np.asarray(list(times), dtype=float)

    if values.size == 0:
        if verbose:
            print("Empty time vector; assuming 1 second per frame.")
        return fallback

    if values.size != n_frames:
        if allow_truncate and values.size > n_frames:
            if verbose:
                print(f"Using the first {n_frames} of {values.size} time points "
                      f"({values[0]:g} to {values[n_frames - 1]:g}).")
            return values[:n_frames]
        print(
            f"Warning: got {values.size} time points for {n_frames} frames. "
            "Assuming 1 second per frame."
            + (" If you loaded a subset with max_frames, pass truncate_times=True."
               if values.size > n_frames else "")
        )
        return fallback

    if verbose:
        print(f"Loaded {values.size} time points ({values[0]:g} to {values[-1]:g}).")
    return values


# --------------------------------------------------------------------------------
# Pre-processing
# --------------------------------------------------------------------------------

def gaussian_denoise(volume: np.ndarray, sigma: float, verbose: bool = True) -> np.ndarray:
    """Apply a 2D Gaussian blur to every frame. ``sigma=0`` returns the input unchanged.

    Reproduces MATLAB's ``imgaussfilt(im, sigma)`` exactly, so the Python and MATLAB
    implementations denoise identically - see the comment in the loop below.
    """
    if sigma is None or sigma <= 0:
        return volume
    volume = np.asarray(volume)
    if volume.dtype != np.float32:
        # scipy would promote integer input to float64, which wastes memory here.
        volume = volume.astype(np.float32)
    out = np.empty_like(volume, dtype=np.float32)
    bar = ProgressBar(volume.shape[0], label="Denoising  ", enabled=verbose)
    try:
        for i in range(volume.shape[0]):
            # Matched to MATLAB's imgaussfilt so both implementations denoise the same
            # way. MATLAB truncates the kernel at 2 standard deviations
            # (FilterSize = 2*ceil(2*sigma)+1) and pads by replicating the edge pixel.
            #
            # cv2.GaussianBlur(im, (0, 0), sigma) was used here previously. It picks its
            # own, much wider kernel for float input - 9 px at sigma=1 where MATLAB uses
            # 5 - and reflects at the border rather than replicating. That was worth
            # about 1.3% in mean DSSIM on the bundled example.
            out[i] = gaussian_filter(volume[i], sigma=sigma, truncate=2.0, mode="nearest")
            bar.update(i + 1)
    finally:
        bar.close()
    return out


def _normalise(volume: np.ndarray) -> np.ndarray:
    """Rescale the whole volume to [0, 1]"""
    vol = volume.astype(np.float32, copy=True)
    lo = float(vol.min())
    hi = float(vol.max())
    if hi <= lo:
        return np.zeros_like(vol)
    vol -= lo
    vol /= (hi - lo)
    return vol


# --------------------------------------------------------------------------------
# SSIM / DSSIM core
# --------------------------------------------------------------------------------

def ssim_map(
    im1: np.ndarray,
    im2: np.ndarray,
    radius: float = 3.0,
    exponents: Sequence[float] = (1.0, 1.0, 1.0),
    data_range: float = 1.0,
) -> np.ndarray:
    """Per-pixel SSIM map between two 2D images.

    This had to be rewritten from scratch because I couldn't find a packaged version that supports exponents

    Uses a Gaussian weighting window of standard deviation ``radius``, truncated at
    3 sigma, giving a neighbourhood of ``2 * ceil(3 * radius) + 1`` pixels per side.

    ``exponents`` are ``[alpha, beta, gamma]``, weighting the mean, variance and
    cross correlation terms. With the default ``[1, 1, 1]`` the standard closed form is used.
    """
    im1 = np.asarray(im1, dtype=np.float64)
    im2 = np.asarray(im2, dtype=np.float64)
    if im1.shape != im2.shape:
        raise ValueError(f"Images must have the same shape, got {im1.shape} and {im2.shape}")

    c1 = (K1 * data_range) ** 2
    c2 = (K2 * data_range) ** 2

    def blur(x):
        return gaussian_filter(x, sigma=radius, truncate=_TRUNCATE, mode="nearest")

    mu1 = blur(im1)
    mu2 = blur(im2)
    mu1_sq, mu2_sq, mu1_mu2 = mu1 * mu1, mu2 * mu2, mu1 * mu2

    # Population (not sample) variance, from the Gaussian-weighted window.
    sigma1_sq = np.maximum(blur(im1 * im1) - mu1_sq, 0.0)
    sigma2_sq = np.maximum(blur(im2 * im2) - mu2_sq, 0.0)
    sigma12 = blur(im1 * im2) - mu1_mu2

    alpha, beta, gamma = (float(e) for e in exponents)

    if alpha == 1.0 and beta == 1.0 and gamma == 1.0:
        # Standard SSIM; the C3 = C2/2 substitution lets the variance and
        # cross-correlation terms collapse into a single term.
        numerator = (2 * mu1_mu2 + c1) * (2 * sigma12 + c2)
        denominator = (mu1_sq + mu2_sq + c1) * (sigma1_sq + sigma2_sq + c2)
        return (numerator / denominator).astype(np.float32)

    # General form: separate mean / variance / cross-correlation terms.
    c3 = c2 / 2.0
    sigma1 = np.sqrt(sigma1_sq)
    sigma2 = np.sqrt(sigma2_sq)

    mean_term = (2 * mu1_mu2 + c1) / (mu1_sq + mu2_sq + c1)
    variance_term = (2 * sigma1 * sigma2 + c2) / (sigma1_sq + sigma2_sq + c2)
    xcorr_term = (sigma12 + c3) / (sigma1 * sigma2 + c3)

    # The cross-correlation term can be negative (anti-correlated neighbourhoods).
    # Raising a negative number to a fractional power is undefined, so take the
    # magnitude and restore the sign afterwards.
    def signed_pow(x, e):
        # The cross-correlation term can be negative (anti-correlated neighbourhoods),
        # so the exponent has to be applied with care.
        if e == 1.0:
            return x
        if e == 0.0:
            # Anything to the power 0 is 1, which drops the term entirely.
            return np.ones_like(x)
        if float(e).is_integer():
            # Integer powers are well defined for a negative base, and carry the
            # correct sign themselves - an even power of a negative number is
            # positive. Forcing the sign here would negate that result.
            return x ** int(e)
        # A fractional power of a negative number is not real, so take the magnitude
        # and restore the sign afterwards.
        return np.sign(x) * np.abs(x) ** e

    result = (signed_pow(mean_term, alpha) * signed_pow(variance_term, beta)
              * signed_pow(xcorr_term, gamma))
    return result.astype(np.float32)


def _dssim_volume(volume: np.ndarray, frame_offset: int, radius: float,
                  exponents: Sequence[float], verbose: bool) -> np.ndarray:
    n_out = volume.shape[0] - frame_offset
    out = np.empty((n_out,) + volume.shape[1:], dtype=np.float32)
    bar = ProgressBar(n_out, label="DSSIM      ", enabled=verbose)
    try:
        for i in range(n_out):
            s = ssim_map(volume[i], volume[i + frame_offset], radius=radius,
                         exponents=exponents, data_range=1.0)
            out[i] = (1.0 - s) / 2.0
            bar.update(i + 1)
    finally:
        bar.close()
    return out


# --------------------------------------------------------------------------------
# Display / contrast
# --------------------------------------------------------------------------------

def _contrast_limits(values: np.ndarray, low_pct: float, high_pct: float) -> tuple[float, float]:
    """Intensity limits after discarding the given percentages of extreme values."""
    if low_pct <= 0 and high_pct <= 0:
        return float(values.min()), float(values.max())
    lo, hi = np.percentile(values.reshape(-1), [low_pct, 100.0 - high_pct])
    return float(lo), float(hi)


def _cmap_lut(cmap_name: str, n: int = 256) -> np.ndarray:
    """256-entry uint8 RGB lookup table for a named matplotlib colormap."""
    cmap = matplotlib.colormaps[cmap_name]
    return (np.asarray(cmap(np.linspace(0.0, 1.0, n)))[:, :3] * 255).round().astype(np.uint8)


def _scale_to_uint8_index(frame: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Rescale one frame from [lo, hi] to 0-255 colormap indices."""
    if hi <= lo:
        return np.zeros(frame.shape, dtype=np.uint8)
    scaled = (frame.astype(np.float32, copy=True) - lo) * (255.0 / (hi - lo))
    np.clip(scaled, 0.0, 255.0, out=scaled)
    return scaled.round().astype(np.uint8)


def apply_colormap(volume01: np.ndarray, cmap_name: str = "viridis") -> np.ndarray:
    """Map a [0, 1] volume (or single frame) through a colormap to uint8 RGB.

    Works frame by frame through a 256-entry lookup table, so peak memory is the
    output array plus one frame, not a float64 RGBA copy of the whole volume.
    """
    lut = _cmap_lut(cmap_name)
    vol = np.asarray(volume01)
    single = vol.ndim == 2
    if single:
        vol = vol[None]
    out = np.empty(vol.shape + (3,), dtype=np.uint8)
    for i in range(vol.shape[0]):
        out[i] = lut[_scale_to_uint8_index(vol[i], 0.0, 1.0)]
    return out[0] if single else out


def to_display_rgb(
    volume: np.ndarray,
    contrast_type: str = "constant_contrast",
    low_pct: float = 0.1,
    high_pct: float = 0.1,
    cmap_name: str = "viridis",
) -> np.ndarray:
    """Contrast-adjust a volume and colour-map it to uint8 RGB frames.

    ``contrast_type`` is either:

    - ``"constant_contrast"`` - one intensity scale for the whole video, so pixel brightness
      is comparable between frames. Use this when comparing frames to each other.
    - ``"per_frame_contrast"`` - each frame contrast is applied independently, maximizing
      within-frame detail but making frames non-comparable.

    ``low_pct`` / ``high_pct`` are the percentages of extreme values clipped before
    rescaling. Set both to 0 to disable outlier removal.
    """
    if contrast_type not in ("constant_contrast", "per_frame_contrast"):
        raise ValueError(
            'contrast_type must be "constant_contrast" or "per_frame_contrast", '
            f'got {contrast_type!r}'
        )

    vol = np.asarray(volume)
    lut = _cmap_lut(cmap_name)
    out = np.empty(vol.shape + (3,), dtype=np.uint8)

    if contrast_type == "constant_contrast":
        lo, hi = _contrast_limits(vol, low_pct, high_pct)

    for i in range(vol.shape[0]):
        if contrast_type == "per_frame_contrast":
            lo, hi = _contrast_limits(vol[i], low_pct, high_pct)
        out[i] = lut[_scale_to_uint8_index(vol[i], lo, hi)]
    return out


def _zero_borders(vol: np.ndarray, rgb: np.ndarray, dist: int, cmap_name: str = "viridis"):
    """Zero a border of width ``dist`` on the DSSIM volume and paint it the colormap's
    zero colour in the RGB frames.

    The SSIM window reaches ``dist`` pixels past each edge, so border pixels are
    computed from padded data and are artefacts rather than measurements.
    """
    if dist <= 0:
        return vol, rgb
    h, w = vol.shape[1], vol.shape[2]
    if 2 * dist >= min(h, w):
        raise ValueError(
            f"remove_border_dist={dist} would erase the whole {w}x{h} frame. "
            "Reduce radius, or set remove_border_dist=0."
        )

    vol[:, :dist, :] = 0
    vol[:, -dist:, :] = 0
    vol[:, :, :dist] = 0
    vol[:, :, -dist:] = 0

    zero_colour = _cmap_lut(cmap_name)[0]
    rgb[:, :dist, :, :] = zero_colour
    rgb[:, -dist:, :, :] = zero_colour
    rgb[:, :, :dist, :] = zero_colour
    rgb[:, :, -dist:, :] = zero_colour
    return vol, rgb


# --------------------------------------------------------------------------------
# Result container and main entry point
# --------------------------------------------------------------------------------

@dataclass
class DssimResult:
    """Everything produced by a DSSIM run.

    Attributes
    ----------
    dssim : (n_dssim, h, w) float32
        The raw DSSIM maps. This is the quantitative output.
    rgb_frames : (n_dssim, h, w, 3) uint8
        Contrast-adjusted, colour-mapped frames for display and video export.
    means : (n_dssim,) float
        Mean DSSIM per frame - the structural-change-over-time curve. The blanked
        border counts as zeros, so this value depends on frame size and ``radius``.
        Curves are comparable across runs at fixed settings, but not across different
        radii or crops. Average ``dssim[:, b:-b, b:-b]`` (b = ``remove_border_dist``)
        for an interior-only figure.
    stats : pandas.DataFrame
        Per-DSSIM-frame table: frame indices, source frame numbers, times, and mean value.
    aligned_indices : (n_dssim,) int
        Which original data frames line up with each DSSIM frame, for side-by-side display.
    window_size : int
        Side length in pixels of the SSIM neighbourhood.
    params : dict
        The parameters this run used.
    """
    dssim: np.ndarray
    rgb_frames: np.ndarray
    means: np.ndarray
    stats: pd.DataFrame
    aligned_indices: np.ndarray
    window_size: int
    params: dict = field(default_factory=dict)

    def __repr__(self) -> str:
        n, h, w = self.dssim.shape
        return (
            f"DssimResult({n} frames, {w}x{h} px, "
            f"mean DSSIM {self.means.mean():.4f}, window {self.window_size}px)"
        )


def dssim_analysis(
    data,
    times=None,
    frame_offset: int = 1,
    radius: float = 3.0,
    exponents: Sequence[float] = (1.0, 1.0, 1.0),
    contrast_type: str = "constant_contrast",
    remove_border_dist: int | None = None,
    low_pct: float = 0.1,
    high_pct: float = 0.1,
    cmap_name: str = "viridis",
    truncate_times: bool = False,
    verbose: bool = True,
) -> DssimResult:
    """Run DSSIM analysis on a video, image stack, or array.

    Parameters
    ----------
    data
        Anything ``load_data`` accepts, or an already-loaded (frames, h, w) array.
    times
        Path to a one-column CSV of frame times, a sequence of numbers, or None
        for 1 second per frame.
    frame_offset
        How many frames apart the compared pair is. 1 compares frame n to n+1.
        Larger values look at slower structural change.
    radius
        Standard deviation of the Gaussian neighbourhood, in pixels. Bigger means
        smoother, coarser change maps.
    exponents
        ``[alpha, beta, gamma]`` weights on the mean, variance and cross-correlation
        terms of SSIM. Leave at [1, 1, 1] unless you have a reason not to.
    contrast_type
        ``"constant_contrast"`` or ``"per_frame_contrast"`` - see ``to_display_rgb``.
    remove_border_dist
        Width of the border to blank out. Defaults to ``ceil(radius * 3)``, which is
        exactly how far the SSIM window reaches. Set 0 to keep border pixels.
    low_pct, high_pct
        Percent of extreme values clipped before display scaling. Display only -
        the raw ``dssim`` array is unaffected.
    truncate_times
        Set True when you loaded a subset with ``max_frames`` and your time vector
        covers the full dataset; the extra time points are then trimmed rather than
        the whole vector being rejected.
    """
    if frame_offset < 1:
        raise ValueError(f"frame_offset must be at least 1, got {frame_offset}")
    if radius <= 0:
        raise ValueError(f"radius must be positive, got {radius}")
    if remove_border_dist is None:
        remove_border_dist = _filter_radius(radius)
    remove_border_dist = int(remove_border_dist)

    volume = data if isinstance(data, np.ndarray) and data.ndim == 3 else load_data(data, verbose=verbose)
    n_frames = volume.shape[0]

    if frame_offset >= n_frames:
        raise ValueError(
            f"frame_offset={frame_offset} needs more than {frame_offset} frames, "
            f"but the data has {n_frames}."
        )

    time_values = load_times(times, n_frames, allow_truncate=truncate_times, verbose=verbose)

    if verbose:
        window = 2 * _filter_radius(radius) + 1
        print(f"Computing {n_frames - frame_offset} DSSIM frames "
              f"(offset {frame_offset}, {window}x{window} px window)...")

    vol_dssim = _dssim_volume(_normalise(volume), frame_offset, radius, exponents, verbose)

    rgb = to_display_rgb(vol_dssim, contrast_type=contrast_type,
                         low_pct=low_pct, high_pct=high_pct, cmap_name=cmap_name)

    vol_dssim, rgb = _zero_borders(vol_dssim, rgb, remove_border_dist, cmap_name=cmap_name)

    means = vol_dssim.reshape(vol_dssim.shape[0], -1).mean(axis=1)

    n_dssim = vol_dssim.shape[0]
    frame_1 = np.arange(n_dssim)
    frame_2 = frame_1 + frame_offset
    t1 = time_values[frame_1]
    t2 = time_values[frame_2]

    stats = pd.DataFrame({
        # 0-based throughout, matching Python indexing and result.aligned_indices,
        # so a frame number here indexes straight into result.dssim / your data.
        "dssim_frame_num": frame_1,
        "dssim_mean_value": means,
        "data_frame_1_num": frame_1,
        "data_frame_2_num": frame_2,
        "data_frame_1_time": t1,
        "data_frame_2_time": t2,
        "dssim_frame_mean_time": (t1 + t2) / 2.0,
    })

    # If frame_offset is odd there is no exact middle frame; favour the later one.
    first = math.ceil(frame_offset / 2)
    aligned = np.arange(first, first + n_dssim)

    params = {
        "frame_offset": frame_offset,
        "radius": radius,
        "exponents": list(exponents),
        "contrast_type": contrast_type,
        "remove_border_dist": remove_border_dist,
        "low_pct": low_pct,
        "high_pct": high_pct,
        "cmap_name": cmap_name,
        "n_data_frames": n_frames,
        "n_dssim_frames": n_dssim,
    }

    if verbose:
        print(f"Done. Mean DSSIM across the video: {means.mean():.4f}")

    return DssimResult(
        dssim=vol_dssim,
        rgb_frames=rgb,
        means=means,
        stats=stats,
        aligned_indices=aligned,
        window_size=2 * _filter_radius(radius) + 1,
        params=params,
    )


# --------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------

_FOURCC = {".avi": "XVID", ".mp4": "mp4v", ".m4v": "mp4v", ".mov": "mp4v", ".mkv": "XVID"}


def write_video(frames, output_path: str, fps: int = 10) -> str:
    """Write frames to a video file.

    ``frames`` may be a uint8 array shaped (n, h, w, 3) or (n, h, w), or any iterable
    yielding uint8 RGB frames one at a time. The iterable form lets you export a video
    larger than memory - see ``save_results`` for an example.

    The codec is chosen from the file extension: XVID for ``.avi``, mp4v for ``.mp4``.
    A frame with an odd width or height gains one duplicated edge row or column, because
    most codecs round the frame size down to even and would otherwise discard it.
    Returns the path written.
    """
    if isinstance(frames, np.ndarray):
        if frames.ndim == 3:
            frames = np.repeat(frames[..., None], 3, axis=3)
        if frames.ndim != 4 or frames.shape[-1] != 3:
            raise ValueError(f"Expected frames shaped (n, h, w, 3) or (n, h, w), got {frames.shape}")
        if frames.dtype != np.uint8:
            raise ValueError(
                f"Frames must be uint8, got {frames.dtype}. "
                "Use to_display_rgb() to convert a float DSSIM volume first."
            )

    ext = os.path.splitext(output_path)[1].lower()
    fourcc = cv2.VideoWriter_fourcc(*_FOURCC.get(ext, "XVID"))

    writer = None
    count = 0
    pad_h = pad_w = 0
    try:
        for frame in frames:
            frame = np.asarray(frame)
            if frame.ndim == 2:
                frame = np.repeat(frame[..., None], 3, axis=2)
            if frame.dtype != np.uint8:
                raise ValueError(f"Frames must be uint8, got {frame.dtype}.")
            if writer is None:
                h, w = frame.shape[:2]
                # Most codecs round the frame size down to an even number, which would
                # silently drop the last row or column of the analysis. Repeat the edge
                # pixel instead so nothing is lost.
                pad_h, pad_w = h % 2, w % 2
                writer = cv2.VideoWriter(output_path, fourcc, fps, (w + pad_w, h + pad_h))
                if not writer.isOpened():
                    raise IOError(
                        f"Could not open a video writer for {output_path}. The codec may not "
                        "be available in this OpenCV build; try a .avi extension."
                    )
            if pad_h or pad_w:
                frame = np.pad(frame, ((0, pad_h), (0, pad_w), (0, 0)), mode="edge")
            writer.write(cv2.cvtColor(frame, cv2.COLOR_RGB2BGR))
            count += 1
    finally:
        if writer is not None:
            writer.release()

    if count == 0:
        raise ValueError("No frames to write.")
    if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
        raise IOError(f"Video writing produced no output at {output_path}.")
    return output_path


def side_by_side(left: np.ndarray, right: np.ndarray, gap: int = 0) -> np.ndarray:
    """Join two equal-length uint8 RGB stacks horizontally, optionally with a black gap."""
    left = np.asarray(left)
    right = np.asarray(right)
    if left.ndim == 3:
        left = np.repeat(left[..., None], 3, axis=3)
    if right.ndim == 3:
        right = np.repeat(right[..., None], 3, axis=3)
    if left.shape[0] != right.shape[0]:
        raise ValueError(f"Stacks have different frame counts: {left.shape[0]} and {right.shape[0]}")
    if left.shape[1] != right.shape[1]:
        raise ValueError(f"Stacks have different heights: {left.shape[1]} and {right.shape[1]}")
    if gap > 0:
        spacer = np.zeros((left.shape[0], left.shape[1], gap, 3), dtype=np.uint8)
        return np.concatenate([left, spacer, right], axis=2)
    return np.concatenate([left, right], axis=2)


def save_results(result: DssimResult, run_name: str, output_dir: str = ".",
                 data_frames: np.ndarray | None = None, fps: int = 10,
                 data_contrast_type: str = "per_frame_contrast",
                 verbose: bool = True) -> dict:
    """Write the standard output set for a run.

    Always writes:
      - ``<run_name>_DSSIM_video.avi``     - the color-mapped DSSIM video
      - ``<run_name>_DSSIM_values.csv``    - the per-frame stats table
      - ``<run_name>_DSSIM_run_info.csv``  - the parameters used

    If ``data_frames`` (the original loaded volume) is supplied, also writes:
      - ``<run_name>_side_by_side.avi``    - original data next to the DSSIM map

    Returns a dict of the paths written.
    """
    if data_contrast_type not in ("constant_contrast", "per_frame_contrast"):
        raise ValueError(
            'data_contrast_type must be "constant_contrast" or "per_frame_contrast", '
            f'got {data_contrast_type!r}'
        )

    os.makedirs(output_dir, exist_ok=True)
    written = {}

    def out(name):
        return os.path.join(output_dir, f"{run_name}{name}")

    written["dssim_video"] = write_video(result.rgb_frames, out("_DSSIM_video.avi"), fps=fps)

    result.stats.to_csv(out("_DSSIM_values.csv"), index=False)
    written["values_csv"] = out("_DSSIM_values.csv")

    info = dict(result.params)
    info["run_name"] = run_name
    info["window_size_px"] = result.window_size
    info["data_contrast_type"] = data_contrast_type
    info["exponents"] = str(info["exponents"])
    pd.DataFrame([info]).to_csv(out("_DSSIM_run_info.csv"), index=False)
    written["run_info_csv"] = out("_DSSIM_run_info.csv")

    if data_frames is not None:
        # Streamed one frame at a time: a full-resolution side-by-side stack of a long
        # video is easily several GB, and never needs to exist all at once.
        # A contiguous slice is a view; fancy indexing would copy the whole volume,
        # which for a full-resolution dataset is exactly the allocation this
        # streaming path exists to avoid.
        first_aligned = int(result.aligned_indices[0])
        aligned = data_frames[first_aligned:first_aligned + len(result.aligned_indices)]
        lut_gray = _cmap_lut("gray")
        low_pct = result.params["low_pct"]
        high_pct = result.params["high_pct"]
        if data_contrast_type == "constant_contrast":
            g_lo, g_hi = _contrast_limits(aligned, low_pct, high_pct)

        def _pairs():
            for i in range(len(aligned)):
                if data_contrast_type == "per_frame_contrast":
                    lo, hi = _contrast_limits(aligned[i], low_pct, high_pct)
                else:
                    lo, hi = g_lo, g_hi
                left = lut_gray[_scale_to_uint8_index(aligned[i], lo, hi)]
                yield np.concatenate([left, result.rgb_frames[i]], axis=1)

        written["side_by_side_video"] = write_video(_pairs(), out("_side_by_side.avi"), fps=fps)

    if verbose:
        for label, path in written.items():
            print(f"  wrote {os.path.basename(path)}")
    return written
