"""Self-contained checks for dssim.py.

Run with::

    python test_dssim.py

The important one is the first: it confirms this SSIM implementation agrees with
scikit-image's reference implementation, which is itself validated against Wang et al.
Everything else uses synthetic data, so no test files are needed.

scikit-image is required for this script but not for dssim.py itself.
"""

import math

import numpy as np

import dssim

rng = np.random.default_rng(0)
failures = []


def check(label, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {label}{'  ' + detail if detail else ''}")
    if not condition:
        failures.append(label)


print("1. ssim_map agrees with scikit-image's reference implementation")
try:
    from skimage.metrics import structural_similarity as sk_ssim
    for radius in (1.5, 3.0, 5.0):
        a = rng.random((128, 160)).astype(np.float32)
        b = np.clip(a + rng.normal(0, 0.08, a.shape), 0, 1).astype(np.float32)
        mine = dssim.ssim_map(a, b, radius=radius)
        _, theirs = sk_ssim(a, b, full=True, gaussian_weights=True, sigma=radius,
                            truncate=3.0, use_sample_covariance=False, data_range=1.0)
        pad = math.ceil(radius * 3)          # ignore the border, where padding differs
        diff = np.abs(mine[pad:-pad, pad:-pad] - theirs[pad:-pad, pad:-pad]).max()
        check(f"radius={radius}", diff < 1e-3, f"max|diff| = {diff:.2e}")
except ImportError:
    print("  [SKIP] scikit-image not installed")

print("\n2. Identical frames give SSIM 1 / DSSIM 0")
a = rng.random((64, 64)).astype(np.float32)
s = dssim.ssim_map(a, a, radius=3.0)
check("SSIM == 1 everywhere", np.allclose(s, 1.0), f"min = {s.min():.8f}")

print("\n3. The general-exponent path reduces to the standard form at [1, 1, 1]")
b = np.clip(a + rng.normal(0, 0.1, a.shape), 0, 1).astype(np.float32)
std = dssim.ssim_map(a, b, radius=3.0, exponents=(1, 1, 1))
gen = dssim.ssim_map(a, b, radius=3.0, exponents=(1.0, 1.0, 0.999999))
check("agreement", np.abs(std - gen).max() < 1e-5, f"max|diff| = {np.abs(std - gen).max():.2e}")
for e in [(0.5, 1, 1), (1, 0.5, 1), (1, 1, 0.5), (2, 2, 2)]:
    v = dssim.ssim_map(a, b, radius=3.0, exponents=e)
    check(f"exponents={e} finite", bool(np.isfinite(v).all()), f"range [{v.min():.3f}, {v.max():.3f}]")

print("\n4. A moving square lights up only where it moved")
vol = np.zeros((6, 80, 80), np.float32)
for i in range(6):
    vol[i, 30:50, 10 + i * 8:30 + i * 8] = 1.0
r = dssim.dssim_analysis(vol, frame_offset=1, radius=3, verbose=False)
moved = r.dssim[0][30:50, 10:46].mean()
still = r.dssim[0][5:20, 5:75].mean()
check("shape", r.dssim.shape == (5, 80, 80), str(r.dssim.shape))
check("change region dominates background", moved > 50 * max(still, 1e-9),
      f"moved {moved:.4f} vs static {still:.4f}")

print("\n5. Border blanking width equals ceil(radius * 3)")
for radius in (1.5, 3.0, 4.0):
    r = dssim.dssim_analysis(vol, radius=radius, verbose=False)
    expected = math.ceil(radius * 3)
    ok = r.params["remove_border_dist"] == expected and np.all(r.dssim[:, :expected, :] == 0)
    check(f"radius={radius}", ok, f"border {r.params['remove_border_dist']} px, window {r.window_size} px")

print("\n6. Frame counts and data alignment for various offsets")
for off in (1, 2, 3, 5):
    r = dssim.dssim_analysis(vol, frame_offset=off, radius=1.5, verbose=False)
    ok = r.dssim.shape[0] == 6 - off and r.aligned_indices.max() < 6
    check(f"offset={off}", ok, f"n={r.dssim.shape[0]}, aligned={r.aligned_indices.tolist()}")

print("\n7. The stats table's midpoint time is a real midpoint")
r = dssim.dssim_analysis(vol, times=np.arange(6) * 0.2, frame_offset=2, radius=1.5, verbose=False)
s = r.stats
check("midpoint", bool(np.allclose(s.dssim_frame_mean_time,
                                   (s.data_frame_1_time + s.data_frame_2_time) / 2)))

print("\n8. Loading a folder of TIFFs, a TIFF stack and an array agree")
try:
    import os
    import tempfile

    import tifffile
    with tempfile.TemporaryDirectory() as tmp:
        folder = os.path.join(tmp, "tifs")
        os.makedirs(folder)
        for i, f in enumerate(vol):
            tifffile.imwrite(os.path.join(folder, f"frame_{i:03d}.tif"), f)
        stack = os.path.join(tmp, "stack.tif")
        tifffile.imwrite(stack, vol)
        from_array = dssim.dssim_analysis(vol, radius=1.5, verbose=False).dssim
        from_folder = dssim.dssim_analysis(folder, radius=1.5, verbose=False).dssim
        from_stack = dssim.dssim_analysis(stack, radius=1.5, verbose=False).dssim
    check("folder == array", bool(np.allclose(from_folder, from_array)))
    check("stack == array", bool(np.allclose(from_stack, from_array)))
except ImportError:
    print("  [SKIP] tifffile not installed")

print("\n9. Bad input raises a clear error")
cases = [
    ("single frame", lambda: dssim.load_data(np.zeros((1, 10, 10), np.float32))),
    ("frame_offset too large", lambda: dssim.dssim_analysis(vol, frame_offset=99, verbose=False)),
    ("border would erase frame", lambda: dssim.dssim_analysis(vol, radius=20, verbose=False)),
    ("missing file", lambda: dssim.load_data("/nonexistent/movie.avi")),
    ("bad contrast_type", lambda: dssim.dssim_analysis(vol, contrast_type="bogus", verbose=False)),
]
for label, fn in cases:
    try:
        fn()
        check(label, False, "no error raised")
    except (ValueError, FileNotFoundError, TypeError, IOError) as exc:
        check(label, True, f"{type(exc).__name__}: {str(exc)[:60]}")

print("\n" + "=" * 60)
if failures:
    print(f"{len(failures)} CHECK(S) FAILED: {failures}")
    raise SystemExit(1)
print("ALL CHECKS PASSED")
