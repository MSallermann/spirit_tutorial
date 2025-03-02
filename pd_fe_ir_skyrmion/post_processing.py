import numpy as np
import numpy.typing as npt
from scipy.optimize import curve_fit


def skyrmion_radius(
    positions: npt.NDArray,
    spins: npt.NDArray,
    center: npt.NDArray,
    z_value: float = 0.0,
    ring_radius: float = 1.0,
    ring_radius_growth_factor: float = 1.2,
    min_n_spins_outside: int = 1,
    max_iter: int = 50,
) -> float:
    """
    Computes the skyrmion radius, defined as the radius of the m_z=z_value contour.

    Args:
        positions (npt.NDArray): Cartesian positions of the spins. Shape (N, 3).
        spins (npt.NDArray): Directions of the spins. Shape (N, 3).
        center (npt.NDArray): Cartesian coordinates of the center of the skyrmion.
        z_value (float, optional): z-value defining the radius contour. Defaults to 0.0.
        ring_radius (float, optional): Initial guess for the ring radius in units of the position vector. Defaults to 1.0.
        ring_radius_growth_factor (float, optional): Factor by which to increase the ring radius while searching for the radius contour. Defaults to 1.2.
        min_n_spins_outside (int, optional): Minimum number of spins outside the contour required. Defaults to 1.
        max_iter (int, optional): Maximum number of iterations for increasing the ring radius. Defaults to 50.

    Returns:
        float: The estimated skyrmion radius.
    """
    distances_from_center = np.linalg.norm(positions - center, axis=1)
    spins_z = spins[:, 2]
    spin_z_center = spins_z[np.argmin(distances_from_center)]

    # Define a threshold based on the sign of the central spin.
    threshold = -np.sign(spin_z_center) * z_value

    found = False
    for i in range(max_iter):
        mask_ring = distances_from_center < ring_radius
        spins_z_ring = spins_z[mask_ring]
        n_spins_above_threshold = np.count_nonzero(spins_z_ring > threshold)
        if n_spins_above_threshold >= min_n_spins_outside:
            found = True
            break
        ring_radius *= ring_radius_growth_factor

    if not found:
        raise Exception(
            f"Did not find enough spins with a z_value {'above' if threshold > 0 else 'below'} {z_value}. "
            f"Found {n_spins_above_threshold} spins, need {min_n_spins_outside}."
        )

    distances_from_center_ring = distances_from_center[mask_ring]

    # Separate spins into "inside" (<= threshold) and "outside" (> threshold) groups.
    mask_inside = spins_z_ring <= threshold
    mask_outside = spins_z_ring > threshold

    spins_z_inside = spins_z_ring[mask_inside]
    spins_z_outside = spins_z_ring[mask_outside]

    distances_inside = distances_from_center_ring[mask_inside]
    distances_outside = distances_from_center_ring[mask_outside]

    # For each group, select a fraction (5%) of spins that are closest to z_value.
    # If a group is empty, simply skip it.
    if spins_z_inside.size > 0:
        n_fit_inside = max(int(0.05 * spins_z_inside.size), 1)
        inside_fit_indices = np.argsort(np.abs(spins_z_inside - z_value))[:n_fit_inside]
        distances_inside_fit = distances_inside[inside_fit_indices]
        spins_z_inside_fit = spins_z_inside[inside_fit_indices]
    else:
        distances_inside_fit = np.empty(0)
        spins_z_inside_fit = np.empty(0)

    if spins_z_outside.size > 0:
        n_fit_outside = max(int(0.05 * spins_z_outside.size), 1)
        outside_fit_indices = np.argsort(np.abs(spins_z_outside - z_value))[
            :n_fit_outside
        ]
        distances_outside_fit = distances_outside[outside_fit_indices]
        spins_z_outside_fit = spins_z_outside[outside_fit_indices]
    else:
        distances_outside_fit = np.empty(0)
        spins_z_outside_fit = np.empty(0)

    distances_fit = np.hstack((distances_inside_fit, distances_outside_fit))
    spins_z_fit = np.hstack((spins_z_inside_fit, spins_z_outside_fit))

    if distances_fit.size < 2:
        raise ValueError("Not enough data points for curve fitting.")

    # Define a linear function such that f(b) = z_value; b is the estimated radius.
    def func(x, a, b):
        return a * (x - b) + z_value

    try:
        popt, _ = curve_fit(
            func, distances_fit, spins_z_fit, p0=[1.0, ring_radius / 2.0]
        )
    except Exception as e:
        raise RuntimeError("Curve fitting failed: " + str(e))

    skyrmion_radius_estimated = popt[1]
    return skyrmion_radius_estimated
