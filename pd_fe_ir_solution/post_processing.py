import numpy as np
import numpy.typing as npt
from scipy.optimize import curve_fit
from scipy.interpolate import griddata, interp1d
from dataclasses import dataclass, field


@dataclass
class Sky:
    area: float = 0.0
    circumference: float = 0.0
    circularity: float = 0.0
    average_radius: float = 0.0
    rho_list: list[float] = field(default_factory=list)
    circ_grid: npt.NDArray = field(default_factory=lambda: np.array([]))


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

    distances_from_center = np.linalg.norm(positions - center, axis=1)
    spins_z = spins[:, 2]
    spin_z_center = spins_z[np.argmin(distances_from_center)]

    # We grow the ring_radius until we have more than `min_n_spins_outside` spins with a z-component larger/smaller (depending on the sign of spin_z_center) than spin_z_center
    for i in range(max_iter + 1):
        mask_ring = distances_from_center < ring_radius
        spins_z_ring = spins_z[mask_ring]

        # multiply by -np.sign(spin_z_center) to account for the fact that the z-component of the spins outside the ring needs to be greater than z_value if the spin_z_center is negative and smaller than z_value if spin_z_center is positive
        n_spins_outside = np.count_nonzero(
            spins_z_ring > -np.sign(spin_z_center) * z_value
        )

        # if we have enough spins outside the ring we can exit the loop
        if n_spins_outside >= min_n_spins_outside:
            break
        else:  # else we increase the ring radius
            ring_radius *= ring_radius_growth_factor

        # if we ever reach max_iter, we throw an exception
        if i == max_iter:
            raise Exception(
                f"Did not find enough spins with a z_value above/below {z_value}. Found {n_spins_outside} spins, need {min_n_spins_outside}"
            )

    distances_from_center_ring = distances_from_center[mask_ring]

    # Once we have enough spins, we find the skyrmion radius by fitting a line through the spin_z components
    # find the spin_z components with the smallest values from "inside" and from "outside"
    mask_inside = spins_z_ring <= -np.sign(spin_z_center) * z_value
    spins_z_inside = spins_z_ring[mask_inside]
    n_spins_inside = np.count_nonzero(mask_inside)
    n_spins_fit_inside = max(int(0.05 * n_spins_inside), 1)
    mask_inside_fit = np.argsort(np.abs(spins_z_inside - z_value))[:n_spins_fit_inside]

    mask_outside = spins_z_ring > -np.sign(spin_z_center) * z_value
    spins_z_outside = spins_z_ring[mask_outside]
    n_spins_outside = np.count_nonzero(mask_outside)
    n_spins_fit_outside = max(int(0.05 * n_spins_outside), 1)
    mask_outside_fit = np.argsort(np.abs(spins_z_outside - z_value))[
        :n_spins_fit_outside
    ]

    distances_fit = np.hstack(
        (
            distances_from_center_ring[mask_inside_fit],
            distances_from_center_ring[mask_outside_fit],
        )
    )

    spins_z_fit = np.hstack(
        (spins_z_ring[mask_inside_fit], spins_z_ring[mask_outside_fit])
    )

    # this function is defined such that f(b) = z_value
    def func(x, a, b):
        return a * (x - b) + z_value

    popt, pcov = curve_fit(
        func, distances_fit, spins_z_fit, p0=[1.0, ring_radius / 2.0]
    )

    skyrmion_radius = popt[1]

    return skyrmion_radius
