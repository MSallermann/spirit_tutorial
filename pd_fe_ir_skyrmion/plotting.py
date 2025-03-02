import numpy.typing as npt
import numpy as np
from matplotlib.pyplot import Axes


def get_rgba_colors(
    spins: npt.NDArray,
    opacity: float = 1.0,
    cardinal_a: npt.NDArray = np.array([1, 0, 0]),
    cardinal_b: npt.NDArray = np.array([0, 1, 0]),
    cardinal_c: npt.NDArray = np.array([0, 0, 1]),
) -> npt.NDArray:
    """Returns a colormap in the matplotlib format (an Nx4 array)"""
    import math

    # Annoying OpenGl functions
    def atan2(y, x):
        if x == 0.0:
            return np.sign(y) * np.pi / 2.0
        else:
            return np.arctan2(y, x)

    def fract(x):
        return x - math.floor(x)

    def mix(x, y, a):
        return x * (1.0 - a) + y * a

    def clamp(x, minVal, maxVal):
        return min(max(x, minVal), maxVal)

    def hsv2rgb(c):
        K = [1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0]

        px = abs(fract(c[0] + K[0]) * 6.0 - K[3])
        py = abs(fract(c[0] + K[1]) * 6.0 - K[3])
        pz = abs(fract(c[0] + K[2]) * 6.0 - K[3])
        resx = c[2] * mix(K[0], clamp(px - K[0], 0.0, 1.0), c[1])
        resy = c[2] * mix(K[0], clamp(py - K[0], 0.0, 1.0), c[1])
        resz = c[2] * mix(K[0], clamp(pz - K[0], 0.0, 1.0), c[1])

        return [resx, resy, resz]

    rgba = []

    for spin in spins:
        projection_x = cardinal_a.dot(spin)
        projection_y = cardinal_b.dot(spin)
        projection_z = cardinal_c.dot(spin)
        hue = atan2(projection_x, projection_y) / (2 * np.pi)

        saturation = projection_z

        if saturation > 0.0:
            saturation = 1.0 - saturation
            value = 1.0
        else:
            value = 1.0 + saturation
            saturation = 1.0

        rgba.append([*hsv2rgb([hue, saturation, value]), opacity])

    return rgba


def plot_spins_2d(ax: Axes, positions: npt.NDArray, spins: npt.NDArray, **kwargs):

    x = positions[:, 0]
    y = positions[:, 1]

    u = spins[:, 0]
    v = spins[:, 1]

    c = get_rgba_colors(spins)

    ax.set_aspect("equal")

    ax.scatter(x, y, marker=".", color="black", s=1)
    cb = ax.quiver(x, y, u, v, color=c, **kwargs)

    ax.set_xlabel(r"$x~[\AA]$")
    ax.set_ylabel(r"$y~[\AA]$")
