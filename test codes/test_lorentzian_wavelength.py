import numpy as np
from scipy.optimize import leastsq, curve_fit
import matplotlib.pyplot as plt


def lorentzian(width, central, height, x):
    return height * width / (2 * np.pi) / ((x - central)**2 + width**2 / 4)


def error_func(p, x, y):
    return lorentzian(*p, x) - y


def find_r_squared(f, p, x, y):
    res = y - f(*p, x)
    ss_res = np.sum(res ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - ss_res / ss_tot


def compare_plot(x, y, p):
    fy = lorentzian(*p, x)
    fig, ax = plt.subplots(1)
    ax.plot(x, y)
    ax.plot(x, fy)
    plt.show()


def fit_lorentzian(scattering, wavelength, split=False):
    # remove nans

    to_del = ~np.isnan(scattering)
    scattering = scattering[to_del]
    wavelength = wavelength[to_del]

    # return if not enough points

    if len(scattering) < 5:
        return [np.nan, np.nan, np.nan, np.nan], [0, 0], 0

    # find max and min

    max_sca = np.max(scattering)
    idx_max = np.argmax(scattering)
    idx_min = np.argmin(scattering)

    # init guess and first fit

    init_guess = [100, wavelength[idx_max], max_sca]
    result, cov_x, res_dict, mesg, ier = leastsq(error_func, init_guess, args=(wavelength, scattering),
                                                 full_output=True)
    result[0] = abs(result[0])
    r_squared = find_r_squared(lorentzian, result, wavelength, scattering)

    # if r_squared is too low, split

    if r_squared < 0.9 and split is False:
        wavelength_low = wavelength[:idx_min]
        wavelength_high = wavelength[idx_min:]
        scattering_low = scattering[:idx_min]
        scattering_high = scattering[idx_min:]
        result_low, r_squared_low = fit_lorentzian(scattering_low, wavelength_low, split=True)
        result_high, r_squared_high = fit_lorentzian(scattering_high, wavelength_high, split=True)

        if r_squared_high > r_squared and ~np.isnan(np.sum(result_high)):
            result = result_high
            r_squared = r_squared_high
        if r_squared_low > r_squared and ~np.isnan(np.sum(result_low)):
            result = result_low
            r_squared = r_squared_low

    compare_plot(wavelength, scattering, result)

    return result, r_squared


wavelength = np.arange(570, 740, 10)
params = [50, 700, 100]
scattering = lorentzian(*params, wavelength)

fig, ax = plt.subplots(1)
ax.plot(wavelength, scattering)
plt.show()

result, r_squared = fit_lorentzian(scattering, wavelength)