import numpy as np
from scipy.optimize import leastsq, curve_fit
import matplotlib.pyplot as plt


def lorentzian(p, x):
    return p[0] + p[1] / ((x - p[2]) ** 2 + (0.5 * p[3]) ** 2)


def lorentzian2(p0, p1, p2, p3, x):
    return p0 + p1 / ((x - p2) ** 2 + (0.5 * p3) ** 2)


def lorentzian_wavelength(width, central, height, x):
    return height * width / (2 * np.pi) / ((x - central)**2 + width**2 / 4)


def error_func(p, x, y):
    return lorentzian(p, x) - y


def find_r_squared(f, p, x, y):
    res = y - f(p, x)
    ss_res = np.sum(res ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - ss_res / ss_tot


def compare_plot(x, y, p):
    fy = lorentzian(p, x)
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

    # prep

    wavelength_ev = 1248 / wavelength
    lorentz_ev = np.arange(np.max(wavelength_ev), np.min(wavelength_ev),
                           (np.min(wavelength_ev) - np.max(wavelength_ev)) / (len(wavelength_ev) * 7 - 1))
    lorentz_ev = np.append(lorentz_ev, 1248 / np.max(wavelength))
    lorentz_w1 = 1248 / lorentz_ev
    fit = np.zeros((len(lorentz_w1), 2))

    # find max and min

    max_sca = np.max(scattering[scattering < np.max(scattering)])
    idx_max = np.argmax(scattering[scattering < np.max(scattering)])
    min_sca = np.min(scattering)
    idx_min = np.argmin(scattering)

    # init guess and first fit

    init_1w = abs(2 / (np.pi * max_sca) * np.trapz(scattering, wavelength_ev))
    init_guess = [min_sca, min_sca * init_1w / (2 * np.pi), wavelength_ev[idx_max], init_1w]
    result, cov_x, res_dict, mesg, ier = leastsq(error_func, init_guess, args=(wavelength_ev, scattering),
                                                 full_output=True)
    #result, pcov = curve_fit(lorentzian2, wavelength_ev, scattering, p0=init_guess)
    result[3] = abs(result[3])

    # if fit is bad, retry

    if result[3] < 0.02:
        lb = abs(idx_max - round(len(scattering) / 4, 0)) - 1  # -1 for MATLAB
        hb = abs(round(3 * len(scattering) / 4, 0)) - 1  # -1 for MATLAB
        init_guess = [min_sca, min_sca * init_1w / (2 * np.pi), wavelength_ev[idx_max], init_1w]
        result, cov_x, res_dict, mesg, ier = leastsq(error_func, init_guess,
                                                     args=(wavelength_ev[lb:hb], scattering[lb:hb]))
        #result, pcov = curve_fit(lorentzian2, wavelength_ev[lb:hb], scattering[lb:hb], p0=init_guess)
        result[4] = abs(result[4])

    # calc r-sqaured. If bad, retry

    r_squared = find_r_squared(lambda p, x: lorentzian(p, x), result, wavelength_ev, scattering)

    if r_squared < 0.9:
        result, cov_x, res_dict, mesg, ier = leastsq(error_func, [-10, 100, 1248 / wavelength[idx_max], 0.15],
                                                     args=(wavelength_ev, scattering), full_output=True)
        #result, pcov = curve_fit(lorentzian2, wavelength_ev, scattering,
        #                         p0=[-10, 100, 1248 / wavelength[idx_max], 0.15])
        result[3] = abs(result[3])
        r_squared = find_r_squared(lambda p, x: lorentzian(p, x), result, wavelength_ev, scattering)

    fit[:, 0] = lorentz_w1
    fit[:, 1] = lorentzian(result, lorentz_ev)

    # split

    if r_squared < 0.9 and split is False:
        wavelength_low = wavelength[:idx_min]
        wavelength_high = wavelength[idx_min:]
        scattering_low = scattering[:idx_min]
        scattering_high = scattering[idx_min:]
        result_low, fit_low, r_squared_low = fit_lorentzian(scattering_low, wavelength_low, split=True)
        result_high, fit_high, r_squared_high = fit_lorentzian(scattering_high, wavelength_high, split=True)

        if r_squared_high > r_squared and ~np.isnan(np.sum(result_high)):
            result = result_high
            fit = fit_high
            r_squared = r_squared_high
        if r_squared_low > r_squared and ~np.isnan(np.sum(result_low)):
            result = result_low
            fit = fit_low
            r_squared = r_squared_low

    compare_plot(wavelength_ev, scattering, result)

    return result, fit, r_squared


params = [0, 200000, 640, 100]
wavelength = np.arange(570, 740, 10)
wavelength_ev = 1248 / wavelength
scattering = lorentzian(params, wavelength_ev)

fig, ax = plt.subplots(1)
ax.plot(wavelength_ev, scattering)
plt.show()

params2 = [50, 640, 100]
scattering2 = lorentzian_wavelength(*params2, wavelength)

fig, ax = plt.subplots(1)
ax.plot(wavelength, scattering2)
plt.show()

result, fit_, r_squared = fit_lorentzian(scattering, wavelength)

res_lambda = 1248 / result[2]
res_linewidth = 1000 * result[3]
