import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf
from scipy.stats import norm
import pandas as pd


def error(phi, n, mae, eps=0.1):
    """
    Generates a random error vector of length `n`, that results
    in the mean absolute error (MAE) equal to `mae`.
    The error vector is based on the AR(1) model. The autocorrelation is
    controlled by the parameter `phi`. The first element of the output
    vector is 0 (no error). The higher `phi`, the higher chance of
    getting a vector with some bias.

    Since the output of the AR model is random, the function
    generates sample vectors until it finds a vector with MAE = mae.
    The standard deviation of the noise in the AR model is adapted
    to the expected MAE, so that that the number of 

    :parameter phi: float from 0 to 1, autocorrelation coefficient
    :parameter n: int, output vector length
    :parameter sd: float, standard deviation of the output vector
    :parameter mae: float, mean absolute error
    :parameter eps: float, standard deviation tolerance
    :return: 1D numpy array
    """
    # Initialize `mae_sample` with a high value,
    # so that the condition in `while` is not met
    mae_sample = mae + 2 * eps

    # Initialize the loop counter
    count = 0

    # Analytically derived standard deviation of the noise in AR(1),
    # which minimizes the number of iterations needed
    # to find the error vector for which MAE = mae:
    # (0.674 is the third quartile of the normal distribution)
    sdw = mae / (0.674 * (phi / np.sqrt(1 - phi ** 2) + 1))

    while abs(mae_sample - mae) >= eps:
        count += 1

        # Initialize error vector
        e = np.zeros(n)

        for i in range(2, n, 1):
            # Add next vector element using the AR(1) model
            e[i] = e[i-1] * phi + np.random.normal(loc=0, scale=sdw, size=1)[0]

        # Check tolerance
        mae_sample = np.mean(np.abs(e))

        # TODO: Use logging instead of print
        print("Error generation try number {} >> MAE = {}".format(count, mae_sample))

    return e


if __name__ == "__main__":

    # ==============================================
    # DEMO
    # ==============================================

    # Generate and plot one error vector and its ACF
    # ==============================================
    phi1 = 0.95
    phi2 = 0.05
    n = 100
    mae = 1

    # Results for phi1
    e1 = error(phi1, n, mae, eps=0.01)
    mae1 = np.mean(np.abs(e1))
    ac1 = acf(e1, nlags=n)
    ci1 = norm.ppf((1 + 0.95)/2) / np.sqrt(e1.size)  # 95% confidence interval

    # Results for phi2
    e2 = error(phi2, n, mae, eps=0.01)
    mae2 = np.mean(np.abs(e2))
    ac2 = acf(e2, nlags=n)
    ci2 = norm.ppf((1 + 0.95)/2) / np.sqrt(e2.size)  # 95% confidence interval

    fig, axes = plt.subplots(2, 2, sharex=True, figsize=(10, 7))
    fig.set_dpi(120)

    ax = axes[0][0]
    ax.plot(e1)
    ax.set_title(r"$\phi$ = {}, MAE = {:.2f}".format(phi1, mae1))
    ax.set_ylabel("Absolute error")
    ax.set_ylim(-5 * mae, 5 * mae)
    ax.plot(np.zeros(e1.size), 'k--', lw=1)

    ax = axes[1][0]
    ax.set_ylabel("ACF")
    ax.plot(ac1)
    ax.plot(np.full(e1.size, ci1), 'b--', lw=1)
    ax.plot(np.full(e1.size, -ci1), 'b--', lw=1)
    ax.plot(np.zeros(e1.size), 'k--', lw=1)
    ax.set_xlim(0, n)

    ax = axes[0][1]
    ax.plot(e2)
    ax.set_title(r"$\phi$ = {}, MAE = {:.2f}".format(phi2, mae2))
    ax.set_ylabel("Absolute error")
    ax.set_ylim(-5 * mae, 5 * mae)
    ax.plot(np.zeros(e2.size), 'k--', lw=1)

    ax = axes[1][1]
    ax.set_ylabel("ACF")
    ax.plot(ac2)
    ax.plot(np.full(e2.size, ci2), 'b--', lw=1)
    ax.plot(np.full(e2.size, -ci2), 'b--', lw=1)
    ax.plot(np.zeros(e2.size), 'k--', lw=1)
    ax.set_xlim(0, n)

    plt.show()

    # Test different correlation coefficients (phi)
