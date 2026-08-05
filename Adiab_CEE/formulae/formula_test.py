#!/usr/bin/env python
# coding: utf-8

# In[2]:


from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd


# =============================================================================
# Internal coefficient reader
# =============================================================================

@lru_cache(maxsize=None)
def _read_lambda_ce_3d_coefficients(filename):
    """
    Read and cache one machine-readable 3D coefficient table.
    """
    coefficient_df = pd.read_csv(filename)

    i = coefficient_df[
        "log10_M_Msun_order"
    ].to_numpy(dtype=int)

    j = coefficient_df[
        "log10_R_Rsun_order"
    ].to_numpy(dtype=int)

    k = coefficient_df[
        "log10_eta_order"
    ].to_numpy(dtype=int)

    coefficients = coefficient_df[
        "coefficient"
    ].to_numpy(dtype=float)


    return coefficients, i, j, k


# =============================================================================
# polynomial evaluator
# =============================================================================

def _evaluate_lambda_ce_3d(
    mass,
    radius,
    eta,
    filename,
):
    mass, radius, eta = np.broadcast_arrays(
        np.asarray(mass, dtype=float),
        np.asarray(radius, dtype=float),
        np.asarray(eta, dtype=float),
    )

    coefficients, i_orders, j_orders, k_orders = (
        _read_lambda_ce_3d_coefficients(
            str(Path(filename).resolve())
        )
    )

    log_mass = np.log10(mass)
    log_radius = np.log10(radius)
    log_eta = np.log10(eta)

    log_lambda_ce = np.zeros(
        mass.shape,
        dtype=float,
    )

    for coefficient, i, j, k in zip(
        coefficients,
        i_orders,
        j_orders,
        k_orders,
    ):
        log_lambda_ce += (
            coefficient
            * log_mass**i
            * log_radius**j
            * log_eta**k
        )

    lambda_ce = 10.0**log_lambda_ce

    # Return a Python float for scalar input.
    if lambda_ce.ndim == 0:
        return float(lambda_ce)

    return lambda_ce


# =============================================================================
# FGB
# =============================================================================

def lambda_ce_3d_fgb(
    mass,
    radius,
    eta,
    coefficient_directory=".",
):
    """
    FGB 3D polynomial fit.

    Valid rectangular range:
        0.8 <= mass/Msun <= 10
        1   <= radius/Rsun <= 400
        1e-4 <= eta <= 0.099
    """
    filename = (
        Path(coefficient_directory)
        / "lambda_CE_3D_polycoef_FGB_degree_6.csv"
    )

    return _evaluate_lambda_ce_3d(
        mass=mass,
        radius=radius,
        eta=eta,
        filename=filename,
    )


# =============================================================================
# AGB
# =============================================================================

def lambda_ce_3d_agb(
    mass,
    radius,
    eta,
    coefficient_directory=".",
):
    """
    AGB 3D polynomial fit.

    Valid rectangular range:
        0.8 <= mass/Msun <= 10
        15  <= radius/Rsun <= 200
        1e-4 <= eta <= 0.099
    """
    filename = (
        Path(coefficient_directory)
        / "lambda_CE_3D_polycoef_AGB_degree_10.csv"
    )

    return _evaluate_lambda_ce_3d(
        mass=mass,
        radius=radius,
        eta=eta,
        filename=filename,
    )


# =============================================================================
# HG
# =============================================================================

def lambda_ce_3d_hg(
    mass,
    radius,
    eta,
    coefficient_directory=".",
):
    """
    HG 3D polynomial fit.

    Valid rectangular range:
        1 <= mass/Msun <= 100
        2 <= radius/Rsun <= 2000
        1e-4 <= eta <= 0.099
    """
    filename = (
        Path(coefficient_directory)
        / "lambda_CE_3D_polycoef_HG_degree_5.csv"
    )

    return _evaluate_lambda_ce_3d(
        mass=mass,
        radius=radius,
        eta=eta,
        filename=filename,
    )


# =============================================================================
# Massive stars
# =============================================================================

def lambda_ce_3d_massive(
    mass,
    radius,
    eta,
    coefficient_directory=".",
):
    """
    Massive star 3D polynomial fit.

    Valid range:
        8  <= mass/Msun <= 100
        10 <= radius/Rsun <= 2000
        1e-4 <= eta <= 0.099
    """
    filename = (
        Path(coefficient_directory)
        / "lambda_CE_3D_polycoef_Massive_degree_10.csv"
    )

    return _evaluate_lambda_ce_3d(
        mass=mass,
        radius=radius,
        eta=eta,
        filename=filename,
    )


# you can run the functions with something like this, they can take floats or arrays as input
#coefficient_directory = "REPLACE WITH YOUR PATH"
coefficient_directory = "lambda_CE_formula"

lambda_fgb = lambda_ce_3d_fgb(
    mass=2.0,
    radius=50.0,
    eta=1e-3,
    coefficient_directory=coefficient_directory,
)

lambda_ce_3d_agb(
    mass=2.0,
    radius=100.0,
    eta=1e-2,
    coefficient_directory=coefficient_directory,
)


# for 2D
# =============================================================================
# FGB fits
# Valid range: 0.8 <= M/Msun <= 10, 1 <= R/Rsun <= 400
# =============================================================================

def lambda_ce_fgb_eta_1e3(mass, radius):
    """FGB fit for eta = 1e-3."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        0.049
        + 0.122*y
        + 0.118*x
        - 0.239*y**2
        + 0.526*x*y
        - 0.652*x**2
        + 0.101*y**3
        + 0.480*x*y**2
        - 2.707*x**2*y
        - 0.084*x**3
        - 0.006*y**4
        - 0.087*x*y**3
        - 0.018*x**2*y**2
        + 1.397*x**3*y
        + 0.064*x**4
    )

    return 10.0**loglambda


def lambda_ce_fgb_eta_1e2(mass, radius):
    """FGB fit for eta = 1e-2."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        0.027
        + 0.083*y
        + 0.159*x
        + 0.404*y**2
        - 0.347*x*y
        - 0.409*x**2
        - 0.548*y**3
        + 1.196*x*y**2
        - 1.815*x**2*y
        - 0.607*x**3
        + 0.185*y**4
        - 0.426*x*y**3
        + 0.273*x**2*y**2
        + 0.397*x**3*y
        + 0.634*x**4
    )

    return 10.0**loglambda


# =============================================================================
# AGB fits
# Valid range: 0.8 <= M/Msun <= 10, 15 <= R/Rsun <= 200
# =============================================================================

def lambda_ce_agb_eta_1e3(mass, radius):
    """AGB fit for eta = 1e-3."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        18.499
        - 51.246*y
        + 4.940*x
        + 49.879*y**2
        - 6.052*x*y
        + 4.334*x**2
        - 21.235*y**3
        + 4.012*x*y**2
        - 1.627*x**2*y
        - 8.046*x**3
        + 3.350*y**4
        - 0.632*x*y**3
        - 1.107*x**2*y**2
        + 2.171*x**3*y
        + 2.597*x**4
    )

    return 10.0**loglambda


def lambda_ce_agb_eta_1e2(mass, radius):
    """AGB fit for eta = 1e-2."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        16.581
        - 46.274*y
        + 1.494*x
        + 45.115*y**2
        - 0.047*x*y
        + 4.261*x**2
        - 19.246*y**3
        + 2.676*x*y**2
        - 9.512*x**2*y
        + 0.046*x**3
        + 3.092*y**4
        - 1.314*x*y**3
        + 3.379*x**2*y**2
        - 1.335*x**3*y
        + 1.349*x**4
    )

    return 10.0**loglambda


# =============================================================================
# HG fits
# Valid range: 1 <= M/Msun <= 100, 2 <= R/Rsun <= 2000
# =============================================================================

def lambda_ce_hg_eta_1e3(mass, radius):
    """HG fit for eta = 1e-3."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        -0.238
        - 0.062*y
        + 0.312*x
        - 0.065*y**2
        - 0.284*x*y
        + 0.123*x**2
    )

    return 10.0**loglambda


def lambda_ce_hg_eta_1e2(mass, radius):
    """HG fit for eta = 1e-2."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        -0.267
        + 0.007*y
        + 0.295*x
        - 0.063*y**2
        - 0.317*x*y
        + 0.133*x**2
    )

    return 10.0**loglambda


# =============================================================================
# Massive-star fits
# Valid range: 8 <= M/Msun <= 100, 10 <= R/Rsun <= 2000
# =============================================================================

def lambda_ce_massive_eta_1e3(mass, radius):
    """Massive-star fit for eta = 1e-3."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        -0.793
        + 0.237*y
        + 1.487*x
        - 0.436*y**2
        + 0.201*x*y
        - 1.194*x**2
        + 0.092*y**3
        - 0.112*x*y**2
        - 0.009*x**2*y
        + 0.340*x**3
    )

    return 10.0**loglambda


def lambda_ce_massive_eta_1e2(mass, radius):
    """Massive-star fit for eta = 1e-2."""
    x = np.log10(np.asarray(mass, dtype=float))
    y = np.log10(np.asarray(radius, dtype=float))

    loglambda = (
        -1.227
        + 0.194*y
        + 2.466*x
        - 0.451*y**2
        + 0.340*x*y
        - 1.967*x**2
        + 0.117*y**3
        - 0.197*x*y**2
        + 0.039*x**2*y
        + 0.497*x**3
    )

    return 10.0**loglambda


# In[ ]:




