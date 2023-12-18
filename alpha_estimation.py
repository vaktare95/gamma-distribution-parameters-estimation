from math import pow

import numpy as np

from gamma_func_derivatives import gamma_func_derivatives_central_difference


def alpha_est_newton_raphson_method(
    X_vec: np.ndarray,
    alpha_start: float = 1.0,
    stop_threshold: float = 1e-3,
    derivative_step: float = 1e-6,
    integration_step: float = 1e-2,
    end_integration_threshold: float = 1e-5,
) -> float:
    constant = np.mean(np.log(X_vec)) - np.log(np.mean(X_vec))

    alpha = alpha_start
    obj_func = np.Inf
    while np.abs(obj_func) > stop_threshold:
        gamma_func, gamma_func_deri1, gamma_func_deri2 = gamma_func_derivatives_central_difference(
            alpha,
            derivative_step,
            integration_step,
            end_integration_threshold,
        )

        obj_func = gamma_func_deri1 / gamma_func - np.log(alpha) - constant
        obj_func_deri1 = (gamma_func_deri2 * gamma_func - pow(gamma_func_deri1, 2)) / pow(
            gamma_func, 2
        ) - 1 / alpha
        direction = -obj_func / obj_func_deri1
        alpha += direction

    return alpha
