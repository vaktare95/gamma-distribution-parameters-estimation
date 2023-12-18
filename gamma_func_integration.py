from math import pow
from typing import Dict

import numpy as np


def integrated_function(var: float, alpha: float) -> float:
    try:
        return pow(var, (alpha - 1)) * np.exp(-var)
    except OverflowError:
        raise OverflowError(
            "Too large numbers occured. There was no convergence. Try with higher end_integration_threshold"
        )


def gamma_func_integration_trapezoidal_method(
    integrated_results_dict: Dict[float, float],
    alpha: float,
    integration_step: float = 1e-2,
    end_integration_threshold: float = 1e-5,
) -> None:
    integration_result = 0.0

    var = integration_step
    fun_value = integrated_function(var, alpha)
    while True:
        var_next = var + integration_step
        fun_value_next = integrated_function(var_next, alpha)
        trapezoidal_area = (fun_value + fun_value_next) * integration_step / 2
        integration_result += trapezoidal_area

        if fun_value_next < fun_value < end_integration_threshold:
            break

        var = var_next
        fun_value = fun_value_next

    integrated_results_dict[alpha] = integration_result
