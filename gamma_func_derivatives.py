import multiprocessing as mp
from typing import Tuple

from gamma_func_integration import gamma_func_integration_trapezoidal_method


def gamma_func_derivatives_central_difference(
    alpha: float,
    derivative_step: float = 1e-6,
    integration_step: float = 1e-2,
    end_integration_threshold: float = 1e-5,
) -> Tuple[float, float, float]:
    alphas = [alpha + i * derivative_step for i in range(-2, 3)]

    with mp.Manager() as manager:
        integration_results = manager.dict()
        arguments_list = [
            (integration_results, alpha, integration_step, end_integration_threshold)
            for alpha in alphas
        ]
        with manager.Pool(processes=mp.cpu_count()) as pool:
            pool.starmap(gamma_func_integration_trapezoidal_method, arguments_list)
            integration_results = dict(integration_results)

    d = derivative_step
    derivative_1 = (integration_results[alpha + d] - integration_results[alpha - d]) / (2 * d)
    derivative_1_next = (integration_results[alpha + 2 * d] - integration_results[alpha]) / (2 * d)
    derivative_1_prev = (integration_results[alpha] - integration_results[alpha - 2 * d]) / (2 * d)
    derivative_2 = (derivative_1_next - derivative_1_prev) / (2 * d)
    return integration_results[alpha], derivative_1, derivative_2
