"""Checks the explicit triangular-lattice critical-band formulas."""

from __future__ import annotations

from math import log, pi, sin, sqrt


Q_CRITICAL = 2.0 * sin(pi / 18.0)


def coupling(p: float) -> float:
    if not 0.5 < p < 1.0:
        raise ValueError("p must satisfy 0.5 < p < 1")
    return log(p / (1.0 - p))


def open_probability(p: float, t: float) -> float:
    return p * (1.0 - (1.0 - p) ** t / p**t)


def inverse_time(p: float, q: float) -> float:
    if not 0.0 <= q <= 2.0 * p - 1.0:
        raise ValueError("q must belong to the range of q_p on [0, 1]")
    return -log(1.0 - q / p) / coupling(p)


def beta_critical(p: float) -> float:
    return inverse_time(p, Q_CRITICAL)


def information_time(p: float) -> float:
    return inverse_time(p, (2.0 * p - 1.0) ** 2)


def main() -> None:
    p_sw = (1.0 + Q_CRITICAL) / 2.0
    p_info = (1.0 + sqrt(Q_CRITICAL)) / 2.0
    p_pure = 0.5 + Q_CRITICAL

    beta_at_info = beta_critical(p_info)
    t_chi_at_info = information_time(p_info)

    assert abs(open_probability(p_sw, 1.0) - Q_CRITICAL) < 1e-12
    assert abs(beta_at_info - t_chi_at_info) < 1e-12
    assert abs(2.0 * p_pure - 1.0 - Q_CRITICAL - Q_CRITICAL) < 1e-12

    print(f"q_c       = {Q_CRITICAL:.12f}")
    print(f"p_SW      = {p_sw:.12f}")
    print(f"p_info    = {p_info:.12f}")
    print(f"p_pure    = {p_pure:.12f}")
    print(f"beta_c(p_info) = {beta_at_info:.12f}")
    print(f"t_chi(p_info)  = {t_chi_at_info:.12f}")

    print("\nreference points")
    print("p              beta_c        t_chi")
    for p in (p_info, 0.8358058, p_pure):
        print(f"{p:.10f}   {beta_critical(p):.10f}   {information_time(p):.10f}")


if __name__ == "__main__":
    main()
