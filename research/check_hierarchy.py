"""Exact finite-state checks for the partial Edwards-Sokal cluster sweep.

Only the Python standard library is used. W is defined so that
    pi(s) propto exp(sum(abs(W_e) * 1[W_e*s_i*s_j > 0])).
Thus the ordinary Ising coefficient is W_e / 2.
The enumeration is exact over states/bonds; real arithmetic uses floats.
"""

from itertools import product
from math import exp, fsum


def states(n):
    return list(product((-1, 1), repeat=n))


def energy(s, edges):
    return fsum(abs(w) for i, j, w in edges if w * s[i] * s[j] > 0)


def target(n, edges):
    weights = {s: exp(energy(s, edges)) for s in states(n)}
    z = fsum(weights.values())
    return {s: v / z for s, v in weights.items()}


def bonds_given(s, edges, t):
    choices = [((), 1.0)]
    for k, (i, j, w) in enumerate(edges):
        p = 1.0 - exp(-t * abs(w)) if w * s[i] * s[j] > 0 else 0.0
        new = []
        for a, q in choices:
            if p < 1:
                new.append((a, q * (1.0 - p)))
            if p > 0:
                new.append((a + (k,), q * p))
        choices = new
    return choices


def components(n, edges, a):
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for k in a:
        i, j, _ = edges[k]
        parent[root(j)] = root(i)
    blocks = {}
    for i in range(n):
        blocks.setdefault(root(i), []).append(i)
    return sorted(blocks.values(), key=min)


def flip(s, block):
    out = list(s)
    for i in block:
        out[i] *= -1
    return tuple(out)


def orientation_sweep(s, edges, a, t, *, naive_uniform=False):
    # A deterministic vertex-based ordering is independent of spin values.
    residual = [e for k, e in enumerate(edges) if k not in a]
    law = {s: 1.0}
    for block in components(len(s), edges, a):
        new = {}
        for u, mass in law.items():
            v = flip(u, block)
            p = 0.5 if naive_uniform else 1.0 / (
                1.0 + exp((1.0 - t) * (energy(v, residual) - energy(u, residual)))
            )
            new[u] = new.get(u, 0.0) + mass * p
            new[v] = new.get(v, 0.0) + mass * (1.0 - p)
        law = new
    return law


def kernel(n, edges, t, *, naive_uniform=False):
    k = {}
    for s in states(n):
        row = dict.fromkeys(states(n), 0.0)
        ts = t(s) if callable(t) else t
        for a, p in bonds_given(s, edges, ts):
            for u, q in orientation_sweep(s, edges, a, ts,
                                          naive_uniform=naive_uniform).items():
                row[u] += p * q
        assert abs(fsum(row.values()) - 1.0) < 2e-12
        k[s] = row
    return k


def stationarity_error(pi, k):
    return max(abs(fsum(pi[s] * k[s][u] for s in pi) - pi[u]) for u in pi)


def matrix_error(k, other):
    return max(abs(k[s][u] - other[s][u]) for s in k for u in k[s])


def single_site_sweep(n, edges):
    # Independent reference: assign +/-1 using the local Ising field.
    k = {}
    for s in states(n):
        law = {s: 1.0}
        for i in range(n):
            new = {}
            for u, mass in law.items():
                field = fsum(w * u[b if a == i else a]
                             for a, b, w in edges if i in (a, b))
                plus = 1.0 / (1.0 + exp(-field))
                for spin, p in ((1, plus), (-1, 1.0 - plus)):
                    v = u[:i] + (spin,) + u[i + 1:]
                    new[v] = new.get(v, 0.0) + mass * p
            law = new
        k[s] = law
    return k


def check_conditionals_and_recovery(n, edges, t):
    pi = target(n, edges)
    joint = {}
    for s in pi:
        for a, p in bonds_given(s, edges, t):
            joint.setdefault(a, {})[s] = pi[s] * p
    d = {(i, j): 0.0 for i in range(n) for j in range(n)}
    susceptibility = error = 0.0
    for a, row in joint.items():
        mass = fsum(row.values())
        bayes = {s: row.get(s, 0.0) / mass for s in pi}
        # Independently evaluate the residual law in signed Ising form.
        weights = {
            s: (exp(0.5 * (1.0 - t) * fsum(
                w * s[i] * s[j] for k, (i, j, w) in enumerate(edges) if k not in a
            )) if all(edges[k][2] * s[edges[k][0]] * s[edges[k][1]] > 0
                      for k in a) else 0.0)
            for s in pi
        }
        z = fsum(weights.values())
        error = max(error, max(abs(bayes[s] - weights[s] / z) for s in pi))
        blocks = components(n, edges, a)
        label = {i: k for k, block in enumerate(blocks) for i in block}
        susceptibility += mass * sum(len(block) ** 2 for block in blocks) / n ** 2
        for i, j in d:
            if label[i] != label[j]:
                d[i, j] += mass * fsum(bayes[s] * s[i] * s[j] for s in pi)
    q = fsum(fsum(pi[s] * s[i] * s[j] for s in pi) ** 2 for i, j in d) / n ** 2
    d2 = fsum(value ** 2 for value in d.values()) / n ** 2
    # At fixed observations; averaging this inequality gives the stated bound.
    assert abs(q - d2) <= 2 * susceptibility + 2e-12
    assert error < 2e-12
    return error


def sw_kernel(n, edges):
    k = {}
    for s in states(n):
        row = dict.fromkeys(states(n), 0.0)
        for a, p in bonds_given(s, edges, 1.0):
            blocks = components(n, edges, a)
            for flips in product((False, True), repeat=len(blocks)):
                u = s
                for block, yes in zip(blocks, flips):
                    if yes:
                        u = flip(u, block)
                row[u] += p / (2 ** len(blocks))
        k[s] = row
    return k


def detailed_balance_error(pi, k):
    return max(abs(pi[s] * k[s][u] - pi[u] * k[u][s]) for s in pi for u in pi)


def main():
    cases = [
        ("2 vertices", 2, [(0, 1, 1.1)]),
        ("signed tree", 3, [(0, 1, 0.7), (1, 2, -1.1)]),
        ("frustrated triangle", 3, [(0, 1, 0.7), (1, 2, 1.1), (0, 2, -0.9)]),
        ("signed complete graph", 4, [(0, 1, 0.7), (1, 2, -1.1),
                                    (0, 2, 0.9), (2, 3, 0.4),
                                    (0, 3, -1.3), (1, 3, -0.8)]),
    ]
    for name, n, edges in cases:
        pi = target(n, edges)
        for t in (0.0, 0.37, 1.0):
            k = kernel(n, edges, t)
            error = stationarity_error(pi, k)
            assert error < 2e-12, (name, t, error)
            conditional_error = check_conditionals_and_recovery(n, edges, t)
            if t == 0:
                assert matrix_error(k, single_site_sweep(n, edges)) < 2e-12
            if t == 1:
                assert matrix_error(k, sw_kernel(n, edges)) < 2e-12
            print(f"PASS {name:23s} t={t:4.2f} stationarity={error:.2e}, "
                  f"conditional={conditional_error:.2e}, recovery bound OK")

    edges = cases[2][2]
    pi = target(3, edges)
    naive_error = stationarity_error(pi, kernel(3, edges, 0.37, naive_uniform=True))
    adaptive_error = stationarity_error(pi, kernel(3, edges, lambda s: 1.0 if s[0] == s[1] else 0.0))
    balance_error = detailed_balance_error(pi, kernel(3, edges, 0.37))
    assert naive_error > 1e-5
    assert adaptive_error > 1e-5
    assert balance_error > 1e-5
    print(f"EXPECTED FAILURE: uniform flips at t=.37, stationarity error={naive_error:.6g}")
    print(f"EXPECTED FAILURE: spin-dependent t, stationarity error={adaptive_error:.6g}")
    print(f"NOTE: ordered sweep invariant but not reversible; detailed balance error={balance_error:.6g}")

    # One positive edge, and D_1 records a merger at b <= 1: spins must agree.
    # Holding this hierarchy fixed forces either leaf conditional to a point.
    # Ordinary Glauber instead creates unequal spins with positive probability.
    w = 1.1
    b = 0.4
    pi_edge = target(2, [(0, 1, w)])
    conditional_weights = {
        x: pi_edge[(x, 1)] * (w * exp(-b * w) if x == 1 else 0.0)
        for x in (-1, 1)
    }
    fixed_hierarchy_leaf_flip = conditional_weights[-1] / fsum(conditional_weights.values())
    ordinary_glauber_leaf_flip = 1.0 / (1.0 + exp(w))
    assert ordinary_glauber_leaf_flip > fixed_hierarchy_leaf_flip + 0.1
    print("EXPECTED FAILURE: leaf update with D_1 fixed has flip probability 0; "
          f"ordinary Glauber gives {ordinary_glauber_leaf_flip:.6g}")


if __name__ == "__main__":
    main()
