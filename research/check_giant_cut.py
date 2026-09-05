"""Finite enumeration of a critical cut with its terminal cluster partition.

Python stdlib only. Enumeration is exhaustive; numerical weights use floats.
Posterior: pi(s) proportional to exp(sum(a_e * satisfied_e(s))).
Edge status: 0 absent at 1, 1 early (<=t), 2 late (t,1].
"""

from collections import defaultdict
from itertools import product
from math import exp, expm1, fsum, log, pi as PI, prod, sin, tanh


TOL = 2e-12


def partition(n, edges, opened):
    parent = list(range(n))

    def root(v):
        while v != parent[v]:
            v = parent[v]
        return v

    for k in opened:
        i, j, _ = edges[k]
        parent[root(j)] = root(i)
    blocks = defaultdict(list)
    for v in range(n):
        blocks[root(v)].append(v)
    return tuple(sorted(map(tuple, blocks.values()), key=min))


def joint_law(n, edges, t):
    spins = tuple(product((-1, 1), repeat=n))
    weights = {s: exp(fsum(abs(w) for i, j, w in edges
                           if w * s[i] * s[j] > 0)) for s in spins}
    z = fsum(weights.values())
    pi = {s: v / z for s, v in weights.items()}
    joint = {}
    for s in spins:
        categories = []
        for i, j, w in edges:
            a = abs(w)
            satisfied = int(w * s[i] * s[j] > 0)
            ps = (exp(-a * satisfied),
                  (1 - exp(-t * a)) * satisfied,
                  (exp(-t * a) - exp(-a)) * satisfied)
            assert abs(fsum(ps) - 1) < TOL
            categories.append(tuple((x, p) for x, p in enumerate(ps) if p > 0))
        for choice in product(*categories):
            status = tuple(x for x, _ in choice)
            mass = pi[s]
            for _, p in choice:
                mass *= p
            ac = tuple(k for k, x in enumerate(status) if x == 1)
            a1 = tuple(k for k, x in enumerate(status) if x)
            joint[s, ac, a1] = mass
    assert abs(fsum(joint.values()) - 1) < TOL
    return pi, joint


def check(n, edges, t):
    pi, joint = joint_law(n, edges, t)
    residual = defaultdict(lambda: defaultdict(float))
    grouped = defaultdict(lambda: defaultdict(float))
    marked = defaultdict(lambda: defaultdict(float))
    for (s, ac, a1), mass in joint.items():
        assert set(ac).issubset(a1)
        residual[s, ac][a1] += mass
        grouped[ac, partition(n, edges, a1)][s] += mass
        marked[a1][s] += mass

    # Independent connectivity weight: sum all possible residual open sets B.
    for ac in {ac for ac, _ in grouped}:
        raw = defaultdict(lambda: defaultdict(float))
        rest = tuple(k for k in range(len(edges)) if k not in ac)
        for s in pi:
            satisfied = tuple(w * s[i] * s[j] > 0 for i, j, w in edges)
            if not all(satisfied[k] for k in ac):
                continue
            for flags in product((False, True), repeat=len(rest)):
                b = tuple(k for k, yes in zip(rest, flags) if yes)
                weight = prod(expm1((1 - t) * abs(edges[k][2])) * satisfied[k]
                              for k in b)
                if weight:
                    raw[partition(n, edges, ac + b)][s] += weight
        for r, weights in raw.items():
            row = grouped[ac, r]
            z, z_raw = fsum(row.values()), fsum(weights.values())
            assert z > 0
            assert max(abs(row.get(s, 0) / z - weights.get(s, 0) / z_raw)
                       for s in pi) < TOL
        assert set(raw) == {r for early, r in grouped if early == ac}

    # Conditional continuation from Ac: independent residual clocks.
    for (s, ac), row in residual.items():
        z = fsum(row.values())
        for a1, mass in row.items():
            predicted = 1.0
            for k, (i, j, w) in enumerate(edges):
                if k in ac:
                    continue
                p = ((1 - exp(-(1 - t) * abs(w)))
                     if w * s[i] * s[j] > 0 else 0.0)
                predicted *= p if k in a1 else 1 - p
            assert abs(mass / z - predicted) < TOL

    pairs = tuple(product(range(n), repeat=2))
    m = {ij: fsum(p * s[ij[0]] * s[ij[1]] for s, p in pi.items())
         for ij in pairs}
    tower = dict.fromkeys(pairs, 0.0)
    localized = dict.fromkeys(pairs, 0.0)
    small_mass = 0.0
    for (ac, r), row in grouped.items():
        z = fsum(row.values())
        early = partition(n, edges, ac)
        root_at = {i: c for c in r for i in c}
        early_at = {i: c for c in early for i in c}
        giant = min(r, key=lambda c: (-len(c), min(c)))
        small_mass += z * (sum(len(c) ** 2 for c in early)
                           + sum(len(c) ** 2 for c in r if c != giant)) / n ** 2
        for i, j in pairs:
            cm = fsum(p * s[i] * s[j] for s, p in row.items()) / z
            if root_at[i] != root_at[j]:
                assert abs(cm) < TOL
            tower[i, j] += z * cm
            if i in giant and j in giant and early_at[i] != early_at[j]:
                localized[i, j] += z * cm
    assert max(abs(tower[ij] - m[ij]) for ij in pairs) < TOL
    q = fsum(x * x for x in m.values()) / n ** 2
    q_giant = fsum(x * x for x in localized.values()) / n ** 2
    assert abs(q - q_giant) <= 2 * small_mass + TOL

    # Terminal *marked edges* fix all relative spins inside their components.
    for a1, row in marked.items():
        for block in partition(n, edges, a1):
            for i, j in product(block, repeat=2):
                assert len({s[i] * s[j] for s in row}) == 1
    return q, q_giant, small_mass


def counterexamples():
    # Keeping the terminal root modifies the early-cut Gibbs conditional.
    a, t = 1.1, 0.37
    _, joint = joint_law(2, [(0, 1, a)], t)
    cut = {s: mass for (s, ac, a1), mass in joint.items()
           if not ac and a1}
    z = fsum(cut.values())
    root_correlation = fsum(mass * s[0] * s[1] for s, mass in cut.items()) / z
    ordinary_correlation = tanh((1 - t) * a / 2)
    assert abs(root_correlation - 1) < TOL
    assert ordinary_correlation < 0.5
    print(f"PASS kept root changes conditional: correlation {ordinary_correlation:.6f} -> 1")

    # An unmarked connected root need not fix its internal relative spins.
    edges = [(0, 1, 0.7), (1, 2, 1.1), (0, 2, -0.9)]
    _, joint = joint_law(3, edges, t)
    products = {s[0] * s[1] for (s, ac, a1), mass in joint.items()
                if not ac and len(partition(3, edges, a1)) == 1 and mass > 0}
    assert products == {-1, 1}
    print("PASS unmarked root can retain uncertainty; marked terminal edges freeze it")


def channel_check():
    for p in (0.55, 0.7, 0.81, 0.95):
        for t in (0.0, 0.37, 1.0):
            a = log(p / (1 - p))
            q = p * (1 - exp(-t * a))
            # Channel law given x=+1: (y,b) has these three probabilities.
            channel_plus = {(1, 1): q, (1, 0): p - q, (-1, 0): 1 - p}
            channel_minus = {(-y, b): prob for (y, b), prob in channel_plus.items()}
            eta = 0.0
            for obs in channel_plus.keys() | channel_minus.keys():
                plus = channel_plus.get(obs, 0.0)
                minus = channel_minus.get(obs, 0.0)
                if plus + minus:
                    eta += (plus - minus) ** 2 / (2 * (plus + minus))
            formula = q + (2 * p - 1 - q) ** 2 / (1 - q)
            gain = 4 * q * (1 - p) ** 2 / (1 - q)
            assert abs(eta - formula) < TOL
            assert abs(eta - (2 * p - 1) ** 2 - gain) < TOL
            assert gain >= -TOL
    print("PASS revealing the early cut increases the channel information parameter")


def two_block_check():
    p = 0.81
    a = log(p / (1 - p))
    qc = 2 * sin(PI / 18)
    t = -log(1 - qc / p) / a
    h = (1 - t) * a
    eta = qc + (2 * p - 1 - qc) ** 2 / (1 - qc)
    assert eta > qc
    print(f"p=.81: critical cut={t:.9f}, residual weight={h:.9f}, "
          f"information bound={eta:.9f} > qc={qc:.9f}")
    for positive, negative, expected in ((1, 1, 0), (2, 0, 1), (2, 1, 0.54907042)):
        edges = [(0, 1, a)] * positive + [(0, 1, -a)] * negative
        _, joint = joint_law(2, edges, t)
        masses = defaultdict(float)
        for (s, ac, a1), mass in joint.items():
            if not ac and a1:
                masses[s[0] * s[1]] += mass
        observed = (masses[1] - masses[-1]) / fsum(masses.values())
        plus, minus = expm1(h * positive), expm1(h * negative)
        rho = (plus - minus) / (plus + minus)
        assert abs(observed - rho) < TOL
        assert abs(rho - expected) < 5e-9
        print(f"PASS connected two-block motif ({positive},{negative}): rho={rho:.8f}")


def main():
    cases = [
        ("signed tree", 3, [(0, 1, 0.7), (1, 2, -1.1)]),
        ("frustrated triangle", 3, [(0, 1, 0.7), (1, 2, 1.1), (0, 2, -0.9)]),
        ("signed K4", 4, [(0, 1, 0.7), (1, 2, -1.1), (0, 2, 0.9),
                          (2, 3, 0.4), (0, 3, -1.3), (1, 3, -0.8)]),
    ]
    for name, n, edges in cases:
        for t in (0.0, 0.37, 1.0):
            q, q_giant, error_mass = check(n, edges, t)
            print(f"PASS {name:19} t={t:.2f}: Q={q:.6f}, giant={q_giant:.6f}, "
                  f"error bound={2 * error_mass:.6f}")
    counterexamples()
    channel_check()
    two_block_check()


if __name__ == "__main__":
    main()
