"""Direct measurement of the replicated target D_L^x on the triangular torus.

The module is dependency-free.  It implements, for the homogeneous binary
GSBM on the L x L triangular torus, the first direct estimator of the
replicated cross-block target of note 41,

``D_L^x = (1/n^2) E_O sum_{i,j} d_{O,ij}^2``,
``d_{O,ij} = E_{D|O}[ 1{i,j in R*(D)} 1{A_D(i) != A_D(j)} pi_{O,D}(s_i s_j) ]``,

with an EXACT inner layer and a Monte-Carlo outer layer:

* the posterior mu_O, its correlations <s_i s_j>_O (hence Q_L given O)
  and the full-tree conditional Gibbs pi_{O,D} are computed EXACTLY by
  enumeration of the 2^n spin configurations (n = 16 at L = 4);
* the observation O (Nishimori gauge, planted all-plus) and the
  dendrogram D (posterior replica + exponential clocks, unmarked
  Kruskal) are sampled Monte-Carlo, with documented seeds;
* the square of the D-average is estimated WITHOUT bias by the
  U-statistic (sum^2 - sum of squares) / (K(K-1)) over K independent
  dendrograms sharing the same observation.

Reported quantities per parameter p:

* ``q_full``     -- Q_L = (1/n^2) E_O sum_{ij} <s_i s_j>_O^2 (exact in O);
* ``d_cross``    -- D_L^x as above (average over D BEFORE the square);
* ``j_cross``    -- the single-D Jensen envelope
  (1/n^2) E_O sum_{ij} E_D[X_{ij}^2] (square BEFORE the average);
* ``s_crit``     -- quadratic mass (1/n^2) E sum_{A critical block in R*} |A|^2.

The ratio d_cross / j_cross quantifies, at this small volume, how much
the dendrogram average cancels compared to its Jensen envelope.  No
asymptotic claim is made: L = 4 is far below any scaling regime.
"""

from __future__ import annotations

import itertools
import math
import random
from bisect import bisect_left

Q_C = 2.0 * math.sin(math.pi / 18.0)


def coupling(p: float) -> float:
    if not 0.5 < p < 1.0:
        raise ValueError("p must satisfy 1/2 < p < 1")
    return math.log(p / (1.0 - p))


def beta_c(p: float) -> float:
    """Critical cut level: q_p(beta_c) = q_c (requires p >= (1+q_c)/2)."""

    if p < (1.0 + Q_C) / 2.0:
        raise ValueError("beta_c undefined below p_SW")
    return -math.log(1.0 - Q_C / p) / coupling(p)


def triangular_torus(size: int):
    """Vertices 0..size^2-1 and the 3 s edge directions of the torus."""

    def idx(x: int, y: int) -> int:
        return (x % size) * size + (y % size)

    edges = []
    for x in range(size):
        for y in range(size):
            edges.append((idx(x, y), idx(x + 1, y)))
            edges.append((idx(x, y), idx(x, y + 1)))
            edges.append((idx(x, y), idx(x + 1, y + 1)))
    if size < 3:
        raise ValueError("size >= 3 required to avoid duplicate edges")
    return size * size, edges


class DisjointSet:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, a: int) -> int:
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        self.parent[ra] = rb
        return True


class GsbmExact:
    """Exact inner layer on the L x L triangular torus (2^n enumeration)."""

    def __init__(self, size: int, p: float, rng: random.Random) -> None:
        self.n, self.edges = triangular_torus(size)
        self.p = p
        self.u = coupling(p)
        self.bc = beta_c(p)
        self.rng = rng
        self.nconf = 1 << self.n
        # parity tables: par[e][sigma] = 1 iff sigma_i == sigma_j
        self.par = []
        for (i, j) in self.edges:
            mask_i, mask_j = 1 << i, 1 << j
            tab = bytearray(self.nconf)
            for s in range(self.nconf):
                tab[s] = 1 if (bool(s & mask_i) == bool(s & mask_j)) else 0
            self.par.append(bytes(tab))
        self._flip = bytes(1 - b for b in range(2)) + bytes(254)

    # ---------------- observation and exact posterior ----------------
    def sample_observation(self):
        """Nishimori gauge (planted all-plus): O_e = +1 w.p. p."""

        return [1 if self.rng.random() < self.p else -1 for _ in self.edges]

    def satisfied_tables(self, obs):
        """sat[e][sigma] = 1 iff edge e is satisfied by sigma under obs."""

        return [
            self.par[e] if obs[e] == 1 else self.par[e].translate(self._flip)
            for e in range(len(self.edges))
        ]

    def posterior_weights(self, sat):
        score = [0] * self.nconf
        for tab in sat:
            for s in range(self.nconf):
                score[s] += tab[s]
        u = self.u
        return [math.exp(u * sc) for sc in score]

    def pair_correlation(self, weights, i: int, j: int) -> float:
        mask_i, mask_j = 1 << i, 1 << j
        num = 0.0
        tot = 0.0
        for s, w in enumerate(weights):
            tot += w
            num += w if (bool(s & mask_i) == bool(s & mask_j)) else -w
        return num / tot

    def all_pair_correlations(self, weights):
        n = self.n
        corr = {}
        tot = sum(weights)
        for i in range(n):
            for j in range(i + 1, n):
                mask_i, mask_j = 1 << i, 1 << j
                num = 0.0
                for s, w in enumerate(weights):
                    num += w if (bool(s & mask_i) == bool(s & mask_j)) else -w
                corr[(i, j)] = num / tot
        return corr

    def sample_configuration(self, weights) -> int:
        cum = list(itertools.accumulate(weights))
        x = self.rng.random() * cum[-1]
        return bisect_left(cum, x)

    # ---------------- dendrogram (unmarked Kruskal) ----------------
    def sample_dendrogram(self, sigma: int, sat):
        """Clocks on satisfied edges; returns (merges, blocks, roots).

        merges: list of (beta_v, frozenset cut_edge_indices) for rings <= 1;
        blocks: DisjointSet at level beta_c; roots: DisjointSet at level 1.
        """

        rings = []
        for e, tab in enumerate(sat):
            if tab[sigma]:
                t = self.rng.expovariate(self.u)
                if t <= 1.0:
                    rings.append((t, e))
        rings.sort()
        ds = DisjointSet(self.n)
        components = {v: {v} for v in range(self.n)}
        merges = []
        blocks = DisjointSet(self.n)
        for (t, e) in rings:
            i, j = self.edges[e]
            ri, rj = ds.find(i), ds.find(j)
            if ri == rj:
                continue
            side_a, side_b = components[ri], components[rj]
            cut = frozenset(
                k
                for k, (x, y) in enumerate(self.edges)
                if (x in side_a and y in side_b) or (x in side_b and y in side_a)
            )
            merges.append((t, cut))
            ds.union(ri, rj)
            merged = side_a | side_b
            del components[ri], components[rj]
            components[ds.find(i)] = merged
            if t <= self.bc:
                blocks.union(i, j)
        return merges, blocks, ds

    # ---------------- exact full-tree conditional Gibbs ----------------
    def gibbs_given_dendrogram(self, merges, sat):
        """Return unnormalised nu[sigma] = prod_v K_v exp((1-b_v) u K_v)."""

        nu = [1.0] * self.nconf
        for (t, cut) in merges:
            cut_list = sorted(cut)
            counts = list(sat[cut_list[0]])
            for e in cut_list[1:]:
                tab = sat[e]
                counts = [c + tab[s] for s, c in enumerate(counts)]
            coef = (1.0 - t) * self.u
            factor = [k * math.exp(coef * k) for k in range(len(cut_list) + 1)]
            nu = [nv * factor[c] for nv, c in zip(nu, counts)]
        return nu

    def gibbs_pair_correlation(self, nu, i: int, j: int) -> float:
        return self.pair_correlation(nu, i, j)


def measure(size, p, n_obs, k_dendro, seed):
    """Run the exact-inner / MC-outer measurement; returns a dict."""

    rng = random.Random(seed)
    model = GsbmExact(size, p, rng)
    n = model.n
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    q_acc, d_acc, j_acc, s_acc = 0.0, 0.0, 0.0, 0.0
    for _ in range(n_obs):
        obs = model.sample_observation()
        sat = model.satisfied_tables(obs)
        weights = model.posterior_weights(sat)
        corr = model.all_pair_correlations(weights)
        q_obs = (n + 2.0 * sum(c * c for c in corr.values())) / (n * n)
        draws = {pair: [] for pair in pairs}
        s_local = 0.0
        for _ in range(k_dendro):
            sigma = model.sample_configuration(weights)
            merges, blocks, roots = model.sample_dendrogram(sigma, sat)
            root_sizes = {}
            for v in range(n):
                r = roots.find(v)
                root_sizes[r] = root_sizes.get(r, 0) + 1
            giant = max(root_sizes, key=root_sizes.get)
            block_sizes = {}
            for v in range(n):
                if roots.find(v) == giant:
                    b = blocks.find(v)
                    block_sizes[b] = block_sizes.get(b, 0) + 1
            s_local += sum(sz * sz for sz in block_sizes.values()) / (n * n)
            cross = [
                (i, j)
                for (i, j) in pairs
                if roots.find(i) == giant
                and roots.find(j) == giant
                and blocks.find(i) != blocks.find(j)
            ]
            nu = None
            values = {}
            if cross:
                nu = model.gibbs_given_dendrogram(merges, sat)
                for (i, j) in cross:
                    values[(i, j)] = model.gibbs_pair_correlation(nu, i, j)
            for pair in pairs:
                draws[pair].append(values.get(pair, 0.0))
        d_obs, j_obs = 0.0, 0.0
        for pair in pairs:
            xs = draws[pair]
            k = len(xs)
            total = sum(xs)
            square = sum(x * x for x in xs)
            d_obs += 2.0 * (total * total - square) / (k * (k - 1))
            j_obs += 2.0 * square / k
        q_acc += q_obs
        d_acc += d_obs / (n * n)
        j_acc += j_obs / (n * n)
        s_acc += s_local / k_dendro
    return {
        "size": size,
        "p": p,
        "beta_c": model.bc,
        "n_obs": n_obs,
        "k_dendro": k_dendro,
        "seed": seed,
        "q_full": q_acc / n_obs,
        "d_cross": d_acc / n_obs,
        "j_cross": j_acc / n_obs,
        "s_crit": s_acc / n_obs,
    }


def main() -> None:
    for p in (0.75, 0.81):
        result = measure(size=4, p=p, n_obs=8, k_dendro=6, seed=20260805)
        print(
            "L={size} p={p} beta_c={beta_c:.6f} seed={seed} "
            "n_obs={n_obs} K={k_dendro}".format(**result)
        )
        print(
            "  Q_L={q_full:.6f}  D_x={d_cross:.6f}  "
            "J_x={j_cross:.6f}  S_c={s_crit:.6f}  ratio D/J={ratio:.4f}".format(
                ratio=result["d_cross"] / result["j_cross"]
                if result["j_cross"]
                else float("nan"),
                **result,
            )
        )


if __name__ == "__main__":
    main()
