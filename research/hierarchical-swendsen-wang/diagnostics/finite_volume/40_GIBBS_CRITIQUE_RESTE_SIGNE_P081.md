# Diagnostic fini à $`p=0.81`$ : enveloppe spectrale et reste signé

**Statut : diagnostic reproductible de volume fini ; le single-$D$ énumère
exactement les orientations collectives conditionnellement aux formes
internes fixées, tandis que le double-$D$ énumère les Gibbs conditionnels
complets à $`L=4`$ ; environnements, observations et hiérarchies sont
échantillonnés ; aucune extrapolation thermodynamique et aucune nouvelle
borne de weak recovery.**

Cette note exécute les deux premières portes quantitatives de la
[cible double-géante](../../active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) :

1. tester l'enveloppe plus forte mais plus simple à un dendrogramme ;
2. si elle reste macroscopique, mesurer directement le reste inter-cellules
   **signé** avec deux dendrogrammes indépendants.

Les résultats à $`L=4`$ sont nets :

- l'enveloppe single-$D$ est proche de un et ne ferme pas la preuve ;
- le reste à deux $D$ est exactement décomposé et peut être négatif
  réalisation par réalisation ;
- sa moyenne observée est toutefois positive et macroscopique sur ce très
  petit volume.

Le diagnostic motive donc le passage à la loi jointe de deux corridors,
mais ne fournit encore aucun indice numérique d'une annulation asymptotique.

## 1. Enveloppe spectrale à un dendrogramme

Pour une racine finale $R$, soient $`C_1,\ldots,C_k`$ ses blocs critiques,

```math
M_R(a,b)
=
\mathbb E_{\pi_{D,R}^c(\cdot\mid y)}[z_az_b],
\qquad
W_R
=
\mathrm{diag}(|C_1|,\ldots,|C_k|),
\qquad
A_R
=
W_R^{1/2}M_RW_R^{1/2}.
\qquad\text{(1.1)}
```

Comme $`M_R=\mathbb E[zz^{\mathsf T}]`$, la matrice $`A_R`$ est
semi-définie positive et

```math
\mathrm{tr}(A_R)
=
|R|,
\qquad
\mathrm{tr}(A_R^2)
=
\sum_{a,b}
|C_a||C_b|M_R(a,b)^2.
\qquad\text{(1.2)}
```

Si $`R_\star`$ est la plus grande racine, la quantité testée est

```math
\Lambda_L^{(1)}
=
\frac{
\lambda_{\max}(A_{R_\star})
}{
n_L
}.
\qquad\text{(1.3)}
```

Sa décroissance vers zéro suffirait à faire décroître la persistance
quadratique de la racine géante. Plus précisément, les inégalités
$`(\Lambda_L^{(1)})^2\le
\mathrm{tr}(A_{R_\star}^2)/n_L^2\le\Lambda_L^{(1)}`$ rendent les deux
décroissances équivalentes après élimination des petites racines. La cible
asymptotique du script est la moyenne annealed
$`\mathbb E[\Lambda_L^{(1)}]\to0`$. L'objet reste une enveloppe de Jensen à
hiérarchie commune, plus forte que la cible à deux dendrogrammes.

### 1.1 Protocole

Commande exacte :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/critical_cut_collective_gibbs_diagnostic.py \
  --side 4 \
  --repetitions 256 \
  --p 0.81 \
  --maximum-block-count 16 \
  --seed 20260726
```

Les 256 environnements passent tous le cutoff. L'énumération porte sur toutes
les orientations collectives des blocs critiques de chaque racine.

### 1.2 Résultat

```text
lambda_max(A_R*) / n        = 0.9507358532 +/- 0.0045624262
collective persistence      = 0.9139725034 +/- 0.0070179255
critical diagonal mass      = 0.7516784668 +/- 0.0134804604
off-diagonal persistence    = 0.1622940366 +/- 0.0102785579
used / skipped environments = 256 / 0
```

Les audits algébriques donnent

```text
maximum trace identity error       = 5.33e-15
maximum spectral persistence error = 1.14e-13
```

À ce volume, $`\Lambda_L^{(1)}`$ n'est pas seulement au-dessus de un seuil
de contraction raisonnable : elle est proche de sa borne maximale. La porte
single-$D$ est donc **défavorable à $`L=4`$**. Ce résultat ne décide pas sa
limite lorsque $`L\to\infty`$ ; la masse diagonale critique est elle-même
énorme sur seize sommets.

## 2. Décomposition signée à deux dendrogrammes

Pour deux Gibbs conditionnels exacts, soit

```math
P_{ij}
=
\pi_{O,D^{(1)}}(\sigma_i\sigma_j)
\pi_{O,D^{(2)}}(\sigma_i\sigma_j).
\qquad\text{(2.1)}
```

Dans l'intersection $`G_{12}^\star`$ des deux plus grandes racines, notons
$`\mathcal C_{12}`$ le raffinement commun des deux partitions critiques.
Le code calcule séparément

```math
\mathcal E_{\mathrm{diag},L}
=
\frac1{n_L^2}
\sum_{C\in\mathcal C_{12}}
\sum_{i,j\in C}
P_{ij}
\qquad\text{(2.2)}
```

et le reste signé

```math
\mathcal E_{\mathrm{off},L}
=
\frac1{n_L^2}
\sum_{\substack{
C,C'\in\mathcal C_{12}\\
C\ne C'
}}
\sum_{\substack{
i\in C\\
j\in C'
}}
P_{ij}.
\qquad\text{(2.3)}
```

L'identité auditée réalisation par réalisation est

```math
\mathcal E_{\star,L}^{(2)}
=
\mathcal E_{\mathrm{diag},L}
+
\mathcal E_{\mathrm{off},L}.
\qquad\text{(2.4)}
```

Aucune valeur absolue n'est appliquée à (2.3).

### 2.1 Protocole

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/double_giant_replicated_gibbs_diagnostic.py \
  --side 4 \
  --p 0.81 \
  --observations 4 \
  --replica-pairs 8 \
  --seed 3801
```

Pour chaque observation, la postérieure sur les $`2^{16}`$ spins est
énumérée exactement. Chacune des 32 paires utilise deux références
postérieures et deux hiérarchies indépendantes, puis énumère exactement les
deux Gibbs conditionnels. Les erreurs standards du résumé prennent les
quatre observations indépendantes comme clusters ; avec seulement quatre
clusters, elles sont descriptives.

### 2.2 Résultat

```text
direct posterior persistence       = 0.6872083363 +/- 0.1122135077
independent-hierarchy estimate     = 0.6969441934 +/- 0.1125641845
paired independent minus direct    = 0.0097358571 +/- 0.0051277226
double-largest-root contribution   = 0.6967000528 +/- 0.1126471143
same-cell signed contribution      = 0.4968941343 +/- 0.1083931965
distinct-cell signed remainder     = 0.1998059185 +/- 0.0116327206
```

L'erreur maximale dans (2.4) vaut $`1.11\times10^{-16}`$. Parmi les 32
tirages, trois restes hors-diagonale sont strictement négatifs, avec un
minimum de $`-0.0110581\ldots`$ ; le maximum vaut
$`0.6188823\ldots`$. Les cancellations signées existent donc bien au niveau
d'une réalisation, mais elles ne dominent pas ce petit échantillon.

La proximité des deux premières lignes est cohérente avec l'identité
d'augmentation de données. Elle n'est pas attendue paire par paire et ne
constitue pas un test asymptotique avec quatre observations.

## 3. Décision scientifique

Ce calcul élimine deux raccourcis.

1. **L'enveloppe single-$D$ ne suffit pas à ce volume.** Une preuve qui
   remplace les deux dendrogrammes par un carré quenched commun perd les
   cancellations recherchées.
2. **Le mot “signé” est opérationnel.** Le terme hors-diagonale peut changer
   de signe ; le remplacer par sa valeur absolue détruit une information
   réelle.

Il montre aussi que la simple existence de signes négatifs ne suffit pas :
à $`L=4`$, la moyenne de (2.3) reste grande et positive.

La prochaine expérience utile n'est pas de multiplier les répliques sur le
même $`L=4`$. Il faut augmenter le volume en remplaçant l'énumération globale
par une élimination exacte ou certifiée le long des ports, puis mesurer :

```math
L
\longmapsto
\mathbb E[
\mathcal E_{\mathrm{off},L}
]
\qquad\text{et}\qquad
L
\longmapsto
\mathbb P(
\mathcal E_{\mathrm{off},L}<0
).
\qquad\text{(3.1)}
```

Ces statistiques doivent être clusterisées par observation et accompagnées
de la masse diagonale géométrique. Une décroissance numérique robuste ne
deviendrait une preuve qu'après construction de la Palm jointe de deux
corridors et d'une enveloppe non linéaire.
