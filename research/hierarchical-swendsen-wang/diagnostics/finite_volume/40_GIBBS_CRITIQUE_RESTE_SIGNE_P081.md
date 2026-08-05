# Diagnostic fini à $`p=0.81`$ : enveloppe spectrale et reste signé

**Statut : diagnostic reproductible de volume fini ; le single-$D$ calcule
exactement les orientations collectives conditionnellement à une
représentation interne all-plus fixée, jusqu'à $`L=6`$, tandis que le
double-$D$ énumère les Gibbs conditionnels complets seulement à $`L=4`$ ;
environnements, observations et hiérarchies sont échantillonnés ; aucune
extrapolation thermodynamique et aucune nouvelle borne de weak recovery.**

Cette note exécute les deux premières portes quantitatives de la
[cible double-géante](../../active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) :

1. tester l'enveloppe plus forte mais plus simple à un dendrogramme ;
2. si elle reste macroscopique, mesurer directement le reste inter-cellules
   **signé** avec deux dendrogrammes indépendants.

Les résultats single-$D$ à $`L=4,5,6`$ et double-$D$ à $`L=4`$ sont nets :

- l'enveloppe single-$D$ reste proche de un aux trois tailles et ne ferme
  pas la preuve ;
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

Commandes exactes, avec la même graine aux trois tailles :

```bash
for side in 4 5 6; do
  python3 \
    research/hierarchical-swendsen-wang/computations/critical_cut_collective_gibbs_diagnostic.py \
    --side "$side" --repetitions 256 --p 0.81 \
    --maximum-block-count 16 --seed 20260726
done
```

Pour chaque taille, les 256 environnements passent tous le cutoff.
L'énumération porte sur toutes les orientations collectives des blocs
critiques de chaque racine, conditionnellement à la forme interne all-plus
fixée.

### 1.2 Résultat

| $`L`$ | $`\lambda_{\max}(A_{R_\star})/n`$ | persistance collective | masse diagonale critique | persistance hors diagonale | utilisés / exclus |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.9507358532 ± 0.0045624262 | 0.9139725034 ± 0.0070179255 | 0.7516784668 ± 0.0134804604 | 0.1622940366 ± 0.0102785579 | 256 / 0 |
| 5 | 0.9522258609 ± 0.0035202553 | 0.9124653803 ± 0.0059886145 | 0.7164750000 ± 0.0129157594 | 0.1959903803 ± 0.0103643093 | 256 / 0 |
| 6 | 0.9478030939 ± 0.0034789971 | 0.9038329299 ± 0.0057577803 | 0.6839976370 ± 0.0122504530 | 0.2198352930 ± 0.0107005609 | 256 / 0 |

À $`L=4`$, les audits algébriques donnent

```text
maximum trace identity error       = 5.33e-15
maximum spectral persistence error = 1.14e-13
```

Aux trois volumes, $`\Lambda_L^{(1)}`$ n'est pas seulement au-dessus d'un
seuil de contraction raisonnable : elle reste proche de sa borne maximale.
La masse diagonale critique baisse de $`0.752`$ à $`0.684`$, mais la
persistance hors diagonale augmente et l'enveloppe spectrale ne montre
aucune décroissance. La porte single-$D$ est donc **défavorable jusqu'à
$`L=6`$**. Trois tailles finies ne décident évidemment pas la limite lorsque
$`L\to\infty`$.

### 1.3 Élimination exacte sur le quotient critique

Le module
[`critical_cut_quotient_elimination.py`](../../computations/critical_cut_quotient_elimination.py)
représente la même loi single-$D$ comme un graphe de facteurs binaires sur
les orientations $`z_a`$ des blocs critiques. Pour un bucket hiérarchique
$`u`$, son facteur non constant est

```math
\psi_u(z)
=
N_u(z)\,
\exp\!\left((1-\beta_u)J N_u(z)\right),
```

où $`N_u(z)`$ est le nombre d'arêtes satisfaites du bucket. L'élimination
min-fill calcule chaque corrélation par quatre fonctions de partition
contraintes. Elle est exacte et son coût est exponentiel dans la largeur
induite observée, plutôt que directement dans le nombre total de blocs.

La campagne suivante utilise 64 environnements pour chacune des trois
graines $`41`$, $`3801`$ et $`20260726`$, sans cutoff :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/critical_cut_quotient_elimination.py \
  --side 6 --repetitions 64 --p 0.81 --seed 20260726
```

| $`L`$ | graine | blocs max. dans une racine | portée initiale max. | largeur min-fill max. | temps mural |
|---:|---:|---:|---:|---:|---:|
| 5 | 41 | 9 | 7 | 6 | 0.282 s |
| 5 | 3801 | 14 | 9 | 8 | 1.098 s |
| 5 | 20260726 | 11 | 8 | 7 | 0.561 s |
| 6 | 41 | 10 | 7 | 6 | 0.860 s |
| 6 | 3801 | 12 | 8 | 7 | 1.066 s |
| 6 | 20260726 | 11 | 8 | 7 | 1.262 s |

Les temps sont dépendants de la machine. Les 384 environnements sont tous
éliminés exactement : il n'y a ni exclusion par complexité, ni biais de
sélection. Ces largeurs sont seulement des observations de volume fini, pas
une borne uniforme en $`L`$.

Cette exactitude a une frontière essentielle. La forme interne $`y`$ de
chaque bloc est fixée à la représentation all-plus ; seules ses orientations
collectives sont intégrées. Cette loi de quotient n'est donc **pas** la
Gibbs conditionnelle complète $`\pi_{O,D}`$ sur les $`2^{L^2}`$ spins.
De plus, le champ JSON
`scaled_unconstrained_log_partition` omet les constantes des facteurs
invariants d'orientation et les maxima retirés lors de la normalisation des
tables. Les rapports de fonctions de partition et les corrélations sont
exacts ; cette valeur ne doit pas être lue comme un log-partition absolu.

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
même $`L=4`$. L'énumération directe demanderait $`2^{25}=33\,554\,432`$
états à $`L=5`$, puis
$`2^{36}=68\,719\,476\,736`$ à $`L=6`$, pour la postérieure et pour chaque
Gibbs conditionnelle.

Il faut donc calibrer un **junction tree sur les spins physiques** et y
éliminer exactement les facteurs de la postérieure ainsi que tous les
facteurs hiérarchiques de $`\pi_{O,D}`$. Le moteur de quotient de la
section 1.3 fournit l'algèbre d'élimination et les audits de largeur, mais
pas encore cette marginalisation des formes internes. Toute limite de
largeur ou de portée devra exposer ses exclusions et le biais de sélection
associé. Une MCMC ne remplacerait ce calcul qu'avec diagnostics de mélange
multi-chaînes explicites.

Une fois cette extension réalisée, il faudra mesurer :

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
