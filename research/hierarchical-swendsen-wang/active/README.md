# Programme actif

## Parcours de lecture pédagogique

Pour comprendre le problème avant de travailler dessus :

1. [42 — problème central](../foundations/ancestral/42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md) :
   seuil de fusion d'une paire lointaine et chaîne ancestrale des
   $`\Lambda_v`$ — l'énoncé canonique (P1/P2), les conventions de
   conditionnement et le verrou G1–G3 ;
2. [36 §1–§7 — arbre géant et Gibbs critique](36_ARBRE_GEANT_GIBBS_CRITIQUE.md) :
   calibration exacte de $`\beta_c(p)`$, couplage par rangs, forme exacte
   du facteur ancestral $`\Lambda_{u,p}(y,z)`$ sur le quotient critique ;
3. [08](../foundations/ancestral/08_ANCESTRAL_LAMBDA_CHAIN.md) /
   [10](../foundations/ancestral/10_ANCESTRAL_LAMBDA_ESTIMATION.md) :
   les quatre $`\Lambda_v^{ab}`$ et leur estimation conditionnelle ;
4. [30](30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) +
   [33](33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) : la dissipation
   $`L^2`$ — ce que chaque intégration d'orientations détruit (un seul
   programme en deux notes) ;
5. [35 §6 — transport entropique](35_DISTANCE_ENTROPIE_ERGODICITE.md) :
   le budget d'entropie relative vers la loi inclinée par l'énergie ;
6. [SBM/07 — la dynamique retrouve-t-elle les seuils du SBM ?](../../SBM/07_SEUILS_PAR_LA_DYNAMIQUE.md) :
   la porte de calibration à franchir avant tout retour au GSBM.

## Ordre de priorité de travail

1. [38 — double géante et Gibbs exact répliqué](38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) :
   cible prioritaire — contrôler le reste signé entre cellules de la
   double géante ;
2. [41 — désintégration Palm du reste signé](41_DESINTEGRATION_PALM_RESTE_SIGNE.md) :
   réduction exacte au carré cross-block moyenné en dendrogramme ;
   verrous TRI1-o/TRI2 (embedding des ports, contraction) ;
3. [37 — pilote SBM](37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) et
   [39 — port global du SBM fini](39_PORT_GLOBAL_SBM_RECOVERY.md) :
   calibration broadcast, no-gos, port global exact ;
4. [36](36_ARBRE_GEANT_GIBBS_CRITIQUE.md) : identités exactes de
   référence à un dendrogramme (le test spectral single-$D$ a déjà été
   exécuté — enveloppe $`\approx0.95`$ à $`L=4,5,6`$, saturée par
   construction à ces volumes, cf. la note) ;
5. [35](35_DISTANCE_ENTROPIE_ERGODICITE.md), [30](30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md),
   [33](33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) : outillage
   analytique subordonné à la cible 38/41.

## Notes et diagnostics associés

| note | modules de calcul |
|---|---|
| 37, 38 | [`sbm_broadcast_density_evolution.py`](../computations/sbm_broadcast_density_evolution.py), [`sbm_critical_cut_replica_diagnostic.py`](../computations/sbm_critical_cut_replica_diagnostic.py) |
| 39 | [`sbm_global_port_convolution.py`](../computations/sbm_global_port_convolution.py), [`sbm_recovery_regimes_diagnostic.py`](../computations/sbm_recovery_regimes_diagnostic.py) |
| 38, 41 | [`double_giant_replicated_gibbs_diagnostic.py`](../computations/double_giant_replicated_gibbs_diagnostic.py), [`critical_cut_collective_gibbs_diagnostic.py`](../computations/critical_cut_collective_gibbs_diagnostic.py), [`critical_cut_quotient_elimination.py`](../computations/critical_cut_quotient_elimination.py) |
| 36 | [`giant_component_quotient_diagnostic.py`](../computations/giant_component_quotient_diagnostic.py) |
| 38 (almost/exact triangulaire) | [`triangular_recovery_regimes_diagnostic.py`](../computations/triangular_recovery_regimes_diagnostic.py) |
| 30, 33 | [`two_step_projective_l2_cell.py`](../computations/two_step_projective_l2_cell.py), [`two_step_l2_population_diagnostic.py`](../computations/two_step_l2_population_diagnostic.py), [`nested_projection_l2_diagnostic.py`](../computations/nested_projection_l2_diagnostic.py) |
| SBM/07 | [`sbm_glauber_stability_benchmark.py`](../computations/sbm_glauber_stability_benchmark.py) |

> [!WARNING]
> « Actif » signifie que la piste est jugée plausible et falsifiable. Cela ne
> signifie ni que ses lemmes sont prouvés, ni qu'elle améliore déjà la borne
> $`0.809439`$.

Le [statut scientifique](../CURRENT_STATUS.md) fixe l'ordre de travail et
prévaut sur toute ancienne feuille de route.
