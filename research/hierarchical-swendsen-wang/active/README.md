# Programme actif

La feuille de route canonique est :

1. [38 — double géante et Gibbs exact répliqué](38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md)
   est la cible prioritaire : après deux dendrogrammes entiers et deux coupes
   exactes, contrôler le reste signé entre cellules de la double géante ;
2. [37 — pilote SBM à la coupe critique](37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md)
   calibre le Jacobien broadcast $`d\theta^2`$, contre-audite la coupe
   partagée et isole le port global du SBM fini ;
3. [39 — port global et régimes du SBM fini](39_PORT_GLOBAL_SBM_RECOVERY.md)
   écrit ce port exactement, l'élimine par convolution et sépare weak,
   almost exact et exact recovery ;
4. [36 — arbre géant et Gibbs critique](36_ARBRE_GEANT_GIBBS_CRITIQUE.md)
   fournit d'abord le test spectral plus fort à un dendrogramme fixé ;
5. [35 — distance, entropie et ergodicité](35_DISTANCE_ENTROPIE_ERGODICITE.md)
   fournit un moteur analytique après le test spectral ;
6. [30 — dissipation quadratique](30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md)
   fournit le socle opératoriel exact ;
7. [33 — cellules critiques](33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md)
   développe le bloc local à deux projections.

Les diagnostics associés suivent le même ordre :

1. [`sbm_broadcast_density_evolution.py`](../computations/sbm_broadcast_density_evolution.py)
   calibre exactement le secteur overlap sur le broadcast
   $`\mathrm{PGW}(d)`$ ;
2. [`sbm_critical_cut_replica_diagnostic.py`](../computations/sbm_critical_cut_replica_diagnostic.py)
   calcule exactement l'inflation créée par une coupe partagée ;
3. [`sbm_global_port_convolution.py`](../computations/sbm_global_port_convolution.py)
   contracte exactement le port full-$D$ conditionnellement aux tailles de
   racines fournies ;
4. [`sbm_recovery_regimes_diagnostic.py`](../computations/sbm_recovery_regimes_diagnostic.py)
   calcule l'affinité oracle finie et les trois benchmarks de recovery ;
5. [`triangular_recovery_regimes_diagnostic.py`](../computations/triangular_recovery_regimes_diagnostic.py)
   certifie les obstructions oracle almost/exact à degré six fixé ;
6. [`giant_component_quotient_diagnostic.py`](../computations/giant_component_quotient_diagnostic.py)
   mesure conditionnellement la géométrie postcritique de l'arbre de la
   géante finale et compte les environnements sans paire admissible ;
7. [`critical_cut_collective_gibbs_diagnostic.py`](../computations/critical_cut_collective_gibbs_diagnostic.py)
   énumère le Gibbs collectif et son enveloppe spectrale sur petits tores ;
8. [`double_giant_replicated_gibbs_diagnostic.py`](../computations/double_giant_replicated_gibbs_diagnostic.py)
   audite exactement à $`L=4`$ la décomposition par intersections et le
   reste hors-diagonale signé, sans prétendre fournir une tendance
   asymptotique.

> [!WARNING]
> « Actif » signifie que la piste est jugée plausible et falsifiable. Cela ne
> signifie ni que ses lemmes sont prouvés, ni qu'elle améliore déjà la borne
> $`0.809439`$.

Le [statut scientifique](../CURRENT_STATUS.md) fixe l'ordre de travail et
prévaut sur toute ancienne feuille de route.
