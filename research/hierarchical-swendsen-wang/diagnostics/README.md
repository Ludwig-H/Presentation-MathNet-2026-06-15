# Diagnostics et no-go

Ce dossier rassemble trois types de documents :

- des **benchmarks exacts locaux**, utiles pour calibrer une conjecture ;
- des **expériences de volume fini**, utiles pour choisir le prochain lemme ;
- des **no-go exacts**, qui ferment un raccourci mais pas nécessairement toute
  la stratégie hiérarchique.

## Parcours pour le problème central

Pour le problème « fusion critique d'une paire lointaine + estimation des
$`\Lambda_v`$ ancestraux » ([note 42](../foundations/ancestral/42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md)),
lire dans cet ordre :

1. [09 — oracle de fusion critique](09_CRITICAL_MERGER_ORACLE.md) :
   calibration exacte du bucket au niveau $`\beta_c`$ et, en §7, le seul
   contrôle certifié des quatre $`\Lambda_v^{ab}`$ ancestraux (sandwich
   des comptes) ;
2. [15 — parité au LCA critique](15_CRITICAL_GIANT_PAIR_FLIP.md) :
   probabilité de retournement $(0,0)/(1,1)$ pour la paire lointaine,
   équivalent aigu en fonction de $p$, hypothèses CUT et ANC ;
3. [29 — audit froid](29_AUDIT_FROID_PIVOT_RANGS_REELS.md) : **lecture
   obligatoire**, les deux no-go exacts (criticalisation multiport
   réfutée, état fidèle $`|U|=K`$) qui délimitent ce qu'aucune preuve ne
   peut faire ;
4. [07 — bande critique](07_CRITICAL_BAND_CRITERION.md) et
   [24 — bilan résiduel](24_SIMPLE_RESIDUAL_BALANCE_OBSTRUCTION.md) en
   support (24 porte un bandeau : sa conclusion « corridor $`m=2`$ » est
   périmée) ;
5. [13 — horloges de Nishimori](13_NISHIMORI_HIERARCHICAL_CLOCKS.md) :
   repère entropique exact, latéral au problème central.

Le dossier [`finite_volume/`](finite_volume/) contient les résultats qui ne
doivent jamais être extrapolés sans une preuve supplémentaire. Son
[fichier 40](finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md) contient
le test spectral single-$D$ à $`L=4,5,6`$ et la première mesure directe du
reste signé à deux dendrogrammes à $`p=0.81`$ ; son
[fichier 28](../archive/roadmaps/28_FIRST_CORRIDOR_P0805_RESULTS.md) est un
journal historique dont le résultat principal est réfuté par 29.
L'[index raisonné](../INDEX.md#diagnostics-benchmarks-et-no-go) précise le
rôle de chaque note.
