# Diagnostics et no-go

Ce dossier rassemble trois types de documents :

- des **benchmarks exacts locaux**, utiles pour calibrer une conjecture ;
- des **expériences de volume fini**, utiles pour choisir le prochain lemme ;
- des **no-go exacts**, qui ferment un raccourci mais pas nécessairement toute
  la stratégie hiérarchique.

Le [fichier 29](29_AUDIT_FROID_PIVOT_RANGS_REELS.md) est la lecture
prioritaire : il explique pourquoi la criticalisation multiport et l'état
local fidèle ont été abandonnés.

Le dossier [`finite_volume/`](finite_volume/) contient les résultats qui ne
doivent jamais être extrapolés sans une preuve supplémentaire. Son
[fichier 40](finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md) contient
le premier test spectral single-$D$ et la première mesure directe du reste
signé à deux dendrogrammes à $`p=0.81`$. L'[index
raisonné](../INDEX.md#diagnostics-benchmarks-et-no-go) précise le rôle de
chaque note.
