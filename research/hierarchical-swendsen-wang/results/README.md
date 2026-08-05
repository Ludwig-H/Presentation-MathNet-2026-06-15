# Résultats établis

Un fichier de ce dossier contient un théorème prouvé dans les hypothèses qu'il
annonce. Deux portées différentes sont séparées.

## Résultat quantitatif sur la grille triangulaire

Le dossier [`non_hierarchical/`](non_hierarchical/) contient la dérivation du
canal triangulaire et la preuve canonique

```math
p_{\mathrm{WR}}\ge0.809439.
```

Cette preuve n'utilise pas la dynamique hiérarchique.

## Résultats proprement hiérarchiques

Le dossier [`hierarchical/`](hierarchical/) contient un certificat exact sur
un cactus triangulaire et l'ordre de persistance entre le LCA seul et le
corridor complet. Ces théorèmes ne se transfèrent pas automatiquement à la
grille triangulaire complète.

Voir l'[index raisonné](../INDEX.md#résultats-établis) pour leur portée exacte.

Le certificat P809439 se compose d'un **théorème local** (rationnel,
autonome, vérifiable par script) et d'un **théorème global** conditionnel
à trois références citées (Makur–Polyanskiy, Polyanskiy–Wu, Chayes–Lei).
