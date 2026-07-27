# Recherche

Ce dossier est la porte d'entrée des notes mathématiques et des calculs
reproductibles associés aux présentations.

## Benchmark spécial SBM

Le dossier [SBM](SBM/) donne une présentation autonome et pédagogique du
SBM symétrique binaire : deux répliques, deux dendrogrammes indépendants
coupés au même $\beta_c$, calibration exacte
$\beta_\chi=\beta_c\Longleftrightarrow d\theta^2=1$, puis extensions
almost exact / exact recovery et limite $\beta\downarrow0$ vers Glauber.
Il suit pas à pas l'architecture bayésienne du chapitre 11 et sépare les
identités prouvées des étapes dynamiques encore ouvertes.

## Le résultat à retenir

Pour le GSBM binaire homogène sur le tore triangulaire, le résultat rigoureux
actuel est

```math
p_{\mathrm{WR}}\ge 0.809439.
```

La preuve complète se trouve dans le
[certificat rationnel P809439](hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).
Elle passe par un canal triangulaire multi-état et non par la dynamique
hiérarchique.

## Le projet actif

Le dossier
[`hierarchical-swendsen-wang/`](hierarchical-swendsen-wang/)
cherche à expliquer et renforcer l'obstruction de weak recovery grâce à une
dynamique qui généralise Glauber aux feuilles et Swendsen--Wang aux racines.

Son statut tient en quatre lignes :

| objet | statut actuel |
|---|---|
| borne $`p_{\mathrm{WR}}\ge0.809439`$ | **établie**, voie non hiérarchique |
| mesure jointe et heat baths hiérarchiques | **établis** en volume fini |
| criticalisation uniforme des fusions multiports | **réfutée** par contre-exemple exact |
| réduction au Gibbs inter-blocs d'un arbre géant fixé | **établie** en volume fini, diagnostic oracle |
| dendrogramme commun aux deux répliques | **réfuté comme cible exacte** par le pilote SBM |
| port global du SBM fini | **écrit et éliminé exactement** par convolution ; comparaison au broadcast ouverte |
| réduction à un reste signé inter-cellules sur la double géante | **établie** géométriquement |
| annulation de ce reste à $`p=0.81`$ | **prioritaire**, sans nouveau seuil prouvé |

## Trois parcours de lecture

### 1. Comprendre rapidement où en est la recherche

1. [Vue d'ensemble pédagogique](hierarchical-swendsen-wang/README.md)
2. [Statut scientifique canonique](hierarchical-swendsen-wang/CURRENT_STATUS.md)
3. [Cible exacte sur la double géante](hierarchical-swendsen-wang/active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md)
4. [Calibration broadcast et verrou du SBM fini](hierarchical-swendsen-wang/active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md)
5. [Port global fini et trois régimes SBM](hierarchical-swendsen-wang/active/39_PORT_GLOBAL_SBM_RECOVERY.md)
6. [Premier test spectral et signé à p = 0,81](hierarchical-swendsen-wang/diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md)
7. [Diagnostic à un dendrogramme géant](hierarchical-swendsen-wang/active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md)

### 2. Vérifier la borne rigoureuse

1. [Baseline du chapitre 11](hierarchical-swendsen-wang/foundations/02_CHAPTER_11_BASELINE.md)
2. [Dérivation du canal triangulaire](hierarchical-swendsen-wang/results/non_hierarchical/11_TRIANGLE_BLOCK_SDPI.md)
3. [Certificat canonique à p = 0,809439](hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md)
4. [Scripts et tests](hierarchical-swendsen-wang/computations/README.md)

### 3. Étudier la dynamique hiérarchique

1. [Mesure jointe et dynamique exacte](hierarchical-swendsen-wang/foundations/01_MATHEMATICAL_FRAMEWORK.md)
2. [Critère pairwise de weak recovery](hierarchical-swendsen-wang/foundations/03_HIERARCHICAL_WEAK_RECOVERY.md)
3. [Information des coupes et biais de Palm](hierarchical-swendsen-wang/foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md)
4. [Projections de heat bath](hierarchical-swendsen-wang/foundations/19_FAVORABLE_SWEEP_PROJECTIONS.md)
5. [Projection collapsed du corridor](hierarchical-swendsen-wang/foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md)
6. [Deux no-go qui imposent le pivot](hierarchical-swendsen-wang/diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md)
7. [Pilote SBM](hierarchical-swendsen-wang/active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md)
8. [Port global et régimes du SBM fini](hierarchical-swendsen-wang/active/39_PORT_GLOBAL_SBM_RECOVERY.md)
9. [Double géante et Gibbs répliqué exact](hierarchical-swendsen-wang/active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md)
10. [Diagnostic à un arbre fixé](hierarchical-swendsen-wang/active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md)
11. [Moteur distance–entropie](hierarchical-swendsen-wang/active/35_DISTANCE_ENTROPIE_ERGODICITE.md)

L'[index exhaustif](hierarchical-swendsen-wang/INDEX.md) classe toutes les
notes, y compris les diagnostics et les anciennes feuilles de route.

## Règle de lecture

Les numéros des fichiers indiquent l'ordre historique de création. Le dossier
qui contient le fichier indique son statut scientifique actuel. En cas de
contradiction entre deux anciennes feuilles de route, le
[statut canonique](hierarchical-swendsen-wang/CURRENT_STATUS.md) prévaut.
