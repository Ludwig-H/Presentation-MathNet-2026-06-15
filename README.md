# Community detection on signed graphs

[![Research checks](https://github.com/Ludwig-H/Presentation-MathNet-2026-06-15/actions/workflows/research-checks.yml/badge.svg)](https://github.com/Ludwig-H/Presentation-MathNet-2026-06-15/actions/workflows/research-checks.yml)

Présentations, notes de recherche et calculs reproductibles autour de la
détection de communautés dans les graphes signés, des couplages MCMC et de la
weak recovery.

> [!IMPORTANT]
> La voie de recherche active est la dynamique de Swendsen--Wang
> hiérarchique par horloges exponentielles. Commencer par le
> [programme prioritaire](research/hierarchical-swendsen-wang/00_RESEARCH_PROGRAM.md)
> ou par sa
> [présentation pédagogique](research/hierarchical-swendsen-wang/README.md).

## Axe de recherche prioritaire

L'objectif intermédiaire est de montrer l'impossibilité de weak recovery à

```math
p=\frac45
```

pour le GSBM binaire homogène sur la grille triangulaire. L'expérience
favorable centrale suit deux sommets lointains dont le LCA hiérarchique se
situe juste au seuil de percolation. La preuve cherche ensuite à exploiter le
corridor complet entre les feuilles et ce LCA, plutôt que le seul nœud de
fusion.

Le statut mathématique est explicite dans les notes : la criticalisation est
prouvée à géométrie fixée, le cactus triangulaire fournit un premier
certificat exact, et la domination de la géométrie Palm critique sur la
grille reste ouverte.

## Navigation

| dossier ou fichier | contenu |
|---|---|
| [`research/`](research/) | index des projets de recherche |
| [`research/hierarchical-swendsen-wang/`](research/hierarchical-swendsen-wang/) | programme actif, preuves, audits et calculs |
| [`beamer-presentation/`](beamer-presentation/) | séminaire MathNet, 15 juin 2026 |
| [`beamer-presentation-neo/`](beamer-presentation-neo/) | séminaire NEO, 25 juin 2026 |
| [`beamer-presentation-reunion-2026-07-16/`](beamer-presentation-reunion-2026-07-16/) | réunion du 16 juillet 2026 et slides de dynamique hiérarchique |
| [`ChapII.tex`](ChapII.tex) | source du chapitre de référence |
| [`Manuscrit_de_these.pdf`](Manuscrit_de_thèse.pdf) | manuscrit de thèse |

## Présentations

| date | support | sources | PDF |
|---|---|---|---|
| 15 juin 2026 | Séminaire MathNet | [`beamer-presentation/`](beamer-presentation/) | [`main.pdf`](beamer-presentation/main.pdf) |
| 25 juin 2026 | Séminaire NEO | [`beamer-presentation-neo/`](beamer-presentation-neo/) | [`PresentationNIM_2026-06-25.pdf`](beamer-presentation-neo/PresentationNIM_2026-06-25_LouisHauseux_ABayesianFrameworkForCommunityDetectionOnSignedGraphs.pdf) |
| 16 juillet 2026 | Réunion de recherche | [`beamer-presentation-reunion-2026-07-16/`](beamer-presentation-reunion-2026-07-16/) | [`Presentation_2026-07-16.pdf`](beamer-presentation-reunion-2026-07-16/Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.pdf) |

Chaque présentation se compile depuis son dossier avec :

```bash
make
```

## Reproductibilité de la recherche

Depuis la racine du dépôt :

```bash
python3 .agents/check_math.py
python3 .agents/check_markdown_links.py
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_*.py'
```

Les scripts de recherche n'ont pas de dépendance scientifique externe. Les
résultats numériques ne sont utilisés comme preuves que lorsqu'un certificat
exact, rationnel ou par intervalles est explicitement fourni.

## Convention de statut

Les notes distinguent toujours :

- **établi** : preuve complète dans les hypothèses annoncées ;
- **conditionnel** : implication prouvée sous un lemme nommé ;
- **diagnostic** : calcul fini ou simulation ;
- **conjecture** : cible non utilisée en aval comme un fait.

> [!NOTE]
> Aucune nouvelle borne d'impossibilité à $`p=0.8`$ n'est encore revendiquée.
> Le dépôt documente un programme de preuve audité, ses certificats exacts et
> les verrous restant à fermer.
