# Détection de communautés sur graphes signés

Ce dépôt réunit les présentations de Louis Hauseux, le chapitre de thèse
associé et un dossier de recherche sur les dynamiques de clusters.

## Recherche

Le [dossier research/](research/README.md) contient quatre
notes courtes : [hiérarchie et lois exactes](research/01_HIERARCHIE.md),
[audit mathématique](research/02_AUDIT.md) et
[coupe critique pour la weak recovery](research/03_RECOVERY.md), puis
[audit de la coupe dans la géante](research/04_GEANTE_CRITIQUE.md).

La construction retrouve exactement Glauber à la coupe zéro et
Swendsen–Wang à la coupe un. Une amélioration des bornes de recovery
par cette dynamique reste à démontrer. Les résultats antérieurs sont
référencés dans le dossier.

## Présentations

| date | support | sources | PDF |
|---|---|---|---|
| 15 juin 2026 | Séminaire MathNet | [`beamer-presentation/`](beamer-presentation/) | [`main.pdf`](beamer-presentation/main.pdf) |
| 25 juin 2026 | Séminaire NEO | [`beamer-presentation-neo/`](beamer-presentation-neo/) | [`PresentationNIM_2026-06-25.pdf`](beamer-presentation-neo/PresentationNIM_2026-06-25_LouisHauseux_ABayesianFrameworkForCommunityDetectionOnSignedGraphs.pdf) |
| 16 juillet 2026 | Réunion de recherche | [`beamer-presentation-reunion-2026-07-16/`](beamer-presentation-reunion-2026-07-16/) | [`Presentation_2026-07-16.pdf`](beamer-presentation-reunion-2026-07-16/Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.pdf) |

Chaque présentation est un instantané autonome : son thème, ses images et sa
bibliographie restent dans son propre dossier. Pour la compiler :

```bash
cd beamer-presentation
make
```

## Vérifier le dossier

```bash
python3 .agents/check_math.py
python3 .agents/check_markdown_links.py
python3 research/check_hierarchy.py
python3 research/check_giant_cut.py
```

Les scripts utilisent uniquement la bibliothèque standard Python.
