# Recherche

Ce dossier regroupe les notes mathématiques et les calculs reproductibles
associés aux présentations.

## Projet actif

### Swendsen--Wang hiérarchique et weak recovery

Le projet cherche à adapter le couplage Swendsen--Wang du chapitre 11 à une
dynamique hiérarchique construite par horloges exponentielles. L'objectif est
d'obtenir une obstruction de weak recovery plus forte sur le GSBM
triangulaire. Le point $`p=4/5`$ sert de pré-certificat technique ; la
première cible donnant un seuil strictement supérieur à $`0.8`$ est
$`p_0=0.805`$.

- [Vue d'ensemble pédagogique](hierarchical-swendsen-wang/README.md)
- [Programme de recherche prioritaire](hierarchical-swendsen-wang/00_RESEARCH_PROGRAM.md)
- [Feuille de route technique](hierarchical-swendsen-wang/05_PROOF_ROADMAP.md)
- [Feuille de route vers un seuil strictement supérieur à 0,8](hierarchical-swendsen-wang/26_FEUILLE_DE_ROUTE_PSTAR.md)
- [Calculs reproductibles](hierarchical-swendsen-wang/computations/README.md)

La piste centrale étudie une paire lointaine sur son corridor hiérarchique
réel. À squelette, tailles de coupes et état de bord fixés, elle avance les
canaux tardifs jusqu'au niveau critique, puis mesure la perte de corrélation
sans supposer que le LCA ponctuel est lui-même critique.
