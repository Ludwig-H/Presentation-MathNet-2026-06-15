# Index raisonné des notes

Cet index classe les documents par **statut scientifique actuel**. Les numéros
`00` à `38` indiquent seulement leur ordre de création.

> [!TIP]
> Pour une première lecture, ne parcourez pas les fichiers dans l'ordre
> numérique. Suivez l'un des parcours du [README](README.md#7-parcours-de-lecture).

## Pages de pilotage

| document | rôle |
|---|---|
| [README](README.md) | introduction pédagogique et parcours de lecture |
| [CURRENT_STATUS](CURRENT_STATUS.md) | source de vérité, acquis, verrous et priorités |
| [38 — double géante et Gibbs exact répliqué](active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) | cible prioritaire réduite à un reste signé |
| [37 — pilote SBM](active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) | calibration broadcast, contre-test de coupe partagée et verrou fini |
| [36 — arbre géant à $D$ fixé](active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md) | diagnostic oracle à un dendrogramme |
| [35 — distance, entropie et ergodicité](active/35_DISTANCE_ENTROPIE_ERGODICITE.md) | moteur conditionnel après le test spectral |
| [calculs reproductibles](computations/README.md) | correspondance entre scripts, preuves et diagnostics |
| [littérature](references/LITERATURE.md) | sources primaires et limites de transfert |

## Résultats établis

### Résultat quantitatif non hiérarchique

| note | contenu | portée exacte |
|---|---|---|
| [11 — canal de triangle](results/non_hierarchical/11_TRIANGLE_BLOCK_SDPI.md) | dérivation du canal multi-état et secteurs less-noisy | outil local ; le théorème canonique est 34 |
| [34 — certificat P809439](results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) | certificats rationnels, Sturm et fermeture globale | $`p_{\mathrm{WR}}\ge0.809439`$ sur le tore triangulaire |

### Résultats hiérarchiques dans un domaine restreint

| note | contenu | portée exacte |
|---|---|---|
| [21 — certificat cactus](results/hierarchical/21_CACTUS_COLLAPSED_CERTIFICATE.md) | perte exponentielle du corridor | cactus triangulaire d'articulations, pas grille complète |
| [22 — LCA ou hiérarchie complète](results/hierarchical/22_LCA_VS_FULL_HIERARCHY.md) | ordre des projections et comparaison de persistance | volume fini, à dendrogramme fixé |

## Fondations réutilisables

### Modèle, dynamique et weak recovery

| note | notion principale | lire avant |
|---|---|---|
| [01 — cadre mathématique](foundations/01_MATHEMATICAL_FRAMEWORK.md) | loi jointe, dendrogramme non marqué et heat baths exacts | toute note hiérarchique |
| [02 — baseline du chapitre 11](foundations/02_CHAPTER_11_BASELINE.md) | borne de percolation et point de départ | 03, 11, 34 |
| [03 — critère pairwise](foundations/03_HIERARCHICAL_WEAK_RECOVERY.md) | de la corrélation spin–spin à l'absence de weak recovery | programme actif |
| [04 — GSBM triangulaire](foundations/04_TRIANGULAR_GSBM.md) | géométrie, paramètres et constantes | calculs spécifiques à la grille |

### Coupes, projections et transfert

| note | notion principale | limite à garder en tête |
|---|---|---|
| [18 — transfert répliqué sous Palm](foundations/18_CRITICAL_PALM_REPLICATED_TRANSFER.md) | identités répliquées et pondérations | l'ancienne route spectrale est dépassée |
| [19 — projections de heat bath](foundations/19_FAVORABLE_SWEEP_PROJECTIONS.md) | projections orthogonales et Blackwell mono-bit | la comparaison ne vaut pas uniformément en multiport |
| [20 — corridor collapsed](foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md) | projection du corridor et surrogate factorisé | ne pas identifier le surrogate au corridor réel |
| [25 — information des coupes](foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) | charge, frontière conditionnelle et conventions de Palm | référence canonique pour les facteurs $`mN_\rho`$ et $`N_\rho`$ |

### Trilogie ancestrale

| note | notion principale |
|---|---|
| [08 — chaîne exacte des taux](foundations/ancestral/08_ANCESTRAL_LAMBDA_CHAIN.md) | quatre taux à chaque niveau et incidences ancestrales |
| [10 — estimation des taux](foundations/ancestral/10_ANCESTRAL_LAMBDA_ESTIMATION.md) | moments pondérés et queues |
| [14 — frontière critique](foundations/ancestral/14_CRITICAL_COMPONENT_BOUNDARY.md) | marques résiduelles, majorité et Palm de fusion |

## Programme actif

| note | statut précis | usage actuel |
|---|---|---|
| [30 — dissipation quadratique](active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) | identités exactes + diagnostic fini | socle opératoriel ; accumulation générique abandonnée |
| [33 — cellules critiques](active/33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) | sous-programme conditionnel | géométrie et cellule à deux updates |
| [35 — distance–entropie](active/35_DISTANCE_ENTROPIE_ERGODICITE.md) | moteur conditionnel | seulement après une marge spectrale sur la double géante |
| [36 — arbre géant à $D$ fixé](active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md) | réduction finie exacte, mais oracle | diagnostics de rangs, buckets et Gibbs quenched |
| [37 — pilote SBM](active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) | calibration exacte du broadcast, transfert fini ouvert | sépare équilibre, dynamique, coupe partagée et port global |
| [38 — double géante répliquée](active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) | réduction exacte au reste inter-cellules signé, seuil ouvert | programme prioritaire sur deux Gibbs d'arbres entiers |

## Diagnostics, benchmarks et no-go

Ces fichiers sont utiles pour tester une intuition. Ils ne doivent pas être
cités comme preuves d'un seuil en volume infini.

| note | résultat du diagnostic | décision associée |
|---|---|---|
| [07 — bande critique](diagnostics/07_CRITICAL_BAND_CRITERION.md) | benchmark de localisation et transmission | réutilisable pour définir les fenêtres near-critical |
| [09 — oracle de fusion](diagnostics/09_CRITICAL_MERGER_ORACLE.md) | canal local exact d'une coupe | insuffisant sans message ancestral et géométrie Palm |
| [13 — horloges de Nishimori](diagnostics/13_NISHIMORI_HIERARCHICAL_CLOCKS.md) | calibration entropique de face | identité locale, pas obstruction globale |
| [15 — parité au LCA](diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md) | désintégration locale de la parité | oracle sans fermeture ancestrale |
| [24 — bilan résiduel](diagnostics/24_SIMPLE_RESIDUAL_BALANCE_OBSTRUCTION.md) | test d'une fermeture simple | remplacé par le conditionnement géométrique de 25 |
| [28 — premier cycle P0805](diagnostics/finite_volume/28_FIRST_CORRIDOR_P0805_RESULTS.md) | audits Palm et petites cellules | historique de volume fini |
| [29 — audit aux rangs réels](diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md) | contre-exemple multiport et déficit local nul | lecture obligatoire avant toute nouvelle route |

## Archives

### Oracles locaux abandonnés comme stratégie globale

| note | intérêt historique | raison de l'archivage |
|---|---|---|
| [06 — corrélation au LCA](archive/oracles/06_LCA_SPIN_CORRELATION.md) | calculs locaux détaillés | le LCA seul n'exploite pas la profondeur |
| [16 — chemin descendant](archive/oracles/16_FLIP_PROBABILITIES_DESCENDANT_PATH.md) | formules de flips à tous les niveaux | la factorisation PATH-FAC n'est pas la loi réelle |
| [17 — seuil PATH-FAC](archive/oracles/17_PATH_DECORRELATION_THRESHOLD.md) | seuil du modèle factorisé | oracle non transféré à la grille |

### Feuilles de route remplacées

| note | ancien objectif | document qui la remplace |
|---|---|---|
| [00 — programme prioritaire](archive/roadmaps/00_RESEARCH_PROGRAM.md) | corridor aux rangs réels | [statut courant](CURRENT_STATUS.md) |
| [05 — feuille de route des preuves](archive/roadmaps/05_PROOF_ROADMAP.md) | fermeture hiérarchique initiale | [programme 35](active/35_DISTANCE_ENTROPIE_ERGODICITE.md) |
| [12 — réduction favorable](archive/roadmaps/12_FAVORABLE_HIERARCHICAL_REDUCTION.md) | criticalisation du cas favorable | réfutée en multiport par 29 |
| [23 — obstruction dite optimale](archive/roadmaps/23_OPTIMAL_WEAK_RECOVERY_OBSTRUCTION.md) | accumulation annulaire uniforme | concentration rare observée dans 30–33 |
| [26 — route vers un seuil supérieur à 0,8](archive/roadmaps/26_FEUILLE_DE_ROUTE_PSTAR.md) | première amélioration quantitative | objectif dépassé par 34 |
| [27 — corridor à p = 0,805](archive/roadmaps/27_SUBROADMAP_CORRIDOR_P0805.md) | protocole expérimental initial | audit 29 puis programme 35 |

### Certificats quantitatifs subsumés

| note | borne certifiée | référence canonique actuelle |
|---|---:|---|
| [31 — certificat A0](archive/certificates/31_CERTIFICAT_RATIONNEL_A0.md) | $`0.805`$ | [34](results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) |
| [32 — certificat P809](archive/certificates/32_CERTIFICAT_RATIONNEL_P809.md) | $`0.809`$ | [34](results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) |

## Comment citer une affirmation du dépôt

1. Citer un fichier de [`results/`](results/) pour un théorème final.
2. Citer un fichier de [`foundations/`](foundations/) pour une identité ou un
   cadre réutilisable.
3. Qualifier explicitement de diagnostic tout résultat de
   [`diagnostics/`](diagnostics/).
4. Ne jamais citer [`archive/`](archive/) comme état actuel sans mentionner
   le document qui l'a remplacé.
5. Pour la stratégie en cours, citer
   [`CURRENT_STATUS.md`](CURRENT_STATUS.md) et le
   [programme 38](active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md).
