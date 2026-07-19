# Calculs reproductibles

Ce dossier contient les contre-audits finis des énoncés mathématiques. Une
sortie numérique est toujours étiquetée comme diagnostic tant qu'elle n'est
pas accompagnée d'une preuve ou d'un certificat d'intervalles.

Le contexte, l'ordre de travail et les lemmes servis par ces calculs sont
résumés dans le [programme prioritaire](../00_RESEARCH_PROGRAM.md).

## Validation complète

Depuis la racine du dépôt :

```bash
python3 .agents/check_math.py
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_*.py' -v
python3 -m compileall -q \
  research/hierarchical-swendsen-wang/computations
```

Les scripts n'ont pas de dépendance scientifique externe.

## Voie active

| module | fonction |
|---|---|
| `ancestral_lambda_chain.py` | quatre taux ancestraux et message exact sur un squelette fini |
| `ancestral_lambda_estimation.py` | moments pondérés et certificat de queue des ancêtres |
| `critical_component_boundary.py` | marques de frontière, canal instantané, charge géométrique, taux Palm et critères quatre états |
| `hierarchical_flip_probabilities.py` | probabilités racine, feuille, nœud interne et transfert tordu |
| `joint_hierarchical_sweep.py` | sweep exact top-down/bottom-up sur petits tores |
| `favorable_time_comparison.py` | anti-alignement, Blackwell à taille fixe et incomparabilité cross-size certifiée à $`p=4/5`$ |
| `pair_favorability_diagnostic.py` | comparaison pondérée critique/tardive par classes de paires |
| `collapsed_corridor_transfer.py` | transfert collapsed exact pour un corridor et un prior corrélé |
| `cactus_collapsed_certificate.py` | canal cactus exact, LCA seul contre corridor complet et certificat $`p=0.8`$ |
| `lca_palm_corridor_diagnostic.py` | benchmark snapshot critique, corridor final réel et criticalisation à squelette fixé |
| `triangular_band_collapsed_certificate.py` | premier secteur répliqué E1+ sur une cellule triangulaire neutre à quatre ports |
| `twisted_feynman_kac_composition.py` | composition finie du secteur tordu par un déficit de Feynman--Kac |

Chaque module actif possède un fichier `test_*.py` associé.

Le module `critical_component_boundary.py` contient aussi le contre-audit
des bilans résiduels et les nouveaux calculs conditionnés par une coupe :
moments du vote instantané, charge de Chernoff, fiabilité $`L^2`$ et taux de
fusion $`m u_ps_p(\beta)`$. La fiabilité est recalculée indépendamment à
partir des deux expériences binomiales symétriques dans les tests.

## Calculs auxiliaires conservés

| module | rôle de contre-audit |
|---|---|
| `critical_band_thresholds.py` | constantes triangulaires et inversion des horloges |
| `critical_merger_oracle.py` | canal local critique sans message ancestral |
| `critical_pair_path_geometry.py` | hiérarchie de Kruskal et échantillonnage Palm fini |
| `path_decorrelation_threshold.py` | oracle PATH-FAC et seuils conditionnels |
| `triangle_block_sdpi.py` | canal d'un triangle isolé |
| `nishimori_hierarchical_entropy.py` | identité entropique de face |

Ces modules restent testés, mais ne déterminent plus l'ordre du programme de
recherche.

## Corridor collapsed à $`p=0.8`$

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/collapsed_corridor_transfer.py
```

Le script :

1. énumère
   $`\mathbb E[\mathbb E(F(X)\mid K_1,\ldots,K_h)^2]`$ ;
2. compare les niveaux critiques et tardifs sous un prior uniforme ;
3. répète le calcul avec un prior de chaîne d'Ising corrélé ;
4. affiche la contraction exacte de $`N`$ blocs neutres $`m=2`$.

La sortie de référence est :

```text
uniform: critical=0.232015050844 late=0.047131567858 gap=0.184883482986
correlated: critical=0.426226710965 late=0.221677424071 gap=0.204549286894
neutral m=2 blocks= 5 bound=0.160505443478
neutral m=2 blocks=10 bound=0.025761997386
neutral m=2 blocks=20 bound=0.000663680509319
neutral m=2 blocks=40 bound=4.4047181845e-07
```

Ces nombres valident l'énumération sur un corridor fixé ; ils ne représentent
pas la loi du tore triangulaire.

## Diagnostic LCA-Palm du corridor réel à $`p=0.805`$

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/lca_palm_corridor_diagnostic.py \
  --side 12 --repetitions 50 --p 0.805 \
  --distance-fraction 0.25 --maximum-bucket-size 8 \
  --maximum-charge 1.0 --seed 20260719
```

Le module sépare deux expériences qui ne sont pas interchangeables.

1. Dans le benchmark snapshot à $`q_c`$, les coupes candidates sont
   pondérées par l'intensité pré-saut $`mN_\rho`$.
2. Dans l'arbre final réalisé jusqu'à $`q_1=2p-1`$, chaque nœud est pondéré
   seulement par $`N_\rho`$ ; la course de Kruskal a déjà introduit $`m`$.

Le contre-audit détecte explicitement le faux poids $`m^2N_\rho`$ et vérifie
que la somme des $`N_\rho`$ sur les LCA réalisés est exactement le nombre de
paires ordonnées lointaines connectées. Le benchmark snapshot change le
squelette et ne constitue pas une domination de Blackwell. Seul le corridor
final applique la criticalisation favorable
$`q_v\mapsto\min(q_v,q_c)`$ sur un squelette inchangé.

À $`L=12`$, avec les paramètres ci-dessus, le corridor final contient en
moyenne jackknife

```text
all corridor cuts:       19.002 +/- 0.328
bucket size exactly 2:    2.929 +/- 0.143
favourable proxy G_8,1:   8.687 +/- 0.281
```

Le proxy impose $`2\le m\le8`$ et
$`m h_p(q_v^{\mathrm{fav}})^2\le1`$. Il ne calcule ni le screening, ni les
ports latéraux, ni le potentiel extérieur. Les erreurs sont des jackknives
par environnement de rang, jamais des erreurs i.i.d. par nœud. Le tableau
d'échelle complet et ses limites sont dans le
[fichier 28](../28_FIRST_CORRIDOR_P0805_RESULTS.md).

## Cellule triangulaire répliquée E1+

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/triangular_band_collapsed_certificate.py
```

La cellule possède deux ports gauches, deux ports droits et un triangle. Ses
quatre arêtes sont conditionnées fermées à $`q_c`$. Le programme construit
$`\mathbb E_Z[K_Z\otimes K_Z]`$ avec le même environnement résiduel dans les
deux répliques, puis sépare le bloc de masse et le secteur
$`\chi\otimes\chi`$.

La sortie de référence à $`p=0.805`$ commence par :

```text
scope=E1+ neutral all-closed cell sector test
shared chi-x-chi uniform coefficient=0.293993788340
rational upper bound <0.293993788341 strict=True
independent-environment counterfactual coefficient=0.086432347583
```

Le majorant strictement inférieur à $`0.3`$ est certifié par intervalles
rationnels. Le contre-factuel à environnements indépendants montre pourquoi
les deux répliques doivent partager $`Z`$.

Un champ extérieur non borné donne cependant un no-go exact : le second
moment brut passe de $`0.293993788340`$ à $`B=0`$ à
$`0.998663483928`$ à $`B=8`$, puis tend vers un. Le certificat E1+ ne vaut
donc pas uniformément sur tous les potentiels extérieurs. Ce résultat
n'exclut pas une norme centrée ou annealed avec la polarisation dans l'état.

Cette cellule n'est pas encore E2/T2 : l'arête gagnante, la partition
ouverte, les $`\Lambda`$ ancestraux, les attaches en peigne et la loi Palm
sont absents.

## Composition tordue de Feynman--Kac

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/twisted_feynman_kac_composition.py
```

Pour chaque transfert positif levé déjà placé dans une normalisation
stochastique commune, le module construit

```math
K=\sum_\epsilon T_\epsilon,
\qquad
U=\sum_\epsilon\epsilon T_\epsilon,
\qquad
r=\frac{|U|}{K}\in[0,1].
```

Il certifie en dimension finie l'enveloppe du produit tordu par l'espérance
de $`\prod r`$ sous la chaîne de masse $`K`$. Deux exemples en `Fraction`
comparent exactement la récurrence dynamique à une énumération indépendante
de tous les chemins.

Pour la cellule E1+ à $`p=0.805`$, la sortie contient :

```text
depth= 2 signed=0.0649753038062 FK=0.0738919329503 uniform=0.0864323475826
depth=10 signed=1.89285427006e-07 FK=1.08695758136e-06 uniform=4.82371394009e-06
```

La normalisation de Doob commune à tous les blocs du corridor et
l'identification au transfert LCA-Palm réel restent ouvertes. Le module
certifie le lemme de composition fini, pas ces deux étapes.

## Certificat cactus collapsed

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/cactus_collapsed_certificate.py
```

Le module sépare exactement les fusions directe-première et
latérale-première d'un triangle, puis calcule les coefficients connecté et
pivotal. À $`p=0.8`$ et au rang critique, le début de la sortie est :

```text
p=0.8 q_critical=0.347296355334 beta_critical=0.410716539196
one block: connection=0.426022047760 direct-first|connected=0.564864236889 connected_reliability=0.886752566857 flux_reliability=0.791530736866
blocks=20 lca_only=0.791530736866 full_over_lca=0.101917000003 second_moment=0.0903751613589 conformity=0.545187580679 lca_second=0.0806704381115 lca_conformity=0.540335219056
blocks=40 lca_only=0.791530736866 full_over_lca=0.00921076532048 second_moment=0.00816766979065 conformity=0.504083834895 lca_second=0.00729060386122 lca_conformity=0.503645301931
three path-first blocks: direct=0.334328185717 transfer=0.334328185717 gap=0
```

La preuve fermée, la distinction entre connexion cumulative et densité LCA,
ainsi que les limites du transfert à la grille sont dans le fichier 21. La
comparaison LCA seul/corridor complet est dans le fichier 22. Les tests
comparent la formule à une quadrature, une énumération globale des spins et
marques, un produit local indépendant et des intervalles rationnels.

## Certificat Blackwell lorsque la taille change

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/favorable_time_comparison.py
```

La fin de la sortie donne :

```text
critical m=4 vs late m=2 call gap in [-0.00718430527188, -0.00718430527187]
late m=2 vs critical m=4 call gap in [-0.0445551245997, -0.0445551245997]
```

Les deux bornes sont obtenues avec des `Fraction`, à partir d'encadrements
rationnels de $`q_\triangle`$ et de $`4^{-1/5}`$. Elles prouvent que le
bucket critique de taille quatre et le bucket tardif de taille deux au niveau
$`t=4/5`$ sont incomparables. Les fonctions génériques de comparaison
cross-size utilisent des flottants et restent des diagnostics ; ce certificat
particulier, lui, est une preuve par intervalles exacts.

## Diagnostic HF-S2 sur petits tores

Les trois lignes du fichier 19 se reproduisent par :

```bash
python3 research/hierarchical-swendsen-wang/computations/pair_favorability_diagnostic.py \
  --side 4 --repetitions 200 --sweeps 200 --p-values 0.8 \
  --distance-fraction 0.25 --critical-window 0.05 --seed 20260719
python3 research/hierarchical-swendsen-wang/computations/pair_favorability_diagnostic.py \
  --side 6 --repetitions 120 --sweeps 160 --p-values 0.8 \
  --distance-fraction 0.25 --critical-window 0.05 --seed 20260719
python3 research/hierarchical-swendsen-wang/computations/pair_favorability_diagnostic.py \
  --side 8 --repetitions 60 --sweeps 120 --p-values 0.8 \
  --distance-fraction 0.25 --critical-window 0.05 --seed 20260719
```

La sortie JSON contient les masses de classes, le nombre d'environnements
contributeurs, les deux ordres de sweep, les seconds moments et le contraste
jackknife apparié. Le contraste à $`L=8`$ est compatible avec zéro à environ
une erreur standard : il ne faut pas surinterpréter les six signes positifs.

## Conventions de développement

- Toute probabilité nouvelle doit être calculée de deux façons indépendantes
  lorsque la taille d'état le permet.
- Les tests utilisent des exemples déterministes ou des graines explicites.
- Un estimateur de carré de moyenne doit enlever les termes diagonaux.
- Les deux répliques partagent le même environnement ; seuls leurs aléas de
  heat bath sont indépendants.
- Le bucket d'une fusion contient toutes les arêtes physiques de la coupe.
- L'identité de l'arête gagnante de Kruskal est oubliée dans le dendrogramme
  non marqué.
- Les fichiers de résultats bruts ne sont ajoutés que s'ils sont nécessaires
  à une figure ou à un certificat non reproductible rapidement.

## Prochaine étape

La cellule E1+ a rempli son rôle de test de secteur et a révélé le verrou de
polarisation. Le prochain module doit être une cellule **T2-Kruskal**, et non
une simple bande neutre plus large. Il devra :

1. inclure une fusion réelle et marginaliser correctement l'arête gagnante ;
2. encoder la partition ouverte, au moins trois ports et une attache en
   peigne ;
3. conserver le potentiel extérieur et les quatre
   $`\Lambda_v^{ab}`$ d'un ancêtre ;
4. construire les deux répliques dans le même environnement ;
5. produire un déficit tordu dépendant de l'état, composable par une formule
   de Feynman--Kac ;
6. fournir deux implémentations concordantes avant tout certificat
   d'intervalles.

Le diagnostic Palm doit ensuite enregistrer exactement les ports et états
requis par cette cellule. Il est inutile de lancer une preuve multiscale
d'abondance avant d'avoir identifié ce motif fini. L'ordre détaillé est dans
la [sous-feuille de route](../27_SUBROADMAP_CORRIDOR_P0805.md) et les
[premiers résultats](../28_FIRST_CORRIDOR_P0805_RESULTS.md).
