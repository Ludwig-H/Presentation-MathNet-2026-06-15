# Calculs reproductibles

Ce dossier contient les contre-audits finis des énoncés mathématiques. Une
sortie numérique est toujours étiquetée comme diagnostic tant qu'elle n'est
pas accompagnée d'une preuve ou d'un certificat d'intervalles.

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
| `critical_component_boundary.py` | marques de frontière, biais Palm et critères quatre états |
| `hierarchical_flip_probabilities.py` | probabilités racine, feuille, nœud interne et transfert tordu |
| `joint_hierarchical_sweep.py` | sweep exact top-down/bottom-up sur petits tores |
| `favorable_time_comparison.py` | anti-alignement, Blackwell à taille fixe et incomparabilité cross-size certifiée à $`p=4/5`$ |
| `pair_favorability_diagnostic.py` | comparaison pondérée critique/tardive par classes de paires |
| `collapsed_corridor_transfer.py` | transfert collapsed exact pour un corridor et un prior corrélé |

Chaque module actif possède un fichier `test_*.py` associé.

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

## Prochain module

Le prochain ajout doit être `cactus_collapsed_certificate.py`. Il devra :

1. énumérer un cactus de deux puis trois triangles ;
2. construire l'état de bord minimal ;
3. calculer le noyau répliqué collapsed ;
4. comparer critique et tardif par un programme linéaire fini ;
5. fournir une seconde implémentation indépendante ;
6. préparer une certification par arithmétique d'intervalles à $`p=4/5`$.

Aucune nouvelle simulation de grand tore n'est prioritaire avant ce module.
