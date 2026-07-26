# Calculs reproductibles

Ce dossier contient à la fois des certificats exacts, des briques de calcul
et des diagnostics de volume fini. Leur statut dépend de la sortie produite,
pas du fait qu'ils possèdent un test automatisé.

> [!IMPORTANT]
> Une valeur flottante reste un diagnostic. Elle devient un ingrédient de
> preuve seulement lorsqu'un certificat rationnel, symbolique ou par
> intervalles contrôle toutes les erreurs pertinentes.

Le [statut scientifique](../CURRENT_STATUS.md) fixe l'ordre de travail.
L'[index des notes](../INDEX.md) indique quel énoncé chaque calcul sert.

## Validation complète

Depuis la racine du dépôt :

```bash
python3 .agents/check_math.py
python3 .agents/check_markdown_links.py
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_*.py' -v
python3 -m compileall -q \
  research/hierarchical-swendsen-wang/computations
```

Les scripts n'ont pas de dépendance scientifique externe.

## Choisir le bon calcul

### Certificats utilisés dans des théorèmes

| module | ce qu'il certifie | note canonique |
|---|---|---|
| `rational_a0_less_noisy_certificate.py` | Sturm, dominance diagonale et marge exacte jusqu'à $`p=0.809439`$ | [34](../results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) |
| `cactus_collapsed_certificate.py` | contraction exacte sur le cactus | [21](../results/hierarchical/21_CACTUS_COLLAPSED_CERTIFICATE.md) |
| `ancestral_lambda_chain.py` | quatre taux ancestraux sur un squelette fini | [08](../foundations/ancestral/08_ANCESTRAL_LAMBDA_CHAIN.md) |
| `critical_component_boundary.py` | loi de frontière, charge et taux Palm | [14](../foundations/ancestral/14_CRITICAL_COMPONENT_BOUNDARY.md), [25](../foundations/25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md) |

### Calculs directement utiles au programme actif

| module | question testée | limite actuelle |
|---|---|---|
| `sbm_broadcast_density_evolution.py` | la Gibbs exacte du broadcast possède-t-elle le seuil $`d\theta^2=1`$ ? | le sandwich $`\ell_t\le q_t\le r_t`$ certifie $`\mathrm{PGW}(d)`$ pour toute élimination exacte ; il ne prouve ni un rôle informationnel de $`\beta_c`$, ni un seuil temporel, ni le transfert arbre–graphe |
| `sbm_critical_cut_replica_diagnostic.py` | que coûte le partage de la coupe physique entre deux répliques ? | identité d'arête exacte ; elle distingue le Jacobien marginal $`\theta^2`$ du transfert oracle gonflé |
| `sbm_global_port_convolution.py` | comment contracter exactement balance ou non-arêtes entre les racines d'un full-$D$ SBM ? | convolution exacte conditionnelle aux tailles de racines fournies ; elle ne génère ni $D$, ni sa loi, ni un seuil |
| `sbm_recovery_regimes_diagnostic.py` | quels sont les benchmarks weak, almost exact et exact du SBM classique ? | affinité binomiale finie exacte et constantes de seuil ; oracle seulement, aucune achievability hiérarchique |
| `triangular_recovery_regimes_diagnostic.py` | almost exact ou exact recovery sont-ils possibles à degré six fixé ? | erreur oracle exacte ; nécessités $`p_n\to1`$ et $`n\varepsilon_6(p_n)\to0`$, sans suffisance |
| `double_giant_replicated_gibbs_diagnostic.py` | la décomposition exacte par intersections et cellules ferme-t-elle sur $`L=4`$ ? | le reste hors-diagonale est conservé signé et les Gibbs sont exacts ; observations et hiérarchies restent Monte-Carlo, sans tendance asymptotique |
| `giant_component_quotient_diagnostic.py` | quelle géométrie voit une paire lointaine dans l'arbre de la géante finale, après contraction des blocs critiques ? | diagnostic conditionnel aux environnements ayant une paire admissible, dont le nombre est exposé ; PATH-FAC reste un oracle local factorisé non probant |
| `critical_cut_collective_gibbs_diagnostic.py` | quelle persistance collective et quelle enveloppe spectrale single-$D$ subsistent entre blocs critiques ? | énumération exponentielle sur petits tores ; tout cutoff de blocs est signalé comme biais de sélection |
| `nested_projection_l2_diagnostic.py` | où les projections collapsed dissipent-elles l'énergie ? | énumération exacte à $`L=4`$ |
| `two_step_l2_population_diagnostic.py` | les cellules à deux updates sont-elles enrichies près du critique ? | population finie, cellules recouvrantes |
| `two_step_projective_l2_cell.py` | une cellule possède-t-elle une marge sur les potentiels atteints ? | witness exact, marge uniforme nulle au bord |
| `critical_pair_path_geometry.py` | comment explorer symétriquement le corridor d'une paire ? | brique géométrique, pas grande déviation |

L'ordre de lancement est impératif :

1. utiliser le [pilote SBM](../active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md)
   pour contre-auditer la coupe partagée, sans attribuer $`d\theta^2`$ au
   choix de $`\beta_c`$ ;
2. auditer les trois régimes et le
   [port global fini](../active/39_PORT_GLOBAL_SBM_RECOVERY.md) ;
3. tester à $`p=0.81`$ l'enveloppe spectrale à
   [un dendrogramme fixé](../active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md) ;
4. passer à la
   [double géante triangulaire](../active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md)
   si cette enveloppe reste macroscopique, en conservant les produits signés.

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/sbm_broadcast_density_evolution.py \
  --degree 3 --lambdas 0.8 0.95 1 1.05 1.2 \
  --depth 30 --particles 50000 --batches 8 --seed 20260726
python3 \
  research/hierarchical-swendsen-wang/computations/sbm_critical_cut_replica_diagnostic.py \
  --degree 3 --theta 0.5
python3 \
  research/hierarchical-swendsen-wang/computations/sbm_global_port_convolution.py \
  --root-sizes 3 2 1 --a 4 --b 1
python3 \
  research/hierarchical-swendsen-wang/computations/sbm_recovery_regimes_diagnostic.py \
  --n 100000 --a 30 --b 10 \
  --log-within-coefficient 9 --log-between-coefficient 1
python3 \
  research/hierarchical-swendsen-wang/computations/triangular_recovery_regimes_diagnostic.py \
  --vertices 1000000 --p 0.81
python3 \
  research/hierarchical-swendsen-wang/computations/giant_component_quotient_diagnostic.py \
  --sides 16,32,64 --repetitions 20 --pairs 100 \
  --p 0.809439 --distance-fraction 0.25 --seed 20260726
python3 \
  research/hierarchical-swendsen-wang/computations/critical_cut_collective_gibbs_diagnostic.py \
  --side 4 --repetitions 256 --p 0.81 \
  --maximum-block-count 16 --seed 20260726
python3 \
  research/hierarchical-swendsen-wang/computations/double_giant_replicated_gibbs_diagnostic.py \
  --side 4 --p 0.81 --observations 4 --replica-pairs 8 --seed 3801
```

Le premier JSON affiche $`\widehat q_t`$, l'écart de Nishimori, les bornes
$`\ell_t,r_t`$ et un statut de cohérence empirique. Le seuil du broadcast
est certifié par les bornes déterministes, pas par ce statut Monte-Carlo.
Le deuxième vérifie exactement qu'à $`d=3,\theta=1/2`$ le facteur marginal
vaut $`0.75`$ tandis qu'une coupe commune donne $`1.125`$. Le troisième
contracte le port global full-$D$ et audite la convolution contre les huit
orientations directes des racines de tailles $`(3,2,1)`$ ; les facteurs
internes communs sont omis. Le quatrième calcule l'affinité oracle finie et
refuse explicitement d'en déduire une achievability. Le cinquième calcule
l'erreur oracle triangulaire exacte et les seules conditions nécessaires
qu'elle implique. Le sixième conserve les rangs postcritiques réels, compte
les environnements sans paire admissible et étiquette explicitement PATH-FAC
comme non-preuve. Le septième énumère les orientations collectives exactes,
leur matrice spectrale pondérée et tous les audits de trace ; si le cutoff
exclut une réalisation, ses moyennes sont conditionnelles. Le huitième mesure
le reste signé à deux dendrogrammes et clusterise les erreurs standards du
résumé par observation.

Le [premier protocole à $`p=0.81`$](../diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md)
donne, à $`L=4`$,

```text
single-D lambda_max / n       = 0.9507358532 +/- 0.0045624262
two-D signed off-diagonal     = 0.1998059185 +/- 0.0116327206
```

La première quantité est défavorable et la seconde reste positive en
moyenne. Trois des 32 restes signés sont néanmoins négatifs. Ces nombres
valident les objets calculés, pas leur limite.

### Diagnostics et no-go

| module | verdict principal |
|---|---|
| `multiport_blackwell_counterexample.py` | la criticalisation multiport uniforme est fausse |
| `kruskal_fusion_t2_transfer.py` | l'inversion persiste dans une cellule T2-Kruskal exacte |
| `real_rank_t2_deficit_prototype.py` | l'état fidèle donne $`\lvert U\rvert=K`$ et un déficit nul |
| `last_use_attachment_palm_diagnostic.py` | la dernière incidence libère rarement l'orientation assez tôt |
| `lca_palm_corridor_diagnostic.py` | distingue Palm pré-saut, fusion réalisée et faux poids $`m^2N_\rho`$ |
| `corridor_t2_signature_diagnostic.py` | mesure les petites attaches sans certifier le screening |
| `ancestral_polarization_palm_diagnostic.py` | mesure les messages ancestraux sans les convertir en déficit |
| `joint_real_rank_t2_palm_diagnostic.py` | combine rang réel, attache et message sur les mêmes nœuds |

### Briques et oracles historiques

| module | rôle conservé |
|---|---|
| `ancestral_lambda_estimation.py` | moments pondérés et certificat de queue des ancêtres |
| `hierarchical_flip_probabilities.py` | probabilités de flip racine, feuille et nœud interne |
| `joint_hierarchical_sweep.py` | sweep exact top-down/bottom-up sur petits tores |
| `favorable_time_comparison.py` | comparaison Blackwell à taille fixée et contre-tests cross-size |
| `pair_favorability_diagnostic.py` | comparaison critique/tardive par classes de paires |
| `collapsed_corridor_transfer.py` | surrogate produit mono-bit, distinct du corridor réel |
| `triangular_band_collapsed_certificate.py` | secteur répliqué E1+ d'une cellule neutre |
| `twisted_feynman_kac_composition.py` | algèbre finie des transferts tordus |

Chaque module appelé par la recherche possède un fichier `test_*.py`
associé.

Le module `critical_component_boundary.py` contient aussi le contre-audit
des bilans résiduels et les nouveaux calculs conditionnés par une coupe :
moments du vote instantané, charge de Chernoff, fiabilité $`L^2`$ et taux de
fusion $`m u_ps_p(\beta)`$. La fiabilité est recalculée indépendamment à
partir des deux expériences binomiales symétriques dans les tests.

## Autres calculs auxiliaires conservés

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

## Surrogate factorisé mono-bit à $`p=0.8`$

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

Ces nombres valident l'énumération d'une expérience produit mono-bit. Ils ne
représentent ni le corridor collapsed multiport ni la loi du tore
triangulaire.

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
squelette et ne constitue pas une domination de Blackwell. Le corridor final
calcule aussi le rang tronqué $`q_v\mapsto\min(q_v,q_c)`$, mais celui-ci est
un proxy algébrique, pas une majoration du transfert multiport.

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
[fichier 28](../diagnostics/finite_volume/28_FIRST_CORRIDOR_P0805_RESULTS.md).

## Contre-exemples multiports à $`p=0.805`$

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/multiport_blackwell_counterexample.py
python3 \
  research/hierarchical-swendsen-wang/computations/kruskal_fusion_t2_transfer.py
```

Le premier script certifie par fractions et intervalles

```text
R(-- ) in [-0.071225876442769, -0.071225876442768]
late-minus-critical variance in [0.006261909458020, 0.006261909458020]
```

pour un canal à deux relations et une gagnante marginalisée. Le second
reconstruit une fusion cible puis une attache ancestrale avec les facteurs
$`\Lambda e^{(1-\beta)\Lambda}`$ exacts. Pour le bord $`B=4,J=3`$, il donne

```text
critical=0.735112203 late=0.755637535 gap=+0.020525332
```

Deux voies indépendantes, arête par arête et par comptes groupés, coïncident.
Ces calculs réfutent la domination uniforme ; ils ne réfutent pas une
contraction annealed sous la loi réelle des messages.

## Signatures T2 et polarisation ancestrale sous Palm

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/corridor_t2_signature_diagnostic.py \
  --sides 12 --repetitions 40 --p 0.805 \
  --distance-fraction 0.25 --maximum-bucket-size 6 \
  --maximum-charge 1 --maximum-ports 6 \
  --maximum-attachment-size 4 --seed 21260722
python3 \
  research/hierarchical-swendsen-wang/computations/ancestral_polarization_palm_diagnostic.py \
  --side 12 --repetitions 30 --p 0.805 \
  --distance-fraction 0.25 --maximum-bucket-size 8 \
  --maximum-charge 1 --message-thresholds 1,2,4 --seed 21260725
python3 \
  research/hierarchical-swendsen-wang/computations/joint_real_rank_t2_palm_diagnostic.py \
  --sides 8,12,16 --repetitions 24,10,5 --p 0.805 \
  --distance-fraction 0.25 --maximum-bucket-size 8 \
  --maximum-attachment-size 4 --maximum-charge 1 \
  --message-thresholds 1,2,4 --seed 20260723
```

Aux trois tailles testées, le premier diagnostic observe davantage de petites
attaches en peigne lorsque $`L`$ augmente, tandis qu'un cap fixe sur la
frontière globale les écrase. Cette tendance de volume fini n'est pas une
loi d'échelle. Le second calcule le vrai message externe de tous les ancêtres stricts aux rangs
réalisés. Un message borné n'est pas un certificat de screening. Les tables
et les erreurs par environnement sont dans le
[fichier 29](../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md).

Le troisième diagnostic réunit ces filtres sur les mêmes nœuds et utilise en
premier lieu la charge au **rang réel**. Le nombre moyen de candidats passe
de $`4.476\pm0.205`$ à $`L=8`$ à $`10.213\pm0.693`$ à $`L=16`$ ; après
$`|B|\le2`$, il passe de $`2.275\pm0.151`$ à $`5.763\pm0.637`$. Le
proxy criticalisé donne presque les mêmes comptes, donc le signal ne dépend
pas de l'oracle réfuté. Le module ne calcule encore ni screening, ni déficit
T2, ni borne de weak recovery.

## No-go du déficit local complètement résolu

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/real_rank_t2_deficit_prototype.py \
  --side 8 --p 0.805 --distance-fraction 0.25 \
  --maximum-bucket-size 8 --maximum-attachment-size 4 \
  --seed 20260724
```

Le module extrait une petite attache d'un corridor réel, conserve ses rangs,
les quatre poids du potentiel extérieur et les deux répliques dans le même
environnement. Si l'état cible conserve les configurations de spins
complètes, le twist $`\epsilon`$ est mesurable depuis la sortie. On a donc
exactement $`U=\epsilon K`$, $`|U|=K`$ et un déficit nul sur chaque
transition.

La sortie reproductible sélectionne un bucket $`m=5`$ avec une attache de
taille un au rang $`q=0.3688343\ldots`$ :

```text
faithful_feynman_kac_envelope: 1.0
faithful_logarithmic_attenuation: 0.0
projected_feynman_kac_envelope: 0.9917743479761344
projected_logarithmic_attenuation: 0.008259669371149435
projected_boundary_is_markov_closed: false
composable_t2_deficit_certified: false
```

La projection agrège bien des signes opposés, mais oublie les orientations
nécessaires pour mettre à jour le potentiel ancestral. Son déficit positif
n'est donc pas composable. Une jauge Markov-fermée reste un contre-test ; la
voie active passe aux projections collapsed $`L^2`$, qui n'exigent pas de
fermeture locale de dimension fixe.

## Dernière utilisation des orientations

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/last_use_attachment_palm_diagnostic.py \
  --sides 8,12,16 --repetitions 24 --p 0.805 \
  --distance-fraction 0.25 --maximum-bucket-size 8 \
  --maximum-attachment-size 4 --ancestor-windows 0,1,2,4,8 \
  --seed 20260725
```

Le critère est un certificat structurel conservateur : en prenant la branche
principale comme jauge, une attache peut être marginalisée après le dernier
bucket physique qui lui est incident. Des contributions incidentes peuvent
s'annuler, donc ce rang majore le dernier usage fonctionnel sans
nécessairement lui être égal. Sous la Palm des événements réalisés, pondérée
par $`N_\rho`$, cette dernière incidence a une profondeur moyenne
$`2.153`$, $`3.961`$, puis $`7.551`$ pour $`L=8,12,16`$. Le mécanisme local
n'est donc pas vide.

Il ne marginalise toutefois pas l'orientation globale qui porte le twist de
la paire. La dernière incidence possible de cette orientation est à la
racine dans des fractions $`0.940`$, $`0.939`$, puis $`0.962`$ ; à $`L=16`$,
seulement $`0.106`$ des unions n'ont plus d'incidence après les huit niveaux
suivants.
L'élimination locale de l'attache réduit l'état, mais ne rend pas composable
le déficit projeté. Le diagnostic ne prouve ni fermeture Markov, ni
contraction, ni résultat de weak recovery.

## Dissipation $`L^2`$ globale sur le tore $`L=4`$

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/nested_projection_l2_diagnostic.py \
  --side 4 --repetitions 96 --p 0.805 \
  --distance-fraction 0.5 --seed 20260726
```

Le module énumère les $`2^{16}`$ états de la loi conditionnelle $`\pi_D`$,
prend indépendamment du dendrogramme une paire à distance maximale $`2`$ sur
ce tore, puis intègre des paquets croissants d'orientations sur ses deux
bras. Le choix antérieur $`\texttt{distance-fraction}=0.25`$ ne filtrait en
réalité que les paires distinctes à $`L=4`$ ; il ne doit pas être décrit
comme un test de paire lointaine. Le calcul vérifie
indépendamment à chaque étape

```math
\|M_{k-1}\|_2^2-\|M_k\|_2^2
=
\|M_{k-1}-M_k\|_2^2.
```

À $`p=0.805`$, 94 des 96 paires sont dans une même racine finale. Sur ces
94 environnements, la persistance collapsed moyenne vaut
$`0.82045\pm0.02732`$, contre $`0.79657`$ pour le carré de la moyenne
postérieure. La perte totale vaut $`0.17955`$ en moyenne mais seulement
$`0.04372`$ en médiane. Conditionnellement à une perte non nulle, le nombre
effectif de paquets vaut $`1.547`$ et le paquet dominant porte en moyenne
$`79.8\%`$ de la perte.

Le second paquet strictement pré-LCA a un ratio énergétique agrégé
$`0.02837`$, mais sa perte relative médiane n'est que
$`1.22\times10^{-5}`$. Les $`10\%`$ de cas les plus dissipatifs portent
$`78.0\%`$ de sa perte absolue. Même la prédominance pré-LCA doit être lue
avec prudence : dans 54 cas sur 94, le paquet du LCA n'ajoute aucun générateur
indépendant et sa perte est donc structurellement nulle.

Le pivot $`L^2`$ voit ainsi une cancellation réelle que l'enveloppe fidèle
$`|U|=K`$ perd, mais D2 est un avertissement contre une accumulation
multiscalaire brute. Le seul test encore justifié est d'isoler deux mises à
jour emboîtées et de déterminer si la queue dissipative correspond à une
classe de bords atteignable et certifiable. L'erreur pythagoricienne maximale
vaut $`1.67\times10^{-15}`$ ; aucune extrapolation en $`L`$ n'est faite.

## Cellule D1 exacte à deux projections

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/two_step_projective_l2_cell.py
```

Le module reconstruit un witness réel $`L=4`$, $`p=0.805`$, choisi après un
scan exploratoire. Les deux nœuds strict-arm consécutifs ont des rangs
réalisés $`0.19324`$ et $`0.20258`$ ; leurs trois bits de cluster ont les
tailles $`1,4,1`$. Les pertes exactes sont

```math
\|M_0\|_2^2-\|M_1\|_2^2=0.0308734,
\qquad
\|M_1\|_2^2-\|M_2\|_2^2=0.1322561,
```

soit $`\alpha_2=0.1364694`$. Le nouveau sibling ne contient aucun endpoint :
le projeter seul sur $`f_{ij}`$ donne une perte nulle. La seconde perte est
donc bien calculée sur la fonction $`M_1`$ propagée.

Le calcul conditionne ensuite sur tous les cosets extérieurs positifs. Il
obtient 128 potentiels projectifs réellement atteints et audite la
factorisation du poids postérieur complet à $`1.43\times10^{-14}`$ près. Les
64 potentiels strictement positifs portent $`94.8\%`$ de la masse
postérieure ; leur ratio énergétique agrégé vaut $`0.14420`$ et leur marge
relative minimale $`0.00303`$. Les 64 potentiels de bord font toutefois
tomber la marge globale à zéro.

D1 prouve donc que le mécanisme à deux niveaux n'est pas algébriquement vide.
Il ne prouve pas qu'il soit fréquent : la paire du witness est à distance
$`1`$, la cellule est très précritique et elle a été sélectionnée après
exploration. Aucun poids Palm ni aucune abondance asymptotique n'est inféré.

## Audit de population D1 sans sélection

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/two_step_l2_population_diagnostic.py \
  --repetitions 32 --p 0.805 --distance-fraction 0.5 \
  --critical-rank-window 0.02 --seed 20260729
python3 \
  research/hierarchical-swendsen-wang/computations/two_step_l2_population_diagnostic.py \
  --repetitions 64 --p 0.805 --distance-fraction 0.5 \
  --critical-rank-window 0.02 --seed 20260730
```

Pour chaque dendrogramme, la paire à distance maximale est tirée avant de
regarder ses fusions. Le module énumère ensuite toutes les paires de fusions
strict-arm consécutives et recalcule exactement le postérieur complet, les
potentiels extérieurs atteints et les deux identités pythagoriciennes. Les
cellules se recouvrent et ne sont pas indépendantes ; il s'agit toujours
d'un diagnostic au volume $`L=4`$, pas d'un échantillon Palm asymptotique.

Les deux graines donnent ensemble 96 paires connectées et 302 cellules. Sur
la population entière, la seconde perte est très concentrée : son ratio
énergétique vaut $`0.03722`$, sa médiane relative
$`7.87\times10^{-5}`$, et les $`10\%`$ de cellules les plus dissipatives
portent plus de $`81\%`$ de la perte. Une marge uniforme reste impossible :
$`65.2\%`$ des cellules ont un potentiel atteint de bord à marge nulle.

Le signal nouveau est sa localisation en rang. Dans la fenêtre

```math
|q_{\mathrm{sup}}-q_c|\le0.02,
```

on ne trouve que 14 cellules sur 302. Elles portent $`4.13\%`$ de l'énergie
entrante, mais $`34.12\%`$ de la seconde perte absolue ; leur ratio
énergétique agrégé vaut $`0.30779`$. Pour la fenêtre de largeur $`0.05`$, 42
cellules portent $`58.35\%`$ de la perte, avec un ratio $`0.16897`$.

Les deux graines sont variables : la fenêtre de largeur $`0.02`$ contient
respectivement 10 cellules sur 102 et 4 sur 200, portant $`56.6\%`$ puis
$`17.6\%`$ de la perte. Le résultat ne prouve donc ni une densité positive,
ni une contraction thermodynamique. Il falsifie toutefois l'idée que le rang
critique est sans rapport avec la queue dissipative et remplace le chantier
annulaire générique par un test beaucoup plus étroit : démontrer une
occupation **énergétique** répétée de cellules consécutives dans une fenêtre
near-critical.

## Certificats rationnels less-noisy

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/rational_a0_less_noisy_certificate.py \
  --candidate a0
python3 \
  research/hierarchical-swendsen-wang/computations/rational_a0_less_noisy_certificate.py \
  --candidate p809439
```

Le module conserve aussi les jalons $`p=0.809`$ et $`p=0.8094`$ comme
contre-tests lisibles. Il n'effectue aucun scan flottant. Pour chaque
candidat, il construit
exactement les quatre polynômes nécessaires, compte leurs racines par les
suites de Sturm, vérifie les signes aux extrémités et contrôle l'identité de
dominance diagonale de la matrice polarisée. La sortie renforcée est

```text
candidate: p809439
status: CERTIFIED_PSD
scope: exhaustive
unresolved_regions: 0
variance_gap: 1/50000000
order_slack: 7/500000000
self_dual_slack: 7/500000000
```

Les tests reconstruisent aussi les formes quadratiques du canal physique et
du canal auxiliaire sur l'intérieur et les faces du simplexe. Le certificat
local exact et sa conséquence globale sont détaillés dans les [fichiers
31](../archive/certificates/31_CERTIFICAT_RATIONNEL_A0.md),
[32](../archive/certificates/32_CERTIFICAT_RATIONNEL_P809.md) et
[34](../results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).

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

Pour chaque transfert positif levé, le module construit

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

Pour une suite non normalisée, les fonctions rétrogrades
$`h_{r-1}=K_rh_r`$ construisent automatiquement un transformé de Doob
inhomogène. Les facteurs diagonaux télescopent et les lignes de masse nulle
restent hors support.

Pour la cellule E1+ à $`p=0.805`$, la sortie contient :

```text
depth= 2 signed=0.0649753038062 FK=0.0738919329503 uniform=0.0864323475826
depth=10 signed=1.89285427006e-07 FK=1.08695758136e-06 uniform=4.82371394009e-06
```

L'identification de ces transferts non normalisés au corridor LCA-Palm réel
et le contrôle thermodynamique de la mesure tuée restent ouverts. La
normalisation finie abstraite ne l'est plus.

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

## Diagnostic historique HF-S2 sur petits tores

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
Le contre-exemple multiport montre désormais que ce diagnostic ne teste pas
une domination universelle. Il est conservé comme comparaison
cible-spécifique historique.

## Conventions de développement

- Toute probabilité nouvelle doit être calculée de deux façons indépendantes
  lorsque la taille d'état le permet.
- Les tests utilisent des exemples déterministes ou des graines explicites.
- Un estimateur de carré de moyenne doit enlever les termes diagonaux.
- Les deux répliques partagent la même observation physique, mais leurs
  dendrogrammes auxiliaires et leurs tirages Gibbs sont indépendants
  conditionnellement à cette observation.
- Le bucket d'une fusion contient toutes les arêtes physiques de la coupe.
- L'identité de l'arête gagnante de Kruskal est oubliée dans le dendrogramme
  non marqué.
- Les fichiers de résultats bruts ne sont ajoutés que s'ils sont nécessaires
  à une figure ou à un certificat non reproductible rapidement.

## Prochaine étape

Le [pilote SBM](../active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) sépare trois
objets : dendrogramme figé, coupe partagée et deux Gibbs entiers
indépendants. Seul le troisième donne le Jacobien marginal
$`d\theta^2`$. Le
[port global fini](../active/39_PORT_GLOBAL_SBM_RECOVERY.md) est maintenant
écrit exactement ; sa comparaison au broadcast est ouverte. Le
[transfert GSBM](../active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) porte ensuite
sur la double géante et conserve tous les facteurs postcritiques. Le
[fichier 36](../active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md) reste un diagnostic
à un dendrogramme fixé.

Dans ce transfert, le [programme distance–entropie](../active/35_DISTANCE_ENTROPIE_ERGODICITE.md)
reste un moteur analytique conditionnel : il ne sera assemblé qu'après une
marge spectrale sur l'opérateur overlap répliqué. Le
[fichier 30](../active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md)
fournit le socle opératoriel et le
[fichier 33](../active/33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) la
cellule locale. Une T2 plus riche ou une nouvelle criticalisation uniforme
ne redevient pas prioritaire.
