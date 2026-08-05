# Programme de recherche : un seuil au-delà de 0,809439 par la coupe critique

## 1. Question formelle et objectifs

Sur la suite de tores triangulaires $`\mathbb T_L`$, notons

```math
Q_L(p)
=
\frac1{n_L^2}
\sum_{i,j}
\mathbb E\left[
\langle\sigma_i\sigma_j\rangle_O^2
\right],
\qquad
n_L=L^2,
```

le critère quadratique à deux répliques
([note 04 §6](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md)).
La weak recovery est possible à $p$ si et seulement si
$`\liminf_LQ_L(p)>0`$. Le seuil cherché est

```math
p_c^{\mathrm{WR}}
=
\inf\{p:\ \liminf_LQ_L(p)>0\}.
```

L'état actuel des connaissances est l'encadrement

```math
0{,}809439
\le
p_c^{\mathrm{WR}}
\le
1,
\qquad
p_c^{\mathrm{WR}}
\overset{\text{conj.}}{=}
0{,}835805792367\ldots
```

La borne inférieure est le
[certificat rationnel de la note 34](../hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) ;
la valeur conjecturée est le point multicritique de Nishimori--Ohzeki,
dont l'équation exacte est retrouvée par la
[calibration entropique des horloges](../hierarchical-swendsen-wang/diagnostics/13_NISHIMORI_HIERARCHICAL_CLOCKS.md) :

```math
3h_2(p)-h_2\left(\frac{1+(2p-1)^3}{2}\right)=1.
```

Les objectifs, par ordre de valeur croissante :

- **O1** — améliorer la borne d'impossibilité au-delà de $`0{,}809439`$
  (première cible chiffrée : $`p=0{,}81`$) ;
- **O2** — atteindre un $p$ « difficile » : toute borne d'impossibilité
  $`>0{,}82`$ dépasserait ce que les méthodes par canaux locaux ont donné ;
- **O3** — produire la **première borne d'achievability** $`p_{\mathrm{ach}}<1`$
  du dépôt, sans laquelle aucun énoncé de seuil n'est complet ;
- **O4** — le seuil exact : impossibilité et achievability se rejoignant,
  idéalement au point multicritique.

## 2. État de l'art audité

Toutes les lignes ci-dessous sont vérifiées dans le dépôt (notes citées) ;
$`q_c=2\sin(\pi/18)=0{,}347296\ldots`$

| $p$ | méthode | statut | source |
|---:|---|---|---|
| $`0{,}673648`$ | gel SW : $`2p-1>q_c`$ | rigoureux (volume fini à formaliser) | chapitre 11 ; [note 04 §3](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md) |
| $`0{,}719224`$ | triangles de Chayes--Lei | rigoureux sous hypothèses | [note 04 §3](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md) |
| $`0{,}782432`$ | majorité « vraies tardives contre fausses » | diagnostic (montre l'insuffisance du vote scalaire) | [note 14](../hierarchical-swendsen-wang/foundations/ancestral/14_CRITICAL_COMPONENT_BOUNDARY.md) |
| $`0{,}788675`$ | percolation scalaire de triangles ($`\gamma_2`$) | moins bon que la baseline — écarté | [note 04 §10](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md) |
| $`0{,}794659`$ | information-percolation par arêtes | rigoureux | [note 04 §3](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md) |
| $`0{,}809439`$ | canal triangulaire multi-état, less-noisy, Chayes--Lei | **rigoureux — record actuel** | [note 34](../hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) |
| $`0{,}809909`$ | candidat tangent $`p_\star^{\mathrm{cond}}`$ | conditionnel : marges nulles sur des faces, secteur polarisé ouvert | [note 34](../hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md) ; [note 11 §6/§8](../hierarchical-swendsen-wang/results/non_hierarchical/11_TRIANGLE_BLOCK_SDPI.md) |
| $`0{,}835806`$ | point multicritique Nishimori--Ohzeki | conjecture (première correction de dualité : $`0{,}835985`$) | [note 13](../hierarchical-swendsen-wang/diagnostics/13_NISHIMORI_HIERARCHICAL_CLOCKS.md) |
| $`0{,}847296`$ | bande pure $`(\beta_c,1]`$ seule | trop exigeant — écarté | [note 04 §7](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md) |

Côté achievability : **rien**. C'est l'asymétrie principale du problème.

## 3. Trois faits structurants

Ils délimitent ce que la coupe critique peut et ne peut pas apporter.

### 3.1 La coupe est auto-calibrée au point critique

Sous la loi annealed, la coupe au temps $t$ est une percolation
indépendante de paramètre $`q_p(t)=p(1-e^{-u_pt})`$, et

```math
q_p(\beta_c(p))=q_c
\qquad\text{pour tout }p>0{,}673648.
```

La géométrie des blocs à la coupe est donc **critique pour tout $p$** :
pas de bloc géant (rigoureux), et — sous l'universalité étoile-triangle
de Grimmett--Manolescu pour la percolation **par arêtes**, à importer
théorème par théorème — box-crossing et estimées proche-critiques.
C'est l'échelle exacte demandée par le problème central
([note 42](../hierarchical-swendsen-wang/foundations/ancestral/42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md)) :
la fusion d'une paire lointaine se localise en $`\beta_{ij}\to\beta_c^-`$
sous le conditionnement de gauche — résultat **esquissé sous intrants
standard** (note 42 §4.2), avec le caveat que la masse de l'événement
conditionnant s'annule : aucune borne sur $`Q_L`$ n'en découle sans
réinsérer ce facteur de masse.

### 3.2 Le transfert par arête ne dépend pas de la coupe

La contraction $`\chi^2`$ d'une arête observée vaut $`(2p-1)^2`$ quelle
que soit la position de la coupe, et la calibration temporelle
$`t_\chi(p)>\beta_c(p)\Leftrightarrow(2p-1)^2>q_c`$ redonne exactement la
baseline $`p_{\mathrm{info}}=0{,}794659`$
([note 04 §7](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md)).
Conséquence méthodologique ferme : **couper à $`\beta_c`$ ne crée aucune
information par soi-même**. Sur l'arbre du SBM, la tensorisation par arête
est exacte et cette calibration lit le seuil exact
([SBM/03](../SBM/03_PREUVE_DU_SEUIL_WEAK_RECOVERY.md)) ; sur la grille,
l'écart conjecturé $`0{,}794659\to0{,}8358`$ mesure précisément ce que la
tensorisation par arête perd. Tout gain doit venir du traitement **joint**
d'objets multi-arêtes : c'est déjà la leçon du canal triangulaire
(trois arêtes, $`+0{,}0148`$ sur la baseline), et les coupes $`E_v`$ du
dendrogramme sont des objets multi-arêtes de taille non bornée.

### 3.3 La cible est un carré moyenné en dendrogramme

La chaîne de réduction établie
([note 38](../hierarchical-swendsen-wang/active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md),
[note 41](../hierarchical-swendsen-wang/active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md))
donne, avec $`S_L^c(p)\to0`$ :

```math
Q_L(p)\longrightarrow0
\quad\Longleftrightarrow\quad
\mathcal D_L^\times(p)
=
\frac1{n_L^2}
\mathbb E_O\sum_{i,j}
d_{O,ij}^2
\longrightarrow0,
```

où $`d_{O,ij}`$ moyenne sur le dendrogramme **avant** le carré. Toutes les
enveloppes à un dendrogramme fixé (Jensen) sont saturées aux volumes
accessibles
([note 40](../hierarchical-swendsen-wang/diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md)) :
le seul objet discriminant est la moyenne signée. La
[première mesure directe](03_EXPERIENCE_CIBLE_REPLIQUEE.md) montre que
cette moyenne annule déjà $`81`$ à $`94\,\%`$ de l'enveloppe de Jensen à
$`L=4`$ : le pari quantitatif de la route hiérarchique est que cette
annulation devient totale quand $`L\to\infty`$ pour
$`p<p_c^{\mathrm{WR}}`$.

## 4. Les cinq routes

Chaque route donne : objectif, verrou, porte falsifiable, premier calcul.

### Route A — consolider le tangent $`p_\star=0{,}809909`$ (non hiérarchique)

**Objectif.** Transformer le candidat algébrique
$`p_\star^{\mathrm{cond}}=0{,}8099092892\ldots`$ en certificat rationnel,
soit $`+0{,}00047`$ sur le record.

**Verrou.** La comparaison less-noisy $`P_\star`$ est prouvée pour tout
a priori vérifiant $`\max_x\mu_x\le1/2`$ (corollaire 5.2 de la note 11) ;
sur le secteur polarisé, l'échantillonnage des mineurs principaux n'a
trouvé aucun contre-exemple sans produire de preuve uniforme
([note 11 §6](../hierarchical-swendsen-wang/results/non_hierarchical/11_TRIANGLE_BLOCK_SDPI.md)).
Au point $`p_2=0{,}809439`$, la
[note 34](../hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md)
certifie les deux secteurs ; c'est l'extension jusqu'au point tangent qui
reste ouverte, avec des marges qui s'annulent sur des faces.

**Porte GA.** Produire une décomposition en carrés (SOS) rationnelle des
mineurs sur le secteur polarisé jusqu'au point tangent, ou un
contre-exemple exact. Sans l'un ni l'autre après l'exploration SOS, la
route est gelée.

**Premier calcul.** Le plafond du schéma actuel est déjà établi : la
note 34 §7 démontre que « séparation univariée + dominance diagonale +
sous-criticité stricte » s'arrête strictement avant $`0{,}80944`$ et ne
peut donc pas atteindre $`p_\star`$. Le premier calcul utile est une
certification de la matrice **sans dominance diagonale** : SOS rationnel
des mineurs sur le secteur polarisé $`\max_x\mu_x>1/2`$, par tranches de
masse dominante.

### Route B — annuler la cible répliquée à $`p=0{,}81`$ (hiérarchique, impossibilité)

**Objectif.** Prouver $`\mathcal D_L^\times(0{,}81)\to0`$. Par la
réduction du fait 3.3, cela donnerait $`p_c^{\mathrm{WR}}\ge0{,}81`$ —
première amélioration hiérarchique du record.

**Verrous.** TRI1-o : plonger les espaces de ports dépendant de $D$ dans
un état commun suffisant sans perdre l'embedding physique. TRI2 : la
contraction inhomogène

```math
\left\|
\mathcal L_{L,0}^{(2)}
\cdots
\mathcal L_{L,k-1}^{(2)}
\right\|
\le
Cr^k,
\qquad r<1,
```

sur les produits de transferts moyennés
$`\overline{\mathsf K}_O^\times=\mathbb E_{D\mid O}[\mathsf K_{O,D}^\times]`$
— le carré doit rester **après** la moyenne en $D$
([CURRENT_STATUS §6](../hierarchical-swendsen-wang/CURRENT_STATUS.md)).

**Porte GB.** La décision porte sur $`\mathcal D_L^\times(0{,}81)`$
lui-même : par la réduction du fait 3.3, sa limite est exactement
équivalente à $`Q_L\to0`$. La mesure doit décroître en $L$ sur les
volumes exacts accessibles ($`L=4`$ fait, $`L=5,6`$ par junction tree).
Le ratio $`\mathcal D_L^\times/\mathcal J_L^\times`$ n'est qu'un
diagnostic de cancellation — il peut rester positif alors même que
$`\mathcal D_L^\times\to0`$. Si $`\mathcal D_L^\times`$ ne décroît pas
hors barres d'erreur, deux lectures : ou bien le volume est trop petit
(contrôler $`S_L^c`$, encore macroscopique à $`L=4`$), ou bien — lecture
forte, à ne retenir qu'à des volumes où $`S_L^c`$ est devenu petit — la
weak recovery serait possible à $`0{,}81`$ et la conjecture multicritique
serait fausse. Dans les deux cas la mesure tranche la suite du travail.

**Premier calcul.** Fait : [expérience E1](03_EXPERIENCE_CIBLE_REPLIQUEE.md).
Suivant : E2 (junction tree exact à $`L=5`$).

### Route C — loi du squelette proche-critique sous Palm (structure)

**Objectif.** Fermer le verrou G1 de la
[note 42](../hierarchical-swendsen-wang/foundations/ancestral/42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md) :
la loi groupée du squelette (buckets $`E_v`$, temps $`\beta_v`$, tailles)
au-dessus du nœud de fusion, sous la Palm à deux points, dans la fenêtre
$`\beta\to\beta_c(p)`$.

**Verrou.** Le biais de sélection de Kruskal porte sur la géométrie du
squelette (les noyaux de marques conditionnels, eux, sont exacts,
[note 04 §7](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md)).
Sur la grille, la loi du squelette est un objet de percolation
proche-critique planaire : il faut les fonctions à deux points et les
événements pivotaux au voisinage de $`q_c`$, pour la percolation **par
arêtes** — d'où le pont Grimmett--Manolescu (étoile-triangle,
box-crossing), à importer théorème par théorème.

**Porte GC.** Écrire la loi jointe $`(m,\beta)`$ des fusions dans
$`|\beta-\beta_c|\le\varepsilon`$ conditionnellement à la fusion d'une
paire à distance $`\ell`$, avec bornes RSW explicites, et la valider
contre l'énumération exacte à $`L=4`$. Échec si la loi exige des exposants
d'arms non disponibles pour les arêtes.

**Premier calcul.** Sous la mesure exacte de
[E1](03_EXPERIENCE_CIBLE_REPLIQUEE.md), histogrammer $`(m,\beta_v)`$ des
fusions contenant une paire cross-block donnée et comparer à la
multinomiale de bande $`\mathrm{Mult}(m-1;h_p,\ldots)`$.

### Route D — achievability : la moitié manquante

**Objectif.** La première borne $`p_{\mathrm{ach}}<1`$ : pour
$`p>p_{\mathrm{ach}}`$, la weak recovery est possible.

**Approche.** Sur la ligne de Nishimori, $`Q_L(p)`$ est une fonction de
corrélation du modèle $`\pm J`$ ([note 04 §6](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md)) ;
l'équivalence précise « weak recovery $\Longleftrightarrow$ ordre à
longue portée dans la phase de Nishimori » est l'énoncé de cadrage à
établir en premier — il n'est pas encore consigné dans le dépôt.
Pour $p$ proche de $1$, un argument de contours de Peierls adapté aux
couplages $`\pm J`$ (densité $`1-p`$ de liens frustrés) doit donner un
$`p_{\mathrm{ach}}`$ explicite. Littérature à intégrer : le livre de
Nishimori et les travaux d'Ozeki--Nishimori sur le diagramme de phase
$`\pm J`$ (aucune de ces références n'est encore auditée dans le dépôt).

**Porte GD.** N'importe quel $`p_{\mathrm{ach}}<1`$ prouvé. Sans lui, tout
énoncé de « seuil » reste unilatéral. Bonus hiérarchique ensuite :
l'estimateur naturel « orientation relative des blocs critiques du
dendrogramme » hérite de l'ordre LRO ; son analyse est subordonnée à la
borne de Peierls, pas l'inverse.

**Premier calcul.** Borne de Peierls annealed : espérance du nombre de
contours de longueur $`\ell`$ dont le poids Nishimori dépasse $1$,
sommée sur $`\ell`$ ; premier régime où la somme converge.

### Route E — l'exactitude au point multicritique

**Objectif.** Identifier un fonctionnel hiérarchique dont le changement de
signe est exactement l'équation de Nishimori--Ohzeki, et le promouvoir en
critère de seuil.

**Point d'appui.** La
[note 13](../hierarchical-swendsen-wang/diagnostics/13_NISHIMORI_HIERARCHICAL_CLOCKS.md)
retrouve exactement l'équation $`3h_2(p)-h_2(\cdot)=1`$ au niveau d'une
face : le budget entropique des horloges d'une face égale l'information
qu'elle transmet. La conjecture de travail est que la version
**multi-échelle** de cette balance — entropie des horloges d'un bloc
critique contre information de ses ports — est le critère exact que la
percolation d'information par arêtes approxime par en dessous.

**Porte GE.** Dériver la balance entropie/information pour une cellule
strictement plus grande que la face (bande de largeur deux) et observer si
la racine bouge vers $`0{,}835985`$ (première correction de dualité) ou
reste à $`0{,}835806`$. Si la constante dérive sans loi, l'identification
au seuil est probablement un artefact de petite cellule et la route est
gelée.

**Premier calcul.** Réécrire l'identité de la note 13 pour la bande
$`2\times L`$ périodique par matrices de transfert exactes.

## 5. Expériences

| # | contenu | statut | où |
|---|---|---|---|
| E1 | mesure directe de $`\mathcal D_L^\times`$ à $`L=4`$, intérieur exact | **fait** | [note 03](03_EXPERIENCE_CIBLE_REPLIQUEE.md) |
| E2 | junction tree exact sur les spins physiques, $`L=5,6`$ ; puis $`\mathcal D_L^\times`$ à $`L=5`$ | à faire (priorité) | [CURRENT_STATUS §9](../hierarchical-swendsen-wang/CURRENT_STATUS.md) |
| E3 | produits signés de Jacobiennes sous la Palm cross--cross | à faire | [note 41](../hierarchical-swendsen-wang/active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md) |
| E4 | SDPI multi-état pour la bande de largeur deux (au-dessus du triangle) | à faire | prolonge la [note 11](../hierarchical-swendsen-wang/results/non_hierarchical/11_TRIANGLE_BLOCK_SDPI.md) |

E4 sert à la fois la route A (meilleur canal local) et la route E
(dépendance de la balance en la taille de cellule).

## 6. Risques et obstacles identifiés

1. **Pas de limite locale arborescente.** La fermeture de
   [SBM/08](../SBM/08_PREUVES_COMPLETES_SEUILS.md) (récursion exacte sur
   l'arbre de Poisson--Galton--Watson) n'a aucun analogue : la grille a
   des cycles à toutes les échelles. Le remplacement est le couple
   RSW/corridor planaire — significativement plus dur.
2. **La coupe seule ne crée rien** (fait 3.2). Toute rédaction qui
   attribuerait un gain au seul choix de $`\beta_c`$ serait fausse.
3. **Les enveloppes quenched sont saturées** : plafond structurel
   $`\mathbb E[|R^\star|]/n_L\approx0{,}997`$ à $`L\le6`$
   ([note 36](../hierarchical-swendsen-wang/active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md)).
   Seules les quantités doublement moyennées discriminent.
4. **L'inégalité de transfert est vide aux petits volumes** :
   $`S_L^c\approx0{,}75`$ à $`L=4`$, donc $`2\sqrt{S_L^c}>1`$ — la mesure
   E1 est un diagnostic de cancellation, pas une borne sur $`Q_L`$.
5. **Biais de Kruskal** sur la géométrie du squelette : seules les lois de
   marques conditionnelles sont exactes ; toute loi de squelette doit être
   démontrée, pas supposée.
6. **Signes** : les produits de Jacobiennes ne doivent être ni normalisés
   ni pris en valeur absolue ; les poids Palm restent les seules
   probabilités ([note 41](../hierarchical-swendsen-wang/active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md)).
7. **La conjecture d'identification** (seuil de weak recovery = point
   multicritique) n'est pas un théorème, et la constante elle-même bouge
   sous les corrections de dualité.

## 7. Critères d'arrêt

- **A** : gel si ni SOS ni contre-exemple après l'exploration par
  tranches du secteur polarisé.
- **B** : gel à $`p=0{,}81`$ si le ratio $`\mathcal D/\mathcal J`$ ne
  décroît pas de $`L=4`$ à $`L=6`$ hors barres d'erreur ; reprise à un $p$
  plus proche de $`0{,}8358`$ interdite tant que le cas $`0{,}81`$ n'est
  pas compris.
- **C** : gel si la loi de squelette exige des exposants d'arms non
  transférés par l'universalité par arêtes.
- **D** : pas de critère d'arrêt — cette moitié est nécessaire.
- **E** : gel si la balance de bande $`2\times L`$ dérive sans structure.

## 8. Ordre de travail

1. E2 (junction tree $`L=5,6`$) puis GB : c'est le chemin critique de la
   route B et il recycle directement le module d'E1.
2. GA en parallèle (calcul algébrique indépendant, gain court possible).
3. GD (Peierls annealed) : chantier indépendant, indispensable à O3/O4.
4. GC nourri par les histogrammes d'E1/E2.
5. GE en dernier (E4 d'abord, pour voir si la balance est stable).
