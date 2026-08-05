# Dictionnaire SBM → GSBM : ce qui se transporte de la preuve complète

La [note SBM/08](../SBM/08_PREUVES_COMPLETES_SEUILS.md) démontre depuis
zéro le seuil de weak recovery du SBM classique lu sur la coupe
$`\beta_c`$, puis le seuil d'almost exact recovery en faisant tendre
$`\beta`$ vers $0$. Cette note transporte la preuve ingrédient par
ingrédient vers le tore triangulaire et classe chaque pièce : **passe tel
quel**, **casse**, ou **à reconstruire**. C'est le cahier des charges
précis des routes B, C et D du
[programme](01_PROGRAMME_DE_RECHERCHE.md).

## 1. Tableau de transport

| ingrédient de SBM/08 | SBM classique | tore triangulaire | verdict |
|---|---|---|---|
| identité de Nishimori (vérité $\leftrightarrow$ réplique) | th. I.5 | valable sur tout graphe | **passe tel quel** |
| critère quadratique $`Q_n`$ à deux répliques | th. I.7 | [note 04 §6](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md) | **passe tel quel** |
| invariance de la dynamique coupée (Edwards--Sokal, coupe marginalisée) | th. I.11 | [cadre général](../hierarchical-swendsen-wang/foundations/01_MATHEMATICAL_FRAMEWORK.md) | **passe tel quel** |
| graphe gelé = percolation indépendante exacte | fait I.9 (graphe $`G(n,\frac{a-b}n)`$) | toute coupe $t$ : percolation $`q_p(t)`$, annealed | **passe — et s'améliore** |
| borne de gel Swendsen--Wang (chapitre 11) | sous-criticité du graphe gelé | $`2p-1\le q_c`$ donne $`p_{\mathrm{SW}}=0{,}673648`$ | **passe** (volume fini à formaliser, [note 04 §3](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md)) |
| limite locale arborescente $`\mathrm{PGW}(d)`$ | déf. I.16 + emprunt E1 | la grille a des cycles à toutes les échelles | **casse** |
| fermeture inférieure : second moment du census (EKPS pour $`\lambda=1`$, emprunt E2) | th. I.17 | objet purement arborescent | **casse** |
| fermeture supérieure : SDPI récursive sur l'arbre | th. I.18 | version par arêtes seulement : $`p_{\mathrm{info}}=0{,}794659`$, non fine | **casse (perte de finesse)** |
| calibration $`\beta_\chi=\beta_c^{\mathrm{geom}}\Leftrightarrow`$ seuil | I.15/I.20 : lit le seuil **exact** $`d\theta^2=1`$ | $`t_\chi=\beta_c\Leftrightarrow(2p-1)^2=q_c`$ : lit la **baseline** | **passe formellement, perd l'exactitude** |
| port global (balance / non-arêtes) | I.2/I.10 (présent, hors hiérarchie) ; convolution : [note 39](../hierarchical-swendsen-wang/active/39_PORT_GLOBAL_SBM_RECOVERY.md) | absent : a priori produit uniforme sur le tore | **passe — par disparition** |
| cadran $`\beta\to0`$ vers Glauber | II : $`\beta_{c,n}\to0`$ quand $`\lambda_n\to\infty`$ | $`\beta_c(p)\to0`$ quand $`p\to1`$ | **passe, en changeant de variable** |
| almost exact au cadran $`\beta\to0`$ | II.5--II.6 : seuil $`\lambda_n\to\infty`$ | impossible à $p<1$ fixé : $`\varepsilon_6(p)>0`$ | **casse (degré borné)** |
| exact recovery : $`(\sqrt A-\sqrt B)^2>2`$ | [SBM/05 §4](../SBM/05_ALMOST_EXACT_ET_EXACT_RECOVERY.md), éq. (4.3) — SBM/08 s'arrête à l'exposant Hellinger (II.5--II.6) | nécessite $`1-p_n=o(n^{-1/3})`$ | **à reconstruire dans le régime $`p_n\to1`$** |

## 2. Ce qui passe tel quel

### 2.1 La couche bayésienne complète

L'identité de Nishimori, le critère $`Q_L`$, l'augmentation par
dendrogramme et l'invariance de la dynamique coupée ne font aucune
hypothèse de géométrie : toute la **couche de représentation** de SBM/08
(sa partie I jusqu'au théorème I.11) est disponible sur le tore sans
modification. En particulier, la réduction à la cible répliquée

```math
Q_L(p)\to0
\quad\Longleftrightarrow\quad
\mathcal D_L^\times(p)\to0
```

([note 41](../hierarchical-swendsen-wang/active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md))
est l'analogue exact — et déjà établi côté GSBM — du bookkeeping répliqué
de SBM/08.

### 2.2 Le graphe gelé, en mieux

Sur le SBM, seule la coupe finale $`\beta=1`$ était identifiée à un graphe
d'Erdős--Rényi explicite. Sur le tore, la
[note 04 §2](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md)
donne mieux : **toute** coupe $t$ est, sous la loi annealed, une
percolation indépendante de paramètre $`q_p(t)`$, les coupes étant
couplées de façon monotone par les mêmes horloges. La borne de gel du
chapitre 11 s'applique alors et donne $`p_{\mathrm{SW}}`$ — après la
formalisation du passage en volume fini demandée par la
[note 04 §3](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md).

### 2.3 Le port global disparaît

Le SBM fini imposait de compresser balance ou non-arêtes en un port de
magnétisation ([note 39](../hierarchical-swendsen-wang/active/39_PORT_GLOBAL_SBM_RECOVERY.md)) ;
le no-go des deux géantes anticorrélées vivait là. Sur le tore
triangulaire, l'a priori est produit uniforme et l'observation est portée
par les arêtes : il n'y a **pas de port global**. Seule subsiste la
symétrie de flip $`\pm1`$, absorbée par le travail en corrélations de
paires. C'est une simplification structurelle réelle du GSBM.

## 3. Ce qui casse

### 3.1 La limite locale arborescente

Toute la fermeture de SBM/08 — récursion exacte du broadcast, second
moment du census, SDPI récursive — vit sur l'arbre
$`\mathrm{PGW}(d)`$ obtenu comme limite locale. Le tore triangulaire n'a
pas de limite arborescente : le voisinage de tout sommet contient des
triangles, et les cycles persistent à toutes les échelles. **La fermeture
arborescente (déf. I.16 et théorèmes I.17--I.18) ne se transporte pas.**
Le remplacement doit être
planaire : box-crossing et estimées proche-critiques pour la percolation
par arêtes (pont Grimmett--Manolescu, route C), plus un contrôle de
corridor pour les facteurs postcritiques.

### 3.2 La calibration perd l'exactitude

C'est la perte la plus instructive. Sur l'arbre, la tensorisation par
arête de la SDPI est **fine** : la calibration
$`\beta_\chi=\beta_c^{\mathrm{geom}}`$ lit le seuil exact
$`d\theta^2=1`$. Sur la grille, la même calibration existe
($`t_\chi(p)`$, [note 04 §7](../hierarchical-swendsen-wang/foundations/04_TRIANGULAR_GSBM.md))
mais ne lit que $`p_{\mathrm{info}}=0{,}794659`$, strictement sous la
valeur conjecturée $`0{,}8358`$. L'écart est exactement ce que la
tensorisation par arête perd en présence de cycles ; le canal
triangulaire de la [note 34](../hierarchical-swendsen-wang/results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md)
en récupère une partie ($`+0{,}0148`$) en groupant trois arêtes. La
morale du dictionnaire : sur la grille, **le seuil exact n'est pas une
calibration par arête — c'est un objet multi-arêtes**, et le dendrogramme
est le candidat naturel pour l'exprimer.

### 3.3 L'almost exact à $p$ fixé

Le cadran $`\beta\to0`$ de SBM/08 exigeait $`\lambda_n\to\infty`$ — un
degré divergent. Le tore a un degré fixe $6$ : même avec les six voisins
révélés, l'erreur locale optimale vaut

```math
\varepsilon_6(p)
=
\sum_{k=0}^{2}\binom6kp^k(1-p)^{6-k}
+10p^3(1-p)^3
>0
```

([CURRENT_STATUS §7](../hierarchical-swendsen-wang/CURRENT_STATUS.md)),
soit $`0{,}0505`$ à $`p=0{,}81`$. L'almost exact recovery est donc
**impossible pour tout $`p<1`$ fixé** ; le régime pertinent est
$`p_n\to1`$, et l'exact recovery exige $`1-p_n=o(n^{-1/3})`$.

## 4. Le parallèle des cadrans

Dans SBM/08, le passage weak $\to$ almost exact se lit sur la coupe :
$`\beta_{c,n}^{\mathrm{geom}}\le2/\lambda_n\to0`$, et la dynamique coupée
dégénère en Glauber. Sur le tore, le même cadran existe avec $p$ pour
variable :

```math
\beta_c(p)
=
-\frac{\log\left(1-q_c/p\right)}{u_p}
\xrightarrow[p\to1]{}
0,
\qquad
u_p\to\infty.
```

Quand la fiabilité tend vers $1$, la coupe critique s'écrase vers $0$ et
la hiérarchie dégénère de la même façon vers une dynamique de type
Glauber. La table des régimes se transporte donc ainsi :

| régime | SBM (variable $`\lambda_n`$) | GSBM (variable $`p_n`$) |
|---|---|---|
| weak recovery | $`\lambda>1`$, $`\beta_c^{\mathrm{geom}}`$ fixe | $`p>p_c^{\mathrm{WR}}`$, $`\beta_c(p)`$ fixe |
| almost exact | $`\lambda_n\to\infty`$, $`\beta_{c,n}\to0`$ | $`p_n\to1`$, $`\beta_c(p_n)\to0`$ |
| exact | $`(\sqrt A-\sqrt B)^2>2`$ | $`1-p_n=o(n^{-1/3})`$ (nécessaire) |

La colonne GSBM de la dernière ligne n'a pas de condition suffisante
consignée : c'est un chantier ouvert du régime $`p_n\to1`$, subordonné à
la route D.

## 5. Synthèse

```math
\boxed{
\text{couche de représentation : transportée intégralement}
\ ;\quad
\text{couche de fermeture : à reconstruire en planaire}
\ ;\quad
\text{achievability : à créer.}
}
```

La partie de SBM/08 qui reposait sur des identités bayésiennes et sur la
percolation exacte des coupes est déjà acquise sur le tore — souvent sous
une forme plus propre. La partie qui reposait sur l'arbre (récursions,
second moment, SDPI fine) est exactement ce que les routes B et C doivent
remplacer par de la percolation planaire proche-critique, et la borne
inférieure d'achievability n'a pas encore d'analogue du tout.
