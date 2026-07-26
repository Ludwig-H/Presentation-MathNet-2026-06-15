# Statut scientifique actuel

**Dernière mise à jour : 26 juillet 2026.** Cette page est la source de
vérité du projet. Les anciennes feuilles de route restent consultables dans
[`archive/roadmaps/`](archive/roadmaps/), mais ne fixent plus les priorités.

## 1. Résultat rigoureux et programme ouvert

Le résultat quantitatif actuel est

```math
p_{\mathrm{WR}}
\ge
\frac{809439}{1000000}
=
0.809439.
\qquad\text{(1.1)}
```

Il est démontré par le
[canal triangulaire multi-état](results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md),
pas par la dynamique hiérarchique.

La cible hiérarchique prioritaire est désormais

```math
\text{deux dendrogrammes complets indépendants}
\ \longrightarrow\
\text{deux Gibbs exacts sur tous leurs arbres}
\ \longrightarrow\
\text{deux coupes à }\beta_c(p)
\ \longrightarrow\
\text{overlap sur la double géante}.
\qquad\text{(1.2)}
```

Aucun facteur postcritique n'est supprimé, contracté ou remplacé par un
canal moins informatif. La coupe critique sert de séparateur exact pour
l'élimination du Gibbs entier.

## 2. Ce que le chapitre 11 impose

Le [chapitre 11 du manuscrit](../../ChapII.tex) fournit deux contraintes
méthodologiques.

1. Une transition construite à partir de la vérité ne peut remplacer une
   réplique postérieure que si elle laisse exactement la postérieure
   invariante.
2. La borne de percolation du chapitre utilise davantage que l'invariance :
   elle exige le recoloriage indépendant et uniforme des clusters
   conditionnellement à l'objet gelé.

Le projet actif pousse ces idées plus loin. La mise à jour par le Gibbs
entier de chaque arbre final reste invariante ; après la coupe critique, ses
blocs ne sont toutefois pas indépendants lorsque leurs séparateurs sont
marginalisés. Le théorème de percolation du chapitre ne s'applique donc pas
mécaniquement.

La nouveauté recherchée n'est pas une nouvelle règle locale de recoloriage.
Elle consiste à conserver le dendrogramme complet comme variable auxiliaire,
à calculer exactement le Gibbs joint de chacun de ses arbres, puis à
répliquer indépendamment cette augmentation pour accéder au carré de la
corrélation postérieure. L'objet géométrique qui en résulte est une double
géante munie de deux corridors ancestraux ; il n'apparaît pas dans le
chapitre 11.

## 3. Verdict du pilote broadcast du SBM

Pour le SBM symétrique à deux communautés, posons

```math
d
=
\frac{a+b}{2},
\qquad
\theta
=
\frac{a-b}{a+b},
\qquad
\lambda
=
d\theta^2.
\qquad\text{(3.1)}
```

Le benchmark exact est le broadcast binaire sur
$`\mathrm{PGW}(d)`$. Le seuil exact de reconstruction est

```math
\lambda=1.
\qquad\text{(3.2)}
```

Le [pilote détaillé](active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) donne quatre
conclusions locales et isole un verrou supplémentaire sur le SBM fini.

### 3.1 Un dendrogramme figé échoue sur le broadcast

Sur un arbre local, chaque fusion est portée par une arête unique.
Conditionner par le dendrogramme complet impose alors toutes les parités
internes. Le Gibbs d'une racine finale ne conserve qu'un flip global. Cette
version voit la percolation Swendsen--Wang $`d\theta=1`$, pas
Kesten--Stigum. La contrainte de balance du SBM fini interdit de transposer
ce no-go sans un terme global supplémentaire.

### 3.2 La stratégie impose toujours la coupe à $\beta_c$

Avec

```math
p
=
\frac{1+\theta}{2},
\qquad
u
=
\log\frac p{1-p},
\qquad
q_p(\beta)
=
p(1-e^{-u\beta}),
\qquad\text{(3.3)}
```

la coupe physique critique vérifie

```math
dq_p(\beta_c)=1.
\qquad\text{(3.4)}
```

Elle produit des blocs critiques de degré moyen un et fournit l'échelle de
décomposition demandée. Elle n'est pas sélectionnée par l'information :
après marginalisation exacte, le Jacobien répliqué vaut $`\theta^2`$ pour
tout niveau de coupe.

### 3.3 Partager la coupe donne trop d'information

Si la même coupe critique est révélée aux deux copies, son transfert
quadratique vaut

```math
\eta_{\mathrm{partagée}}
=
\theta^2
+
\frac{
(1/d)(1-\theta)^2
}{
1-1/d
}
>
\theta^2.
\qquad\text{(3.5)}
```

Pour $`d=3`$ et $`\theta=1/2`$,

```math
d\theta^2=0.75,
\qquad
d\eta_{\mathrm{partagée}}=1.125.
\qquad\text{(3.6)}
```

Une hiérarchie commune créerait donc un faux régime supercritique sous le
vrai seuil.

### 3.4 Deux Gibbs entiers indépendants donnent le bon Jacobien à l'équilibre

Tirons deux répliques postérieures et deux dendrogrammes complets
indépendamment conditionnellement à la seule observation. Chaque message
Gibbs transmet $\theta$ au premier ordre ; le message overlap transmet
$`\theta^2`$. Le nombre moyen de branches vaut $d$, donc

```math
\mathcal L_{\mathrm{SBM}}^{(2)}
=
d\theta^2.
\qquad\text{(3.7)}
```

La densité d'évolution exacte du broadcast, complétée par le sandwich global

```math
\left(
\sum_{s=0}^{t}\lambda^{-s}
\right)^{-1}
\le
q_t
\le
r_t,
\qquad
r_{t+1}
=
1-e^{-\lambda r_t},
\qquad\text{(3.8)}
```

retrouve le seuil exact du broadcast, égalité comprise. Cette fermeture
emploie des bornes classiques de reconstruction ; elle n'est pas produite par
$`\beta_c`$.

Ce résultat est une validation du bookkeeping à l'équilibre de la stratégie
répliquée. Il
ne prouve pas encore :

- que le seuil temporel d'un nombre fixé de sweeps hiérarchiques est
  Kesten--Stigum ;
- le transfert du broadcast au graphe SBM fini ;
- un seuil nouveau pour le GSBM triangulaire.

Sur le SBM fini, les deux formulations usuelles ont chacune un port global.
Dans le planted bisection, la contrainte $`\sum_iX_i=0`$ couple les
orientations des racines finales. Avec des labels i.i.d., les non-arêtes
créent un potentiel proportionnel à $`(\sum_iX_i)^2`$. La
[note 39](active/39_PORT_GLOBAL_SBM_RECOVERY.md) compresse exactement ces
deux couplages en un port de magnétisation et l'élimine par convolution.
L'indépendance par arbres est fausse sur le graphe ; la comparaison
asymptotique de ce port avec le broadcast reste à prouver.

## 4. Gibbs exact sur un arbre complet

Pour une observation $O$, la mesure augmentée est

```math
\nu_O(d\sigma,dD)
=
\mu_O(d\sigma)R_O(dD\mid\sigma).
\qquad\text{(4.1)}
```

À $D$ fixé et sous l'a priori produit uniforme du GSBM triangulaire ou du
broadcast edge-only, la Gibbs se factorise entre racines finales. Cette
phrase n'inclut pas le port global du SBM fini décrit plus haut. Dans une
racine, tous les facteurs du dendrogramme restent présents. Si
$`A\in\Pi_{\beta_c}(D)`$ et
$`\partial_D^+A`$ est son ensemble de ports postcritiques, on élimine
exactement son intérieur pour produire un message

```math
Z_{D,A}(s_A)
=
\sum_{\substack{
\sigma_A\\
\sigma_{\partial_D^+A}=s_A
}}
\prod_{\substack{
u\subseteq A\\
\beta_u\le\beta_c
}}
F_{u,p}^{D}(\sigma_A).
\qquad\text{(4.2)}
```

Les messages internes sont ensuite multipliés par **tous** les facteurs
$`F_{u,p}^{D}`$ tels que $`\beta_u>\beta_c`$, puis les ports sont tirés
conjointement. Les intérieurs de blocs peuvent utiliser des hasards
indépendants seulement conditionnellement à tous leurs ports.

## 5. Décomposition exacte et cible de la double géante

Conditionnellement à $O$, tirons

```math
(\sigma^{(1)},D^{(1)}),
(\sigma^{(2)},D^{(2)})
\overset{\mathrm{i.i.d.}}{\sim}
\nu_O.
\qquad\text{(5.1)}
```

Pour $`f_{ij}(\sigma)=\sigma_i\sigma_j`$,

```math
\mu_O(f_{ij})^2
=
\mathbb E
\left[
\pi_{O,D^{(1)}}(f_{ij})
\pi_{O,D^{(2)}}(f_{ij})
\mid O
\right].
\qquad\text{(5.2)}
```

Soit $`R_{r,\star}`$ la racine géante de la réplique $r$. La seule
intersection susceptible d'être macroscopique est

```math
G_{12}^\star
=
R_{1,\star}\cap R_{2,\star}.
\qquad\text{(5.3)}
```

À la coupe critique, les cellules sont

```math
C_{A_1,A_2}
=
A_1\cap A_2\cap G_{12}^\star.
\qquad\text{(5.4)}
```

Leur contribution diagonale vérifie

```math
\sum_{A_1,A_2}
|C_{A_1,A_2}|^2
\le
\min
\left\{
\sum_{A_1}|A_1|^2,
\sum_{A_2}|A_2|^2
\right\}.
\qquad\text{(5.5)}
```

Les intersections hors double géante vérifient aussi l'inégalité
déterministe

```math
\sum_{(a,b)\ne(\star,\star)}
|R_{1,a}\cap R_{2,b}|^2
\le
\sum_{a\ne\star}|R_{1,a}|^2
+
\sum_{b\ne\star}|R_{2,b}|^2.
\qquad\text{(5.6)}
```

Les faits marginaux de percolation supercritique et critique font tendre les
membres droits de (5.5) et (5.6), divisés par $`n_L^2`$, vers zéro. La
[note prioritaire](active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) en déduit la
réduction exacte

```math
Q_L(p)
=
\mathcal E_{\mathrm{off},L}^{(2),\star}(p)
+
o(1),
\qquad\text{(5.7)}
```

où le reste est le produit **signé** des deux corrélations Gibbs entre
cellules critiques distinctes de la double géante. Ainsi, améliorer la borne
jusqu'à $`p=0.81`$ équivaut maintenant à montrer
$`\mathcal E_{\mathrm{off},L}^{(2),\star}(0.81)\to0`$.

La [note 36](active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md) reste un diagnostic
exact à un dendrogramme fixé ; sa quantité quenched est une enveloppe de
Jensen plus forte, non la cible finale.

## 6. Stratégie analytique pour la géante

Le programme ne cherche pas la loi globale de tout l'arbre géant. Il tire
une paire dans $`G_{12}^\star`$ avec le biais quadratique exact et conserve
seulement les deux corridors ancestraux complets entre ses endpoints.

L'état local contient :

- les deux partitions critiques ;
- leurs cellules d'intersection ;
- les deux systèmes de ports ;
- les rangs réels et les buckets multiports ;
- les messages extérieurs issus de tous les facteurs postcritiques.

Les deux corridors ne sont pas indépendants sous la loi annealed : ils sont
couplés par l'observation commune. L'opérateur linéarisé candidat est encore
une notation schématique,

```math
(\mathcal L_p^{(2)}h)(s_1,s_2)
=
\mathbb E_{\mathrm{Palm}}^{(2),\star}
\left[
\sum_c
J_c^{(1)}J_c^{(2)}
h(s_{1,c},s_{2,c})
\ \middle|\
s_1,s_2
\right].
\qquad\text{(6.1)}
```

Il reste à définir sa Palm jointe sans supposer déjà l'overlap petit, son
espace de messages extérieurs, sa normalisation et le produit **signé** des
Jacobiennes. Une valeur absolue après révélation des corridors détruirait les
cancellations qui réparent le no-go SBM.

La première porte numérique est donc l'enveloppe single-$D$ plus forte mais
déjà définie dans la note 36 :

```math
\mathbb E
\left[
\frac1{n_L}
\lambda_{\max}
\left(
W_{R_L^\star}^{1/2}
M_{R_L^\star}^c
W_{R_L^\star}^{1/2}
\right)
\right]
\longrightarrow0
\quad\text{à }p=0.81.
\qquad\text{(6.2)}
```

Si elle échoue, la route à deux dendrogrammes devient nécessaire. Après
construction non circulaire de $`\mathcal L_p^{(2)}`$, la question suivante
sera

```math
\rho(\mathcal L_{0.81}^{(2)})<1.
\qquad\text{(6.3)}
```

Une réponse négative arrête cette famille de blocs. Une réponse positive
avec marge demande encore une enveloppe non linéaire, puis la fermeture du
reste signé (5.7).

Le [premier diagnostic](diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md)
donne $`0.9507359\pm0.0045624`$ à $`L=4,p=0.81`$ sur 256 environnements
sans exclusion. L'enveloppe est donc très macroscopique à cette taille. Ce
résultat rend la route à deux dendrogrammes nécessaire pour le diagnostic ;
il ne détermine pas la limite de (6.2).

## 7. Almost exact et exact recovery

Sur le tore triangulaire de degré six, almost exact et exact recovery sont
impossibles pour tout $`p<1`$ fixé. Même avec les six labels voisins révélés,
l'erreur locale optimale vaut

```math
\varepsilon_6(p)
=
\sum_{k=0}^{2}
\binom6k p^k(1-p)^{6-k}
+
10p^3(1-p)^3
>
0.
\qquad\text{(7.1)}
```

Elle vaut $`0.0505275\ldots`$ à $`p=0.81`$. Une version triangulaire
non triviale de ces objectifs doit faire tendre $p$ vers un, augmenter le
degré ou répéter les observations.

Avec $`\delta=1-p`$,

```math
\varepsilon_6(1-\delta)
=
10\delta^3-15\delta^4+6\delta^5
\sim10\delta^3.
\qquad\text{(7.2)}
```

L'almost exact exige donc au minimum $`p_n\to1`$. Un packing linéaire
d'étoiles disjointes montre que l'exact recovery exige
$`n\varepsilon_6(p_n)\to0`$, donc $`1-p_n=o(n^{-1/3})`$. Ces conditions sont
seulement nécessaires. Le module
[`triangular_recovery_regimes_diagnostic.py`](computations/triangular_recovery_regimes_diagnostic.py)
audite exactement l'identité polynomiale et ces échelles.

Pour le SBM divergent, la même géométrie peut organiser les calculs mais le
lift probabiliste change :

| régime | fonctionnel pilote | cible SBM |
|---|---|---|
| weak | $`\mathbb E[M^2]`$ | $`d\theta^2=1`$ |
| almost exact | $`\mathbb E[e^{-L/2}]`$ | exposant local tendant vers l'infini |
| exact | $`n\,\mathbb P(LX\le0)`$ | $`(\sqrt A-\sqrt B)^2=2`$ dans le régime logarithmique |

Le lift Hellinger marginalise séparément les deux fonctions de partition
hiérarchiques avant leur moyenne géométrique. Ces seuils SBM sont des
benchmarks ultérieurs ; la contraction quadratique ne les démontre pas.

## 8. Portes go/no-go

| porte | test | verdict requis |
|---|---|---|
| SBM0 | Gibbs d'un dendrogramme figé sur le broadcast | reproduire le no-go Swendsen--Wang |
| SBM1 | deux Gibbs entiers et deux coupes marginalisées séparément sur le broadcast | obtenir $`d\theta^2`$ pour toute coupe |
| SBM2 | fermeture non linéaire sur le broadcast | retrouver (3.8) sans l'attribuer à la coupe physique |
| SBM-F0 | écrire et éliminer exactement balance ou non-arêtes comme port global | fermé par la note 39 |
| SBM-F1 | comparer l'overlap avec port au broadcast | encore ouvert |
| TRI0 | réduire à l'énergie inter-cellules signée | établi par (5.5)–(5.7) |
| TRI0b | tester l'enveloppe single-$D$ (6.2) | $`L=4`$ défavorable ; route double nécessaire pour ce diagnostic fini |
| TRI1 | construire la Palm jointe de deux corridors | conserver tous les ports, facteurs postcritiques et cancellations signées |
| TRI2 | certifier (6.3) | abandonner cette famille si le rayon dépasse un |
| TRI3 | fermer le régime non linéaire et la moyenne pairwise | n'annoncer un seuil qu'après cette étape |

## 9. Priorités immédiates

1. comparer au broadcast le
   [port global exact](active/39_PORT_GLOBAL_SBM_RECOVERY.md) du SBM fini ;
2. prolonger en volume le test spectral single-$D$ : à $`L=4,p=0.81`$, sa
   valeur normalisée $`0.9507\ldots`$ est défavorable ;
3. prolonger au-delà de $`L=4`$ la décomposition signée exacte : sur ce
   volume, le reste moyen est encore positif et macroscopique ;
4. mesurer conjointement le reste signé (5.7), sa fréquence de changement
   de signe et la masse diagonale du raffinement commun ;
5. définir la Palm jointe des deux corridors sans hypothèse produit ni
   circularité ;
6. estimer les produits signés de Jacobiennes sous le bon biais de paire ;
7. ne lancer la fermeture multiscalaire que si le test spectral **répliqué**
   possède une marge robuste et une enveloppe non linéaire compatible.

## 10. Références de navigation

- [cible prioritaire à deux Gibbs entiers](active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md) ;
- [pilote SBM et contre-test de coupe partagée](active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md) ;
- [port global exact et régimes de recovery du SBM fini](active/39_PORT_GLOBAL_SBM_RECOVERY.md) ;
- [premier test spectral et signé à $`p=0.81`$](diagnostics/finite_volume/40_GIBBS_CRITIQUE_RESTE_SIGNE_P081.md) ;
- [diagnostic à un dendrogramme géant fixé](active/36_ARBRE_GEANT_GIBBS_CRITIQUE.md) ;
- [cadre mathématique de la mesure jointe](foundations/01_MATHEMATICAL_FRAMEWORK.md) ;
- [critère pairwise](foundations/03_HIERARCHICAL_WEAK_RECOVERY.md) ;
- [moteur distance–entropie, désormais subordonné au test spectral](active/35_DISTANCE_ENTROPIE_ERGODICITE.md) ;
- [certificat rigoureux à $`0.809439`$](results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).
