# Certificat collapsed sur un cactus triangulaire critique

Cette note réalise le premier calcul fini annoncé dans les fichiers 19--20.
Elle traite une chaîne de triangles partageant seulement des sommets
d'articulation et répond exactement aux deux conditionnements suivants :

1. les deux extrémités sont déjà dans la même composante au rang de
   percolation $`q`$ ;
2. leur LCA de Kruskal fusionne précisément au rang $`q`$, au sens de la
   densité de Palm du flux de fusions.

Dans le second cas, au rang critique $`q=q_\triangle`$, on obtient un
équivalent exact de la probabilité de préserver la relation de la réplique
de référence de Nishimori après un heat bath collapsed complet. À
$`p=0.8`$, cette probabilité tend exponentiellement vers $`1/2`$ lorsque le
nombre de triangles séparant les extrémités tend vers l'infini.

Ce résultat est un théorème sur le cactus. Il ne constitue pas encore une
borne pour la grille triangulaire : le cactus supprime précisément les
cycles chevauchants qui créent l'état de bord difficile.

Le mot « critique » désigne ici le rang $`q_\triangle`$ importé de la grille
triangulaire. La chaîne de cactus est quasi unidimensionnelle et son propre
seuil de percolation infinie vaut $`1`$. Elle ne possède donc pas de
composante géante à $`q_\triangle`$ : le calcul porte sur le canal
conditionnel favorable, pas sur la géométrie typique de la géante critique.

## 1. Modèle et hiérarchie non marquée

Le cactus $`G_h`$ possède les sommets d'articulation

```math
a_0,a_1,\ldots,a_h
```

et un sommet latéral $`b_r`$ dans chaque bloc $`1\le r\le h`$. Le triangle
$`r`$ a pour arêtes

```math
(a_{r-1},a_r),
\qquad
(a_{r-1},b_r),
\qquad
(b_r,a_r).
\tag{1.1}
```

Deux blocs consécutifs ne partagent que $`a_r`$. En particulier, tout chemin
de $`a_0`$ à $`a_h`$ doit connecter $`a_{r-1}`$ à $`a_r`$ dans chacun des
$`h`$ triangles.

Sous la loi jointe annealed de Nishimori, utilisons la coordonnée de rang du
lemme 4.3 du fichier 20. Chaque arête reçoit un rang indépendant
$`R_e\sim\mathrm{Unif}[0,1]`$ tel que :

- elle est conforme à la réplique plantée si $`R_e\le p`$ ;
- elle est ouverte au rang $`q\le2p-1`$ si $`R_e\le q`$.

La correspondance avec le niveau d'horloge demandé dans les slides est

```math
q=q_p(\beta)=p(1-e^{-u_p\beta}),
\qquad
u_p=\log\frac p{1-p},
\tag{1.2}
```

donc

```math
\beta_c(p)
=
-\frac1{u_p}\log\left(1-\frac{q_\triangle}{p}\right).
\tag{1.3}
```

Ce niveau appartient à $`[0,1]`$ exactement lorsque
$`q_\triangle\le2p-1`$, soit
$`p\ge(1+q_\triangle)/2=0.673648\ldots`$. Dans toute la note, comparer les
rangs $`q`$ revient à comparer les niveaux $`\beta`$, car $`q_p`$ est
strictement croissante.

La « conformité » de cette note a un sens précis. On part d'une réplique
$`\sigma^{(0)}\sim\mu_O`$, on construit $`D`$ depuis cette réplique, puis le
heat bath collapsed produit $`\sigma^{(1)}`$. On mesure

```math
\sigma^{(1)}_{a_0}\sigma^{(1)}_{a_h}
=
\sigma^{(0)}_{a_0}\sigma^{(0)}_{a_h}.
\tag{1.4}
```

Par Nishimori, $`\sigma^{(0)}`$ peut servir de vérité plantée dans cette loi
jointe augmentée. Mais (1.4) n'est pas la probabilité quenched qu'un
estimateur, à $`O`$ fixé, retrouve la vérité originale. Son rôle exact est
de calculer la persistance répliquée qui entre dans le critère $`L^2`$ de
weak recovery du fichier 20.

La hiérarchie $`D`$ est le dendrogramme de partitions **non marqué** : elle
retient les deux composantes fusionnées et le rang de fusion, mais pas
l'identité de l'arête gagnante dans un bucket de taille deux.

Notons $`\mathcal C_h(q)`$ l'événement

```math
a_0\longleftrightarrow a_h
\quad\text{au rang }q.
```

Par la structure d'articulation,

```math
\mathcal C_h(q)
=
\bigcap_{r=1}^h
\{a_{r-1}\longleftrightarrow a_r
\text{ dans le triangle }r\}.
\tag{1.5}
```

Les événements de droite portent sur des ensembles d'arêtes disjoints et
sont donc indépendants.

En particulier,

```math
\mathbb P(\mathcal C_h(q))=c(q)^h,
\tag{1.6}
```

avec $`c(q)<1`$ pour $`q<1`$. Le conditionnement étudié devient donc
exponentiellement rare sur le cactus. C'est volontaire : il isole le canal
le plus favorable demandé, tandis que la probabilité et le biais d'apparition
d'une paire dans la géante de la grille restent un verrou séparé.

## 2. Décomposition exacte d'un triangle

Fixons un triangle, ses articulations gauche et droite, et appelons
**directe** l'arête qui les relie. Les deux autres arêtes forment le chemin
latéral de longueur deux.

La probabilité de connexion des articulations au rang $`q`$ vaut

```math
c(q)
=
q+q^2-q^3.
\tag{2.1}
```

En effet, elles sont connectées si l'arête directe est ouverte ou si les
deux arêtes latérales le sont.

### 2.1 L'arête directe fusionne en premier

Si l'arête directe a le premier rang du triangle, le premier bucket est de
taille un. Le facteur hiérarchique

```math
\Lambda e^{(1-\beta)\Lambda}
```

s'annule après inversion de la parité des articulations. Le dendrogramme
révèle donc parfaitement leur relation.

La masse des histoires directes-premières qui connectent avant $`q`$ est

```math
d(q)
=
\int_0^q(1-r)^2\,dr
=
q-q^2+\frac{q^3}{3}.
\tag{2.2}
```

### 2.2 Une arête latérale fusionne en premier

Supposons que $`(a_{r-1},b_r)`$ soit la première arête. Elle force la
relation entre ces deux spins par son bucket de taille un. La seconde fusion
utilise le bucket non marqué

```math
\{(a_{r-1},a_r),(b_r,a_r)\},
```

de taille deux. Si son rang est $`r`$, la densité non normalisée de cette
histoire, après sommation des deux choix de première arête et des deux choix
de gagnante du second bucket, est

```math
g(r)=4r(1-r),
\qquad 0<r<1.
\tag{2.3}
```

Sa masse avant $`q`$ vaut donc

```math
a(q)
=
\int_0^qg(r)\,dr
=
2q^2-\frac{4q^3}{3}.
\tag{2.4}
```

On vérifie exactement $`c(q)=d(q)+a(q)`$.

### Lemme 2.1 — le second bucket est un effacement, statut : établi

Conditionnellement à son rang $`r`$, l'arête non gagnante du bucket de
taille deux est conforme avec probabilité

```math
s_p(r)=\frac{p-r}{1-r}.
\tag{2.5}
```

Après marginalisation de l'identité de la gagnante :

- si cette arête est conforme, les comptes des deux parités sont $`(2,0)`$
  et la relation correcte est forcée ;
- si elle est non conforme, les comptes sont $`(1,1)`$ et les deux parités
  ont exactement le même poids.

Le transfert de la parité des articulations est donc un effacement de
fiabilité $`s_p(r)`$.

#### Preuve

Écrivons

```math
\phi_t(k)
=
u_pk\exp((1-t)u_pk),
\qquad
u_p=\log\frac p{1-p}.
\tag{2.6}
```

Le bucket de taille un fixe la relation portée par la première arête. Dans
le bucket suivant, si les deux observations sont conformes, la parité vraie
donne $`\phi_t(2)>0`$ et la parité inversée $`\phi_t(0)=0`$. Si une seule
observation est conforme, l'inversion échange les deux arêtes et donne
$`\phi_t(1)`$ dans les deux états. Enfin, conditionnellement à
$`R_e>r`$,

```math
\mathbb P(R_e\le p\mid R_e>r)
=
\frac{p-r}{1-r}.
```

L'oubli de l'identité gagnante est indispensable : avec une gagnante
marquée, l'état qui échange gagnante et résiduelle serait artificiellement
interdit.

## 3. Transfert répliqué d'un bloc connecté

Pour une réalisation du dendrogramme et des marques, la moyenne collapsed de
la parité vaut ici soit $`1`$, soit $`0`$. Son premier moment annealed et son
second moment annealed coïncident donc.

### Théorème 3.1 — coefficient connecté, statut : établi

Conditionnellement au fait que les deux articulations sont connectées avant
$`q`$, leur fiabilité répliquée collapsed est

```math
\boxed{
\kappa_{\rm conn}(p,q)
=
\frac{1+(2p-1)q-q^2}{1+q-q^2}.
}
\tag{3.1}
```

#### Preuve

La masse informative est la somme de la partie directe, parfaite, et de la
partie latérale non effacée :

```math
\begin{aligned}
n(q)
&=
d(q)+\int_0^q4r(1-r)s_p(r)\,dr\\
&=
d(q)+4\int_0^qr(p-r)\,dr\\
&=
q+(2p-1)q^2-q^3.
\end{aligned}
\tag{3.2}
```

Le quotient $`n(q)/c(q)`$ donne (3.1).

Sur l'état répliqué
$`\delta=\sigma^{(1)}_{a_{r-1}}\sigma^{(1)}_{a_r}
\sigma^{(2)}_{a_{r-1}}\sigma^{(2)}_{a_r}`$, le noyau vaut exactement

```math
K_{\rm conn}(p,q)
=
\frac12
\begin{pmatrix}
1+\kappa_{\rm conn}&1-\kappa_{\rm conn}\\
1-\kappa_{\rm conn}&1+\kappa_{\rm conn}
\end{pmatrix}.
\tag{3.3}
```

Le mode constant a valeur propre $`1`$ et le mode de parité valeur propre
$`\kappa_{\rm conn}`$.

### Contre-audit : le seul $`s_c`$ est trop optimiste

Conditionner artificiellement chaque triangle à être de type latéral-premier
donnerait le coefficient $`s_p(r)`$. Mais, sous le vrai biais de connexion,
une proportion

```math
\frac{d(q)}{c(q)}
```

des triangles révèle directement la relation et possède coefficient un. Le
coefficient physique correct est donc $`\kappa_{\rm conn}`$, strictement plus
grand que la seule fiabilité d'un bucket $`m=2`$.

## 4. Conditionnement exact du LCA au rang critique

L'événement « connecté à $`q`$ » est cumulatif. Pour formaliser
« le LCA fusionne au rang $`q`$ », il faut utiliser la densité du flux de
fusions. La désintégration au niveau $`\beta`$ donne le même quotient : le
jacobien strictement positif $`q_p'(\beta)`$ multiplie simultanément les
densités informative et totale, puis s'annule.

On a

```math
c'(q)=1+2q-3q^2,
\qquad
n'(q)=1+2(2p-1)q-3q^2.
\tag{4.1}
```

### Théorème 4.1 — coefficient pivotal, statut : établi

Sous la loi de Palm qui fixe le rang de fusion des deux articulations à
$`q`$, le coefficient du bloc pivotal vaut

```math
\boxed{
\kappa_{\rm flux}(p,q)
=
\frac{1+(4p-2)q-3q^2}{1+2q-3q^2}.
}
\tag{4.2}
```

#### Preuve

La densité directe-première vaut $`(1-q)^2`$. La densité
latérale-première vaut $`4q(1-q)`$, dont la partie informative vaut

```math
4q(1-q)s_p(q)=4q(p-q).
```

Le quotient de la somme informative par
$`(1-q)^2+4q(1-q)=c'(q)`$ donne (4.2). Équivalemment,
$`\kappa_{\rm flux}=n'/c'`$.

### Théorème 4.2 — paire distante sur $`h`$ triangles, statut : établi

Sous la densité de Palm

```math
q_{\rm LCA}(a_0,a_h)=q,
```

le second moment collapsed exact est

```math
\boxed{
A_h^{\rm LCA}(p,q)
=
\kappa_{\rm flux}(p,q)
\kappa_{\rm conn}(p,q)^{h-1}.
}
\tag{4.3}
```

La probabilité moyenne que le heat bath collapsed préserve la relation de la
réplique de référence est

```math
\boxed{
P_h^{\rm conf}(p,q)
=
\frac12
\left[
1+A_h^{\rm LCA}(p,q)
\right].
}
\tag{4.4}
```

#### Preuve

Le rang du LCA global est le maximum des $`h`$ rangs de connexion locaux.
Sa fonction de répartition est $`c(q)^h`$ et sa densité

```math
h c(q)^{h-1}c'(q).
```

La masse informative cumulative est $`n(q)^h`$ et sa densité

```math
h n(q)^{h-1}n'(q).
```

Le quotient des deux densités est précisément (4.3). Il y a presque sûrement
un unique triangle pivotal au rang maximal ; les $`h-1`$ autres sont
conditionnés seulement à être déjà connectés.

Puisque la moyenne conditionnelle de la parité relative vaut $`0`$ ou $`1`$
sur chaque réalisation, son premier et son second moment sont égaux. La
probabilité de préserver une relation de moyenne $`H`$ vaut $`(1+H)/2`$,
d'où (4.4). L'interprétation par rapport à la vérité est exactement celle de
(1.4), et non une probabilité quenched de succès d'estimation.

### Corollaire 4.3 — perte sur le cactus, statut : établi

Pour tout $`p<1`$ et tout $`q\in(0,2p-1]`$,

```math
1-\kappa_{\rm conn}(p,q)
=
\frac{2(1-p)q}{1+q-q^2}
>0.
\tag{4.5}
```

Ainsi

```math
A_h^{\rm LCA}(p,q)\longrightarrow0,
\qquad
P_h^{\rm conf}(p,q)\longrightarrow\frac12
\tag{4.6}
```

exponentiellement lorsque $`h\to\infty`$.

Plus précisément, l'équivalent demandé est ici une identité exacte :

```math
P_h^{\rm conf}(p,q)-\frac12
=
\frac{\kappa_{\rm flux}(p,q)}{2\kappa_{\rm conn}(p,q)}
\exp\{-h\gamma(p,q)\},
\qquad
\gamma(p,q)=-\log\kappa_{\rm conn}(p,q)>0.
\tag{4.7}
```

Ce corollaire porte sur la dynamique hiérarchique collapsed conditionnée par
son dendrogramme non marqué. Il ne remplace pas le théorème géométrique encore
requis sur la grille bidimensionnelle.

## 5. Le seuil critique est bien le cas favorable sur le cactus

Les deux coefficients s'écrivent

```math
\kappa_{\rm conn}(p,q)
=
1-\frac{2(1-p)q}{1+q-q^2},
\tag{5.1}
```

```math
\kappa_{\rm flux}(p,q)
=
1-\frac{4(1-p)q}{1+2q-3q^2}.
\tag{5.2}
```

Leurs dérivées en $`q`$ sont strictement négatives :

```math
\partial_q\kappa_{\rm conn}
=
-\frac{2(1-p)(1+q^2)}{(1+q-q^2)^2},
\tag{5.3}
```

```math
\partial_q\kappa_{\rm flux}
=
-\frac{4(1-p)(1+3q^2)}{(1+2q-3q^2)^2}.
\tag{5.4}
```

Par conséquent, pour $`q_\triangle\le q_1\le q_2\le2p-1`$,

```math
A_h^{\rm LCA}(p,q_2)
\le
A_h^{\rm LCA}(p,q_1)
\le
A_h^{\rm LCA}(p,q_\triangle).
\tag{5.5}
```

Le cas où la paire fusionne au seuil de percolation est donc exactement le
cas postcritique le plus favorable à la conservation de la relation de
référence sur ce cactus.

Au niveau répliqué, cette comparaison est une dégradation explicite. Pour
$`\alpha\in\{\mathrm{conn},\mathrm{flux}\}`$ et $`q_1\le q_2`$,

```math
K_{\alpha}(p,q_2)
=
K_{\alpha}(p,q_1)
K_{\rm BSC}
\left(
\frac{\kappa_{\alpha}(p,q_2)}
{\kappa_{\alpha}(p,q_1)}
\right).
\tag{5.6}
```

Ici $`K_{\rm flux}`$ désigne la matrice (3.3) avec
$`\kappa_{\rm conn}`$ remplacé par $`\kappa_{\rm flux}`$.

Cette domination intègre ici le changement aléatoire de forme du
dendrogramme. Elle est plus forte que la domination à taille fixée du fichier
20, mais repose entièrement sur la géométrie d'articulation du cactus.

## 6. Valeurs certifiées à $`p=0.8`$

Au rang critique triangulaire,

```math
q_\triangle=0.347296355333860\ldots,
```

et, pour $`p=0.8`$,

```math
\beta_c(0.8)=0.410716539196\ldots.
```

On obtient

```math
\begin{aligned}
c(q_\triangle)&=0.426022047760\ldots,\\
\frac{d(q_\triangle)}{c(q_\triangle)}
&=0.564864236889\ldots,\\
s_c(0.8)&=0.693582222752\ldots,\\
\kappa_{\rm conn}(0.8,q_\triangle)
&=0.886752566857\ldots,\\
\kappa_{\rm flux}(0.8,q_\triangle)
&=0.791530736866\ldots.
\end{aligned}
\tag{6.1}
```

Le taux et le préfacteur de (4.7) valent alors

```math
\gamma(0.8,q_\triangle)=0.120189290653\ldots,
\qquad
\frac{\kappa_{\rm flux}}{2\kappa_{\rm conn}}
=0.446308680939\ldots.
\tag{6.2}
```

La différence entre $`s_c`$ et $`\kappa_{\rm conn}`$ mesure exactement le
prix des histoires où l'arête directe révèle parfaitement la relation.

| $`h`$ | connecté à $`q_\triangle`$ : $`\kappa_{\rm conn}^h`$ | LCA à $`q_\triangle`$ : $`A_h^{\rm LCA}`$ | $`P_h^{\rm conf}`$ sous LCA |
|---:|---:|---:|---:|
| 2 | $`0.786330114827`$ | $`0.701891912662`$ | $`0.850945956331`$ |
| 3 | $`0.697280247720`$ | $`0.622404455209`$ | $`0.811202227605`$ |
| 5 | $`0.548292457256`$ | $`0.489415366734`$ | $`0.744707683367`$ |
| 10 | $`0.300624618684`$ | $`0.268342754045`$ | $`0.634171377023`$ |
| 20 | $`0.090375161359`$ | $`0.080670438112`$ | $`0.540335219056`$ |
| 40 | $`0.008167669791`$ | $`0.007290603861`$ | $`0.503645301931`$ |

Le script encadre $`q_\triangle`$ par deux fractions décimales, vérifie le
changement de signe de $`q^3-3q+1`$, puis propage exactement ces fractions
dans (3.1) et (4.3). Les signes et les intervalles annoncés ne reposent donc
pas sur un arrondi flottant.

## 7. Deux calculs indépendants

Le fichier
[`cactus_collapsed_certificate.py`](computations/cactus_collapsed_certificate.py)
contient trois voies de calcul.

1. **Énumération globale.** Pour un cactus latéral-premier fixé de $`h`$
   triangles, il somme les $`2^{2h}`$ configurations de spins après fixation
   du flip global et les $`3^h`$ configurations de marques résiduelles. Le
   poids utilisé est le produit exact des facteurs
   $`\Lambda e^{(1-\beta)\Lambda}`$.
2. **Transfert local.** Il multiplie les coefficients
   $`s_p(r_1)\cdots s_p(r_h)`$ des buckets de taille deux.
3. **Intégration des formes.** Il recalcule (3.1) par quadrature de la densité
   $`4r(1-r)`$, indépendamment de la forme simplifiée.

Pour trois blocs latéraux-premiers aux rangs
$`(q_\triangle-0.002,q_\triangle-0.001,q_\triangle)`$, les deux premiers
calculs donnent

```text
direct=0.334328185717
transfer=0.334328185717
```

Les tests vérifient en outre que les niveaux des buckets de taille un
disparaissent bien du transfert de la paire, comme l'impose la preuve.

## 8. État de bord et contre-audit final

Sur le cactus libre, le seul état de transfert nécessaire est la parité
répliquée $`\delta\in\{-1,+1\}`$. Si l'événement d'effacement reçoit un
message extérieur $`B`$ vérifiant $`|B|\le b`$, le coefficient local devient

```math
\kappa(b)
=
\kappa(0)
+[1-\kappa(0)]\tanh^2(b/2)
<1
\tag{8.1}
```

pour tout $`b<\infty`$. Cette borne ne peut être multipliée sur la grille
sans contrôler la dépendance des messages de bord.

Les six limites exactes du résultat sont les suivantes.

1. Les triangles du cactus ne partagent aucune arête et leur graphe de blocs
   est un arbre.
2. Le conditionnement LCA utilise une densité de Palm ; l'événement
   $`q_{\rm LCA}=q`$ a probabilité nulle avant désintégration.
3. Le dendrogramme est non marqué. Conserver la gagnante détruirait le canal
   d'effacement du lemme 2.1.
4. Sur la grille triangulaire, une coupe de Kruskal peut agréger plusieurs
   cycles et le message de bord n'est plus un bit.
5. Le rang $`q_\triangle`$ est critique pour la grille, pas pour le cactus ;
   la connexion des extrémités a la masse rare $`c(q_\triangle)^h`$.
6. Le calcul porte sur le heat bath collapsed pair-spécifique. Il fournit le
   couplage suffisant du théorème 2.2 du fichier 20, mais ne donne pas la loi
   après un unique sweep top-down ou bottom-up.

## 9. Conséquence pour l'arbre de recherche

Le cactus valide les trois mécanismes recherchés :

- le cas critique est effectivement le cas postcritique le plus favorable ;
- même sous ce conditionnement favorable, la conformité Nishimori tend vers
  $`1/2`$ ;
- le transfert répliqué et son état de bord se calculent exactement.

Le prochain modèle doit être une bande triangulaire de largeur deux. Son
état de bord doit conserver la partition des sommets de coupe et les deux
parités répliquées. L'objectif n'est plus de deviner un coefficient local,
mais de construire une matrice finie $`\mathscr U_{p,2}`$ et de certifier

```math
r(\mathscr U_{0.8,2})<1
```

ainsi que la dégradation critique/postcritique. C'est le premier endroit où
les cycles se chevauchent et où le mécanisme du cactus peut réellement
échouer. Une bande de largeur fixée reste toutefois quasi unidimensionnelle :
elle teste le transfert de bord, pas encore la loi de la géante critique
bidimensionnelle.

## 10. Audit synthétique

| affirmation | statut | limite |
|---|---|---|
| Le bucket latéral de taille deux est un effacement | Établi exactement | dendrogramme non marqué |
| La connexion avant $`q`$ se factorise sur la chaîne | Établi | propriété d'articulation |
| Les extrémités sont typiquement dans une géante à $`q_\triangle`$ | Faux sur le cactus | masse conditionnante $`c(q_\triangle)^h`$ |
| Le LCA global à $`q`$ donne un bloc pivotal et $`h-1`$ blocs connectés | Établi par les densités | égalités de rang de probabilité nulle |
| Le seuil critique est le cas favorable parmi les rangs postcritiques | Établi sur le cactus | dérivées (5.3)--(5.4) |
| À $`p=0.8`$, la conformité Nishimori tend vers $`1/2`$ | Établi sur le cactus | persistance du couplage, pas succès quenched d'un estimateur |
| Le même équivalent vaut après un sweep séquentiel | Non démontré | le collapsed est une projection différente |
| Le coefficient d'un triangle est seulement $`s_c`$ | Faux | les fusions directes-premières sont parfaites |
| Ce calcul prouve $`p=0.8`$ sur la grille triangulaire | Faux | état de bord et cycles chevauchants ouverts |

## 11. LCA seul contre corridor complet

Le [fichier 22](22_LCA_VS_FULL_HIERARCHY.md) répond à la question de
profondeur laissée ouverte ici. Sous le même conditionnement Palm,

```math
A_h^{\rm LCA\ only}(p,q)=\kappa_{\rm flux}(p,q),
\qquad
\frac{A_h^{\rm full}(p,q)}{A_h^{\rm LCA\ only}(p,q)}
=
\kappa_{\rm conn}(p,q)^{h-1}.
```

Le bloc pivotal critique est donc un point de départ de bonne qualité, mais
il ne voit pas la distance. La perte asymptotique utilise nécessairement les
deux bras descendants jusqu'aux feuilles.
