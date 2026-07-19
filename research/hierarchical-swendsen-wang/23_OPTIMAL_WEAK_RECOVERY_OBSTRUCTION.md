# Stratégie optimale pour une obstruction de weak recovery

Cette note répond à la question suivante : parmi toutes les façons d'utiliser
la dynamique hiérarchique, quelle architecture de preuve a le plus de chances
de transformer la fusion critique d'une paire lointaine en une impossibilité
de weak recovery ?

La réponse courte est :

1. utiliser le **corridor descendant collapsed** propre à la paire, et non le
   seul LCA ;
2. rendre chaque canal postcritique artificiellement aussi informatif qu'au
   seuil, **sans changer le squelette ni les tailles de buckets** ;
3. décomposer le corridor en blocs annulaires screenés ;
4. prouver une contraction $`\chi^2`$ uniforme du transfert répliqué sur
   chaque bon bloc ;
5. montrer sous la loi Palm de fusion que le nombre de bons blocs diverge.

Cette stratégie conserve exactement l'idée du cas le plus favorable, mais
elle la formule sous une forme que l'ordre de Blackwell permet réellement de
prouver. La comparaison plus forte avec un véritable corridor critique reste
un lemme géométrique ouvert.

## 1. Quantité maître

Sur le tore triangulaire $`G_L`$, soient $`I_L,J_L`$ uniformes indépendants et

```math
f_{ij}(\sigma)=\sigma_i\sigma_j.
```

Conditionnellement aux observations $`O`$ et au dendrogramme non marqué
$`D`$, notons $`P_{ij}^{\downarrow}`$ le heat bath conjoint de toutes les
orientations sur les deux bras qui vont de $`i,j`$ à leur LCA. Les facteurs
ancestraux $`\Lambda_v`$ restent tous présents dans la loi conditionnelle.
Posons

```math
A_{ij}(O,D)
=
\|P_{ij}^{\downarrow}f_{ij}\|_{L^2(\pi_D)}^2,
\qquad 0\le A_{ij}\le1.
\tag{1.1}
```

Le théorème 2.2 du fichier 20 donne déjà

```math
Q_L
\le
\mathbb E[A_{I_LJ_L}(O,D)],
\tag{1.2}
```

où $`Q_L`$ est le second moment moyen des corrélations postérieures. Il suffit
donc de montrer

```math
\mathbb E[A_{I_LJ_L}]\longrightarrow0.
\tag{1.3}
```

Le caractère pair-spécifique du corridor n'est pas une difficulté : (1.2)
est obtenu par Jensen séparément pour chaque paire. En revanche, étendre le
bloc à tout le dendrogramme ne serait plus une simplification ; cela
reformulerait directement le rééchantillonnage de la postérieure complète.

## 2. Décomposition par rang de fusion

Il est préférable de travailler dans la coordonnée de percolation

```math
q_p(t)=p(1-e^{-u_pt}),
\qquad q_c=q_\triangle=2\sin(\pi/18),
\tag{2.1}
```

et de noter $`T_{ij}=q_p(\beta_{ij})`$ le rang du LCA. Pour une paire à
distance macroscopique et tout $`\delta>0`$ fixé,

```math
\begin{aligned}
\mathbb E[A_{I_LJ_L}]
={}&
\mathbb E[A;T<q_c-\delta]
+\mathbb E[A;|T-q_c|\le\delta]\\
&+\mathbb E[A;T>q_c+\delta,\ T\le q_p(1)]
+\mathbb E[A;T>q_p(1)].
\tag{2.2}
\end{aligned}
```

Les termes extrêmes sont déjà contrôlés :

- le premier est $`o_L(1)`$ par décroissance sous-critique, puisque
  $`A\le1`$ ;
- le dernier est exactement nul, car deux racines finales sont recolorées
  indépendamment.

Il reste la fenêtre critique et les fusions postcritiques. L'ordre des
limites doit être $`L\to\infty`$ à $`\delta`$ fixé, puis
$`\delta\downarrow0`$. L'événement $`T_{ij}=q_c`$ a probabilité nulle en
volume fini ; toute égalité au seuil signifie une densité de Palm du flux de
fusions ou une limite de fenêtres.

### Audit du mot « géante »

À $`q_c`$, la percolation bidimensionnelle infinie ne possède pas de
composante infinie de densité positive. Sur un tore fini, « la composante
géante au seuil » doit donc désigner une composante critique macroscopique ou
traversante. L'événement que deux sommets uniformes appartiennent à une même
telle composante est rare asymptotiquement. On ne peut pas simplement
conditionner dessus et oublier son complément : il faut montrer que cette
expérience rare est une **enveloppe favorable** des autres paires.

## 3. Ce qui est vraiment favorable au seuil

Il faut séparer la qualité des marques de la géométrie du corridor.

### Proposition 3.1 — criticalisation à squelette fixé, statut : établi

Fixons un corridor, ses incidences, son état de bord et ses tailles
$`m_r`$. Remplaçons chaque niveau $`t_r>\beta_c`$ par $`\beta_c`$ et laissons
les niveaux plus précoces inchangés :

```math
t_r^{\rm fav}=\min(t_r,\beta_c).
\tag{3.1}
```

Alors, pour toute loi corrélée des parités latentes et toute cible $`F`$,

```math
\mathscr R((m_r,t_r)_r)
\le
\mathscr R((m_r,t_r^{\rm fav})_r).
\tag{3.2}
```

#### Justification

À taille fixée, le canal du compte du bucket au niveau le plus précoce
Blackwell-domine celui du niveau tardif. Les noyaux de dégradation se
tensorisent conditionnellement au vecteur complet des parités, même si ce
vecteur a un prior corrélé. C'est le théorème 4.1 du fichier 20.

Cette proposition est exactement la version rigoureuse de « prendre les
liens les plus nombreux et de meilleure qualité ». Elle ne modifie aucune
taille et ne confond pas les arêtes internes avec les arêtes de la coupe.

### Conjecture 3.2 — enveloppe géométrique critique, statut : ouvert

Notons $`\overline A_L(q)`$ l'espérance de Palm de (1.1) par rapport à la
mesure de flux des LCA de rang $`q`$. La formulation forte de l'idée serait

```math
\sup_{q\ge q_c}\overline A_L(q)
\le
\sup_{|r-q_c|\le\delta_L}\overline A_L(r)+o(1)
\qquad(\delta_L\downarrow0).
\tag{3.3}
```

Elle est exacte sur le cactus triangulaire du fichier 21. Elle n'est pas
démontrée sur la grille : changer $`q`$ change les composantes de Kruskal,
les coupes et les tailles $`m_r`$. Or deux expériences de tailles différentes
peuvent être incomparables au sens de Blackwell.

### Choix méthodologique

La meilleure preuve ne doit pas dépendre entièrement de (3.3). Elle doit :

- utiliser immédiatement la proposition 3.1 sur chaque corridor réel ;
- prouver le mécanisme de contraction uniformément pour les lois de
  corridors de rang $`q\ge q_c`$ ;
- conserver (3.3) comme raccourci possible si un couplage préservant les
  tailles est découvert.

Ainsi, le cas critique est l'oracle le plus favorable pour les **canaux** ;
sa domination de toutes les **géométries** reste à établir, et n'est plus un
axiome caché de la preuve.

## 4. Le moteur multiscale recommandé

Prenons des annuli séparés

```math
\mathcal A_k
=
B(i,\rho^{k+1})\setminus B(i,\rho^k),
\qquad 1\le k\le K_L\asymp\log L,
\tag{4.1}
```

et, symétriquement, des annuli autour de $`j`$. Un bloc peut aussi être
centré sur une portion intrinsèque du corridor plutôt que sur un endpoint.
Le choix final doit garantir que les blocs utilisent des ensembles d'arêtes
disjoints.

### 4.1 État de bord exact

Un bloc ne peut pas être résumé par une seule majorité. Son état minimal doit
contenir au moins :

```math
Z_k=(\Pi_k,R_k^{(1)},R_k^{(2)},P_k),
\tag{4.2}
```

où $`\Pi_k`$ est la partition de connectivité des ports, les deux
$`R_k^{(a)}` sont les parités dans les deux répliques et $`P_k`$ indique le
statut de fusion/pivotalité. Si le nombre de ports n'est pas borné, l'espace
d'états n'est pas fini. Un bon bloc doit précisément fournir un événement de
screening qui le borne.

### 4.2 Bon bloc annulaire

Une définition utilisable doit imposer simultanément :

1. **ports bornés** : au plus $`b`$ interfaces macroscopiques traversent le
   bloc ;
2. **screening latéral** : des séparateurs duaux empêchent le bord extérieur
   de contourner le canal retenu ;
3. **ambiguïté hiérarchique** : au moins deux routes ou arêtes candidates
   participent à une coupe de Kruskal pertinente ;
4. **contraction uniforme** : après sommation exacte des configurations
   internes, l'opérateur répliqué sur le secteur impair vérifie

```math
\|\mathscr U_{p,k}g\|_2^2
\le
\kappa(p)\|g\|_2^2,
\qquad \kappa(p)<1,
\tag{4.3}
```

   pour chaque état de bord admis.

Le screening peut être formulé par une borne $`|B_k|\le B_0`$ sur le
log-rapport de vraisemblance extérieur. Pour un bucket critique $`m=2`$, le
calcul local donne alors

```math
\kappa_2(B_0;p)
=
s_c(p)+(1-s_c(p))\tanh^2(B_0/2)<1,
\qquad
s_c(p)=\frac{p-q_c}{1-q_c}.
\tag{4.4}
```

La formule (4.4) est une calibration, pas un remplacement du transfert
complet : les routes parallèles et la partition de bord doivent être
intégrées dans $`\mathscr U_{p,k}`$.

### 4.3 Lemme de canal à viser

Soit $`G_k`$ l'événement que le bloc $`k`$ est bon et
$`N_L=\sum_{k\le K_L}\mathbf1_{G_k}`$. Après découpage exact du corridor,
montrer

```math
A_{ij}
\le
\kappa(p)^{N_L}+\varepsilon_L^{\rm tr},
\qquad
\mathbb E[\varepsilon_L^{\rm tr}]=o(1).
\tag{4.5}
```

Cette inégalité doit être conditionnelle aux géométries et états de bord
retenus. Elle ne requiert alors aucune indépendance entre les blocs : la
composition suit de la sous-multiplicativité des normes d'opérateurs.

Une variante plus réaliste remplace $`\kappa^{N_L}`$ par le rayon spectral
du produit de matrices positives normalisées. Le certificat de bande de
largeur deux doit déterminer laquelle des deux formulations est exacte.

### 4.4 Lemme géométrique à viser

Sous la loi Palm de connexion/fusion d'une paire lointaine, trouver des
annuli espacés tels que, pour une filtration révélant successivement les
échelles,

```math
\mathbb P(G_k\mid\mathcal F_{k-1},\mathrm{Palm})
\ge a(p)>0.
\tag{4.6}
```

Uniformément dans le rang postcritique pertinent, (4.6) implique une
abondance linéaire en $`K_L`$ des bons blocs. En combinant (4.5)--(4.6), on
obtient soit

```math
\mathbb E[A_{ij}\mid\mathrm{Palm}]
\le
\bigl(1-a(p)(1-\kappa(p))\bigr)^{K_L}+o(1),
\tag{4.7}
```

soit la borne équivalente

```math
\mathbb E[A_{ij}\mid\mathrm{Palm}]
\le
\mathbb P(N_L<c\log L)+\kappa(p)^{c\log L}+o(1).
\tag{4.8}
```

Dans les deux cas, le membre de droite tend polynomialement vers zéro.

## 5. Pourquoi la géométrie planaire peut suffire

L'objectif (4.6) ne nécessite pas, en première intention, les exposants
critiques exacts. Les outils robustes à viser sont :

- RSW pour obtenir des traversées et circuits à probabilité non dégénérée ;
- quasi-multiplicativité et séparation des interfaces pour conserver ces
  bornes sous un conditionnement de connexion lointaine ;
- sprinkling proche-critique pour comparer des fenêtres de rang ;
- finite energy pour insérer un motif ambigu local une fois les interfaces
  séparées.

Le conditionnement est toutefois une loi Palm **à deux points** portant sur
tout le backbone. RSW non conditionné ne donne pas automatiquement (4.6).
Il faut un vrai lemme de rapport de probabilités, uniforme par rapport aux
configurations intérieure et extérieure compatibles.

Les travaux sur les mesures pivotales suggèrent la bonne filtration
multiscale, mais un pivot isolé ne fournit pas une contraction : il correspond
précisément à un bucket $`m=1`$, donc à un canal parfait. Les bons objets sont
des blobs ambigus autour des goulots, pas les goulots seuls.

## 6. Théorème conditionnel de clôture

### Théorème 6.1 — critère annulaire, statut : conditionnel

Fixons $`p>1/2`. Supposons :

1. la localisation sous-critique de la section 2 ;
2. l'annulation exacte pour les racines distinctes ;
3. la criticalisation à squelette fixé de la proposition 3.1 ;
4. un découpage de transfert satisfaisant (4.5) ;
5. le lemme géométrique (4.6), soit sous toutes les lois Palm
   postcritiques, soit au seuil complété par (3.3).

Alors

```math
\mathbb E[A_{I_LJ_L}]\longrightarrow0,
```

et la weak recovery est impossible à ce $`p`$.

#### Preuve

Les classes précoce et racines distinctes disparaissent par 1--2. La
proposition 3.1 majore les canaux tardifs sur leur squelette réel par leurs
versions criticalisées. Les hypothèses 4--5 et (4.7) font tendre vers zéro
la contribution critique/postcritique. L'équation (1.2) donne $`Q_L\to0`$,
ce qui interdit la weak recovery.

Si ce théorème est obtenu à $`p_0=4/5`$, la dégradation BSC des observations
étend l'impossibilité à tout $`p\le p_0`$.

## 7. Jalon quantitatif à $`p=0.8`$

Au seuil triangulaire,

```math
s_c(0.8)=0.693582222752\ldots.
\tag{7.1}
```

Le cactus fournit déjà un bloc physique exact de coefficient

```math
\kappa_{\rm conn}(0.8,q_c)
=0.886752566857\ldots<1.
\tag{7.2}
```

Il démontre que la marge locale est confortable. Le jalon suivant n'est pas
un nouveau calcul scalaire : il faut certifier, sur une bande triangulaire de
largeur deux, un opérateur avec partition de bord et deux répliques dont le
rayon spectral est strictement inférieur à un. Ce calcul doit aussi extraire
une liste finie de configurations de ports qui puissent devenir les bons
blocs de (4.6).

La bande reste un certificat de **canal**. Elle ne possède pas la géométrie
critique bidimensionnelle et ne prouve pas l'abondance annulaire.

## 8. Audit et contre-audit

| affirmation | verdict | conséquence |
|---|---|---|
| Le LCA seul est la meilleure dynamique d'obstruction | Faux | il est plus persistant ; le corridor collapsed est plus contractant |
| Les états $`(0,0)`$ et $`(1,1)`$ décorrèlent $`i,j`$ | Faux | ils conservent exactement $`\sigma_i\sigma_j`$ |
| Mettre un canal tardif au niveau critique aide la récupération | Établi à squelette et taille fixés | domination de Blackwell, donc oracle favorable valide |
| Une vraie géométrie critique domine toute géométrie tardive | Ouvert | les tailles et incidences changent ; contre-exemples cross-size |
| Une paire dans la composante critique est un événement typique | Faux | au seuil, il n'y a pas de géante de densité positive |
| On peut ignorer le complément du conditionnement critique | Faux | il faut (3.3) ou un lemme uniforme postcritique |
| Un pivot fournit un bloc contractant | Faux | un bucket $`m=1`$ transmet parfaitement |
| Une majorité locale stricte suffit | Faux | le message de bord peut l'écraser |
| Un coefficient local $`<1`$ suffit | Faux | il faut une abondance divergente et contrôler les contournements |
| RSW ordinaire donne directement (4.6) | Faux | le biais Palm à deux points exige quasi-multiplicativité et séparation |
| Les exposants exacts de la percolation de sites triangulaire s'appliquent automatiquement ici | Faux | le modèle géométrique présent est une percolation par arêtes ; éviter tout transfert d'universalité non démontré |
| Un sweep top-down est équivalent au collapsed | Faux | les projections ne commutent pas ; seule leur limite cyclique donne le bloc conjoint |
| La preuve à $`p=0.8`$ est terminée | Faux | les lemmes (4.5), (4.6) et la porte postcritique restent ouverts |

### Contre-audit structurel

Une preuve qui donnerait (4.6) avec screening uniforme pour tout $`p<1`$
conduirait abusivement à une impossibilité jusqu'à $`p=1`$. Le point où une
vraie transition peut apparaître est précisément l'échec du screening : les
messages latéraux ou les routes parallèles deviennent capables de conserver
une information macroscopique. Le seuil recherché doit donc sortir du
spectre du transfert complet, pas de la seule inégalité locale
$`\kappa_2(B_0;p)<1`$.

## 9. Ordre de travail recommandé

1. **Test simple.** Chercher d'abord sous Palm un nombre divergent de buckets
   bornés, en priorité $`m=2`$, dont le message ancestral et les
   contournements sont screenés. Le fichier 24 donne le théorème conditionnel
   exact correspondant.
2. **Certificat fini si nécessaire.** Si le test simple ne permet pas de
   fermer l'état de bord, construire le transfert exact de largeur deux,
   doublement répliqué, et identifier le plus petit motif de ports réellement
   contractant à $`p=0.8`$.
3. **Lemme annulaire critique.** Traduire ce motif en événement planaire et
   prouver (4.6) sous la loi Palm à deux points au seuil.
4. **Porte postcritique.** D'abord tenter la version uniforme de (4.6) après
   criticalisation à squelette fixé ; ne chercher (3.3) que si cette voie
   échoue.
5. **Clôture pairwise.** Insérer (4.7) dans (2.2), puis dans (1.2), et
   appliquer la dégradation BSC.

Cette organisation sépare trois difficultés qui doivent rester distinctes :
la qualité des marques, le transfert hiérarchique avec tous les
$`\Lambda_v`$, et la géométrie Palm du corridor.

## Conclusion

La meilleure stratégie n'est ni de s'arrêter au nœud de fusion, ni de suivre
un unique chemin gagnant de la MSF. Il faut utiliser toute la hiérarchie
descendante de la paire, mais la regrouper en blocs planaires dont le transfert
répliqué est calculable. La fusion au seuil sert d'**oracle favorable de
canal**. Le vrai théorème manquant est qu'une paire lointaine, même sous ce
biais favorable, traverse logarithmiquement beaucoup de blobs hiérarchiques
screenés et strictement contractants.
