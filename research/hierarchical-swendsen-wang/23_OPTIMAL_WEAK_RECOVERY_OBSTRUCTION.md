# Stratégie optimale pour une obstruction de weak recovery

> [!NOTE]
> La version synthétique et l'ordre de travail courant sont dans le
> [programme prioritaire](00_RESEARCH_PROGRAM.md). Le
> [contre-audit multiport](29_AUDIT_FROID_PIVOT_RANGS_REELS.md) invalide
> l'étape historique de criticalisation ci-dessous. L'architecture annulaire
> reste pertinente, mais le [pivot $`L^2`$](30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md)
> remplace désormais le transfert local borné par des projections collapsed
> et une occupation pondérée par l'énergie.

Cette note répond à la question suivante : parmi toutes les façons d'utiliser
la dynamique hiérarchique, quelle architecture de preuve a le plus de chances
de transformer la fusion critique d'une paire lointaine en une impossibilité
de weak recovery ?

La réponse courte est :

1. utiliser le **corridor descendant collapsed** propre à la paire, et non le
   seul LCA ;
2. conserver chaque update postcritique à son rang réalisé et intégrer des
   blocs collapsed imbriqués ;
3. décomposer le corridor en blocs annulaires screenés ;
4. minorer la variance perdue sur la fonction effectivement propagée, pas une
   norme uniforme qui contient les constantes ;
5. montrer sous la loi Palm que l'énergie portée par les bons blocs est
   abondante.

L'ordre de Blackwell reste utile pour un bucket mono-bit, mais ne fournit pas
d'enveloppe du corridor collapsed. La comparaison avec un corridor critique
est désormais un benchmark, pas une réduction.

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
\qquad\text{(1.1)}
```

Le théorème 2.2 du fichier 20 donne déjà

```math
Q_L
\le
\mathbb E[A_{I_LJ_L}(O,D)],
\qquad\text{(1.2)}
```

où $`Q_L`$ est le second moment moyen des corrélations postérieures. Il suffit
donc de montrer

```math
\mathbb E[A_{I_LJ_L}]\longrightarrow0.
\qquad\text{(1.3)}
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
\qquad\text{(2.1)}
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
\qquad\text{(2.2)}
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

### Proposition 3.1 — surrogate mono-bit, statut : établi sous factorisation

Fixons une expérience dans laquelle chaque bucket observe un unique bit
latent et les observations sont indépendantes conditionnellement au vecteur
de bits. Remplaçons chaque niveau $`t_r>\beta_c`$ par $`\beta_c`$ et laissons
les niveaux plus précoces inchangés :

```math
t_r^{\rm fav}=\min(t_r,\beta_c).
\qquad\text{(3.1)}
```

Alors, pour toute loi corrélée des parités latentes et toute cible $`F`$,

```math
\mathscr R((m_r,t_r)_r)
\le
\mathscr R((m_r,t_r^{\rm fav})_r).
\qquad\text{(3.2)}
```

#### Justification

À taille fixée, le canal mono-bit du compte au niveau le plus précoce
Blackwell-domine celui du niveau tardif. Les noyaux de dégradation se
tensorisent conditionnellement au vecteur complet des parités, même si ce
vecteur a un prior corrélé. C'est le théorème 4.1 du fichier 20.

Cette proposition ne s'applique pas au corridor multiport réel. Deux groupes
d'incidence variant séparément suffisent à la réfuter, sans changer la taille
ni le squelette ; voir le fichier 29.

### Conjecture 3.2 — enveloppe géométrique critique, statut : ouvert

Notons $`\overline A_L(q)`$ l'espérance de Palm de (1.1) par rapport à la
mesure de flux des LCA de rang $`q`$. La formulation forte de l'idée serait

```math
\sup_{q\ge q_c}\overline A_L(q)
\le
\sup_{|r-q_c|\le\delta_L}\overline A_L(r)+o(1)
\qquad(\delta_L\downarrow0).
\qquad\text{(3.3)}
```

Elle est exacte sur le cactus triangulaire du fichier 21. Elle n'est pas
démontrée sur la grille : changer $`q`$ change les composantes de Kruskal,
les coupes et les tailles $`m_r`$. Or deux expériences de tailles différentes
peuvent être incomparables au sens de Blackwell.

### Choix méthodologique

La meilleure preuve ne doit pas dépendre entièrement de (3.3). Elle doit :

- utiliser la proposition 3.1 seulement à l'intérieur d'une vraie coordonnée
  mono-bit, jamais comme remplacement du corridor multiport ;
- construire le mécanisme de contraction directement aux rangs réalisés,
  sous la loi de bord annealed ;
- conserver (3.3) comme raccourci possible si un couplage préservant les
  tailles est découvert.

Ainsi, le cas critique reste un benchmark favorable sur le cactus et pour un
canal mono-bit. Il n'est pas un oracle supérieur pour le canal multiport
réel. La conjecture géométrique (3.3) est un éventuel raccourci séparé, pas
un axiome du programme actif.

## 4. Le moteur multiscale recommandé

Prenons des annuli séparés

```math
\mathcal A_k
=
B(i,\rho^{k+1})\setminus B(i,\rho^k),
\qquad 1\le k\le K_L\asymp\log L,
\qquad\text{(4.1)}
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
\qquad\text{(4.2)}
```

où $`\Pi_k`$ est la partition de connectivité des ports, les deux
$`R_k^{(a)}`$ sont les parités dans les deux répliques et $`P_k`$ indique le
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
\qquad\text{(4.3)}
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
\qquad\text{(4.4)}
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
\qquad\text{(4.5)}
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
\qquad\text{(4.6)}
```

Uniformément dans le rang postcritique pertinent, (4.6) implique une
abondance linéaire en $`K_L`$ des bons blocs. En combinant (4.5)--(4.6), on
obtient soit

```math
\mathbb E[A_{ij}\mid\mathrm{Palm}]
\le
\bigl(1-a(p)(1-\kappa(p))\bigr)^{K_L}+o(1),
\qquad\text{(4.7)}
```

soit la borne équivalente

```math
\mathbb E[A_{ij}\mid\mathrm{Palm}]
\le
\mathbb P(N_L<c\log L)+\kappa(p)^{c\log L}+o(1).
\qquad\text{(4.8)}
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

Fixons $`p>1/2`$. Supposons :

1. la localisation sous-critique de la section 2 ;
2. l'annulation exacte pour les racines distinctes ;
3. une domination cible-spécifique démontrée sous la vraie loi de bord, ou
   une borne directe de Feynman--Kac sur les transferts aux rangs réalisés ;
4. un découpage de transfert satisfaisant (4.5) ;
5. le lemme géométrique (4.6), soit sous toutes les lois Palm
   postcritiques, soit au seuil complété par (3.3).

Alors

```math
\mathbb E[A_{I_LJ_L}]\longrightarrow0,
```

et la weak recovery est impossible à ce $`p`$.

#### Preuve

Les classes précoce et racines distinctes disparaissent par 1--2.
L'hypothèse 3 traite directement les canaux postcritiques réels. Les
hypothèses 4--5 et (4.7) font tendre vers zéro leur contribution. L'équation
(1.2) donne $`Q_L\to0`$, ce qui interdit la weak recovery.

Si ce théorème est obtenu à $`p_0=4/5`$, la dégradation BSC des observations
étend l'impossibilité à tout $`p\le p_0`$.

## 7. Jalon quantitatif à $`p=0.8`$

Au seuil triangulaire,

```math
s_c(0.8)=0.693582222752\ldots.
\qquad\text{(7.1)}
```

Le cactus fournit déjà un bloc physique exact de coefficient

```math
\kappa_{\rm conn}(0.8,q_c)
=0.886752566857\ldots<1.
\qquad\text{(7.2)}
```

Il démontre que la marge locale est confortable dans cette cellule. Le jalon
suivant n'est ni un nouveau calcul scalaire, ni directement une bande de
largeur deux : il faut d'abord construire un quotient de ports Markov-fermé
où une orientation portant le twist est éliminée après sa dernière
interaction. Conserver le micro-état complet donne le no-go $`|U|=K`$ du
fichier 29.

La bande reste un certificat de **canal**. Elle ne possède pas la géométrie
critique bidimensionnelle et ne prouve pas l'abondance annulaire.

## 8. Audit et contre-audit

| affirmation | verdict | conséquence |
|---|---|---|
| Le LCA seul est la meilleure dynamique d'obstruction | Faux | il est plus persistant ; le corridor collapsed est plus contractant |
| Les états $`(0,0)`$ et $`(1,1)`$ décorrèlent $`i,j`$ | Faux | ils conservent exactement $`\sigma_i\sigma_j`$ |
| Mettre un canal tardif mono-bit au niveau critique aide la récupération | Établi pour le bucket scalaire ; faux universellement en multiport | le corridor réel exige un transfert direct |
| Une vraie géométrie critique domine toute géométrie tardive | Ouvert | les tailles et incidences changent ; contre-exemples cross-size |
| Une paire dans la composante critique est un événement typique | Faux | au seuil, il n'y a pas de géante de densité positive |
| On peut ignorer le complément du conditionnement critique | Faux | traiter directement la Palm d'événement aux rangs réalisés |
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

1. **Test de dernière utilisation.** Mesurer sous la Palm d'événement quand
   une orientation interne cesse d'affecter les buckets ancestraux.
2. **Certificat fini.** Construire une jauge de ports Markov-fermée ou un
   bloc multi-update qui élimine cette orientation ; le micro-état complet
   doit rendre automatiquement un déficit nul.
3. **Jointure réelle.** Sur les mêmes cellules, joindre rang, message,
   signature de ports et déficit composable à $`p=0.805`$.
4. **Lemme annulaire marqué.** Traduire le motif en outlet T2 protégé et
   contrôler la transformée de Laplace du nombre de visites avec potentiel
   modéré, sans supposer les annuli indépendants.
5. **Clôture pairwise.** Insérer (4.7) dans (2.2), puis dans (1.2), et
   appliquer la dégradation BSC.

Cette organisation sépare trois difficultés qui doivent rester distinctes :
la qualité des marques, le transfert hiérarchique avec tous les
$`\Lambda_v`$, et la géométrie Palm du corridor.

## Conclusion

La meilleure stratégie n'est ni de s'arrêter au nœud de fusion, ni de suivre
un unique chemin gagnant de la MSF. Il faut utiliser toute la hiérarchie
descendante de la paire, mais la regrouper en blocs planaires dont le transfert
répliqué est calculable. La fusion au seuil sert seulement de benchmark sur
le cactus et dans une coordonnée mono-bit. Le vrai théorème manquant est
qu'une paire lointaine sous la Palm d'événement réelle traverse assez
d'outlets T2 protégés où le quotient de frontière revient dans une zone
modérée et accumule un déficit composable.
