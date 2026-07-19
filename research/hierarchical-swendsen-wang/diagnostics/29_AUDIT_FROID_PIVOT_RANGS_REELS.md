# Audit à froid et pivot vers les rangs réels

> [!IMPORTANT]
> La feuille de route fondée sur la **criticalisation globale du corridor
> collapsed est réfutée**. La dynamique hiérarchique elle-même reste une voie
> possible, mais le transfert Feynman--Kac local à séparateur borné est lui
> aussi une impasse probable. Cette note motive le socle de
> [dissipation quadratique](../active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) ;
> son assemblage le plus récent est le
> [programme distance–entropie](../active/35_DISTANCE_ENTROPIE_ERGODICITE.md).
> Aucun résultat de cette note ne donne une nouvelle borne de weak recovery.

Cette note répond à la question stratégique : sommes-nous dans une bonne
voie, ou dans une impasse ? Le verdict est double.

- La voie écrite dans les fichiers 20, 23, 26 et 27 est une impasse logique :
  son oracle favorable principal est faux sur le corridor réel.
- Le programme hiérarchique réparé n'est pas réfuté. Sur les tailles finies
  testées, les diagnostics indiquent deux ressources quantitatives : des
  petites attaches en peigne et des nœuds dont le message ancestral réel
  reste modéré.
- Un premier transfert fidèle révèle toutefois un second no-go : conserver
  dans l'état toute l'orientation répliquée rend le twist mesurable et force
  un déficit local nul. La prochaine inconnue n'est donc plus seulement la
  taille des ports, mais leur élimination dynamique exacte.
- Il faut arrêter de chercher une domination uniforme. La quantité active est
  la perte relative $`L^2`$ du caractère de paire le long d'une filtration
  collapsed ; le déficit additif local ne subsiste que conditionnellement à
  une compression spéciale encore inconnue.

## 1. Le point de rupture exact

### 1.1 Ce qui reste vrai pour un bucket scalaire

Si toutes les arêtes d'un bucket codent un même bit latent
$`X\in\{-1,+1\}`$, le compte satisfait $`K`$ vérifie

```math
K\mid X=+1
\overset{\mathrm d}=1+\mathrm{Bin}(m-1,s),
\qquad
K\mid X=-1
\overset{\mathrm d}=\mathrm{Bin}(m-1,1-s).
\qquad\text{(1.1)}
```

Pour $`1/2\le s_2\le s_1\le1`$, le canal de paramètre $`s_1`$
Blackwell-domine celui de paramètre $`s_2`$. Le théorème 4.2 du
[fichier 19](../foundations/19_FAVORABLE_SWEEP_PROJECTIONS.md) reste donc valide.

Cette situation correspond à un update local dont les descendants sont
figés : toutes les relations physiques de la coupe se complémentent ensemble
quand on retourne le bit relatif des deux enfants.

### 1.2 Pourquoi le corridor collapsed n'est pas une collection de ces bits

Dans un heat bath conjoint de plusieurs niveaux, les flips descendants font
varier séparément plusieurs groupes d'incidence d'un bucket ancestral. Le
facteur réel est

```math
F_v(\sigma)
=
\Lambda_v(\sigma)
\exp\!\left((1-\beta_v)\Lambda_v(\sigma)\right),
\qquad\text{(1.2)}
```

et $`\Lambda_v`$ dépend alors d'un vecteur de relations, pas d'un unique bit
$`X_v`$. L'identité de l'arête gagnante est marginalisée dans le dendrogramme
non marqué. C'est précisément ce mélange qui détruit l'ordre de Blackwell.

Le théorème 4.1 du [fichier 20](../foundations/20_COLLAPSED_CORRIDOR_BLACKWELL.md) reste un
théorème abstrait correct pour une expérience artificielle factorisée, avec
un bit global par bucket. Son identification au **vrai** corridor collapsed
multiport était incorrecte.

## 2. Contre-exemple minimal certifié

### 2.1 Canal à deux incidences

Prenons $`x=(x_1,x_2)\in\{-1,+1\}^2`$. Un bucket contient deux arêtes,
portant les deux relations $`x_1`$ et $`x_2`$. On tire une gagnante
$`G`$ uniformément dans $`\{1,2\}`$, on impose $`y_G=x_G`$, puis l'autre
marque est correcte avec probabilité $`s`$. La gagnante n'est pas observée.

Le canal est additif sur $`H=\{-1,+1\}^2`$ :

```math
W_s(y\mid x)=P_s(yx),
\qquad
P_s(++ )=s,
\quad
P_s(+-)=P_s(-+)=\frac{1-s}{2},
\quad
P_s(--)=0.
\qquad\text{(2.1)}
```

Ses coefficients de Fourier non triviaux sont

```math
\widehat P_s(\chi_1)=\widehat P_s(\chi_2)=s,
\qquad
\widehat P_s(\chi_{12})=2s-1.
\qquad\text{(2.2)}
```

Supposons qu'un canal critique $`W_{s_c}`$ se dégrade en un canal tardif
$`W_{s_l}`$. En moyennant les translatés de n'importe quel noyau de
dégradation, on obtiendrait un noyau covariant. Il serait la convolution par
une probabilité $`R`$. Comme tous les coefficients de Fourier critiques sont
non nuls, cette probabilité serait unique et vérifierait

```math
\widehat R(\chi_1)=\widehat R(\chi_2)=\frac{s_l}{s_c},
\qquad
\widehat R(\chi_{12})=\frac{2s_l-1}{2s_c-1}.
\qquad\text{(2.3)}
```

L'inversion de Fourier donne

```math
R(--)
=
\frac14
\left(
1-2\frac{s_l}{s_c}
+\frac{2s_l-1}{2s_c-1}
\right)
=
\frac{(1-h_c)(h_l-h_c)}{4h_c(1+h_c)}<0,
\qquad
h_a=2s_a-1.
\qquad\text{(2.4)}
```

Il n'existe donc aucune dégradation de Blackwell dès que
$`1/2<s_l<s_c<1`$.

### 2.2 Inversion pour la cible de parité

Au point

```math
p=\frac{161}{200},
\qquad
q_c=2\sin\!\left(\frac{\pi}{18}\right),
\qquad
q_l=\frac{11}{20},
\qquad\text{(2.5)}
```

les probabilités résiduelles sont

```math
s_c=\frac{p-q_c}{1-q_c}=0.701242667183598\ldots,
\qquad
s_l=\frac{p-q_l}{1-q_l}=\frac{17}{30}.
\qquad\text{(2.6)}
```

Le certificat rationnel et par intervalles obtient

```math
R(--)
\in
[-0.071225876442769,-0.071225876442768].
\qquad\text{(2.7)}
```

Prenons, dans l'ordre $`(++,+-,-+,--)`$, le prior exact

```math
\mu
=
\left(
\frac{19}{20},
\frac1{50},
\frac1{40},
\frac1{200}
\right)
\qquad\text{(2.8)}
```

et la cible $`f(x)=x_1x_2`$. Alors

```math
\mathbb E_\mu f=\frac{91}{100},
\qquad
\mathrm{Var}_\mu(f)=\frac{1719}{10000}.
\qquad\text{(2.9)}
```

La variance expliquée par la sortie vaut

```math
\begin{aligned}
\mathrm{Var}\!\left(\mathbb E[f\mid Y_c]\right)
&=0.022535603548554\ldots,\\
\mathrm{Var}\!\left(\mathbb E[f\mid Y_l]\right)
&=0.028797513006574\ldots.
\end{aligned}
\qquad\text{(2.10)}
```

Le canal tardif conserve donc **davantage** la cible, avec un écart certifié

```math
0.006261909458020\ldots>0.
\qquad\text{(2.11)}
```

Le calcul reproductible est
[`multiport_blackwell_counterexample.py`](../computations/multiport_blackwell_counterexample.py).

## 3. Contre-audit dans une cellule quotient T2-Kruskal

Le contre-exemple abstrait ne contient pas encore les autres facteurs de la
hiérarchie. La cellule quotient T2 suivante vérifie que l'inversion persiste
en présence de deux facteurs de fusion exacts. Elle n'est pas embarquée dans
une géométrie Palm réelle et ne prouve pas que le bord polarisé choisi y
possède une masse positive.

- La fusion cible possède deux arêtes portant
  $`\chi=x_1x_2`$ et survient exactement à $`\beta_c`$.
- Une fusion ancestrale possède deux groupes d'une arête, portant
  respectivement $`x_1`$ et $`x_2`$.
- La gagnante de chaque fusion est marginalisée.
- Les deux facteurs $`\Lambda e^{(1-\beta)\Lambda}`$ sont conservés.
- Le prior extérieur est
  $`\pi_{B,J}(x)\propto e^{B(x_1+x_2)+Jx_1x_2}`$.
- Les deux répliques partagent le même environnement observé.

Deux constructions indépendantes coïncident : énumération arête par arête,
puis élimination analytique de la gagnante par comptes groupés.

À $`p=0.805`$, avec $`\beta_l=0.8`$, la fiabilité normalisée de
$`\chi`$ est :

| bord $`(B,J)`$ | ancêtre critique | ancêtre tardif | tardif moins critique |
|---:|---:|---:|---:|
| $`(0,0)`$ | $`0.749639695`$ | $`0.707168938`$ | $`-0.042470757`$ |
| $`(2,1)`$ | $`0.736267288`$ | $`0.755666309`$ | $`+0.019399022`$ |
| $`(4,3)`$ | $`0.735112203`$ | $`0.755637535`$ | $`+0.020525332`$ |

L'ordre attendu tient au bord neutre et s'inverse au bord polarisé. Il ne
reste donc pas de domination uniforme à sauver par une description plus fine
des ports. Le calcul est
[`kruskal_fusion_t2_transfer.py`](../computations/kruskal_fusion_t2_transfer.py).

## 4. Ce qui survit intégralement

Le contre-exemple n'affecte pas la réduction de la weak recovery à la
corrélation spin--spin répliquée :

```math
Q_L
=
\frac1{|V_L|^2}
\sum_{i,j}
\mathbb E\!\left[
\langle\sigma_i\sigma_j\rangle^2
\right].
\qquad\text{(4.1)}
```

Il n'affecte pas non plus le théorème H pair-spécifique : pour tout heat bath
exact du corridor,

```math
Q_L
\le
\mathbb E\!\left[
\|P_{I_LJ_L}f_{I_LJ_L}\|_2^2
\right].
\qquad\text{(4.2)}
```

Les éléments suivants survivent également :

- l'optimalité $`L^2`$ de la projection collapsed parmi les sweeps des mêmes
  coordonnées ;
- l'annulation exacte lorsque les deux sommets appartiennent à deux racines
  finales distinctes ;
- la loi de Palm des **événements de fusion réalisés**, pondérée par
  $`N_\rho`$ seulement : la taille de coupe a déjà biaisé la course de
  Kruskal et ne doit pas être multipliée une seconde fois ;
- le lemme de Blackwell scalaire du fichier 19, avec sa portée restreinte ;
- la composition positive/tordue de Feynman--Kac.

La dynamique hiérarchique reste donc le bon langage. C'est l'oracle temporel
qui disparaît.

## 5. La composition finie n'a plus besoin d'une normalisation commune

Soient des transferts de masse positifs $`K_r`$ et des transferts tordus
$`U_r`$ aux **rangs réels**, satisfaisant l'inégalité entrée par entrée

```math
|U_r(z,z')|\le K_r(z,z').
\qquad\text{(5.1)}
```

Posons

```math
r_r(z,z')
=
\begin{cases}
|U_r(z,z')|/K_r(z,z'),&K_r(z,z')>0,\\
0,&K_r(z,z')=0,
\end{cases}
\qquad
d_r=
\begin{cases}
-\log r_r,&r_r>0,\\
+\infty,&r_r=0.
\end{cases}
\qquad\text{(5.2)}
```

On utilise la convention $`e^{-d_r}=0`$ lorsque $`d_r=+\infty`$.

Pour un horizon fini $`N`$, choisissons une fonction terminale non négative
$`h_N`$, avec $`|g|\le h_N`$, et définissons rétrogradement

```math
h_{r-1}=K_rh_r.
\qquad\text{(5.3)}
```

Sur le support $`h_{r-1}(z)>0`$, le noyau

```math
P_r(z,z')
=
\frac{K_r(z,z')h_r(z')}{h_{r-1}(z)}
\qquad\text{(5.4)}
```

est stochastique. Si $`|g|\le h_N`$, l'expansion des chemins et la
téléscopie des facteurs diagonaux donnent

```math
\frac{|U_1\cdots U_Ng|(z_0)}{h_0(z_0)}
\le
\mathbb E_{z_0}^{P_1,\ldots,P_N}
\left[
\exp\!\left(-\sum_{r=1}^N d_r(Z_{r-1},Z_r)\right)
\frac{|g(Z_N)|}{h_N(Z_N)}
\right].
\qquad\text{(5.5)}
```

Les états tels que $`h_r=0`$ sont simplement hors support ; aucune ligne
stochastique artificielle ne doit leur être ajoutée.

Ce théorème est fini, exact et inhomogène. Il ferme le problème abstrait de
normalisation commune signalé dans le fichier 28. Il ne ferme pas
l'identification de $`K_r,U_r`$ au corridor réel, ni le contrôle
thermodynamique du membre droit de (5.5).

L'implémentation rationnelle et ses audits sont dans
[`twisted_feynman_kac_composition.py`](../computations/twisted_feynman_kac_composition.py).

## 6. Diagnostics de la loi réelle à $`p=0.805`$

### 6.1 Les petites attaches en peigne croissent aux tailles testées

Sur le corridor final réel, on retient les nœuds stricts des deux bras tels
que

```math
2\le m\le6,
\qquad
J^{\mathrm{fav}}\le1,
\qquad
\text{taille de l'attache orientée}\le4.
\qquad\text{(6.1)}
```

Une paire uniforme est tirée dans chaque événement LCA réalisé et les
événements sont pondérés par $`N_\rho`$. Les erreurs sont des jackknifes
delete-one-environnement.

| $`L`$ | environnements | longueur du corridor | petites attaches | avec au plus 6 ports globaux |
|---:|---:|---:|---:|---:|
| 8 | 40 | $`12.645\pm0.265`$ | $`4.734\pm0.241`$ | $`0.748\pm0.052`$ |
| 12 | 40 | $`18.528\pm0.364`$ | $`7.316\pm0.268`$ | $`0.566\pm0.038`$ |
| 16 | 40 | $`26.085\pm0.513`$ | $`10.598\pm0.376`$ | $`0.381\pm0.029`$ |
| 24 | 20 | $`39.531\pm0.842`$ | $`17.338\pm0.472`$ | $`0.301\pm0.026`$ |

Le signal est net : les petites attaches croissent, mais le filtre
« frontière totale avec un nombre fixé de ports » les écrase. Même le cap
$`12`$ donne seulement $`3.177,3.087,2.726,2.298`$ pour
$`L=8,12,16,24`$.

La cellule T2 ne doit donc pas exiger une frontière globale bornée. Elle doit
soit localiser les ports dans une fenêtre, soit transporter l'extérieur dans
un potentiel projectif. Le calcul et tous les audits de reconstruction sont
dans
[`corridor_t2_signature_diagnostic.py`](../computations/corridor_t2_signature_diagnostic.py).

### 6.2 Les messages ancestraux modérés ne disparaissent pas

Pour chaque candidat $`2\le m\le8`$ de charge au plus $`1`$, le diagnostic
calcule le vrai log-odds externe $`B`$ produit par tous ses ancêtres stricts,
avec les facteurs aux rangs réalisés. Les colonnes suivantes comptent les
candidats par paire Palm.

| $`L`$ | environnements | candidats | $`\lvert B\rvert\le1`$ | $`\lvert B\rvert\le2`$ | $`\lvert B\rvert\le4`$ | médiane de $`\lvert B\rvert`$ fini |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 40 | $`5.361\pm0.225`$ | $`1.476\pm0.141`$ | $`2.338\pm0.183`$ | $`3.696\pm0.217`$ | $`2.217`$ |
| 12 | 30 | $`8.364\pm0.374`$ | $`2.723\pm0.250`$ | $`4.134\pm0.291`$ | $`6.110\pm0.319`$ | $`1.925`$ |
| 16 | 20 | $`11.533\pm0.563`$ | $`3.969\pm0.234`$ | $`5.850\pm0.359`$ | $`8.556\pm0.481`$ | $`1.764`$ |
| 24 | 8 | $`19.088\pm0.636`$ | $`8.053\pm0.416`$ | $`11.322\pm0.355`$ | $`15.026\pm0.577`$ | $`1.239`$ |

La taille $`L=24`$ n'a que huit environnements et doit rester exploratoire.
Sur $`L\le24`$, les données n'observent donc pas le scénario simple où tous
les candidats seraient progressivement polarisés par leurs ancêtres. Elles
ne prouvent ni une tendance asymptotique, ni screening latéral, ni
indépendance, ni déficit positif du transfert T2.

Le calcul exact du message est dans
[`ancestral_polarization_palm_diagnostic.py`](../computations/ancestral_polarization_palm_diagnostic.py).

### 6.3 Le filtre conjoint survit aux rangs réels

Les deux diagnostics précédents ne portaient pas exactement sur les mêmes
nœuds. Un troisième calcul repart donc des **mêmes environnements, mêmes
événements LCA, mêmes paires et mêmes petites attaches**. Son filtre
principal utilise le rang réalisé, sans écrêtage :

```math
2\le m\le8,
\qquad
\text{taille de l'attache orientée}\le4,
\qquad
m h_p(q_v)^2\le1.
\qquad\text{(6.2)}
```

Il intersecte ensuite ce filtre avec le vrai message ancestral $`B`$. Les
erreurs sont de nouveau des jackknifes delete-one-environnement.

| $`L`$ | environnements | charge au rang réel | puis $`\lvert B\rvert\le1`$ | puis $`\lvert B\rvert\le2`$ | puis $`\lvert B\rvert\le4`$ |
|---:|---:|---:|---:|---:|---:|
| 8 | 24 | $`4.476\pm0.205`$ | $`1.387\pm0.118`$ | $`2.275\pm0.151`$ | $`3.425\pm0.200`$ |
| 12 | 10 | $`6.491\pm0.347`$ | $`2.289\pm0.275`$ | $`3.447\pm0.342`$ | $`5.267\pm0.322`$ |
| 16 | 5 | $`10.213\pm0.693`$ | $`4.164\pm0.744`$ | $`5.763\pm0.637`$ | $`7.906\pm0.600`$ |

Le proxy qui remplace $`q_v`$ par $`\min(q_v,q_c)`$ donne respectivement
$`4.456`$, $`6.467`$ et $`10.141`$ candidats avant le filtre sur $`B`$.
L'abondance observée ne provient donc pas de la criticalisation invalide.
En revanche, le nombre géométrique moyen de ports globaux par nœud du
corridor croît de $`10.17`$ à $`15.91`$, puis $`21.10`$. C'est une raison
supplémentaire de transporter un état frontière ou un twist extérieur, pas
un résultat de screening.

La taille $`L=16`$ ne comporte que cinq environnements et reste
exploratoire. Sur les trois tailles testées, ce calcul observe des comptes
croissants de **candidats** ; il ne démontre ni leur abondance asymptotique,
ni $`K_r,U_r`$, ni le déficit $`d_r`$, ni une borne de weak recovery. Les
audits conjoints sont dans
[`joint_real_rank_t2_palm_diagnostic.py`](../computations/joint_real_rank_t2_palm_diagnostic.py).

### 6.4 Un update fidèle complètement résolu a un déficit nul

Soit un transfert positif levé $`T(z,z',\epsilon)`$, avec
$`\epsilon\in\{-1,+1\}`$, et posons

```math
K(z,z')
=
\sum_{\epsilon}T(z,z',\epsilon),
\qquad
U(z,z')
=
\sum_{\epsilon}\epsilon T(z,z',\epsilon).
\qquad\text{(6.3)}
```

Si le signe est mesurable depuis la transition conservée, c'est-à-dire s'il
existe $`g(z,z')\in\{-1,+1\}`$ tel que
$`T(z,z',\epsilon)>0`$ implique $`\epsilon=g(z,z')`$, alors

```math
U(z,z')=g(z,z')K(z,z'),
\qquad
|U(z,z')|=K(z,z'),
\qquad
d(z,z')=0
\quad\text{sur le support}.
\qquad\text{(6.4)}
```

C'est exactement ce qui arrive pour un heat bath d'une petite attache lorsque
l'état cible conserve les deux configurations de spins complètes. Le twist
répliqué $`\epsilon=s_1s_2`$ se lit dans les deux orientations de la branche
portant l'endpoint. L'enveloppe de Feynman--Kac vaut donc **exactement un**,
quel que soit le rang, le bucket ou le potentiel extérieur.

Sur une cellule réelle sélectionnée à $`p=0.805`$, de rang
$`q=0.3688343\ldots`$, bucket $`m=5`$ et attache de taille un, projeter sur
les seules orientations relatives donne

```math
\sum_{z'}|U_{\mathrm{proj}}(z,z')|
=
0.991774347976\ldots,
\qquad
-\log(0.991774347976\ldots)
=
0.008259669371\ldots.
\qquad\text{(6.5)}
```

Cette cancellation est exacte en un pas mais **non composable** : la
projection oublie les orientations qui transforment le potentiel ancestral
$`\Psi`$ au pas suivant. Le nombre (6.5) n'est donc pas un déficit T2
certifié.

La conséquence structurelle est nette. Un déficit composable ne peut
apparaître qu'au moment où plusieurs valeurs de $`\epsilon`$ sont sommées
dans une même transition coarse-grained, tout en conservant assez d'état pour
le futur. Deux constructions restent légitimes :

1. une jauge de ports dont la transition de $`\Psi`$ est exactement fermée ;
2. un bloc multi-update qui élimine une orientation seulement après sa
   dernière interaction avec tous les facteurs futurs.

Le prototype et ses contre-audits sont dans
[`real_rank_t2_deficit_prototype.py`](../computations/real_rank_t2_deficit_prototype.py).

### 6.5 La dernière incidence locale ne libère pas le twist global

Pour une fusion du corridor dont un enfant est une petite attache, prenons
l'autre enfant, situé sur le bras principal, comme jauge. Le dernier ancêtre
dont le bucket physique possède une arête incidente à l'attache donne une
**borne supérieure structurelle** sur son dernier usage : après ce rang,
tous les facteurs futurs sont exactement invariants sous le flip relatif de
l'attache. La réciproque est fausse, car des contributions incidentes peuvent
s'annuler dans les quatre $`\Lambda_v^{ab}`$.

Le diagnostic
[`last_use_attachment_palm_diagnostic.py`](../computations/last_use_attachment_palm_diagnostic.py)
reconstruit cette borne à partir des LCA de toutes les arêtes physiques. Il
utilise les rangs réalisés postcritiques, pondère chaque événement par
$`N_\rho`$ seulement, tire une paire lointaine uniforme par événement et
calcule les erreurs par suppression d'un environnement entier. Avec
$`p=0.805`$, $`m\le8`$, une attache de taille au plus quatre et 24
environnements par taille, on obtient :

| $`L`$ | candidats par paire | aucune incidence future de l'attache | profondeur moyenne de la dernière incidence de l'attache | dernière incidence de l'union dans huit niveaux | dernière incidence de l'union à la racine |
|---:|---:|---:|---:|---:|---:|
| 8 | $`0.667\pm0.122`$ | $`0.354`$ | $`2.153\pm0.452`$ | $`0.938`$ | $`0.940`$ |
| 12 | $`1.246\pm0.184`$ | $`0.402`$ | $`3.961\pm0.578`$ | $`0.443`$ | $`0.939`$ |
| 16 | $`2.825\pm0.270`$ | $`0.382`$ | $`7.551\pm0.809`$ | $`0.106`$ | $`0.962`$ |

L'élimination terminale des petites attaches est donc réellement utile pour
réduire l'état. Mais elle ne marginalise pas l'orientation commune de
l'union, qui porte le twist de la paire : le seul certificat par absence
d'incidence la garde presque toujours jusqu'à la racine aux tailles testées.
Des cancellations exactes peuvent avancer son dernier usage fonctionnel,
mais elles doivent être incorporées dans une jauge Markov-fermée ; les
supposer reviendrait à répéter la projection non composable de (6.5).

Ces données de volume fini ne prouvent aucune loi asymptotique. Combinées au
no-go exact $`|U|=K`$ sur l'état fidèle, elles ne soutiennent toutefois plus
la version simple de R1 dans laquelle une succession d'attaches terminales
produirait directement des déficits locaux additifs. Sans une nouvelle jauge
de ports ou un opérateur global, cette sous-route est une impasse probable
aux tailles testées, pas un théorème d'impossibilité universel.

## 7. Quantité maître conditionnelle à une jauge locale

Si R1 fournit une jauge Markov-fermée de complexité contrôlable, il faut
abandonner la persistance d'un corridor artificiellement criticalisé. À une
trajectoire de la chaîne de masse réelle, associons alors

```math
D_N
=
\sum_{r=1}^N d_r(Z_{r-1},Z_r).
\qquad\text{(7.1)}
```

Le théorème fini montre que la cible pertinente est

```math
\mathbb E^{P_1,\ldots,P_N}
\left[e^{-D_N}\right]
\longrightarrow0.
\qquad\text{(7.2)}
```

Si l'on remplace les fonctions rétrogrades exactes $`h_r`$ par une jauge
locale ou stationnaire calculable, un cocycle de changement de jauge
$`G_N`$ apparaît. Le critère devient alors

```math
D_N-G_N\longrightarrow+\infty
\quad\text{en probabilité, avec intégrabilité uniforme suffisante.}
\qquad\text{(7.3)}
```

Cette formulation accepte les passages très polarisés où $`d_r\simeq0`$.
Elle demande seulement assez de visites régénératives à déficit positif.

## 8. Dernier contre-test de la voie locale

Les étapes R1--R5 ci-dessous sont désormais **conditionnelles**. Elles ne
constituent plus la voie prioritaire : R1 est conservée comme dernier
contre-test fini d'une compression spéciale, puis R2--R5 restent suspendues
à son éventuel succès. Le pivot actif par projections collapsed et
dissipation du secteur impair est détaillé dans le
[fichier 30](../active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md).

### R0 — corrigendum logique

Retirer de toutes les notes canoniques la criticalisation multiport comme
théorème. Conserver explicitement le lemme scalaire et marquer les courbes
criticalisées des anciens diagnostics comme des proxies sans ordre.

**Sortie :** aucune implication de weak recovery ne dépend d'un oracle faux.

### R1 — transfert T2 aux rangs réels

Le nouvel état ne doit pas conserver le twist lui-même. Construire un
séparateur de frontière minimal :

```math
z=(\Pi_{\mathrm{loc}},\Psi,\xi_{\mathrm{front}},\mathfrak p),
\qquad\text{(8.1)}
```

où $`\Pi_{\mathrm{loc}}`$ est la partition signée des ports dans une fenêtre,
$`\Psi`$ le potentiel projectif extérieur, $`\xi_{\mathrm{front}}`$ les
orientations encore nécessaires aux facteurs futurs et $`\mathfrak p`$ le
type de fusion ou d'attache. Garder le rang $`q_r`$ ou le temps $`\beta_r`$
observé, sans écrêtage. Une orientation interne n'est sommée qu'après sa
dernière incidence avec le séparateur.

**Sortie :** une transition quotient Markov-fermée de $`\Psi`$ avec un
déficit positif, ou un no-go. Le scan de dernière incidence peut éliminer les
attaches hors-spine, mais il ne suffit pas pour l'orientation globale ; un
bloc qui la conserve jusqu'à la racine ne fournit aucune densité locale de
déficit. Les tables $`K_r,U_r,r_r,d_r`$ doivent être contre-auditées par
élimination directe et par comptes groupés. Si le twist reste mesurable dans
l'état cible, le test doit échouer automatiquement avec $`d_r=0`$.

### R2 — joindre géométrie, message et déficit

La jointure géométrie--rang--message et le scan de dernière incidence sont
maintenant implémentés pour la même petite attache. **R2 est suspendue à la
réussite de R1.** Il resterait alors à leur ajouter :

- la géométrie locale des ports ;
- la borne de dernière incidence et les cancellations exactes du quotient ;
- le déficit T2 **composable** calculé par la table de R1.

Pour des tentatives spatiales $`\tau_1,\ldots,\tau_{k_L}`$, le diagnostic doit
compter les cellules simultanément protégées, contractantes et visitées avec
un potentiel modéré. La statistique falsifiable est

```math
S_{L,M}
=
\sum_{\ell=1}^{k_L}
\mathbf 1_{G_\ell}
\mathbf 1_{\{|\Psi_{\tau_\ell^-}|\le M\}}.
\qquad\text{(8.2)}
```

**Porte go/no-go :** poursuivre seulement si la transformée de Laplace de
$`S_{L,M}`$ décroît avec $`L`$ pour au moins un $`M`$ robuste, et si R1 donne
sur ces mêmes visites un déficit de bloc strictement positif. Une médiane du
déficit seule ne contrôle pas les hiérarchies rares qui dominent le second
moment annealed.

### R3 — blocs régénératifs spatiaux

Ne chercher ni indépendance le long de la profondeur brute, ni propriété de
Markov du MST. Revenir aux labels d'arêtes i.i.d. avant Kruskal et explorer
des anneaux spatialement séparés. Le motif candidat est un **outlet épaissi
protégé** : deux bras ouverts, deux bras duaux, au moins deux arêtes locales
concurrentes entre les enfants et une petite attache orientée. Un pivot
isolé de bucket $`m=1`$ est parfait et ne donne aucun déficit.

Le gadget ne doit pas demander de bras macroscopiques supplémentaires par
rapport à l'événement pivotal déjà conditionné. S'il exige typiquement un
événement à six bras, sa fréquence par échelle peut être sommable et la
route devient une impasse.

**Outils plausibles :** RSW, séparation de bras, finite energy et
quasi-multiplicativité marquée, appliqués aux labels bruts puis transportés
vers l'événement T2. Les annuli restent dépendants après conditionnement par
la connexion de la paire, la gagnante de Kruskal et la mesure de Doob. Une
estimation de Laplace suffit ; il n'est pas nécessaire de déterminer la loi
limite complète du MST ou de la Palm.

### R4 — théorème annealed de déficit

La troncature par $`\max_r|\Psi_r|\le M`$ est une mauvaise cible : pour une
chaîne longue à potentiel non borné, sa probabilité peut tendre vers zéro et
rendre l'énoncé vrai de façon vide. Il faut permettre les excursions
polarisées et contrôler le nombre de retours modérés.

Prouver, sous la Palm réelle et la mesure de Doob $`P^h`$, qu'il existe
$`k_L\ge c\log L`$, $`\lambda>0`$ et, pour un $`M`$ fixé,

```math
\mathbb E_{\mathrm{ann},P^h}
\left[
e^{-\lambda S_{L,M}}
\,\middle|\,
I_L\leftrightarrow J_L\text{ dans }\Pi_1
\right]
\le
C_M\theta_M^{k_L}+o_L(1),
\qquad
\theta_M<1.
\qquad\text{(8.3)}
```

Si R1 certifie sur chaque visite comptée

```math
\Delta D_\ell\ge\delta_M>0,
\qquad\text{(8.4)}
```

alors $`D_{N_L}\ge\delta_M S_{L,M}`$ et (8.3) ferme directement (7.2).
Il suffit donc d'une récurrence ou occupation annealed de la zone modérée,
pas d'une contraction uniforme en $`\Psi`$.

### R5 — certificat à $`p_0=161/200`$

Rationaliser seulement les cellules effectivement visitées, certifier les
intervalles de $`d_r`$, puis insérer (7.2) dans le théorème H. Tester
$`p=0.81`$ uniquement après fermeture de $`0.805`$.

## 9. Porte de décision honnête

### Continuer la voie du déficit local si

- une jauge de ports Markov-fermée agrège effectivement des signes opposés
  sans conserver l'orientation globale jusqu'à la racine ;
- sa complexité de séparateur reste contrôlable le long du corridor réel ;
- les petites attaches conservent une densité positive dans des fenêtres
  locales séparables ;
- la transformée de Laplace de $`S_{L,M}`$ décroît pour une zone modérée ;
- le déficit T2 reste positif sur cette zone malgré les ports extérieurs ;
- les orientations portant le twist quittent le séparateur dans un nombre
  divergent de blocs ;
- la masse des potentiels extrêmes ne compense pas ce déficit ;
- la normalisation rétrograde peut être localisée avec un cocycle sous-linéaire.

### Arrêter la voie du déficit local si

- $`D_{N_L}`$ reste tendu sous la chaîne de masse réelle ;
- les seuls déficits visibles exigent un bornage gratuit du bord ;
- tout état Markov-fermé de complexité contrôlable conserve le twist jusqu'à
  la racine et rend $`d=0`$ localement ;
- les outlets T2 protégés exigent un événement de bras trop rare ;
- les retours de $`\Psi`$ dans toute zone modérée restent tendus ;
- le coût de jauge $`G_N`$ croît au même rythme que le déficit ou plus vite.

Le test de fermeture R1 doit précéder R2. Il ne doit toutefois plus retarder
le diagnostic $`L^2`$ global du fichier 30, qui n'exige pas de quotient
Markov de dimension fixe.

## 10. Verdict stratégique

La route **inchangée** doit être abandonnée, et la réparation naïve par
attaches terminales est elle aussi une impasse probable. Deux faits exacts se
renforcent : l'état fidèle donne $`d=0`$, tandis que la projection qui donne
$`d>0`$ n'est pas Markov-fermée. Le scan géométrique ajoute que la seule
absence d'incidence ne libère presque jamais l'orientation globale avant la
racine aux tailles testées.

Une jauge exacte existe formellement en transportant l'orbite complète de
$`\Psi`$. Mais avec $`b`$ ports, un potentiel positif projectif général a
$`2^{b-1}-1`$ degrés de liberté, et le nombre de ports observé augmente avec
la taille. C'est donc du bookkeeping à frontière croissante, pas la fermeture
locale bornée recherchée. Cela ne réfute pas une compression spéciale des
potentiels effectivement atteignables ; R1 reste son dernier contre-test.

La voie prioritaire devient la dissipation $`L^2`$ du caractère de paire sous
des blocs collapsed imbriqués. Elle possède une identité de Pythagore exacte
et conserve les cancellations perdues par l'enveloppe $`|U|`$. Son verrou est
une inégalité annealed **pondérée par l'énergie du signal** sur des annuli, et
non une norme d'opérateur uniforme, impossible à cause des constantes. La
[feuille de pivot](../active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) donne le
théorème suffisant, le diagnostic fini et ses portes d'arrêt.

Le premier diagnostic exact $`L=4`$ est mixte mais non trivial. À
$`p=0.805`$, la perte logarithmique collapsed moyenne vaut $`0.47603`$ sur
les paires connectées : $`0.34691`$ avant le LCA et $`0.12912`$ au LCA. La
cancellation globale existe donc réellement. En revanche, malgré $`3.15`$
niveaux à perte positive, les pertes absolues ne représentent que $`1.209`$
niveaux effectifs. Ce volume fini justifie le test à deux updates du fichier
30 ; il ne prouve aucune accumulation avec l'échelle.

Une borne hiérarchique à $`p=0.805`$ est donc **spéculative**, non
« plausible » au sens probatoire. Le seuil exact de Nishimori demanderait
presque certainement une idée supplémentaire. La voie locale A0 par triangle
reste un chantier parallèle pour une amélioration numérique, sans être
confondue avec la preuve hiérarchique recherchée ici.
