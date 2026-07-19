# Information portée par les coupes conditionnées par la géométrie

> [!NOTE]
> Cette note fournit les lemmes de coupe utilisés par le
> [programme prioritaire](00_RESEARCH_PROGRAM.md). La feuille de route
> technique complète est le [fichier 05](05_PROOF_ROADMAP.md).

Cette note corrige le niveau de conditionnement de la voie « épuisement du
signal ». Un bilan sur toutes les arêtes du graphe, ou même sur toutes les
arêtes internes à une composante, ne mesure pas le taux $`\Lambda_u`$ d'un
nœud hiérarchique. Le bon objet est une coupe physique entre deux blocs de la
partition ouverte juste avant leur fusion.

Le calcul exact donne toutefois une conclusion plus nuancée que l'intuition
initiale.

> **Conclusion principale.** Conditionnellement à la partition complète au
> temps $`\beta`$, les marques des arêtes de frontière sont encore
> indépendantes et ont toutes le même biais résiduel $`h_p(\beta)`$. La
> géométrie intervient par la taille $`m`$ de la coupe et par sa sélection
> Palm. Une coupe instantanée perd son information lorsque
> $`m h_p(\beta)^2\to0`$, non lorsqu'un bilan global change de signe.

Une fusion réelle contient en outre une arête gagnante conforme. Sans message
extérieur, cette correction Palm laisse une fiabilité $`1/m`$ à
$`\beta=1`$ et rend une coupe de taille un parfaitement informative. Enfin,
le LCA d'une paire lointaine
ne voit pas une coupe typique : sa loi est repondérée à la fois par $`m`$ et
par le nombre de paires lointaines séparées par la coupe.

Ces faits établissent une partie exacte du programme géométrique. Ils ne
prouvent pas encore une obstruction globale à la weak recovery, car les
messages ancestraux et les routes latérales doivent aussi être contrôlés.

## 1. Cadre et trois expériences à ne pas confondre

Dans la jauge plantée,

```math
O_e=\Sigma_x\Sigma_y Z_e,
\qquad
\mathbb P(Z_e=+1)=p,
\qquad
u_p=\log\frac{p}{1-p}.
\tag{1.1}
```

Une arête conforme $`Z_e=+1`$ reçoit une horloge exponentielle de taux
$`u_p`$ ; une arête fausse a une horloge infinie. Posons

```math
A_e(\beta)
=
\mathbf 1_{\{Z_e=+1,\ \xi_e\le\beta\}},
\qquad
\Pi_\beta
=
\mathrm{cc}\bigl(V,\{e:A_e(\beta)=1\}\bigr).
\tag{1.2}
```

Pour deux blocs distincts $`C_1,C_2\in\Pi_\beta`$, définissons

```math
E(C_1,C_2)
=
\{\{x,y\}\in E:x\in C_1,\ y\in C_2\},
\qquad
m(C_1,C_2)=|E(C_1,C_2)|.
\tag{1.3}
```

Il faut séparer trois expériences.

1. **Coupe instantanée.** On observe deux blocs présents à un temps
   déterministe $`\beta`$, sans conditionner sur leur prochaine fusion.
2. **Coupe de fusion.** On conditionne sur le fait que cette coupe est
   activée à $`\beta`$. Une arête gagnante conforme est alors forcée.
3. **Coupe LCA-Palm.** On choisit la fusion qui est le LCA d'une paire
   lointaine. Aux deux conditionnements précédents s'ajoute un biais par le
   nombre de paires séparées par les deux enfants.

Le premier objet isole le rôle de la taille géométrique. Le deuxième est le
bucket réellement utilisé dans la dynamique hiérarchique. Le troisième est
celui qui intervient dans le critère pairwise de weak recovery.

## 2. La géométrie interne ne contamine pas la loi d'une frontière fixée

Pour une partition $`\pi`$, notons

```math
\partial\pi
=
\{e=\{x,y\}\in E:x\text{ et }y
\text{ appartiennent à deux blocs distincts de }\pi\}.
\tag{2.1}
```

Définissons la probabilité résiduelle d'une arête conforme par

```math
s_p(\beta)
=
\frac{pe^{-u_p\beta}}{1-p+pe^{-u_p\beta}},
\qquad
h_p(\beta)
=
2s_p(\beta)-1
=
\tanh\left(\frac{u_p(1-\beta)}2\right).
\tag{2.2}
```

### Théorème 2.1 — factorisation conditionnelle de coupe, établi

Fixons une partition réalisable $`\pi`$ et supposons
$`\mathbb P(\Pi_\beta=\pi)>0`$. Conditionnellement à
$`\Pi_\beta=\pi`$, les catégories des arêtes de $`\partial\pi`$ sont
indépendantes et ont la loi commune

```math
(\text{vraie tardive},\text{vraie censurée},\text{fausse})
\sim
\left(
h_p(\beta),
\frac{1-h_p(\beta)}2,
\frac{1-h_p(\beta)}2
\right).
\tag{2.3}
```

En particulier, chaque arête de frontière est conforme avec probabilité
$`s_p(\beta)`$.

#### Preuve

L'événement $`\{\Pi_\beta=\pi\}`$ se factorise en deux parties portant sur
des familles disjointes d'arêtes :

```math
\{\Pi_\beta=\pi\}
=
\left(
\bigcap_{C\in\pi}
\{G_\beta[C]\text{ est connexe}\}
\right)
\cap
\left(
\bigcap_{e\in\partial\pi}\{A_e(\beta)=0\}
\right).
\tag{2.4}
```

La première parenthèse ne dépend que des arêtes internes aux blocs. La
seconde est un produit de contraintes à une arête. L'indépendance initiale
des marques montre donc que les arêtes de frontière restent indépendantes,
chacune étant seulement conditionnée à ne pas avoir sonné avant
$`\beta`$.

Avant normalisation, les trois masses résiduelles sont

```math
pe^{-u_p\beta}-(1-p),
\qquad
1-p,
\qquad
1-p.
\tag{2.5}
```

Leur somme vaut $`1-p+pe^{-u_p\beta}`$. La normalisation donne
(2.3). Enfin $`e^{-u_p}=(1-p)/p`$, ce qui donne (2.2).

### Corollaire 2.2 — sélection mesurable par la partition, établi

Soit $`(C_1,C_2)`$ une paire de blocs choisie par une règle mesurable par
rapport à $`\Pi_\beta`$ et à une géométrie extérieure indépendante des
marques résiduelles. Conditionnellement à la partition et à cette paire, les
$`m(C_1,C_2)`$ marques de la coupe suivent encore (2.3).

Le corollaire s'applique notamment au bloc contenant un sommet fixé, au plus
grand bloc après une règle de départage déterministe, ou à deux blocs choisis
à partir de leurs positions. Il ne s'applique pas directement à la coupe qui
fusionnera ensuite : ce choix dépend d'une horloge future et requiert le
conditionnement Palm de la section 4.

En particulier, conditionner $`i,j`$ à appartenir à une composante géante,
ou sélectionner la plus grande composante, change la loi de la géométrie
$`m(C_1,C_2)`$ mais pas la loi des marques une fois la partition et la taille
fixées. Après mélange sur les partitions ayant la même taille $`m`$, le
nombre d'arêtes conformes d'une coupe mesurable par la partition reste
$`\mathrm{Bin}(m,s_p(\beta))`$. L'hypothèse de composante géante agit donc
sur les tailles, les formes et les états de bord, non sur le paramètre
unitaire $`s_p(\beta)`$.

### Ce que l'intuition géométrique avait correctement identifié

Toutes les arêtes déjà ouvertes sont conformes et internes à un bloc de
$`\Pi_\beta`$. La loi des arêtes internes est donc fortement biaisée par la
contrainte de connexité et ne factorise pas. Elles ne doivent pas être mises
dans le vote du nœud courant.

### La surprise du calcul conditionnel

Cette concentration interne des arêtes ouvertes ne retire pas une quantité
supplémentaire de vraies arêtes **fermées** aux frontières. Après
conditionnement par la partition complète, toute arête de frontière a
exactement la loi résiduelle (2.3). La géométrie n'altère donc pas son biais
unitaire ; elle détermine combien de telles arêtes sont présentes et comment
la coupe est sélectionnée.

## 3. Canal exact d'une coupe instantanée

Fixons une coupe de taille $`m\ge1`$ conditionnellement à la partition. Soit
$`X\in\{-1,+1\}`$ la parité relative des deux blocs, avec $`X=+1`$ pour la
parité conforme à la ground truth. Notons $`K`$ le nombre d'arêtes de coupe
satisfaites par cette parité.

Sous les deux hypothèses,

```math
K\mid X=+1
\sim
\mathrm{Bin}(m,s_p(\beta)),
\qquad
K\mid X=-1
\sim
\mathrm{Bin}(m,1-s_p(\beta)).
\tag{3.1}
```

Posons

```math
a_p(\beta)
:=
\log\frac{s_p(\beta)}{1-s_p(\beta)}
=
u_p(1-\beta).
\tag{3.2}
```

### Lemme 3.1 — log-rapport et moments, établi

Pour $`0\le k\le m`$, le log-rapport de vraisemblance est

```math
L^{\mathrm{snap}}_{m,k}
=
a_p(\beta)(2k-m).
\tag{3.3}
```

Si $`V=2K-m`$ sous $`X=+1`$, alors

```math
\mathbb E[V]
=
m h_p(\beta),
\qquad
\mathrm{Var}(V)
=
m\bigl(1-h_p(\beta)^2\bigr).
\tag{3.4}
```

Le rapport signal sur bruit du vote vaut donc exactement

```math
\boxed{
\mathrm{SNR}^{\mathrm{snap}}_{m,\beta}
=
\frac{m h_p(\beta)^2}{1-h_p(\beta)^2}.
}
\tag{3.5}
```

La preuve est la division des deux masses binomiales, puis le calcul des
moments d'une binomiale.

### Charge informationnelle exacte

Le coefficient de Bhattacharyya des deux expériences à une arête est
$`\sqrt{1-h_p(\beta)^2}`$. L'information de Chernoff de la coupe est donc

```math
\boxed{
\mathcal I^{\mathrm{snap}}_{m,\beta}
=
-\frac m2\log\bigl(1-h_p(\beta)^2\bigr).
}
\tag{3.6}
```

Dans le régime de faible biais,

```math
\mathcal I^{\mathrm{snap}}_{m,\beta}
\sim
\frac12\mathcal J_{m,\beta},
\qquad
\mathcal J_{m,\beta}
:=
m h_p(\beta)^2.
\tag{3.7}
```

La variable $`\mathcal J`$ est le premier résumé géométrique pertinent. Elle
ne remplace pas le canal exact lorsque le biais n'est pas petit.

### Fiabilité $`L^2`$

Sous un prior uniforme et sans message extérieur, définissons

```math
\widetilde\Gamma_m(\beta;p)
=
\mathbb E_{+}\left[
\tanh^2\left(
\frac{a_p(\beta)(2K-m)}2
\right)
\right].
\tag{3.8}
```

Par symétrie, cette quantité est aussi

```math
\widetilde\Gamma_m(\beta;p)
=
\frac12
\sum_{k=0}^m
\frac{(P_+(k)-P_-(k))^2}{P_+(k)+P_-(k)}.
\tag{3.9}
```

### Théorème 3.2 — diagramme asymptotique géométrique, établi

Soient $`m_L\ge1`$, $`\beta_L\in[0,1]`$ et
$`h_L=h_p(\beta_L)`$.

1. Si $`m_Lh_L^2\to0`$, alors
   $`\widetilde\Gamma_{m_L}(\beta_L;p)\to0`$.
2. Si $`m_Lh_L^2\to\infty`$, alors
   $`\widetilde\Gamma_{m_L}(\beta_L;p)\to1`$.
3. Si en outre $`m_L\to\infty`$, $`h_L\to0`$ et
   $`m_Lh_L^2\to\lambda\in(0,\infty)`$, alors, pour
   $`Z\sim\mathcal N(0,1)`$,

```math
\widetilde\Gamma_{m_L}(\beta_L;p)
\longrightarrow
\mathbb E\left[
\tanh^2\bigl(\sqrt\lambda Z+\lambda\bigr)
\right].
\tag{3.10}
```

#### Preuve du régime sans information

Écrivons $`a_L=2\,\mathrm{artanh}(h_L)`$ et $`L_L=a_LV_L`$. Si
$`m_Lh_L^2\to0`$, alors $`h_L\to0`$ et $`a_L\sim2h_L`$. Par (3.4),

```math
\mathbb E[L_L^2]
=
a_L^2
\left[
m_L(1-h_L^2)+m_L^2h_L^2
\right]
\longrightarrow0.
\tag{3.11}
```

Comme $`\tanh^2(x/2)\le x^2/4`$, la première conclusion suit.

#### Preuve du régime informatif

Sous $`m_Lh_L^2\to\infty`$,

```math
\frac{V_L-m_Lh_L}{m_Lh_L}
\longrightarrow0
\quad\text{en probabilité}.
\tag{3.12}
```

De plus $`a_Lm_Lh_L\to\infty`$ : si $`h_L\to0`$, ce produit est
asymptotique à $`2m_Lh_L^2`$ ; sinon il diverge linéairement en $`m_L`$.
Ainsi $`L_L\to+\infty`$ en probabilité sous l'hypothèse vraie. La convergence
bornée donne la deuxième conclusion.

#### Preuve de la fenêtre critique

Le théorème central limite donne

```math
\frac{V_L}{\sqrt{m_L}}
\Longrightarrow
Z+\sqrt\lambda,
\qquad
a_L\sqrt{m_L}
\longrightarrow
2\sqrt\lambda.
\tag{3.13}
```

Donc $`L_L/2\Longrightarrow\sqrt\lambda Z+\lambda`$. La fonction
$`\tanh^2`$ est continue et bornée, ce qui prouve (3.10).

### Corollaire 3.3 — l'instant de perte dépend de la taille de coupe

Au voisinage de la censure,

```math
h_p(\beta)
=
\frac{u_p}{2}(1-\beta)
+O\bigl((1-\beta)^3\bigr).
\tag{3.14}
```

Une coupe instantanée de taille $`m`$ entre donc dans sa fenêtre de perte
autour de

```math
\boxed{
1-\beta
\asymp
\frac{2}{u_p\sqrt m}.
}
\tag{3.15}
```

Si la taille elle-même dépend du temps, le critère devient

```math
m(\beta)(1-\beta)^2\asymp1.
\tag{3.16}
```

Il n'existe donc aucun seuil déterministe en $`\beta`$ sans théorème sur la
géométrie de $`m(\beta)`$. Une croissance de la coupe peut compenser, voire
surcompenser, la baisse du biais par arête.

## 4. Une fusion réelle est une coupe Palm, pas une coupe instantanée

### Lemme 4.1 — taux de fusion direct d'une coupe, établi

Conditionnellement à $`\Pi_\beta=\pi`$, une arête de frontière est conforme
avec probabilité $`s_p(\beta)`$. Conditionnellement à être conforme et à ne
pas avoir sonné, la propriété sans mémoire donne un taux instantané $`u_p`$.
Par conséquent, une coupe $`E(C_1,C_2)`$ de taille $`m`$ a le taux direct

```math
\boxed{
r_\beta(C_1,C_2)
=
m(C_1,C_2)u_ps_p(\beta).
}
\tag{4.1}
```

À partition et temps fixés, la prochaine coupe est donc sélectionnée avec un
biais linéaire en $`m`$. C'est un premier effet géométrique exact : les
grandes coupes sont surreprésentées parmi les fusions observées.

### Formule de Campbell pour le LCA d'une paire lointaine

Sur un tore de côté $`L`$, fixons $`\rho>0`$ et posons

```math
N_\rho(A,B)
=
\#\{(x,y)\in A\times B:d_L(x,y)\ge\rho L\}
+
\#\{(x,y)\in B\times A:d_L(x,y)\ge\rho L\}.
\tag{4.2}
```

Pour toute fonction prévisible bornée $`F`$, la compensation des sauts du
processus de composantes donne

```math
\begin{aligned}
&\mathbb E\left[
\sum_{u:\,\beta_u\le1}
N_\rho(C_{u,1},C_{u,2})
F(\Pi_{\beta_u-},C_{u,1},C_{u,2},\beta_u)
\right]
\\
&\quad=
\int_0^1
\mathbb E\left[
\sum_{\{A,B\}\subset\Pi_\beta}
u_ps_p(\beta)m(A,B)N_\rho(A,B)
F(\Pi_\beta,A,B,\beta)
\right]d\beta.
\end{aligned}
\tag{4.3}
```

À un temps déterministe, $`\Pi_{\beta-}=\Pi_\beta`$ presque sûrement, ce qui
justifie la version à droite. L'identité (4.3) est simplement la formule
d'intensité (4.1), sommée sur les coupes, puis repondérée par le nombre de
paires dont le LCA est créé par le saut.

Ainsi, à niveau fixé, la loi LCA-Palm d'une paire lointaine est repondérée par

```math
\boxed{
m(A,B)N_\rho(A,B).
}
\tag{4.4}
```

Ce facteur est absent d'un calcul sur une coupe typique. Il favorise à la
fois les interfaces riches en arêtes et les fusions de deux enfants portant
beaucoup de paires macroscopiquement séparées.

### Lemme 4.2 — canal exact de la coupe de fusion, établi

Conditionnellement à la géométrie $`(C_1,C_2,m)`$ et au niveau de fusion
$`\beta`$, l'arête gagnante est conforme et uniforme parmi les $`m`$ arêtes
de coupe. Les $`m-1`$ marques non gagnantes conservent la loi (2.3). Par
conséquent,

```math
K\mid X=+1
\sim
1+\mathrm{Bin}(m-1,s_p(\beta)),
\qquad
K\mid X=-1
\sim
\mathrm{Bin}(m-1,1-s_p(\beta)).
\tag{4.5}
```

Le log-rapport local vaut, pour $`1\le k\le m-1`$,

```math
L^{\mathrm{merge}}_{m,k}
=
\log\frac{k}{m-k}
+
a_p(\beta)(2k-m),
\tag{4.6}
```

avec les conventions infinies aux deux extrémités. Le premier terme est
exactement la correction Palm de l'arête gagnante.

Sous la parité vraie, le vote $`V=2K-m`$ vérifie

```math
\mathbb E[V]
=
1+(m-1)h_p(\beta),
\qquad
\mathrm{Var}(V)
=
(m-1)\bigl(1-h_p(\beta)^2\bigr).
\tag{4.7}
```

Pour $`m\ge2`$, son rapport signal sur bruit est

```math
\mathrm{SNR}^{\mathrm{merge}}_{m,\beta}
=
\frac{[1+(m-1)h_p(\beta)]^2}
{(m-1)[1-h_p(\beta)^2]}.
\tag{4.8}
```

### Théorème 4.3 — ce qui subsiste à la censure, établi

Si $`\Gamma_m^{\mathrm{merge}}(\beta;p)`$ désigne la fiabilité locale
définie à partir de (4.6), alors

```math
\boxed{
\Gamma_m^{\mathrm{merge}}(1;p)=\frac1m.
}
\tag{4.9}
```

En particulier, un bucket $`m=1`$ est parfait à tout niveau. Pour $`m`$ fixé,
une fusion ne devient donc jamais exactement non informative. Si en revanche
$`m_L\to\infty`$ et $`h_p(\beta_L)\to0`$, les trois régimes du théorème 3.2
restent valables pour la fusion : le terme gagnant est négligeable à
l'échelle $`\sqrt m`$, et la limite dans la fenêtre
$`m_Lh_p(\beta_L)^2\to\lambda`$ est encore (3.10).

Pour le vérifier, écrivons encore $`V_L=2K_L-m_L`$. Dans la fenêtre finie,

```math
\frac{V_L}{\sqrt{m_L}}
\Longrightarrow
Z+\sqrt\lambda,
\qquad
\log\frac{K_L}{m_L-K_L}
=
2\,\mathrm{artanh}\left(\frac{V_L}{m_L}\right)
=
o_{\mathbb P}(1).
\tag{4.9a}
```

Le terme $`a_p(\beta_L)V_L/2`$ a donc la même limite que dans (3.13).
Lorsque $`m_Lh_p(\beta_L)^2\to0`$, les deux termes du log-rapport tendent
vers zéro en probabilité et dans $`L^2`$. Lorsqu'il tend vers l'infini, le
terme résiduel diverge avec le bon signe en probabilité. Cela prouve les trois
régimes annoncés sans assimiler les deux expériences à taille finie.

L'identité (4.9) vient de

```math
\tanh\left(
\frac12\log\frac K{m-K}
\right)
=
\frac{2K-m}{m}
\tag{4.10}
```

et de $`\mathbb E[(2K-m)^2]=m`$ lorsque
$`K=1+\mathrm{Bin}(m-1,1/2)`$.

## 5. Conséquence au temps critique et à $`p=0.8`$

Sur la grille triangulaire, posons

```math
q_c=2\sin(\pi/18),
\qquad
\beta_c(p)
=
-\frac1{u_p}\log\left(1-\frac{q_c}{p}\right).
\tag{5.1}
```

Lorsque $`\beta_c\le1`$,

```math
h_c(p)
:=
h_p(\beta_c)
=
\frac{2p-1-q_c}{1-q_c}.
\tag{5.2}
```

À $`p=4/5`$,

```math
u_p=\log4,
\qquad
\beta_c=0.410716539196\ldots,
\qquad
h_c=0.387164445505\ldots.
\tag{5.3}
```

La charge d'une coupe instantanée critique est donc

```math
\mathcal J_{m,\beta_c}
=
0.149896\ldots\,m.
\tag{5.4}
```

Deux conclusions rigoureuses en découlent.

1. Une coupe critique dont la taille tend vers l'infini devient presque
   parfaitement informative. Elle ne fournit pas l'obstruction recherchée.
2. Le conditionnement LCA-Palm favorise justement les grandes coupes par le
   facteur $`mN_\rho`$ de (4.4). Montrer une perte d'information au LCA
   critique est donc plus difficile que sur une coupe choisie uniformément.

Le scénario « $`i,j`$ fusionnent dès $`\beta_c`$ » reste favorable par la
qualité unitaire maximale parmi les niveaux postcritiques. Mais son caractère
globalement le plus favorable ne suit pas de cette seule monotonie : les
lois de $`m`$, de $`N_\rho`$ et des états de bord changent avec le niveau.

## 6. Où la géométrie peut réellement produire une obstruction

Soit $`u`$ le LCA de $`i,j`$. Pour un parcours complet de la hiérarchie, le
bucket de $`u`$ n'est qu'un facteur. Il faut suivre les coupes du corridor
descendant et les coupes ancestrales qui contribuent aux quatre poids du
heat bath.

Pour chaque coupe pertinente $`v`$, enregistrons au minimum

```math
\mathscr G_v
=
(m_v,\beta_v,Z_v,B_v),
\tag{6.1}
```

où $`Z_v`$ est l'état géométrique des ports latéraux et $`B_v`$ le message
extérieur produit par les autres facteurs. Le diagnostic scalaire est

```math
\mathcal J_v
=
m_vh_p(\beta_v)^2.
\tag{6.2}
```

### Cas screené

Si $`B_v=o(1)`$ et si les routes latérales sont conditionnellement coupées,
alors $`\mathcal J_v\to0`$ implique que la coupe instantanée ne transmet plus
de parité. Pour une fusion de taille croissante, la même conclusion vaut.

### Cas non screené

Même lorsque $`\mathcal J_v\to0`$, un message $`B_v`$ non nul laisse la
fiabilité $`\tanh^2(B_v/2)`$. Une route latérale peut aussi contourner la
coupe. L'absence d'information dans le vote local ne suffit donc pas à
annuler la corrélation de $`i,j`$.

### Ancêtres stricts de $`u`$

Pour $`v\succ u`$, les deux flips descendants découpent la frontière de
$`v`$ en plusieurs groupes d'incidence. Le facteur exact est

```math
F_v(\Lambda_v^{ab})
=
\Lambda_v^{ab}
\exp\bigl((1-\beta_v)\Lambda_v^{ab}\bigr),
\qquad
a,b\in\{0,1\}.
\tag{6.3}
```

La factorisation du théorème 2.1 donne la loi des marques sur la frontière de
$`v`$ une fois la partition fixée. Elle ne permet pas de remplacer les
quatre $`\Lambda_v^{ab}`$ par une seule majorité : le groupe invariant reste
à l'intérieur de la fonction non linéaire $`F_v`$. Il faut donc estimer la
loi jointe des tailles de groupes et de l'état de bord sous la Palm de la
paire.

### Critère suffisant à viser

Une obstruction hiérarchique exploitable prendrait la forme suivante. Sous
la Palm favorable où le LCA de la paire lointaine appartient à la fenêtre
critique, extraire des coupes $`v_1,\ldots,v_{N_L}`$ telles que :

1. leurs transferts sont screenés conditionnellement par les états de bord ;
2. leur coefficient répliqué exact vérifie
   $`\eta_r\le\phi(\mathcal J_{v_r},Z_{v_r})<1`$ ;
3. les coefficients se composent ;
4. $`\sum_{r=1}^{N_L}-\log\eta_r\to\infty`$ en probabilité.

Alors le second moment de la parité de $`i,j`$ tend vers zéro. La quantité
$`\mathcal J_v`$ sert à localiser les coupes candidates ; le certificat final
doit porter sur le transfert répliqué complet.

## 7. Nouveau verrou géométrique précis

La prochaine question n'est plus « quelle proportion globale d'arêtes est
fausse ? », mais :

```math
\boxed{
\text{Quelle est, sous la loi LCA-Palm critique, la loi jointe de }
(m_v,\beta_v,Z_v,B_v)
\text{ le long du corridor de }i,j\,?
}
\tag{7.1}
```

Pour tester l'intuition d'une perte tardive, il faut en particulier estimer

```math
\#\left\{
v:\ m_vh_p(\beta_v)^2\le M,
\ |B_v|\le B_0,
\ v\text{ screené}
\right\}.
\tag{7.2}
```

Une version forte montrerait que ce nombre diverge pour des constantes
$`M,B_0`$ et que les transferts correspondants contractent uniformément. Une
version plus fine utiliserait directement la somme des déficits
$`-\log\eta_v`$ sans imposer une taille de coupe bornée.

La formule Palm (4.3) fournit le point de départ exact pour cette étude. Elle
montre aussi le danger d'un échantillonnage naïf : compter les coupes de la
partition uniformément ne reproduit pas la loi du LCA d'une paire lointaine.

## 8. Audit et contre-audit

| affirmation | verdict | justification |
|---|---|---|
| Les arêtes internes déjà ouvertes sont plus informatives | Vrai | elles sont conformes et sélectionnées par la connexité |
| Elles votent dans le $`\Lambda_u`$ de la fusion courante | Faux | seules les arêtes de $`E(C_1,C_2)`$ changent avec la parité relative |
| La frontière est appauvrie arbitrairement en vraies arêtes par la formation des clusters | Faux conditionnellement à $`\Pi_\beta`$ | factorisation exacte du théorème 2.1 |
| Une arête de frontière est encore plus souvent vraie que fausse avant $`1`$ | Vrai | $`h_p(\beta)>0`$ pour $`\beta<1`$ |
| Il existe un temps universel où toute coupe cesse d'être informative | Faux | la charge dépend de $`m h_p(\beta)^2`$ |
| Une grande coupe à temps fixé $`\beta<1`$ perd l'information | Faux | sa fiabilité tend vers un |
| Une coupe instantanée perd l'information si $`m h_p(\beta)^2\to0`$ | Établi | théorème 3.2 |
| Une fusion réelle a exactement la même loi qu'une coupe instantanée | Faux | l'arête gagnante ajoute le terme $`\log(k/(m-k))`$ |
| Une grande fusion terminale perd l'information | Établi | fiabilité $`1/m\to0`$ |
| Une fusion terminale de taille bornée est sans information | Faux | fiabilité $`1/m`$ ; $`m=1`$ est parfait |
| Le LCA voit une coupe typique de la partition | Faux | repondération $`mN_\rho`$ dans (4.4) |
| Un vote local non informatif suffit pour exclure la weak recovery | Faux | les messages ancestraux et routes latérales peuvent porter la parité |
| À $`p=0.8`$, une grande coupe au LCA critique donne l'obstruction | Faux | $`h_c>0`$ fixé et $`\mathcal J_c\asymp m`$ |

## 9. Résultats établis et lemmes encore ouverts

### Établi exactement

- factorisation des marques de toute frontière après conditionnement par la
  partition complète ;
- canal binomial d'une coupe instantanée ;
- transition gouvernée par $`m h_p(\beta)^2`$ ;
- taux de fusion direct $`m u_p s_p(\beta)`$ ;
- repondération LCA-Palm par $`mN_\rho`$ ;
- correction de l'arête gagnante et canal de fusion décalé ;
- reliquat terminal exact $`1/m`$.

### Ouvert sur le GSBM triangulaire

- loi asymptotique de $`m_v`$ sous la Palm d'une paire lointaine critique ;
- contrôle conjoint de $`m_v`$ et $`\beta_v`$ sur les descendants et les
  ancêtres ;
- screening des messages $`B_v`$ et des ports latéraux ;
- composition des contractions dans le sweep top-down ou bottom-up ;
- domination géométrique des paires postcritiques par l'expérience critique
  favorable ;
- déduction d'une nouvelle borne d'impossibilité, notamment à $`p=0.8`$.

## Conclusion

L'objection géométrique est juste : il faut raisonner coupe par coupe, jamais
avec la proportion globale d'arêtes. Après ce changement de conditionnement,
la perte d'information possède une formulation nette :

```math
\boxed{
\text{faible information locale}
\quad\Longleftrightarrow\quad
m(C_1,C_2)h_p(\beta)^2\ll1
}
\tag{9.1}
```

pour une coupe instantanée screenée, et asymptotiquement pour une fusion de
taille croissante. Le seuil correspondant dépend donc de la croissance
géométrique de la coupe. Sous le conditionnement LCA d'une paire lointaine,
la loi pertinente est en outre biaisée par $`mN_\rho`$.

La voie la plus rigoureuse vers une obstruction à $`p=0.8`$ est maintenant
de mesurer cette charge le long du corridor hiérarchique favorable, puis de
certifier la contraction des coupes de faible charge après screening. Le
bucket critique seul est généralement trop informatif ; la distance entre
$`i`$ et $`j`$ doit être exploitée par l'accumulation des coupes et par le
parcours complet de la hiérarchie.
