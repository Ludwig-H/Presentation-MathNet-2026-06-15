# Majorité conforme postcritique et chaîne des $`\Lambda_v`$

Cette note formalise l'idée suivante, dans la dynamique hiérarchique des
slides 31--33 : pour tester la weak recovery dans le cas le plus favorable,
on suit deux sommets lointains $`i,j`$ appartenant au même arbre et dont le
LCA $`u`$ naît au seuil de percolation. Le nœud $`u`$ voit alors le résidu le
plus informatif parmi les fusions postcritiques. Il faut déterminer si les
liens de ce résidu ont une majorité stricte conforme à la ground truth, puis
vérifier si cette majorité survit dans les quatre taux
$`\Lambda_v^{ab}`$ de **chaque** ancêtre $`v\succeq u`$.

La conclusion exacte est double.

1. Au niveau d'un lien encore fermé au temps $t$, le biais conforme est
   exactement la masse conditionnelle des horloges conformes dans
   $`(t,1]`$. Il est maximal à $`t=\beta_c`$ parmi tous les temps
   postcritiques.
2. Cette propriété scalaire ne suffit pas à contrôler le heat bath au LCA.
   Il faut des majorités **groupées** le long de la chaîne ancestrale et,
   finalement, une inégalité entre les quatre poids du heat bath.

Tous les énoncés probabilistes ci-dessous sont sous la loi annealed de
Nishimori, conditionnellement au squelette de Kruskal non marqué lorsque ce
conditionnement est indiqué. Aucune conclusion de weak recovery n'est
déduite sans le lemme global de domination HF du fichier 12.

## 1. Modèle, horloges et domaine du temps critique

Dans le GSBM homogène binaire, écrivons

```math
O_{xy}=\Sigma_x\Sigma_y Z_{xy},
\qquad
\mathbb P(Z_{xy}=+1)=p,
\qquad
\mathbb P(Z_{xy}=-1)=1-p.
```

Un lien est **conforme** à la ground truth si $`Z_e=+1`$ et **faux** sinon.
Posons

```math
q:=1-p,
\qquad
u_p:=\log\frac pq.
```

Dans la configuration de référence $`\sigma=\Sigma`$, un lien conforme est
satisfait et reçoit une horloge $`\xi_e\sim\mathrm{Exp}(u_p)`$ ; un lien faux
est insatisfait et reçoit $`\xi_e=+\infty`$. La probabilité d'ouverture avant
$t$ vaut donc

```math
q_p(t)=p(1-e^{-u_pt}).
```

Sur la grille triangulaire, notons

```math
q_c=2\sin(\pi/18).
```

Pour tout $`p>1/2`$, l'équation $`q_p(\beta_c)=q_c`$ possède la solution
positive

```math
\boxed{
\beta_c(p)
=
-\frac1{u_p}\log\left(1-\frac{q_c}{p}\right).
}
\tag{1.1}
```

Mais le dendrogramme utilisé par la dynamique est censuré à $t=1$. Ainsi,

```math
\beta_c(p)\le1
\quad\Longleftrightarrow\quad
2p-1\ge q_c
\quad\Longleftrightarrow\quad
p\ge p_{\mathrm{SW}}:=\frac{1+q_c}{2}.
\tag{1.2}
```

Pour $`1/2<p<p_{\mathrm{SW}}`$, le temps critique existe bien sur l'axe
non censuré, mais l'événement « fusion critique dans le même arbre » est vide
pour la dynamique arrêtée à $1$. Toute utilisation du cas favorable dans ce
dossier suppose donc $`p\ge p_{\mathrm{SW}}`$.

## 2. Lemme des quatre catégories postcritiques

Fixons $`0\le t\le1`$ et un lien encore fermé juste avant $t$. Avant
normalisation, il appartient à exactement une des trois catégories
résiduelles suivantes :

| catégorie | masse annealed | conforme ? | compté dans $`\Lambda`$ ? |
|---|---:|:---:|:---:|
| faux, $`Z_e=-1`$ | $`q`$ | non | non |
| vrai tardif, $`Z_e=+1`$, $`t<\xi_e\le1`$ | $`pe^{-u_pt}-q`$ | oui | oui |
| vrai censuré, $`Z_e=+1`$, $`\xi_e>1`$ | $`pe^{-u_p}=q`$ | oui | oui |

La quatrième catégorie, déjà retirée du résidu, est celle des vrais liens
précoces $`\xi_e\le t`$, de masse $`p(1-e^{-u_pt})`$.

La ligne « vrai censuré » est essentielle : la phrase « après $t$, il reste
les faux liens et les vrais liens d'horloge dans $`(t,1]`$ » est exacte
seulement **après annulation des masses symétriques**. Le paquet complet
$`E_v`$ contient aussi les vrais liens censurés après $1$, et
$`\Lambda_v`$ les compte puisqu'ils sont satisfaits.

### Lemme 2.1 — décomposition exacte, statut : établi

Conditionnellement à $`\xi_e>t`$, posons

```math
h_p(t)
:=
\frac{pe^{-u_pt}-q}{q+pe^{-u_pt}}
=
\tanh\left(\frac{u_p(1-t)}2\right),
\tag{2.1}
```

et

```math
s_p(t)
:=
\mathbb P(Z_e=+1\mid\xi_e>t)
=
\frac{pe^{-u_pt}}{q+pe^{-u_pt}}
=
\frac{1+h_p(t)}2.
\tag{2.2}
```

Alors les probabilités conditionnelles des catégories

```math
(\text{vrai tardif},\text{vrai censuré},\text{faux})
```

valent exactement

```math
\boxed{
\left(
h_p(t),
\frac{1-h_p(t)}2,
\frac{1-h_p(t)}2
\right).
}
\tag{2.3}
```

En particulier,

```math
\boxed{
2s_p(t)-1=h_p(t).
}
\tag{2.4}
```

Autrement dit, les vrais liens censurés et les faux liens se compensent
exactement ; tout l'excès de majorité conforme vient des vrais liens dont
l'horloge appartient à $`(t,1]`$.

### Preuve

On utilise $`e^{-u_p}=q/p`$. La masse totale encore fermée est
$`q+pe^{-u_pt}`$, la masse des vrais liens tardifs est
$`p(e^{-u_pt}-e^{-u_p})=pe^{-u_pt}-q`$, et les deux autres masses valent
$q$. La normalisation donne (2.3). La somme des deux catégories conformes
donne (2.2), puis (2.4).

## 3. Spécialisation au seuil de percolation

Supposons désormais $`p\ge p_{\mathrm{SW}}`$. L'identité
$`q_p(\beta_c)=q_c`$ donne

```math
pe^{-u_p\beta_c}=p-q_c.
\tag{3.1}
```

Les quatre masses **non conditionnelles** au temps critique sont donc

```math
\boxed{
\begin{array}{c|c}
\text{catégorie}&\text{masse}\cr
\hline
\text{vrai précoce }(\xi_e\le\beta_c)&q_c\cr
\text{vrai tardif }(\beta_c<\xi_e\le1)&2p-1-q_c\cr
\text{vrai censuré }(\xi_e>1)&1-p\cr
\text{faux}&1-p
\end{array}
}
\tag{3.2}
```

Elles somment à $1$. Parmi les liens encore fermés à $`\beta_c`$,

```math
s_c(p)
=
\frac{p-q_c}{1-q_c},
\qquad
h_c(p)
=
\frac{2p-1-q_c}{1-q_c}
=
\frac{2(p-p_{\mathrm{SW}})}{1-q_c}.
\tag{3.3}
```

### Proposition 3.1 — optimalité scalaire critique, statut : établi

Pour $`t\in[\beta_c,1]`$,

```math
\frac{d}{dt}h_p(t)
=
-\frac{u_p}{2}
\left[1-h_p(t)^2\right]
<0.
\tag{3.4}
```

Ainsi $`h_p(t)\le h_c(p)`$ et $`s_p(t)\le s_c(p)`$. Pour un ensemble fixé
de liens encore fermés, le temps critique fournit donc simultanément :

- la plus grande proportion conditionnelle de liens conformes ;
- la plus grande masse conditionnelle de vrais liens activables avant $1$ ;
- le meilleur rapport signal sur bruit par lien.

C'est le sens rigoureux dans lequel la fusion au seuil est le cas
postcritique le plus favorable. Cette proposition ne compare pas les tailles
aléatoires des coupes $`|E_v|`$ : leur géométrie dépend du nœud sélectionné.

## 4. Deux notions de majorité qu'il ne faut pas confondre

L'idée initiale isole le sous-pool

```math
\{\text{faux}\}
\cup
\{\text{vrais d'horloge dans }(\beta_c,1]\}.
```

Sa fraction conforme vaut

```math
\rho_{\mathrm{late}}(p)
=
\frac{2p-1-q_c}{p-q_c}.
\tag{4.1}
```

Il possède une majorité conforme stricte si et seulement si

```math
\boxed{
p>p_{\mathrm{late}}
:=
\frac{2+q_c}{3}
=0.7824321184\ldots
}
\tag{4.2}
```

Le vrai paquet $`\Lambda_v`$ ne jette toutefois pas les vrais liens censurés.
Dans le résidu complet, sa fraction conforme est $`s_c(p)`$ et elle dépasse
$`1/2`$ si et seulement si

```math
\boxed{p>p_{\mathrm{SW}}.}
\tag{4.3}
```

Le seuil (4.2) est donc un certificat volontairement conservateur : il exige
que les vrais liens tardifs battent à eux seuls tous les faux liens. Le seuil
(4.3) est celui de la majorité brute réellement vue par $`\Lambda`$.

## 5. Loi exacte d'un bucket sélectionné par Kruskal

Fixons le squelette non marqué, un bucket $`E_u`$ de taille $`m\ge1`$ et son
temps de fusion $`\beta_u=t`$. L'identité de l'arête gagnante est
marginalisée. Sous la ground truth, cette gagnante est conforme. Pour les
$`m-1`$ autres arêtes, notons

- $R$ le nombre de vrais liens tardifs $`t<\xi_e\le1`$ ;
- $S$ le nombre de vrais liens censurés $`\xi_e>1`$ ;
- $U$ le nombre de faux liens.

### Théorème 5.1 — loi conditionnelle du vote, statut : établi

```math
\boxed{
(R,S,U)
\sim
\mathrm{Mult}
\left(
m-1;
h_p(t),
\frac{1-h_p(t)}2,
\frac{1-h_p(t)}2
\right).
}
\tag{5.1}
```

Le nombre $K_u$ de liens conformes, donc satisfaits par la ground truth,
vaut

```math
\boxed{
K_u=1+R+S
\stackrel d=
1+\mathrm{Bin}(m-1,s_p(t)).
}
\tag{5.2}
```

Le vote signé $`V_u:=2K_u-m`$ vérifie exactement

```math
\mathbb E[V_u\mid m,t]
=
1+(m-1)h_p(t),
\tag{5.3}
```

```math
\mathrm{Var}(V_u\mid m,t)
=
(m-1)\left[1-h_p(t)^2\right].
\tag{5.4}
```

Enfin,

```math
\boxed{
\mathbb P(V_u>0\mid m,t)
=
\sum_{r=\lfloor m/2\rfloor}^{m-1}
\binom{m-1}{r}
s_p(t)^r[1-s_p(t)]^{m-1-r}.
}
\tag{5.5}
```

Pour $`h_p(t)>0`$ fixé, cette probabilité tend vers $1$ lorsque
$`m\to\infty`$. La fenêtre de transition est gouvernée par

```math
m h_p(t)^2.
\tag{5.6}
```

Près de la censure,

```math
h_p(t)
=
\frac{u_p}{2}(1-t)+O((1-t)^3),
\tag{5.7}
```

donc un vote fortement majoritaire demande $`m(1-t)^2\to\infty`$.

### Audit du conditionnement de Kruskal

La loi (5.1) est exacte **conditionnellement au squelette complet** : la
gagnante est uniforme dans le bucket, elle est conforme, et les autres liens
sont indépendamment conditionnés par $`\xi_e>t`$. En revanche, le choix d'une
paire lointaine dont le LCA est critique biaise la loi de

```math
(m_u,E_u,\beta_u).
```

Il est donc licite d'utiliser (5.1) après avoir conditionné par ces variables,
mais pas de remplacer leur loi sélectionnée par celle d'une coupe uniforme.

## 6. Tous les ancêtres : majorités groupées et quatre taux

Soit $`v\succ u`$. Avec les notations du fichier 08, décomposons

```math
E_v=E_v^{(0)}\mathbin{\dot\cup}E_v^{(1)}\mathbin{\dot\cup}E_v^{(2)},
```

où les groupes $1$ et $2$ touchent respectivement les deux fils $`C_1,C_2`$
de $u$, tandis que le groupe $0$ ne touche aucun d'eux. Dans le cas homogène,
notons

```math
m_{v,r}:=|E_v^{(r)}|,
\qquad
K_{v,r}:=\#\{e\in E_v^{(r)}:Z_e=+1\},
\qquad
M_{v,r}:=2K_{v,r}-m_{v,r}.
```

Conditionnellement au squelette et à la catégorie gagnante
$`G_v\in\{0,1,2\}`$,

```math
\mathbb P(G_v=r\mid\mathscr D)
=
\frac{m_{v,r}}{m_{v,0}+m_{v,1}+m_{v,2}},
\tag{6.1}
```

et, conditionnellement à $`G_v`$,

```math
\boxed{
K_{v,r}
=
\mathbf1_{\{G_v=r\}}
+
\mathrm{Bin}
\left(
m_{v,r}-\mathbf1_{\{G_v=r\}},
s_p(\beta_v)
\right).
}
\tag{6.2}
```

Les quatre taux sous les flips $`(a,b)\in\{0,1\}^2`$ sont

```math
\boxed{
\frac{\Lambda_v^{ab}}{u_p}
=
K_{v,0}
+
\begin{cases}
K_{v,1},&a=0,\cr
m_{v,1}-K_{v,1},&a=1,
\end{cases}
+
\begin{cases}
K_{v,2},&b=0,\cr
m_{v,2}-K_{v,2},&b=1.
\end{cases}
}
\tag{6.3}
```

Par conséquent,

```math
\Lambda_v^{00}-\Lambda_v^{10}=u_pM_{v,1},
\qquad
\Lambda_v^{00}-\Lambda_v^{01}=u_pM_{v,2},
\tag{6.4}
```

```math
\Lambda_v^{00}-\Lambda_v^{11}
=
u_p(M_{v,1}+M_{v,2}).
\tag{6.5}
```

Une majorité sur tout $`E_v`$ ne suffit donc pas. Pour que le taux de la
ground truth $`\Lambda_v^{00}`$ domine séparément les trois concurrents, il
suffit d'avoir des majorités strictes dans **les deux groupes affectés** :

```math
M_{v,1}>0,
\qquad
M_{v,2}>0.
\tag{6.6}
```

Pour chaque groupe, l'échelle de concentration pertinente est

```math
m_{v,r}h_p(\beta_v)^2.
\tag{6.7}
```

Or $`\beta_v>\beta_u\simeq\beta_c`$, donc

```math
h_p(\beta_v)<h_c(p).
\tag{6.8}
```

Attribuer la qualité critique à tous les $`\Lambda_v`$ surestime donc
systématiquement l'information ancestrale.

## 7. De la majorité des taux au heat bath de parité

Les slides 31--33 donnent, pour $`a,b\in\{0,1\}`$,

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
F_v(\Lambda_v^{ab}),
\qquad
F_v(x):=xe^{(1-\beta_v)x}.
\tag{7.1}
```

Sous l'a priori uniforme, la probabilité que le heat bath conserve la parité
conforme entre $`i`$ et $`j`$ est

```math
P_u^{\mathrm{keep}}
=
\frac{q_u^{00}+q_u^{11}}
{q_u^{00}+q_u^{01}+q_u^{10}+q_u^{11}}.
\tag{7.2}
```

Posons

```math
L_u
:=
\log\frac{q_u^{00}+q_u^{11}}
{q_u^{10}+q_u^{01}}.
\tag{7.3}
```

Alors

```math
P_u^{\mathrm{keep}}>\frac12
\quad\Longleftrightarrow\quad
L_u>0,
\tag{7.4}
```

et le biais signé du heat bath vaut $`\tanh(L_u/2)`$.

Lorsque les quatre poids sont strictement positifs, définissons

```math
D_{ab}:=\log\frac{q_u^{00}}{q_u^{ab}}.
```

Le critère exact devient

```math
\boxed{
1+e^{-D_{11}}
>
e^{-D_{10}}+e^{-D_{01}}.
}
\tag{7.5}
```

Une condition suffisante simple est

```math
\min(D_{10},D_{01})>\log2.
\tag{7.6}
```

Les cas où un taux est nul doivent être traités directement dans (7.2), sans
prendre de logarithme.

### Théorème 7.1 — certificat de majorité hiérarchique, statut : établi

Supposons l'a priori i.i.d. uniforme. Dans le cas pondéré, posons

```math
X_{v,r}:=2\lambda_{v,r}-T_{v,r}.
```

Si

```math
2\Lambda_u-T_u\ge0
\tag{7.7}
```

et si, pour tout ancêtre strict $`v\succ u`$,

```math
X_{v,1}\ge0,
\qquad
X_{v,2}\ge0,
\tag{7.8}
```

alors

```math
\boxed{
q_u^{00}+q_u^{11}
\ge
q_u^{10}+q_u^{01}.
}
\tag{7.9}
```

Si la majorité locale (7.7) est stricte, alors (7.9) est stricte. Dans le cas
homogène, (7.8) est exactement $`M_{v,1},M_{v,2}\ge0`$.

### Preuve

Pour $`v\succ u`$, posons

```math
f_v(a,b):=F_v(\Lambda_v^{ab}),
\qquad
F_v(x)=xe^{(1-\beta_v)x}.
```

Comme $`0\le\beta_v\le1`$, la fonction $`F_v`$ est croissante et convexe sur
$`[0,+\infty)`$ :

```math
F_v'(x)
=
e^{(1-\beta_v)x}\left[1+(1-\beta_v)x\right]>0,
```

```math
F_v''(x)
=
(1-\beta_v)e^{(1-\beta_v)x}
\left[2+(1-\beta_v)x\right]
\ge0.
```

Écrivons les quatre coefficients de Walsh de $`f_v`$ :

```math
\widehat f_v(\varnothing)
=
\frac14(f_v^{00}+f_v^{10}+f_v^{01}+f_v^{11}),
```

```math
\widehat f_v(1)
=
\frac14(f_v^{00}+f_v^{01}-f_v^{10}-f_v^{11}),
```

```math
\widehat f_v(2)
=
\frac14(f_v^{00}+f_v^{10}-f_v^{01}-f_v^{11}),
```

```math
\widehat f_v(12)
=
\frac14(f_v^{00}+f_v^{11}-f_v^{10}-f_v^{01}).
\tag{7.10}
```

Le coefficient constant est positif. Sous (7.8), les deux coefficients
linéaires sont non négatifs par croissance. Pour le dernier, les quatre
arguments de $`F_v`$ forment le
rectangle additif

```math
L,\quad L-X_{v,1},\quad L-X_{v,2},
\quad L-X_{v,1}-X_{v,2}.
```

La convexité donne

```math
F_v(L)+F_v(L-X_{v,1}-X_{v,2})
\ge
F_v(L-X_{v,1})+F_v(L-X_{v,2}),
```

donc $`\widehat f_v(12)\ge0`$.

Au nœud $u$, le facteur local prend une valeur sur les états pairs et une
autre sur les états impairs. Sous (7.7), ses seuls coefficients de Walsh non
nuls, $`\widehat f_u(\varnothing)`$ et $`\widehat f_u(12)`$, sont non
négatifs.

Enfin, la transformée de Walsh d'un produit pointwise est la convolution des
transformées. Le cône « quatre coefficients de Walsh non négatifs » est donc
stable par produit. Comme l'a priori uniforme est constant,

```math
\frac14
\left(
q_u^{00}+q_u^{11}-q_u^{10}-q_u^{01}
\right)
\ge0.
```

Si (7.7) est stricte, le coefficient local $`\widehat f_u(12)`$ est strictement
positif et le coefficient constant de chaque autre facteur est positif ; la
convolution donne la stricte positivité finale.

> **Portée.** Le théorème est un certificat suffisant, pas une équivalence.
> Le heat bath peut préférer la parité conforme même si un groupe ancestral
> a une majorité négative, grâce aux autres facteurs. L'a priori uniforme et
> la structure factorisée sont essentiels : pour quatre poids positifs
> arbitraires, la domination séparée de $`q^{00}`$ sur $`q^{10},q^{01}`$ ne
> suffit pas.

### Corollaire 7.2 — probabilité conditionnelle du certificat

Conditionnons par le squelette et par toutes les catégories gagnantes.
Pour $`n\ge0`$, $`g\in\{0,1\}`$ et $`s\in[0,1]`$, définissons les deux
queues binomiales exactes

```math
\mathcal A_{\ge}(n,g,s)
:=
\sum_{b=\lceil(n-g)/2\rceil}^{n}
\binom nb s^b(1-s)^{n-b},
\tag{7.11}
```

et

```math
\mathcal A_{>}(n,g,s)
:=
\sum_{b=\lfloor(n-g)/2\rfloor+1}^{n}
\binom nb s^b(1-s)^{n-b}.
\tag{7.12}
```

Pour le nœud local, posons $`n_u=m_u-1`$. Pour un groupe ancestral,

```math
g_{v,r}:=\mathbf1_{\{G_v=r\}},
\qquad
n_{v,r}:=m_{v,r}-g_{v,r}.
```

L'indépendance conditionnelle des comptes groupés donne exactement la
probabilité du certificat du théorème 7.1 :

```math
\mathcal C_u(\mathscr D,G)
=
\mathcal A_{>}\left(
n_u,1,s_p(\beta_u)
\right)
\prod_{v\succ u}\prod_{r=1}^2
\mathcal A_{\ge}\left(
n_{v,r},g_{v,r},s_p(\beta_v)
\right).
\tag{7.13}
```

Par conséquent,

```math
\boxed{
\mathbb P\left(
P_u^{\mathrm{keep}}>\frac12
\;\middle|\;
\mathscr D,(G_v)_{v\succ u}
\right)
\ge
\mathcal C_u(\mathscr D,G).
}
\tag{7.14}
```

Cette borne est exacte pour l'événement suffisant, mais pas nécessairement
pour l'événement de préférence lui-même.

Pour obtenir un certificat plus lisible, posons aussi

```math
\varepsilon(n,g,h)
:=
\begin{cases}
0,&n=0,\cr
\exp\left[-\dfrac{(nh+g)^2}{2n}\right],&n\ge1.
\end{cases}
\tag{7.15}
```

Alors

```math
\mathbb P\left(
P_u^{\mathrm{keep}}>\frac12
\;\middle|\;
\mathscr D,(G_v)_{v\succ u}
\right)
\ge
\left[
1
-\varepsilon(n_u,1,h_p(\beta_u))
-\sum_{v\succ u}\sum_{r=1}^2
\varepsilon(n_{v,r},g_{v,r},h_p(\beta_v))
\right]_+.
\tag{7.16}
```

Ici $`[x]_+:=\max(x,0)`$.

En effet, un groupe de taille $`m=n+g`$ a une marge

```math
M=g+\sum_{\ell=1}^n X_\ell,
\qquad
\mathbb E[X_\ell]=h,
\qquad
X_\ell\in\{-1,+1\}.
```

Pour $`n\ge1`$, l'inégalité de Hoeffding borne l'événement $`M\le0`$ par
(7.15) ; pour $`n=0`$, chacun des certificats pertinents est déterministe et
valide. Une union bound sur la majorité locale stricte et les deux majorités
ancestrales non négatives, suivie du théorème 7.1, donne (7.16).

Ce corollaire rend le verrou quantitatif explicite. Pour une chaîne de $H$
ancêtres, un régime suffisant est

```math
\min\left\{
(m_u-1)h_p(\beta_u)^2,\,
\min_{\substack{v\succ u\\r\in\{1,2\}}}
m_{v,r}h_p(\beta_v)^2
\right\}
\gg
\log H.
\tag{7.17}
```

Il n'est pas nécessaire : le critère exact (7.5) peut rester favorable
lorsque le certificat de majorité échoue.

### Conséquence exacte au seul nœud $u$

Au nœud de fusion,

```math
\Lambda_u^{00}=\Lambda_u^{11}=u_pK_u,
\qquad
\Lambda_u^{10}=\Lambda_u^{01}=u_p(m-K_u).
```

Comme $`F_u`$ est strictement croissante, si l'a priori et les ancêtres sont
neutralisés, le heat bath préfère la parité conforme si et seulement si
$`K_u>m/2`$. Avec les ancêtres présents, cette équivalence locale disparaît :
seul (7.5) décide.

## 8. Conséquence pour une borne de weak recovery

Les seuils scalaires pertinents s'ordonnent comme suit :

```math
\begin{array}{c|c|c}
\text{seuil}&\text{valeur}&\text{interprétation}\cr
\hline
p_{\mathrm{SW}}&(1+q_c)/2=0.673648\ldots
&\beta_c\le1\text{ et majorité du paquet }\Lambda\cr
p_{\mathrm{late}}&(2+q_c)/3=0.782432\ldots
&\text{vrais tardifs majoritaires face aux faux}\cr
p_{\mathrm{info}}&(1+\sqrt{q_c})/2=0.794659\ldots
&\text{baseline d'impossibilité information--percolation}\cr
p_{\mathrm N}^{(0)}&0.835805\ldots
&\text{calibration conjecturale de face}
\end{array}
\tag{8.1}
```

Deux conclusions rigoureuses en résultent.

1. Le seul test « vrais tardifs contre faux » ne peut pas améliorer la
   baseline $`p_{\mathrm{info}}`$, puisque
   $`p_{\mathrm{late}}<p_{\mathrm{info}}`$. Il devient déjà favorable avant
   la zone encore ouverte.
2. Une majorité conforme locale, même avec probabilité tendant vers $1$,
   n'est pas une preuve de weak recovery. Le score LCA est un majorant de la
   corrélation réelle ; la persistance d'un oracle favorable ne construit ni
   estimateur non oracle ni borne inférieure sur l'overlap.

Pour obtenir une meilleure borne par la dynamique hiérarchique, il faut donc
exploiter davantage que le signe de la majorité : les amplitudes des quatre
poids, les tailles des trois groupes, la dégradation
$`h_p(\beta_v)`$ le long des ancêtres, et la masse des paires qui voient ce
squelette. Le candidat exact à calculer est (7.3), pas seulement (4.1).

Pour une preuve d'impossibilité, la logique favorable reste valide sous HF :
si même l'expérience critique, dans le même arbre et avec toute la chaîne
ancestrale, a un score macroscopique nul, alors les paires moins favorables ne
peuvent pas sauver la weak recovery. Le succès de l'oracle n'implique pas la
réciproque.

## 9. Audit et contre-audit

| affirmation | audit | contre-audit |
|---|---|---|
| La fusion critique donne les meilleurs liens postcritiques. | Vrai à coupe fixée : (3.4). | La taille et la forme de la coupe sélectionnée ne sont pas ordonnées. |
| Le résidu est « faux + vrais tardifs ». | Vrai pour l'excès après compensation. | Le paquet $`\Lambda`$ contient aussi les vrais censurés, de masse exactement égale aux faux. |
| Une majorité conforme dans $`E_u`$ favorise la vérité. | Vrai pour le facteur local isolé. | Les ancêtres modifient les quatre poids ; il faut (7.5). |
| Une majorité dans chaque ancêtre suffit. | Une majorité globale ne suffit pas. | La majorité locale et les deux majorités groupées de chaque ancêtre suffisent sous a priori uniforme, par le théorème 7.1. |
| Tous les ancêtres ont la qualité critique. | Faux. | Leurs temps sont plus grands et leur biais est strictement plus faible : (6.8). |
| Conditionner sur une paire critique donne une borne globale. | Faux sans comparaison. | Le conditionnement biaise le squelette ; HF reste à prouver. |
| Le seuil de majorité est un seuil de weak recovery. | Faux. | Les seuils de majorité sont locaux/oracles et sont déjà sous $`p_{\mathrm{info}}`$. |

## 10. Lemme manquant désormais bien posé

Pour une paire lointaine critique, il faut déterminer la loi asymptotique du
processus marqué

```math
\left(
\beta_v,
m_{v,0},m_{v,1},m_{v,2},
G_v
\right)_{v\succeq u}
```

sous le biais $`|C_{u,1}||C_{u,2}|`$. Conditionnellement à ce processus, les
comptes $`K_{v,r}`$ sont donnés exactement par (6.2), les quatre taux par
(6.3), et le bit gardé par (7.2). Le verrou n'est donc plus la loi des marques,
mais la géométrie du squelette critique et la sommabilité des ancêtres dont

```math
m_{v,r}h_p(\beta_v)^2
```

reste d'ordre un ou plus.

Le premier calcul certifiable doit être fait sur un cactus de triangles,
puis sur une bande triangulaire de largeur fixée. Il doit comparer, sous une
même loi :

1. le critère exact (7.5) ;
2. le certificat groupé du théorème 7.1 ;
3. le test conservateur (4.2).

Le module
[`computations/postcritical_ground_truth_majority.py`](computations/postcritical_ground_truth_majority.py)
vérifie toutes les identités scalaires, la loi binomiale du bucket et les
quatre taux groupés. Son test unitaire effectue les mêmes contrôles par des
formules indépendantes et par énumération finie.
