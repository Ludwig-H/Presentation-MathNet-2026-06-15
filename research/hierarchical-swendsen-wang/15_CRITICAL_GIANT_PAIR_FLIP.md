# Paire lointaine critique : probabilité de parité paire

Cette note répond à la question suivante dans le GSBM homogène binaire :
partant de la ground truth, quelle est la probabilité que le heat bath au
nœud

```math
u:C_u=C_1\mathbin{\dot\cup}C_2
```

choisisse l'un des deux états pairs $`(0,0)`$ ou $`(1,1)`$, lorsque deux
points lointains $`i\in C_1`$ et $`j\in C_2`$ coalescent au temps critique ?

Le verdict mathématique est volontairement séparé en trois niveaux.

1. **Établi, volume fini.** Conditionnellement au graphe ouvert à un temps
   $`\beta`$, la loi des marques des arêtes internes encore fermées est
   explicite et n'est pas modifiée par un conditionnement de composante
   géante mesurable par ce graphe.
2. **Établi, oracle local critique.** Pour une coupe critique de taille $`m`$
   et sans message ancestral, la probabilité moyenne de parité paire est
   exactement $`(1+\Gamma_m^c(p))/2`$. Elle vaut
   $`1/2+1/(2m)`$ au bord $`p_{\mathrm{SW}}`$ et tend vers $`1`$ pour tout
   $`p>p_{\mathrm{SW}}`$ fixé lorsque $`m\to\infty`$. Son déficit possède un
   équivalent complet, avec exposant, préfacteur $`m^{-1/2}`$ et correction
   de parité explicites.
3. **Conditionnel, dynamique hiérarchique complète.** La même limite vaut si
   la coupe du LCA croît et si le message ancestral est sous-linéaire devant
   sa taille. Ni la distance des points ni le mot « géante » n'impliquent à
   eux seuls ces deux propriétés.

La probabilité locale est donc calculée exactement. Un équivalent ne
dépendant que de $`p`$ pour la grille entière reste conditionnel à deux lemmes
géométriques et ancestraux clairement identifiés ci-dessous.

## 1. Modèle et convention critique

Posons

```math
q:=1-p,
\qquad
u_p:=\log\frac pq,
\qquad
A_e(\beta)
:=
\mathbf1_{\{Z_e=+1,\ \xi_e\le\beta\}}.
```

Les $`A_e(\beta)`$ sont i.i.d. de paramètre

```math
a_\beta
:=
p(1-e^{-u_p\beta}).
```

Sur la grille triangulaire, écrivons

```math
q_c:=2\sin(\pi/18),
\qquad
p_{\mathrm{SW}}:=\frac{1+q_c}{2},
```

et, pour $`p\ge p_{\mathrm{SW}}`$,

```math
\beta_c(p)
=
-\frac1{u_p}\log\left(1-\frac{q_c}{p}\right).
```

Il faut distinguer les paramètres : $`q_c`$ est le seuil de percolation en
paramètre d'ouverture, tandis que $`p`$ est la probabilité qu'un lien soit
conforme à la ground truth. La condition pour atteindre la percolation avant
la censure $`\beta=1`$ est $`p\ge p_{\mathrm{SW}}`$, et non simplement
$`p>q_c`$.

Le paramètre de percolation au temps $`\beta_c`$ est toujours exactement
$`q_c`$, même lorsque $`p>p_{\mathrm{SW}}`$. En particulier, sur la grille
plane infinie, il n'existe pas de composante infinie de densité positive au
temps critique. Dans cette note, « composante géante critique » désigne donc
l'un des objets finis suivants, à annoncer explicitement :

- la plus grande composante d'une exhaustion finie ;
- une composante de diamètre macroscopique sélectionnée par une paire
  lointaine ;
- la composante géante de $`\Pi_1`$ dont deux branches ont un LCA dans une
  fenêtre critique.

Comme les temps sont continus,
$`\mathbb P(\beta_{ij}=\beta_c)=0`$ à volume fini. Une affirmation « au
temps critique » signifie une désintégration régulière au temps
$`\beta_c`$, une mesure de Palm du flux de fusions, ou la limite d'une
fenêtre, par exemple

```math
\mathcal F_{L,\rho,\varepsilon}
=
\left\{
d_L(I_L,J_L)\ge\rho L,
\ \beta_c-\varepsilon\le\beta_{I_LJ_L}\le\beta_c
\right\},
```

avec l'ordre $`L\to\infty`$, puis $`\varepsilon\downarrow0`$.

## 2. Arêtes internes d'un cluster à temps fixé

Fixons $`0\le\beta\le1`$ et posons

```math
b_\beta:=pe^{-u_p\beta},
\qquad
1-a_\beta=q+b_\beta.
```

Les quatre catégories non conditionnelles d'une arête sont

| catégorie | masse |
|---|---:|
| vraie ouverte, $`\xi_e\le\beta`$ | $`a_\beta`$ |
| vraie future, $`\beta<\xi_e\le1`$ | $`b_\beta-q`$ |
| vraie censurée, $`\xi_e>1`$ | $`q`$ |
| fausse | $`q`$ |

### Lemme 2.1 — marquage résiduel, statut : établi

Conditionnellement au graphe ouvert complet $`A(\beta)`$, les catégories des
arêtes fermées sont indépendantes et leur loi commune est

```math
\boxed{
(\text{vraie future},\text{vraie censurée},\text{fausse})
\sim
\left(
\frac{b_\beta-q}{q+b_\beta},
\frac q{q+b_\beta},
\frac q{q+b_\beta}
\right).
}
```

Le résultat subsiste après tout conditionnement mesurable par rapport à
$`A(\beta)`$, notamment « le cluster sélectionné est la plus grande
composante » ou « les deux points appartiennent à la même composante ».

#### Preuve

Pour chaque arête, $`A_e(\beta)=0`$ est la réunion disjointe des trois
catégories du membre droit, de masses respectives
$`b_\beta-q,q,q`$. Le conditionnement se factorise arête par arête. Un
événement mesurable par $`A(\beta)`$ ne lit aucune des marques résiduelles.

### Lemme 2.2 — proportions internes, statut : établi

Soit $`C`$ un cluster choisi par une règle mesurable par $`A(\beta)`$. Posons

```math
M_C:=|E(C)|,
\qquad
O_C:=\sum_{e\in E(C)}A_e(\beta),
\qquad
N_C:=M_C-O_C.
```

Si $`F_C`$ est le nombre d'arêtes internes fausses et $`U_C`$ le nombre
d'arêtes internes vraies dont l'horloge dépasse $`\beta`$, alors,
conditionnellement à $`A(\beta)`$,

```math
\boxed{
F_C\sim\mathrm{Bin}\left(N_C,\frac q{1-a_\beta}\right),
\qquad
U_C=N_C-F_C.
}
```

Ainsi, si $`N_C\to\infty`$ et $`O_C/M_C\to\alpha`$ sous n'importe lequel
des conditionnements précédents,

```math
\boxed{
\frac{F_C}{M_C}
\longrightarrow
(1-\alpha)\frac q{1-a_\beta},
\qquad
\frac{U_C}{M_C}
\longrightarrow
(1-\alpha)\frac {b_\beta}{1-a_\beta}.
}
```

Le conditionnement de géante ne change donc pas la proportion
fausse/vraie-tardive **parmi les arêtes fermées**. Il peut changer la
proportion parmi toutes les arêtes internes uniquement par la densité
ouverte $`\alpha`$.

Une borne de Hoeffding donne, lorsque $`N_C>0`$,

```math
\mathbb P\left(
\left|
\frac{F_C}{N_C}-\frac q{1-a_\beta}
\right|>\delta
\ \middle|\ A(\beta)
\right)
\le
2e^{-2\delta^2N_C}.
```

En particulier, pour l'événement $`\mathcal G_{ij}`$ disant que $`i,j`$
appartiennent à la composante sélectionnée au temps $`\beta`$, dès que cet
événement et la règle de sélection ne lisent que $`A(\beta)`$,

```math
\mathcal L(F_C\mid A(\beta),\mathcal G_{ij})
=
\mathrm{Bin}\left(N_C,\frac q{1-a_\beta}\right).
```

L'hypothèse de composante géante n'a donc **aucun effet supplémentaire** sur
la qualité des marques fermées une fois $`A(\beta)`$ connu. Elle a un effet
géométrique, potentiellement majeur, sur la loi de $`(M_C,O_C,N_C)`$, puis
sur la coupe du LCA et sur $`B_u`$.

### Lemme 2.3 — biais exact de connexité, statut : établi

Supposons $`0<a_\beta<1`$. Fixons l'ensemble de sommets $`C`$ avec
$`M:=|E(C)|\ge1`$, et supposons que $`G[C]`$ admette au moins un sous-graphe
connexe. Soit

```math
R_C(a)
:=
\mathbb P_a\bigl(G[C]\text{ est connexe}\bigr).
```

Conditionnellement au fait que $`C`$ soit exactement une composante de
$`\Pi_\beta`$, les arêtes internes sont des Bernoulli $`(a_\beta)`$
conditionnées à rendre $`G[C]`$ connexe et

```math
\boxed{
\mathbb E[O_C\mid C\in\Pi_\beta]
=
a_\beta M
+a_\beta(1-a_\beta)
\frac{d}{da}\log R_C(a)\bigg|_{a=a_\beta}.
}
```

Par conséquent,

```math
\boxed{
\mathbb E\left[\frac{F_C}{M}\middle|C\in\Pi_\beta\right]
=
q\left[
1-\frac{a_\beta}{M}
\frac{d}{da}\log R_C(a)\bigg|_{a=a_\beta}
\right],
}
```

et la même formule vaut pour les vraies non activées en remplaçant $`q`$
par $`b_\beta`$. La dérivée est non négative puisque la connexité est un
événement croissant. Le conditionnement de composante ne peut donc augmenter
ces deux proportions moyennes au-dessus de leurs masses non conditionnelles.

#### Contre-audit

Si $`G[C]`$ est un arbre, sa connexité force ses $`M=|C|-1`$ arêtes à être
ouvertes. Il n'existe alors aucune arête potentielle interne fausse ou
tardive. Aucune proportion interne parmi **toutes** les arêtes ne peut donc
être une fonction universelle de $`(p,\beta)`$.

Ces lemmes sont une calibration. Les arêtes internes à un enfant de $`u`$
n'entrent dans aucun $`\Lambda_v`$ pour $`v\succeq u`$ : elles restent dans
le même enfant et s'annulent sous les flips globaux. Les taux hiérarchiques
utilisent les frontières successives, pas ces réservoirs internes.

## 3. Probabilité paire exacte au nœud de fusion

Pour $`a,b\in\{0,1\}`$, les poids exacts des slides sont

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
```

La probabilité point par point de conserver la relation entre $`i`$ et
$`j`$, c'est-à-dire de choisir $`(0,0)`$ ou $`(1,1)`$, est

```math
\boxed{
P_u^{\mathrm{pair}}
=
\frac{q_u^{00}+q_u^{11}}
{q_u^{00}+q_u^{01}+q_u^{10}+q_u^{11}}
=
\frac1{1+e^{-L_u}},
}
```

où

```math
L_u
:=
\log\frac{q_u^{00}+q_u^{11}}{q_u^{10}+q_u^{01}}.
```

Avec $`m_u:=|E_u|`$, $`K_u`$ le nombre d'arêtes de la coupe satisfaites par
la ground truth et

```math
B_u
:=
\mathrm{LSE}(\Phi_u^{00},\Phi_u^{11})
-
\mathrm{LSE}(\Phi_u^{01},\Phi_u^{10}),
```

la décomposition exacte au temps critique est

```math
\boxed{
L_u
=
B_u
+\ell_{m_u,K_u}^c(p),
}
```

```math
\ell_{m,k}^c(p)
=
\log\frac{k}{m-k}
+a_c(p)(2k-m),
```

avec les conventions usuelles aux coins, et

```math
a_c(p)
=
u_p(1-\beta_c)
=
\log\frac{p-q_c}{1-p}.
```

Conditionnellement au squelette non marqué de la coupe,

```math
\boxed{
K_u
\stackrel d=
1+\mathrm{Bin}(m_u-1,s_c(p)),
\qquad
s_c(p)=\frac{p-q_c}{1-q_c}.
}
```

## 4. Oracle local critique : formule et équivalents

Dans cette section uniquement, posons $`B_u=0`$. Il s'agit d'un sous-modèle
local soluble, et non d'une approximation automatiquement valable pour la
dynamique complète.

Pour une réalisation $`K=k`$, posons

```math
\pi_{m,k}^c(p)
:=
\frac1{1+e^{-\ell_{m,k}^c(p)}}.
```

### Proposition 4.1 — probabilité moyenne exacte, statut : établi

Sous la parité vraie,

```math
\boxed{
\overline P_m^c(p)
:=
\mathbb E\left[\pi_{m,K}^c(p)\right]
=
\frac{1+\Gamma_m^c(p)}2,
}
```

où

```math
\Gamma_m^c(p)
=
\mathbb E\left[
\tanh^2\left(\frac{\ell_{m,K}^c(p)}2\right)
\right]
```

est la fiabilité du fichier 09.

#### Preuve

Soient $`P_+`$ et $`P_-`$ les lois du compte sous les deux parités. Elles
satisfont $`P_-(k)=P_+(m-k)`$. Le biais postérieur

```math
r(k)
=
\frac{P_+(k)-P_-(k)}{P_+(k)+P_-(k)}
=
\tanh\left(\frac{\ell_{m,k}^c(p)}2\right)
```

vérifie $`r(m-k)=-r(k)`$. En groupant $`k`$ et $`m-k`$, on obtient

```math
\mathbb E_{P_+}[r(K)]
=
\mathbb E_{P_+}[r(K)^2]
=
\Gamma_m^c(p).
```

Comme $`\pi_{m,k}^c=(1+r(k))/2`$, la formule suit.

### Corollaire 4.2 — bord géométrique, statut : établi

Au point

```math
p=p_{\mathrm{SW}},
\qquad
\beta_c=1,
\qquad
s_c=\frac12,
\qquad
a_c=0,
```

on a exactement, pour tout $`m\ge1`$,

```math
\boxed{
\overline P_m^c(p_{\mathrm{SW}})
=
\frac12+\frac1{2m}.
}
```

### Théorème 4.3 — grande coupe, statut : établi

Pour tout $`p>p_{\mathrm{SW}}`$ fixé,

```math
\boxed{
\overline P_m^c(p)\longrightarrow1
\qquad(m\to\infty).
}
```

Plus précisément, avec

```math
h_c(p)
:=
2s_c(p)-1
=
\frac{2p-1-q_c}{1-q_c},
```

on a l'équivalent typique, conditionnel au compte,

```math
\boxed{
-\frac1m\log\left(1-\pi_{m,K}^c(p)\right)
\longrightarrow
a_c(p)h_c(p)
\quad\text{en probabilité}.
}
```

Pour la moyenne, un véritable équivalent demande de distinguer les tailles
paires et impaires. Posons

```math
\varepsilon_0:=0,
\qquad
\varepsilon_1:=\frac12,
```

et, pour $`r\in\{0,1\}`$,

```math
\boxed{
C_r(p)
=
\frac1{2s_c(p)\sqrt{2\pi}}
\sum_{j\in\mathbb Z}
\frac1{\cosh(a_c(p)(j+\varepsilon_r))}.
}
```

La série converge absolument dès que $`p>p_{\mathrm{SW}}`$. Alors

```math
\boxed{
1-\overline P_m^c(p)
\sim
\frac{C_{m\bmod 2}(p)}{\sqrt m}
e^{-mI_c(p)}.
}
```

En particulier, l'équivalent logarithmique est

```math
\boxed{
1-\overline P_m^c(p)
=
\exp\left[-mI_c(p)+o(m)\right],
}
```

où

```math
\boxed{
I_c(p)
=
D\left(\frac12\middle\|s_c(p)\right)
=
-\frac12\log\left(1-h_c(p)^2\right)
=
\log\cosh\left(\frac{a_c(p)}2\right).
}
```

Ici, $`D(x\|s)`$ désigne la divergence binaire

```math
D(x\|s)
=
x\log\frac xs
+(1-x)\log\frac{1-x}{1-s}.
```

#### Preuve du taux moyen

On a exactement

```math
1-\overline P_m^c
=
\sum_k
\frac{P_+(k)P_-(k)}{P_+(k)+P_-(k)}.
```

Pour $`x,y\ge0`$,

```math
\frac12\min(x,y)
\le
\frac{xy}{x+y}
\le
\min(x,y).
```

Le déficit possède donc le même exposant que l'erreur de Bayes entre les
deux expériences symétriques. Le rapport de vraisemblance a le signe de
$`2K-m`$ ; cette erreur est, à un facteur de bord près, la queue
$`\mathbb P_{P_+}(K\le m/2)`$. Le théorème de Cramér pour
$`K=1+\mathrm{Bin}(m-1,s_c)`$ donne l'exposant
$`D(1/2\|s_c)`$. Les deux dernières identités suivent de
$`h_c=2s_c-1=\tanh(a_c/2)`$.

#### Preuve du préfacteur

Écrivons $`s=s_c`$, $`a=a_c`$ et

```math
H_m(k)
:=
\frac{P_+(k)P_-(k)}{P_+(k)+P_-(k)}.
```

Si $`m=2n`$, prenons $`A_{2n}:=P_+(n)`$. Pour chaque
$`j\in\mathbb Z`$ fixé,

```math
\frac{H_{2n}(n+j)}{A_{2n}}
\longrightarrow
\frac1{2\cosh(aj)}.
```

En prolongeant $`H_{2n}`$ par zéro hors de son support, l'inégalité
$`H_m\le\min(P_+,P_-)`$, la maximalité du coefficient binomial central et
$`s/(1-s)=e^a`$ donnent un majorant uniforme
$`C(a)e^{-a|j|}`$. Il est sommable sur $`\mathbb Z`$ ; la convergence dominée
s'applique donc à la somme sur $`j`$. Enfin, Stirling donne

```math
A_{2n}
\sim
\frac{e^{-2nI_c}}{s\sqrt{4\pi n}}
=
\frac{e^{-mI_c}}{s\sqrt{2\pi m}}.
```

Si $`m=2n+1`$, prenons $`A_{2n+1}:=P_+(n+1)`$. Pour chaque $`j`$ fixé,

```math
\frac{H_{2n+1}(n+1+j)}{A_{2n+1}}
\longrightarrow
\frac1{2e^{a/2}\cosh(a(j+1/2))},
```

et le même argument de domination s'applique. Cette fois,

```math
A_{2n+1}
\sim
\frac{e^{-mI_c}}{\sqrt{2\pi s(1-s)m}}.
```

Comme $`e^{a/2}=\sqrt{s/(1-s)}`$, les deux calculs donnent exactement les
constantes $`C_0`$ et $`C_1`$ annoncées. Ils montrent aussi pourquoi un
préfacteur unique, indépendant de la parité de $`m`$, serait en général
faux.

### Corollaire 4.4 — fenêtre $`m^{-1/2}`$, statut : établi

Si

```math
p_m
=
p_{\mathrm{SW}}
+\frac{(1-q_c)\alpha}{2\sqrt m},
```

alors

```math
\boxed{
\overline P_m^c(p_m)
\longrightarrow
\frac12\left[
1+\mathbb E\tanh^2(\alpha Z+\alpha^2)
\right],
\qquad
Z\sim\mathcal N(0,1).
}
```

## 5. De la distance des points à la taille de la coupe

Sous la loi de paire critique choisie, posons

```math
M_L:=|E_{u_{I_LJ_L}}|,
\qquad
B_L:=B_{u_{I_LJ_L}}.
```

Dans l'oracle local $`B_L=0`$, le conditionnement par $`M_L=m`$ donne
exactement $`\overline P_m^c(p)`$. Ainsi, pour toute suite déterministe
$`m_L\to\infty`$,

```math
\boxed{
1-
\mathbb P\bigl((a,b)\text{ pair}\mid M_L=m_L, B_L=0\bigr)
\sim
\frac{C_{m_L\bmod2}(p)}{\sqrt{m_L}}
e^{-m_LI_c(p)}.
}
```

Si $`M_L`$ est aléatoire, la formule exacte devient

```math
1-\overline P_L^{\mathrm{loc}}(p)
=
\mathbb E^\star\left[1-\overline P_{M_L}^c(p)\right].
```

Sa valeur dépend alors des petites tailles de coupe et de leur parité. Le
conditionnement « $`i,j`$ sont dans la composante critique sélectionnée »
agit précisément à travers cette loi géométrique ; il ne peut être effacé au
profit de $`p`$ seul.

Les deux propriétés dont on aurait besoin sont les suivantes.

### Hypothèse CUT — croissance de la coupe critique

```math
M_L\longrightarrow\infty
\qquad\text{en probabilité sous la loi de paire critique.}
```

### Hypothèse ANC — absence d'annulation ancestrale extensive

```math
\frac{B_L}{M_L}\longrightarrow0
\qquad\text{en probabilité sous la même loi.}
```

### Théorème 5.1 — limite conditionnelle complète, statut : établi sous CUT et ANC

Pour tout $`p>p_{\mathrm{SW}}`$ fixé, sous CUT et ANC,

```math
\boxed{
P_{u_{I_LJ_L}}^{\mathrm{pair}}
\longrightarrow1
\qquad\text{en probabilité}.
}
```

La probabilité moyenne conditionnelle tend également vers $`1`$.

#### Preuve

Conditionnellement à $`M_L`$ et au squelette, la loi du bucket donne

```math
\frac{2K_L-M_L}{M_L}
\longrightarrow
h_c(p)>0
```

en probabilité dès que CUT tient. Le terme logarithmique divisé par $`M_L`$
tend vers zéro, donc

```math
\frac{\ell_{M_L,K_L}^c(p)}{M_L}
\longrightarrow
a_c(p)h_c(p)>0.
```

ANC donne alors $`L_L/M_L\to a_ch_c>0`$, d'où
$`P_L^{\mathrm{pair}}=(1+e^{-L_L})^{-1}\to1`$. La convergence des moyennes
suit puisque les probabilités sont bornées par $`1`$.

### Variante 5.2 — message ancestral extensif

Plus généralement, si

```math
\frac{B_L}{M_L}\longrightarrow b(p),
```

alors

```math
\frac{L_L}{M_L}
\longrightarrow
a_c(p)h_c(p)+b(p).
```

La probabilité paire tend vers $`1`$ si cette limite est positive, vers
$`0`$ si elle est négative, et aucun équivalent universel ne suit lorsqu'elle
est nulle.

## 6. Pourquoi CUT et ANC ne sont pas automatiques

### Contre-audit 1 — distance $`\not\Rightarrow`$ grande coupe

Deux amas de diamètre macroscopique peuvent être reliés par une unique arête
pivotale. Le cas $`M_L=1`$ est même parfaitement informatif localement : le
taux de parité impaire est nul et $`P_u^{\mathrm{pair}}=1`$. La distance
entre $`i`$ et $`j`$ ne détermine donc ni $`M_L`$ ni son régime asymptotique.

### Contre-audit 2 — la géante ne fixe pas la coupe du LCA

Conditionner les points à appartenir à la plus grande composante biaise la
loi des formes et des interfaces, mais ne fournit aucune concentration de
$`|E_u|`$ sans théorème supplémentaire. Au temps critique, la composante
sélectionnée est en outre macroscopique en diamètre et non géante en densité.

### Contre-audit 3 — $`B_u=0`$ n'est pas une marginalisation

Le message

```math
B_u
=
\mathrm{LSE}(\Phi_u^{00},\Phi_u^{11})
-
\mathrm{LSE}(\Phi_u^{01},\Phi_u^{10})
```

utilise tous les quatre taux $`\Lambda_v^{ab}`$ de chaque ancêtre. Il peut
renforcer ou annuler le message local. Le supprimer ne donne en général ni
un majorant ni un minorant de $`P_u^{\mathrm{pair}}`$.

### Contre-audit 4 — limite en probabilité et moyenne exponentielle

ANC suffit à obtenir $`P_u^{\mathrm{pair}}\to1`$, mais ne suffit pas à
transporter l'exposant $`I_c(p)`$. Des événements ancestraux rares peuvent
dominer $`\mathbb E[1-P_u^{\mathrm{pair}}]`. Un taux moyen exige une borne de
grandes déviations jointe pour $`(K_L,B_L,M_L)`$.

Même dans l'oracle local, le préfacteur précis oscille entre $`C_0(p)`$ et
$`C_1(p)`$. Écrire un équivalent avec une constante unique sans fixer la
parité de $`M_L`$ serait faux ; pour une coupe aléatoire, il faut également
connaître sa loi de parité.

### Contre-audit 5 — probabilité point par point contre moyenne stationnaire

Pour une réalisation initiale fixée,
$`P_u^{\mathrm{pair}}=(1+\tanh(L_u/2))/2`$. Après moyenne sur l'orientation
initiale stationnaire, conditionnellement à une sélection mesurable par le
quotient du heat bath,

```math
\boxed{
\mathbb E[P_u^{\mathrm{pair}}]
=
\frac12\left[
1+
\mathbb E\tanh^2\left(\frac{L_u}{2}\right)
\right].
}
```

Remplacer ponctuellement le biais signé par son carré serait faux.

## 7. Quantité exacte restant à estimer

Pour la loi de Palm de la paire critique, la réponse hiérarchique complète
est exactement

```math
\boxed{
\overline P_{L,\rho,\varepsilon}^{\mathrm{hier}}(p)
=
\frac12\left[
1+
\mathbb E_{L,\rho,\varepsilon}^{\star}
\tanh^2\left(
\frac{\ell_{M_L,K_L}^c(p)+B_L}{2}
\right)
\right].
}
```

Cette formule est fonction de $`p`$ **et** de la loi du squelette ancestral
vu depuis la paire. La réduire à une fonction de $`p`$ seul demande :

1. la limite de la loi de $`M_L`$ ;
2. la limite, ou une domination suffisante, de $`B_L/M_L`$ ;
3. pour un taux exponentiel moyen, une grande déviation jointe ;
4. pour une conclusion de weak recovery, la masse non conditionnelle des
   paires concernées ou le lemme de domination HF du fichier 12.

## 8. Prochaine étape rigoureuse

Le premier objet à mesurer sur tores triangulaires et à calculer exactement
sur cactus et bandes est la loi jointe

```math
\left(
M_L,
\left(
\beta_v,m_{v,0},m_{v,1},m_{v,2}
\right)_{v\succ u},
B_L
\right)
```

sous le biais de paires lointaines de la fenêtre critique. Les tests doivent
enregistrer séparément :

- la probabilité de $`M_L=1`$ et les queues de $`M_L`$ ;
- $`B_L/M_L`$ et les annulations $`B_L\simeq-\ell_{M_L,K_L}^c`$ ;
- la probabilité paire signée et la fiabilité carrée ;
- le poids non conditionnel de l'événement de paire.

Une limite numérique proche de $`1`$ sans ces quatre diagnostics ne constitue
pas un équivalent mathématique.
