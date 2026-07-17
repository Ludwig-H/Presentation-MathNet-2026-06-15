# Probabilités de flip à tous les niveaux et lecture par descendants

Cette note répond à deux questions distinctes.

1. Comment calculer explicitement les probabilités de retournement à une
   racine, une feuille et un nœud interne du dendrogramme ?
2. Peut-on remplacer le message des ancêtres par une lecture des descendants
   ou du chemin entre deux sommets dans l'arbre couvrant de Kruskal ?

Le verdict est le suivant.

- **Racine finale : établi.** Sous a priori binaire uniforme, son orientation
  globale est rééchantillonnée avec probabilités $`1/2,1/2`$.
- **Feuille : établi.** Le noyau exact de la construction est le heat bath
  mono-site de Glauber. Metropolis--Hastings est un noyau alternatif de même
  cible, mais sa probabilité d'acceptation n'est pas la probabilité heat bath.
- **Nœud interne : établi.** Les quatre probabilités sont des fonctions
  explicites de trois coefficients $`(h_1,h_2,J)`$. La probabilité paire et
  les deux probabilités marginales de flip sont données ci-dessous.
- **Descendants : séparation nécessaire.** Leurs facteurs s'annulent dans le
  heat bath d'un nœud fixé. Ils redeviennent pertinents pour un balayage de
  plusieurs nœuds ou pour un heat bath collapsed. Le produit de canaux le
  long du chemin est exact dans un oracle factorisé clairement défini, mais
  pas dans la dynamique hiérarchique complète sans théorème d'indépendance.
- **Arbre couvrant marqué : établi mais trop informatif.** Le produit des
  signes sur le chemin donne exactement la relation des spins de la réplique
  qui a généré les horloges. Révéler les arêtes gagnantes change cependant la
  variable auxiliaire et donc les probabilités de flip.

## 1. Loi conditionnelle commune à tous les niveaux

Fixons le dendrogramme non marqué $`D`$ et posons

```math
\Lambda_v(\sigma)
:=
\sum_{e=\{x,y\}\in E_v}
|W_e|\mathbf1_{\{W_e\sigma_x\sigma_y>0\}},
```

puis

```math
F_v(x):=x e^{(1-\beta_v)x},
```

avec $`F_v(0)=0`$. La cible conditionnelle des heat baths est

```math
\boxed{
\nu_O(\sigma\mid D)
\propto
\mu_0(\sigma)
\prod_{v\in D}F_v(\Lambda_v(\sigma)).
}
```

Chaque probabilité de flip est donc obtenue en restreignant cette mesure à
l'orbite du mouvement considéré. Les trois cas ci-dessous diffèrent par la
taille de cette orbite.

## 2. Orientation globale d'une racine finale

Soit $`R`$ une composante finale de $`\Pi_1`$ et $`\sigma^R`$ la
configuration obtenue en retournant tous les spins de $`R`$.

### Lemme 2.1 — heat bath de racine, statut : établi

On a

```math
\boxed{
\mathbb P(R\text{ retournée}\mid\sigma,D)
=
\frac{\mu_0(\sigma^R)}
{\mu_0(\sigma)+\mu_0(\sigma^R)}.
}
```

Sous l'a priori binaire i.i.d. uniforme,

```math
\boxed{
\mathbb P(R\text{ conservée}\mid\sigma,D)
=
\mathbb P(R\text{ retournée}\mid\sigma,D)
=
\frac12.
}
```

#### Preuve

Pour chaque fusion $`v`$ située dans $`R`$, les deux extrémités de toute
arête de $`E_v`$ sont retournées ensemble. Sa satisfaction, donc
$`\Lambda_v`$, est inchangée. Les facteurs des autres racines sont également
inchangés. Seul le rapport d'a priori subsiste.

### Contre-audit — racine finale contre nœud de fusion supérieur

L'orientation globale de $`R`$ est une orbite à deux états. Le dernier nœud
de fusion $`u:C_1\mathbin{\dot\cup}C_2=R`$ possède, lui, quatre états : les
deux enfants peuvent être retournés séparément. Même sous a priori uniforme,
ces quatre états ne sont généralement pas équiprobables, car le facteur de
la coupe $`E_u`$ distingue les parités paire et impaire.

## 3. Heat bath d'une feuille et Metropolis--Hastings

Soit $`i`$ une feuille. Notons $`\mathcal A(i)`$ la chaîne des nœuds de
fusion qui contiennent $`i`$, et $`\sigma^i`$ la configuration où seul
$`\sigma_i`$ est retourné.

### Lemme 3.1 — rapport mono-site exact, statut : établi

Posons

```math
R_i(\sigma,D)
:=
\frac{\mu_0(\sigma^i)}{\mu_0(\sigma)}
\prod_{v\in\mathcal A(i)}
\frac{F_v(\Lambda_v(\sigma^i))}
{F_v(\Lambda_v(\sigma))}.
```

Lorsque les taux concernés sont positifs, son logarithme vaut

```math
\Delta_i
=
\log\frac{\mu_0(\sigma^i)}{\mu_0(\sigma)}
+
\sum_{v\in\mathcal A(i)}
\left[
\log\frac{\Lambda_v(\sigma^i)}{\Lambda_v(\sigma)}
+(1-\beta_v)
(\Lambda_v(\sigma^i)-\Lambda_v(\sigma))
\right].
```

La probabilité du heat bath de retourner la feuille est

```math
\boxed{
P_i^{\mathrm{HB}}
=
\frac{R_i}{1+R_i}
=
\frac1{1+e^{-\Delta_i}}.
}
```

Avec la proposition déterministe « retourner $`i`$ », l'acceptation de
Metropolis--Hastings est en revanche

```math
\boxed{
P_i^{\mathrm{MH}}
=
\min(1,R_i)
=
\min(1,e^{\Delta_i}).
}
```

Les deux noyaux vérifient la balance détaillée pour la même loi
$`\nu_O(\cdot\mid D)`$, mais

```math
P_i^{\mathrm{HB}}\ne P_i^{\mathrm{MH}}
```

en général. Le premier est le heat bath de Glauber, ou règle de Barker pour
une proposition de flip ; le second est Metropolis--Hastings.

#### Preuve

La restriction de $`\nu_O(\cdot\mid D)`$ à l'orbite
$`\{\sigma,\sigma^i\}`$ possède le rapport de masses $`R_i`$. Sa
normalisation donne la règle logistique. Pour la proposition déterministe de
flip, le rapport des acceptations Metropolis dans les deux directions vaut

```math
\frac{\min(1,R_i)}{\min(1,R_i^{-1})}=R_i,
```

ce qui donne la balance détaillée.

### Spécialisation homogène

Dans le GSBM homogène, $`|W_e|=u_p`$. Pour $`v\in\mathcal A(i)`$, notons
$`m_{v,i}`$ le nombre d'arêtes de $`E_v`$ incidentes à $`i`$ et $`k_{v,i}`$
le nombre de ces arêtes satisfaites. Si

```math
K_v:=\frac{\Lambda_v(\sigma)}{u_p},
```

alors

```math
\boxed{
\frac{\Lambda_v(\sigma^i)}{u_p}
=
K_v+m_{v,i}-2k_{v,i}.
}
```

Ainsi, même à une feuille, le calcul exact lit toute sa chaîne ancestrale.
L'homogénéité remplace des sommes pondérées par des comptes entiers ; elle ne
supprime pas ces comptes.

## 4. Les quatre probabilités d'un nœud interne

Soit

```math
u:C=C_1\mathbin{\dot\cup}C_2.
```

Pour $`a,b\in\{0,1\}`$, retournons $`C_1`$ si $`a=1`$ et $`C_2`$ si
$`b=1`$. Les poids exacts sont

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}F_v(\Lambda_v(\sigma^{ab})),
```

et

```math
\boxed{
p_u^{ab}
:=
\mathbb P((a,b)\mid\sigma,D,u)
=
\frac{q_u^{ab}}
{q_u^{00}+q_u^{01}+q_u^{10}+q_u^{11}}.
}
```

Cette formule avec les poids non normalisés reste la définition sûre lorsque
certains taux sont nuls.

### Lemme 4.1 — paramétrisation de Walsh, statut : établi

Supposons les quatre poids strictement positifs et posons

```math
g_{ab}:=\log q_u^{ab},
\qquad
x:=(-1)^a,
\qquad
y:=(-1)^b.
```

Il existe une constante $`C`$ telle que

```math
g_{ab}=C+h_1x+h_2y+Jxy,
```

avec

```math
h_1
=
\frac{g_{00}+g_{01}-g_{10}-g_{11}}4,
```

```math
h_2
=
\frac{g_{00}+g_{10}-g_{01}-g_{11}}4,
```

```math
J
=
\frac{g_{00}+g_{11}-g_{01}-g_{10}}4.
```

Écrivons

```math
\mathcal Z
:=
e^J\cosh(h_1+h_2)
+e^{-J}\cosh(h_1-h_2).
```

Alors les quatre probabilités sont explicitement

```math
\boxed{
\begin{aligned}
p_u^{00}&=\frac{e^{h_1+h_2+J}}{2\mathcal Z},
&
p_u^{01}&=\frac{e^{h_1-h_2-J}}{2\mathcal Z},\\
p_u^{10}&=\frac{e^{-h_1+h_2-J}}{2\mathcal Z},
&
p_u^{11}&=\frac{e^{-h_1-h_2+J}}{2\mathcal Z}.
\end{aligned}
}
```

On en déduit les marginales

```math
\boxed{
\mathbb P(a=1)
=
\frac{e^{-h_1}\cosh(h_2-J)}{\mathcal Z},
\qquad
\mathbb P(b=1)
=
\frac{e^{-h_2}\cosh(h_1-J)}{\mathcal Z}.
}
```

La probabilité de conserver l'orientation relative des deux fils est

```math
\boxed{
P_u^{\mathrm{pair}}
=
p_u^{00}+p_u^{11}
=
\frac{e^J\cosh(h_1+h_2)}
{e^J\cosh(h_1+h_2)+e^{-J}\cosh(h_1-h_2)}.
}
```

Ses log-odds sont donc

```math
\boxed{
L_u
=
2J
+\log\cosh(h_1+h_2)
-\log\cosh(h_1-h_2).
}
```

La parité et l'orientation commune se séparent en outre exactement :

```math
\boxed{
\begin{aligned}
\mathbb P((0,0)\mid\text{pair})
&=\frac1{1+e^{-2(h_1+h_2)}},
&
\mathbb P((1,1)\mid\text{pair})
&=\frac1{1+e^{2(h_1+h_2)}},\\
\mathbb P((0,1)\mid\text{impair})
&=\frac1{1+e^{-2(h_1-h_2)}},
&
\mathbb P((1,0)\mid\text{impair})
&=\frac1{1+e^{2(h_1-h_2)}}.
\end{aligned}
}
```

Ainsi $`J`$ porte le couplage direct de parité, tandis que $`h_1+h_2`$
partage explicitement les états $`(0,0)`$ et $`(1,1)`$ une fois la parité
paire connue. Les champs influencent aussi les odds de parité par la
correction en différences de $`\log\cosh`$ affichée plus haut.

#### Preuve

Les quatre fonctions $`1,x,y,xy`$ forment la base de Walsh de
$`\{-1,+1\}^2`$. L'inversion donne les trois coefficients affichés et une
constante commune qui disparaît à la normalisation. La somme des états pairs
vaut $`2e^J\cosh(h_1+h_2)`$ ; celle des états impairs vaut
$`2e^{-J}\cosh(h_1-h_2)`$. Les formules marginales s'obtiennent en sommant
respectivement sur $`y`$ et sur $`x`$.

### Lemme 4.2 — contribution locale contre message extérieur

Supposons d'abord $`0<\Lambda_u<T_u`$ ; les cas de bord s'obtiennent dans
les réels étendus, où

```math
T_u:=\sum_{e\in E_u}|W_e|.
```

Posons

```math
\ell_u
:=
\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u).
```

Le facteur du nœud $`u`$ ajoute exactement $`\ell_u/2`$ au couplage
$`J`$ et n'ajoute rien aux champs $`h_1,h_2`$. L'a priori et les ancêtres
fournissent $`h_1,h_2,J_{\mathrm{ext}}`$, d'où

```math
J=J_{\mathrm{ext}}+\frac{\ell_u}{2}.
```

Cette formule est équivalente à $`L_u=B_u+\ell_u`$ du fichier 08, mais elle
donne en plus chacune des quatre probabilités et les deux marginales de flip.

#### Preuve

Le facteur local prend la même valeur $`F_u(\Lambda_u)`$ dans les états
pairs et $`F_u(T_u-\Lambda_u)`$ dans les états impairs. Ses deux
coefficients de champ sont donc nuls. Son coefficient de Walsh $`J`$ vaut
la moitié du log-rapport pair/impair, soit

```math
\frac12\log\frac{F_u(\Lambda_u)}{F_u(T_u-\Lambda_u)}
=
\frac{\ell_u}{2}.
```

### Corollaire 4.3 — nœud supérieur sans message extérieur

Sous a priori uniforme, si $`u`$ n'a pas d'ancêtre strict, alors

```math
h_1=h_2=J_{\mathrm{ext}}=0.
```

Par conséquent,

```math
\boxed{
p_u^{00}=p_u^{11}=\frac12\frac1{1+e^{-\ell_u}},
}
```

```math
\boxed{
p_u^{01}=p_u^{10}=\frac12\frac1{1+e^{\ell_u}}.
}
```

Il ne faut donc pas confondre ce nœud supérieur à quatre états avec le flip
global de sa racine, qui reste uniforme à deux états.

## 5. Bucket homogène à un niveau arbitraire

Cette section donne la formule locale explicite pour tous les descendants,
pas seulement au seuil critique.

Comme les horloges sont continues, le conditionnement par un temps déterminé
$`t`$ désigne la désintégration de la densité de fusion ; les formules à
$`t=1`$ sont leurs limites au bord.

On se place sous le couplage de Nishimori qui identifie la réplique générant
les horloges à la vérité, et l'on conditionne par la filtration juste avant
la fusion, les deux composantes et leur coupe $`E_u`$. Les arêtes de cette
coupe ont alors survécu jusqu'à $`t`$ ; l'indépendance résiduelle des horloges
donne l'expérience ci-dessous.

Posons

```math
u_p:=\log\frac p{1-p},
\qquad
a_p(t):=u_p(1-t),
```

```math
s_p(t)
:=
\frac{pe^{-u_pt}}{1-p+pe^{-u_pt}}
=
\frac1{1+e^{-a_p(t)}},
```

```math
h_p(t):=2s_p(t)-1=\tanh\frac{a_p(t)}2.
```

### Lemme 5.1 — probabilité locale exacte, statut : établi

Conditionnellement à un bucket non marqué de taille
$`m:=|E_u|`$ fusionnant au temps $`t`$, le nombre $`K`$ d'arêtes satisfaites
sous la parité vraie suit

```math
\boxed{
K\stackrel d=1+\mathrm{Bin}(m-1,s_p(t)).
}
```

Conditionnellement à $`K=k`$ et après neutralisation du message extérieur,

```math
\boxed{
\ell_{m,k}(t;p)
=
\log\frac{k}{m-k}
+a_p(t)(2k-m),
}
```

avec les conventions $`\ell_{m,0}=-\infty`$ et
$`\ell_{m,m}=+\infty`$. Ainsi

```math
\boxed{
\pi_{m,k}(t;p)
:=
\mathbb P((a,b)\text{ pair}\mid K=k,B=0)
=
\frac1{1+e^{-\ell_{m,k}(t;p)}}.
}
```

La correction $`\log(k/(m-k))`$ vient de l'arête gagnante satisfaite dont
l'identité a été marginalisée.

#### Preuve

Pour une arête,

```math
\mathbb P(T_e>t)=1-p+pe^{-u_pt}.
```

Conditionnellement à $`T_e>t`$, sa probabilité d'être satisfaite est
$`s_p(t)`$. Dans un bucket de taille $`m`$, l'arête qui atteint le minimum en
$`t`$ est uniforme et satisfaite ; les $`m-1`$ autres arêtes sont
indépendantes conditionnellement à leur survie au-delà de $`t`$. Cela donne
la loi de $`K`$.

Sous les deux parités, les lois du compte sont

```math
P_+(k)
=
\binom{m-1}{k-1}s^{k-1}(1-s)^{m-k},
```

```math
P_-(k)=P_+(m-k).
```

Pour $`1\le k\le m-1`$,

```math
\log\frac{P_+(k)}{P_-(k)}
=
\log\frac{k}{m-k}
+(2k-m)\log\frac{s}{1-s}.
```

Comme $`\log(s_p(t)/(1-s_p(t)))=a_p(t)`$, ce rapport est exactement
$`\ell_{m,k}(t;p)`$.

### Lemme 5.2 — moyenne sous la vérité, statut : établi

Définissons

```math
\Gamma_m(t;p)
:=
\mathbb E\left[
\tanh^2\left(\frac{\ell_{m,K}(t;p)}2\right)
\right].
```

Alors

```math
\boxed{
\overline P_m(t;p)
:=
\mathbb E[\pi_{m,K}(t;p)]
=
\frac{1+\Gamma_m(t;p)}2.
}
```

#### Preuve de l'identité moyenne

Le biais postérieur

```math
r(k)
:=
\frac{P_+(k)-P_-(k)}{P_+(k)+P_-(k)}
=
\tanh\frac{\ell_{m,k}(t;p)}2
```

est antisymétrique sous $`k\leftrightarrow m-k`$. En regroupant les deux
valeurs,

```math
\mathbb E_{P_+}[r(K)]
=
\mathbb E_{P_+}[r(K)^2]
=
\Gamma_m(t;p).
```

Or $`\pi_{m,k}=(1+r(k))/2`$, ce qui prouve la formule.

Au niveau terminal $`t=1`$,

```math
\boxed{
\Gamma_m(1;p)=\frac1m,
\qquad
\overline P_m(1;p)=\frac12+\frac1{2m}.
}
```

À $`t=1`$, $`a_p(1)=0`$ et

```math
\tanh\left(\frac12\log\frac K{m-K}\right)
=
\frac{2K-m}{m}.
```

Avec $`K=1+\mathrm{Bin}(m-1,1/2)`$, le numérateur a un second moment égal à
$`m`$, d'où $`\Gamma_m(1;p)=1/m`$.

Pour tout $`t<1`$ fixé,

```math
\overline P_m(t;p)\longrightarrow1
\qquad(m\to\infty).
```

Plus précisément, avec

```math
I(t;p)
:=
D\left(\frac12\middle\|s_p(t)\right)
=
-\frac12\log(1-h_p(t)^2)
=
\log\cosh\frac{a_p(t)}2,
```

le même calcul de point selle que dans le fichier 15 donne

```math
1-\overline P_m(t;p)
\sim
\frac{C_{m\bmod2}(t;p)}{\sqrt m}
e^{-mI(t;p)},
```

où, pour $`r\in\{0,1\}`$, $`\varepsilon_0=0`$,
$`\varepsilon_1=1/2`$ et

```math
C_r(t;p)
=
\frac1{2s_p(t)\sqrt{2\pi}}
\sum_{j\in\mathbb Z}
\frac1{\cosh(a_p(t)(j+\varepsilon_r))}.
```

La spécialisation $`t=\beta_c(p)`$ redonne exactement le fichier 15.

L'équivalent précis suit mot pour mot de la preuve du fichier 15 après les
substitutions $`s_c\mapsto s_p(t)`$ et $`a_c\mapsto a_p(t)`$ : le déficit
est la somme harmonique

```math
\sum_k\frac{P_+(k)P_-(k)}{P_+(k)+P_-(k)},
```

les valeurs $`k-m/2=O(1)`$ donnent le préfacteur en série, et les queues sont
dominées géométriquement.

### Corollaire 5.3 — probabilités moyennes des quatre états

Dans l'oracle local sans message extérieur, les champs sont nuls. Sous la
parité vraie,

```math
\boxed{
\begin{aligned}
\mathbb E[p_u^{00}]
&=
\mathbb E[p_u^{11}]
=
\frac{1+\Gamma_m(t;p)}4,\\
\mathbb E[p_u^{01}]
&=
\mathbb E[p_u^{10}]
=
\frac{1-\Gamma_m(t;p)}4.
\end{aligned}
}
```

En particulier, au bord terminal,

```math
\mathbb E[p_u^{00}]
=
\mathbb E[p_u^{11}]
=
\frac14+\frac1{4m},
```

tandis que, pour tout $`t<1`$ fixé,

```math
\mathbb E[p_u^{00}],\mathbb E[p_u^{11}]\longrightarrow\frac12,
\qquad
\mathbb E[p_u^{01}],\mathbb E[p_u^{10}]\longrightarrow0.
```

Le déficit total des deux états conformes possède l'équivalent du lemme 5.2.
Avec un message ancestral, les égalités $`00=11`$ et $`01=10`$ peuvent
disparaître ; les quatre formules du lemme 4.1 restent alors la réponse
exacte.

#### Preuve

Lorsque $`h_1=h_2=0`$, les deux états de chaque parité ont le même poids.
Le lemme 5.2 donne la masse moyenne de la parité paire ; la division par
$`2`$ donne les quatre identités.

### Audit de l'intuition « les descendants sont de meilleure qualité »

Si $`0\le t_1<t_2<1`$, alors

```math
s_p(t_1)>s_p(t_2),
\qquad
a_p(t_1)>a_p(t_2),
\qquad
I(t_1;p)>I(t_2;p).
```

La probabilité d'une majorité vraie stricte dans le bucket est donc plus
grande au niveau $`t_1`$ par couplage monotone des binomiales. Pour une grande
coupe, l'erreur moyenne possède également un meilleur exposant au niveau le
plus précoce.

Point par point, il faut cependant conserver le signe : pour $`k>m/2`$, une
diminution de $`t`$ augmente $`\pi_{m,k}`$, tandis que pour $`k<m/2`$ elle
la diminue, car elle rend le mauvais verdict plus confiant. La monotonie
utile est une affirmation sous la loi vraie ou sur l'exposant, pas une
inégalité uniforme en $`k`$.

## 6. Pourquoi les descendants disparaissent du heat bath de $`u`$

### Lemme 6.1 — annulation exacte, statut : établi

Soit $`w`$ un descendant strict de $`u`$. Toutes les arêtes de $`E_w`$ ont
leurs deux extrémités dans le même enfant $`C_1`$ ou $`C_2`$ de $`u`$.
Pour tout $`a,b\in\{0,1\}`$,

```math
\boxed{
\Lambda_w(\sigma^{ab})=\Lambda_w(\sigma).
}
```

Le facteur $`F_w(\Lambda_w)`$ est donc commun aux quatre poids et disparaît
de leur normalisation. Les branches disjointes de $`u`$ s'annulent de la même
façon. C'est pourquoi le heat bath exact de $`u`$ ne contient que $`u`$, ses
ancêtres et l'a priori.

#### Preuve

Un flip global de $`C_1`$ ou $`C_2`$ multiplie par $`-1`$ les deux spins de
toute arête interne à cet enfant. Leur produit ne change pas. C'est le cas de
toutes les arêtes de chaque coupe descendante $`E_w`$.

Les descendants restent présents indirectement dans le motif de spins figé
à l'intérieur de $`C_1,C_2`$, lequel détermine les comptes satisfaits sur les
frontières de $`u`$ et de ses ancêtres. Leurs seuls temps de fusion ne
déterminent pas ces comptes.

### Contre-exemple minimal aux seuls niveaux

Prenons deux sommets $`x_1,x_2\in C_1`$ et $`y\in C_2`$, avec deux arêtes de
même poids entre $`C_1`$ et $`C_2`$. Les mêmes niveaux de fusion internes à
$`C_1`$ sont compatibles avec deux arêtes de frontière satisfaites, une
seule, ou aucune selon leurs signes et le motif interne. On obtient alors
des valeurs différentes de $`\Lambda_u`$ et de $`p_u^{ab}`$ avec les mêmes
niveaux descendants. L'égalité des poids ne supprime donc pas l'information
de frontière.

## 7. Identité exacte le long d'un balayage descendant

La lecture par chemin devient pertinente si l'on applique plusieurs
mouvements, et non un seul heat bath au LCA.

Considérons une suite de flips de clusters $`S_1,\ldots,S_N`$. Soit
$`A_r\in\{0,1\}`$ l'indicatrice que $`S_r`$ est effectivement retourné et

```math
\chi_r(i,j)
:=
\mathbf1_{\{i\in S_r\}}
\mathbin{\oplus}
\mathbf1_{\{j\in S_r\}}.
```

Un tirage à quatre états d'un nœud interne est représenté par les deux
clusters enfants placés consécutivement dans cette liste. Leurs indicatrices
gardent leur loi jointe ; aucune indépendance n'est introduite, et l'ordre
est sans effet puisque les deux flips commutent.

### Lemme 7.1 — parité de chemin, statut : établi

Après le balayage,

```math
\boxed{
\frac{\sigma_i'\sigma_j'}{\sigma_i\sigma_j}
=
(-1)^{\sum_{r=1}^N A_r\chi_r(i,j)}
=
\prod_{r:\chi_r(i,j)=1}(-1)^{A_r}.
}
```

Par conséquent,

```math
\boxed{
\mathbb P(\text{relation conservée}\mid O,D)
=
\frac12\left[
1+
\mathbb E\left[
\prod_{r:\chi_r(i,j)=1}(-1)^{A_r}
\middle|O,D
\right]
\right].
}
```

#### Preuve

Le flip de $`S_r`$ multiplie $`\sigma_i\sigma_j`$ par $`-1`$ si et seulement
si $`S_r`$ contient exactement un des deux sommets. La multiplication de ces
effets élémentaires donne la première identité ; l'indicatrice de l'événement
« produit final positif » vaut $`(1+\prod_r(-1)^{A_r\chi_r})/2`$.

Dans l'arbre de fusions, les clusters qui contiennent exactement un des deux
sommets se trouvent sur les deux bras allant de $`i`$ et $`j`$ à leur LCA.
Les ancêtres stricts du LCA contiennent les deux sommets et ne figurent pas
dans le produit déterministe.

Si la réplique initiale est identifiée à la ground truth, cette formule est
exactement la probabilité que $`i`$ et $`j`$ soient soit tous deux conformes,
soit tous deux inversés après le balayage.

### Lemme 7.2 — quatre statuts après recoloration globale équitable

Supposons $`i,j`$ dans la même racine et terminons par sa recoloration
uniforme. Posons

```math
R_i:=\sigma_i'\Sigma_i,
\qquad
R_j:=\sigma_j'\Sigma_j,
```

et

```math
c_{ij}^{\mathrm{path}}
:=
\mathbb E[R_iR_j].
```

La symétrie globale impose alors

```math
\boxed{
\begin{aligned}
\mathbb P(R_i=+1,R_j=+1)
&=
\mathbb P(R_i=-1,R_j=-1)
=
\frac{1+c_{ij}^{\mathrm{path}}}{4},\\
\mathbb P(R_i=+1,R_j=-1)
&=
\mathbb P(R_i=-1,R_j=+1)
=
\frac{1-c_{ij}^{\mathrm{path}}}{4}.
\end{aligned}
}
```

Ainsi, calculer la corrélation signée du produit de chemin donne chacune des
quatre probabilités de conformité/inversion après moyenne sur l'orientation
globale. Sans cette moyenne, les deux probabilités à l'intérieur d'une même
parité ne sont pas nécessairement égales ; les champs $`h_1,h_2`$ du lemme
4.1 les distinguent. La recoloration équitable de racine détruit
l'orientation absolue, mais pas la relation de paire pertinente modulo le
flip global pour la weak recovery.

#### Preuve

La recoloration uniforme rend la loi de $`(R_i,R_j)`$ invariante sous
$`(R_i,R_j)\mapsto(-R_i,-R_j)`$. Les deux marginales signées sont donc
nulles. La transformée de Walsh d'une loi sur deux signes donne alors

```math
\mathbb P(R_i=r,R_j=s)
=
\frac14(1+rs\,c_{ij}^{\mathrm{path}}),
```

pour $`r,s\in\{-1,+1\}`$.

### Corollaire 7.3 — entrée exacte dans l'obstruction de weak recovery

Pour un même balayage hiérarchique $`S`$ appliqué à toutes les paires,
posons $`\zeta_x=\sigma_x\sigma_x'`$. Alors

```math
\boxed{
H_S(i,j)
:=
\mathbb E[\zeta_i\zeta_j\mid O,\sigma,D]
=
\mathbb E\left[
\prod_{r:\chi_r(i,j)=1}(-1)^{A_r}
\middle|O,\sigma,D
\right].
}
```

En effet,
$`\zeta_i\zeta_j=(\sigma_i'\sigma_j')/(\sigma_i\sigma_j)`$, puis le lemme
7.1 s'applique. La matrice $`H_S`$ est celle du théorème d'obstruction du
fichier 03. La récursion exacte de la proposition 10.2 calcule donc ses
entrées lorsque l'état de frontière est contrôlable.

Un produit PATH-FAC petit suggère une meilleure borne d'impossibilité que la
seule taille de composante, mais ne peut être inséré à la place de $`H_S`$
sans une domination matricielle ou une borne directe sur
$`\lambda_{\max}(H_S)`$.

### Contre-audit — les ancêtres restent dans la loi

Bien qu'ils ne figurent pas dans le produit de signes, les ancêtres entrent
dans la probabilité conditionnelle de chaque $`A_r`$. Les heat baths de
nœuds emboîtés sont dépendants et une mise à jour change les taux utilisés
par les suivantes. En général,

```math
\mathbb E\left[\prod_r(-1)^{A_r}\middle|O,D\right]
\ne
\prod_r\mathbb E\left[(-1)^{A_r}\middle|O,D\right].
```

L'indépendance conditionnelle des **marques de buckets** ne donne pas
l'indépendance des **décisions de heat bath**.

Le contre-exemple abstrait minimal est déjà décisif : prenons $`A_1`$
équitable et $`A_2=A_1`$. Alors

```math
\mathbb E[(-1)^{A_1}]
=
\mathbb E[(-1)^{A_2}]
=0,
```

mais

```math
\mathbb E[(-1)^{A_1+A_2}]=1.
```

Le produit des biais marginaux prédit donc $`0`$ alors que la parité jointe
est parfaite. Le test numérique associé est inclus dans
`test_hierarchical_flip_probabilities.py`.

## 8. Oracle factorisé du chemin

Il est néanmoins utile de définir un sous-modèle soluble qui formalise
l'intuition descendante.

### Définition 8.1 — oracle PATH-FAC

Pour chaque nœud $`w`$ des deux bras entre $`i,j`$ et leur LCA :

1. conserver seulement son bucket non marqué $`(m_w,t_w,K_w)`$ ;
2. poser son message extérieur à zéro ;
3. tirer sa parité avec le heat bath local ;
4. rendre ces parités indépendantes entre nœuds.

Ce modèle est une approximation factorisée ou un oracle de calcul. Il n'est
ni un majorant ni un minorant automatique de la vraie dynamique.

### Proposition 8.2 — formule produit, statut : établi dans PATH-FAC

Sous la vérité, le biais signé moyen du canal au nœud $`w`$ vaut
$`\Gamma_{m_w}(t_w;p)`$. Par indépendance,

```math
\boxed{
P_{ij}^{\mathrm{PATH-FAC}}
=
\frac12\left[
1+
\prod_{w\in\mathcal P(i,j)}
\Gamma_{m_w}(t_w;p)
\right].
}
```

Après recoloration globale équitable, PATH-FAC donne donc explicitement

```math
\mathbb P(\text{tous deux conformes})
=
\mathbb P(\text{tous deux inversés})
=
\frac14\left[
1+\prod_{w\in\mathcal P(i,j)}\Gamma_{m_w}(t_w;p)
\right].
```

#### Preuve

Sous la parité vraie, le signe moyen d'un heat bath local vaut
$`2\overline P_m(t;p)-1=\Gamma_m(t;p)`$. L'indépendance imposée dans
PATH-FAC transforme l'espérance du produit de signes en produit de leurs
espérances. Le lemme 7.2 donne ensuite les quatre statuts.

Ainsi, des niveaux individuellement bons ne suffisent pas : leur nombre et
la sommabilité de leurs pertes comptent.

### Corollaire 8.3 — critère de produit infini

Pour une suite infinie de canaux avec $`0<\Gamma_w\le1`$,

```math
\prod_w\Gamma_w>0
\quad\Longleftrightarrow\quad
\sum_w(1-\Gamma_w)<\infty.
```

Cette équivalence est le critère classique des produits infinis : pour
$`x`$ proche de $`1`$, $`-\log x`$ est comparable à $`1-x`$ ; les facteurs
qui restent loin de $`1`$ ne peuvent être qu'en nombre fini lorsque le
produit est positif.

Pour une famille de chemins $`\mathcal P_L`$ :

- si $`\sum_{w\in\mathcal P_L}(1-\Gamma_w)\to0`$, alors
  $`P_{ij}^{\mathrm{PATH-FAC}}\to1`$ ;
- s'il existe $`\varepsilon,c>0`$ tels qu'au moins
  $`c|\mathcal P_L|`$ nœuds vérifient $`\Gamma_w\le1-\varepsilon`$, alors
  $`P_{ij}^{\mathrm{PATH-FAC}}\to1/2`$.

### Corollaire 8.4 — tailles de coupe suffisantes

Fixons $`\beta<1`$ et supposons $`t_w\le\beta`$. Posons

```math
a_\beta:=u_p(1-\beta),
\qquad
h_\beta:=\tanh(a_\beta/2),
```

```math
c_\beta
:=
\min\left(
\frac{h_\beta^2}{8},
\frac{a_\beta h_\beta}{2}
\right)>0.
```

L'argument de Hoeffding du fichier 09 donne uniformément

```math
1-\Gamma_m(t;p)
\le
5e^{-(m-1)c_\beta},
\qquad
t\le\beta.
```

Si $`H_L:=|\mathcal P_L|\to\infty`$ et

```math
\min_{w\in\mathcal P_L}m_w
\ge
\left(\frac1{c_\beta}+\varepsilon\right)\log H_L
```

pour un $`\varepsilon>0`$, alors PATH-FAC donne
$`P_{ij}^{\mathrm{PATH-FAC}}\to1`$.

En effet,

```math
\sum_{w\in\mathcal P_L}(1-\Gamma_w)
\le
5H_L\exp[-(m_{\min,L}-1)c_\beta]
\le
5e^{c_\beta}H_L^{-\varepsilon c_\beta}
\longrightarrow0.
```

Inversement, pour $`p<1`$ et un entier $`M`$ fixé, il existe
$`\delta(p,\beta,M)>0`$ tel que

```math
\Gamma_m(t;p)\le1-\delta(p,\beta,M)
```

pour $`0\le t\le\beta`$ et $`2\le m\le M`$. Une proportion positive de
telles petites coupes force donc la limite $`1/2`$ dans PATH-FAC. Les coupes
$`m=1`$ sont exclues : l'arête gagnante y révèle parfaitement la parité
locale et $`\Gamma_1=1`$.

Pour justifier l'inverse, pour chaque $`2\le m\le M`$ la fonction
$`t\mapsto\Gamma_m(t;p)`$ est continue sur le compact $`[0,\beta]`$ et
strictement inférieure à $`1`$ : une valeur intérieure de $`K`$ a une masse
strictement positive et un log-odds fini. Le maximum sur la réunion finie de
ces compacts est donc strictement inférieur à $`1`$.

### Corollaire 8.5 — spécialisation au GSBM triangulaire critique

Pour la percolation de liens sur la grille triangulaire, posons

```math
q_\triangle=2\sin(\pi/18),
\qquad
p_{\mathrm{SW}}=\frac{1+q_\triangle}{2}.
```

Si $`p>p_{\mathrm{SW}}`$, le niveau critique accessible par les horloges est

```math
\beta_c(p)
=
-\frac1{u_p}\log\left(1-\frac{q_\triangle}{p}\right)<1.
```

Les buckets des deux bras sous le LCA critique ont presque sûrement des
niveaux $`t_w<\beta_c(p)`$. La borne précédente s'applique donc avec

```math
a_c(p)
=
\log\frac{p-q_\triangle}{1-p},
\qquad
h_c(p)
=
\frac{2p-1-q_\triangle}{1-q_\triangle},
```

```math
c_c(p)
:=
\min\left(\frac{h_c(p)^2}{8},
\frac{a_c(p)h_c(p)}2\right).
```

Ainsi, dans PATH-FAC, un chemin de longueur $`H_L`$ dont toutes les tailles
de coupe sont au moins

```math
\left(\frac1{c_c(p)}+\varepsilon\right)\log H_L
```

conserve asymptotiquement la relation vraie avec probabilité $`1-o(1)`$.
Une densité positive de coupes de tailles comprises entre $`2`$ et une
constante fixe donne au contraire la limite $`1/2`$. Ce résultat est exact
pour PATH-FAC ; le transporter à la dynamique hiérarchique sur la grille
demande toujours une comparaison des décisions jointes.

## 9. Le chemin dans l'arbre couvrant marqué

Soit $`T`$ la forêt couvrante minimale **marquée**, qui conserve l'identité
de chaque arête gagnante. Pour $`i,j`$ dans le même arbre, notons
$`P_T(i,j)`$ son chemin physique.

### Lemme 9.1 — identité télescopique, statut : établi

Chaque arête $`e=\{x,y\}`$ de $`T`$ possède une horloge finie seulement si
elle est satisfaite par la réplique $`\sigma`$ ayant généré les horloges.
Ainsi

```math
\mathrm{sign}(W_e)\sigma_x\sigma_y=1.
```

En multipliant le long du chemin,

```math
\boxed{
\sigma_i\sigma_j
=
\prod_{e\in P_T(i,j)}\mathrm{sign}(W_e).
}
```

Le principe minimax de Kruskal donne aussi

```math
\max_{e\in P_T(i,j)}\xi_e=\beta_{ij};
```

tous les autres niveaux du chemin sont inférieurs à $`\beta_{ij}`$ presque
sûrement.

Si l'on identifie la réplique initiale à la ground truth par Nishimori, le
chemin marqué révèle donc exactement $`\Sigma_i\Sigma_j`$.

#### Preuve

Chaque arête retenue par Kruskal a une horloge finie ; par construction des
horloges, elle est donc satisfaite par $`\sigma`$. Sur un chemin
$`i=x_0,x_1,\ldots,x_r=j`$, la multiplication des identités

```math
\mathrm{sign}(W_{x_{k-1}x_k})
=
\sigma_{x_{k-1}}\sigma_{x_k}
```

fait apparaître chaque spin intérieur deux fois et laisse
$`\sigma_i\sigma_j`$. Enfin, la propriété minimax d'un arbre couvrant
minimal dit que le plus grand poids du chemin entre $`i`$ et $`j`$ est le
plus petit seuil auquel ils sont connectés. C'est précisément
$`\beta_{ij}`$. La continuité des horloges rend ce maximum unique presque
sûrement.

### Contre-audit décisif — ce n'est pas le dendrogramme utilisé

Le dossier emploie le dendrogramme de partitions **non marqué**. Il oublie
l'identité de l'arête gagnante et conserve seulement la somme de taux
$`\Lambda_u`$. Si l'arête gagnante $`e_u`$ est révélée, la densité contient
son indicatrice de satisfaction. Tout état qui la rend insatisfaite reçoit
un poids nul. Sur un chemin marqué, les relations deviennent donc des
contraintes dures et les flips internes qui les cassent sont interdits.

Révéler $`T`$ produit un oracle parfait mais une autre variable auxiliaire et
une autre dynamique. Cela ne permet ni de calculer les $`q_u^{ab}`$ non
marqués, ni de conclure à la weak recovery à partir de l'observation $`O`$
seule.

L'égalité des poids $`|W_e|=u_p`$ rend l'arête gagnante uniforme parmi les
candidates d'une coupe. Elle ne révèle pas son identité et ne transforme pas
les niveaux $`\xi_e`$ en signes observables.

## 10. Une vraie récursion par descendants : message de frontière

Une représentation descendante exacte existe après marginalisation des
configurations internes, mais son état n'est généralement pas scalaire.

Pour un cluster $`C`$, notons

```math
\partial_E C
:=
\{x\in C:\exists y\notin C,\ \{x,y\}\in E\}.
```

Après sommation des spins de $`C\setminus\partial_E C`$, le message exact du
sous-arbre enraciné en $`C`$ est une fonction de la configuration
$`\sigma_{\partial_E C}`$. Si $`\partial_E C`$ est non vide, il possède
donc jusqu'à

```math
2^{|\partial_E C|-1}
```

états sous symétrie globale, après quotient par le flip global. Les facteurs
des ancêtres lisent ce
motif de frontière à travers leurs comptes satisfaits.

### Lemme 10.1 — homogénéité sans fermeture scalaire

Dans le GSBM homogène,

```math
\Lambda_v=u_pK_v.
```

Cette identité remplace chaque poids satisfait par un compte, mais deux
motifs de frontière de même taille peuvent donner des comptes $`K_v`$
différents, car les extrémités et les signes observés sont différents. Aucun
message scalaire dépendant seulement de $`(|C|,\beta_C)`$ ne ferme donc la
récursion sur la grille triangulaire.

### Proposition 10.2 — récursion de transfert tordue, statut : établi

Ordonnons les $`H`$ mises à jour pertinentes du chemin. Soit $`X_r`$ un état
suffisant juste avant la mise à jour $`r`$ ; il peut contenir le motif de
frontière, les comptes ancestraux courants et les décisions antérieures.
Conditionnellement à $`O,D`$, écrivons

```math
Q_r(x,a,dx')
=
\mathbb P(A_r=a,X_{r+1}\in dx'\mid X_r=x,O,D).
```

Pour une fonction test $`f`$, définissons l'opérateur tordu

```math
(\mathcal T_r f)(x)
:=
\sum_{a=0}^1
(-1)^{a\chi_r(i,j)}
\int f(x')Q_r(x,a,dx').
```

Si la loi initiale de $`X_1`$ est $`\lambda`$, alors la corrélation exacte du
chemin est

```math
\boxed{
\mathbb E\left[
\prod_{r=1}^H(-1)^{A_r\chi_r(i,j)}
\middle|O,D
\right]
=
\lambda\mathcal T_1\mathcal T_2\cdots\mathcal T_H\mathbf1.
}
```

#### Preuve

Développer le membre de droite somme successivement sur
$`(A_1,X_2),\ldots,(A_H,X_{H+1})`$. Le produit des noyaux $`Q_r`$ est la loi
jointe donnée par la règle de la chaîne, tandis que les facteurs tordus se
multiplient en $`(-1)^{\sum_rA_r\chi_r}`$. On retrouve le lemme 7.1.

Lorsque les séparateurs ont taille bornée, les $`\mathcal T_r`$ sont des
matrices finies et cette formule est un algorithme exact. PATH-FAC correspond
à l'écrasement non justifié de chaque état $`X_r`$ en un singleton et au
remplacement de $`\mathcal T_r`$ par le scalaire
$`\Gamma_{m_r}(t_r;p)`$.

La stratégie descendante devient néanmoins calculable lorsque les
séparateurs restent petits :

- sur un cactus, par messages finis attachés aux sommets d'articulation ;
- sur une bande de largeur $`W`$, par matrice de transfert exponentielle en
  $`W`$ ;
- sur la grille entière, seulement après un nouveau théorème de compression
  ou de concentration des messages de frontière.

Pour une paire $`i,j`$, on peut intégrer les branches latérales de bas en
haut, puis effectuer un passage forward--backward sur les deux bras vers le
LCA. C'est la version rigoureuse de l'idée « regarder le chemin », mais son
état doit conserver les interfaces avec les branches latérales.

## 11. Audit et contre-audit des conclusions

| Affirmation | Statut | Raison |
|---|---|---|
| Une racine finale est retournée avec probabilité $`1/2`$. | Établi sous a priori uniforme | tous les facteurs internes sont invariants |
| Une feuille redonne exactement Metropolis--Hastings. | Faux littéralement | le noyau construit est Glauber/Barker ; MH est une substitution valide |
| Les quatre probabilités d'un nœud sont explicites. | Établi | réduction exacte à $`(h_1,h_2,J)`$ |
| Les descendants contribuent directement au $`q_u^{ab}`$ d'un flip global des deux fils. | Faux | leurs facteurs sont communs aux quatre états |
| Les niveaux descendants suffisent à calculer le flip de $`u`$. | Faux | les comptes de frontière et le message ancestral manquent |
| La relation après un balayage est un produit de signes sur les deux bras. | Établi | identité déterministe de parité |
| L'espérance de ce produit factorise. | Non établi pour la dynamique complète | décisions emboîtées et messages ancestraux dépendants |
| La corrélation jointe se calcule par opérateurs tordus. | Établi | règle de la chaîne sur un état de frontière suffisant |
| La formule PATH-FAC est exacte. | Établi dans l'oracle défini | indépendance imposée explicitement |
| Le chemin physique marqué révèle la relation de la réplique génératrice. | Établi | toutes les arêtes gagnantes sont satisfaites |
| Ce chemin marqué prouve la weak recovery. | Faux | oracle dépendant de la réplique et variable auxiliaire différente |
| Les poids homogènes donnent un message descendant scalaire. | Faux en général | la frontière géométrique reste de dimension croissante |

## 12. Programme de calcul prioritaire

La voie descendante peut maintenant être testée sans ambiguïté.

1. Sur cactus triangulaires, calculer exactement le message de frontière et
   comparer le véritable balayage au produit PATH-FAC.
2. Sur bandes, mesurer l'écart entre la corrélation jointe
   $`\mathbb E\prod_w(-1)^{A_w}`$ et le produit des marginales.
3. Sur les tores triangulaires, enregistrer pour les deux bras d'une paire
   lointaine

   ```math
   (t_w,m_w,K_w,B_w)_{w\in\mathcal P(i,j)}
   ```

   ainsi que les tailles de frontière des branches latérales.
4. Comparer explicitement trois auxiliaires : dendrogramme non marqué,
   gagnants marqués et oracle PATH-FAC. Une amélioration observée seulement
   après marquage ne concerne pas la dynamique des slides.
5. Chercher une borne de contraction conditionnelle qui remplace
   l'indépendance dans la formule produit. C'est le lemme manquant pour
   transformer le chemin descendant en résultat de weak recovery.

Le calcul par descendants est donc une piste réelle, mais son objet correct
est une corrélation de décisions le long d'un balayage ou un message de
frontière collapsed. Il ne remplace pas algébriquement le message ancestral
du heat bath à un nœud unique.

Le [fichier 17](17_PATH_DECORRELATION_THRESHOLD.md) poursuit cette voie : il
caractérise exactement la limite $`1/2`$ par l'atténuation cumulée du chemin,
calcule les fenêtres en $`p`$ pour interfaces bornées ou logarithmiques et
remplace la factorisation, dans la dynamique jointe, par un critère de normes
$L^2$ des opérateurs tordus.
