# Bilan simple des arêtes résiduelles et obstruction hiérarchique

Cette note formalise l'intuition suivante : les premières arêtes activées
sont informatives, puis la qualité des arêtes restant entre les clusters se
dégrade avec le niveau $`\beta`$. Elle cherche volontairement le mécanisme le
plus simple qui puisse encore conduire à une obstruction de weak recovery.

> **Portée corrigée.** Les bilans globaux des sections 1--6 sont seulement
> des préliminaires. Ils ne décrivent pas la coupe sélectionnée par la
> géométrie. Le [fichier 25](25_GEOMETRY_CONDITIONED_CUT_INFORMATION.md)
> conditionne par la partition complète, prouve le critère
> $`m h_p(\beta)^2`$, puis dérive le biais LCA-Palm exact
> $`mN_\rho`$. Il constitue désormais la référence pour toute affirmation
> selon laquelle une coupe « cesse d'être informative ».

Le bilan exact donne une conclusion en deux parties.

1. L'intuition qualitative est correcte : les clusters précoces sont formés
   par des arêtes vraies dans la jauge de Nishimori, et le biais des arêtes de
   frontière décroît strictement avec $`\beta`$.
2. Il n'existe toutefois aucun temps $`\beta<1`$ auquel les **fausses** arêtes
   deviennent majoritaires face à **toutes** les vraies arêtes non encore
   activées. L'égalité n'a lieu qu'à $`\beta=1`$.

Le mécanisme simple qui subsiste est plus fin :

- une grande coupe devient presque non informative dans la fenêtre terminale
  $`1-\beta\asymp m^{-1/2}`$ ;
- un bucket de taille deux est exactement un canal d'effacement de coefficient
  $`s_p(\beta)<1`$ ;
- un nombre divergent de buckets bornés, ambigus et screenés suffit donc à
  faire disparaître la corrélation.

Le dernier point est un théorème sur un corridor factorisé et un théorème
conditionnel sur la grille. Le verrou restant est géométrique : obtenir ces
buckets sous la loi Palm d'une paire lointaine dont le LCA est critique, tout
en contrôlant les messages ancestraux.

## 1. Jauge plantée et quatre catégories d'arêtes

Écrivons l'observation d'une arête $`e=\{x,y\}`$ sous la forme

```math
Y_e=\Sigma_x\Sigma_y Z_e,
\qquad
\mathbb P(Z_e=+1)=p,
\qquad
\mathbb P(Z_e=-1)=1-p.
\qquad\text{(1.1)}
```

Sous le couplage annealed de Nishimori, on peut prendre la réplique qui génère
les horloges comme vérité de référence. Dans cette jauge,

```math
T_e\sim\mathrm{Exp}(u_p)
\quad\text{si }Z_e=+1,
\qquad
T_e=+\infty
\quad\text{si }Z_e=-1,
\qquad\text{(1.2)}
```

avec

```math
u_p=\log\frac p{1-p}.
\qquad\text{(1.3)}
```

Ainsi, toute arête ouverte est vraie dans cette jauge. Au temps
$`\beta\in[0,1]`$, les masses non conditionnelles se décomposent exactement en

| catégorie | événement | masse |
|---|---|---:|
| vraie précoce | $`T_e\le\beta`$ | $`p(1-e^{-u_p\beta})`$ |
| vraie future | $`\beta<T_e\le1`$ | $`p(e^{-u_p\beta}-e^{-u_p})`$ |
| vraie censurée | $`T_e>1`$ | $`pe^{-u_p}=1-p`$ |
| fausse | $`Z_e=-1`$ | $`1-p`$ |

Les quatre masses somment à un. La probabilité d'ouverture vaut

```math
q_p(\beta)=p(1-e^{-u_p\beta}).
\qquad\text{(1.4)}
```

Cette table ne doit être appliquée qu'aux arêtes d'une coupe ou d'une
frontière conditionnée à être encore fermée. Les marques des arêtes internes
à un cluster sont biaisées par la contrainte de connexité ; surtout, elles ne
figurent pas dans le $`\Lambda_u`$ de la coupe courante.

## 2. Premier bilan : toutes les vraies non activées contre les fausses

La masse des vraies arêtes dont l'horloge n'a pas encore sonné à
$`\beta`$ est

```math
M_{\rm true}^{>\beta}
=
pe^{-u_p\beta}.
\qquad\text{(2.1)}
```

La masse des fausses arêtes est

```math
M_{\rm false}=1-p=pe^{-u_p}.
\qquad\text{(2.2)}
```

### Lemme 2.1 — pas de majorité fausse avant la censure, statut : établi

Pour tout $`p>1/2`$,

```math
M_{\rm true}^{>\beta}-M_{\rm false}
=
p(e^{-u_p\beta}-e^{-u_p})
\begin{cases}
>0,&0\le\beta<1,\\
=0,&\beta=1.
\end{cases}
\qquad\text{(2.3)}
```

#### Preuve

Comme $`u_p>0`$, la fonction $`\beta\mapsto e^{-u_p\beta}`$ est strictement
décroissante. Pour $`\beta<1`$,

```math
e^{-u_p\beta}>e^{-u_p}=\frac{1-p}{p}.
```

La conclusion suit après multiplication par $`p`$.

### Conséquence

L'énoncé brut « à partir d'un temps intérieur, il y a davantage de fausses
arêtes que de vraies arêtes non encore sonnées » est faux pour cette dynamique.
Le temps d'égalité exact est le bord terminal $`\beta=1`$.

## 3. Deuxième bilan : vraies qui sonneront encore avant 1 contre fausses

Si l'on exclut les vraies arêtes censurées après $`1`$, la masse du signal
encore activable vaut

```math
M_{\rm future}(\beta)
=
p(e^{-u_p\beta}-e^{-u_p})
=
pe^{-u_p\beta}-(1-p).
\qquad\text{(3.1)}
```

Comparer cette masse à $`1-p`$ donne un croisement plus précoce.

### Lemme 3.1 — croisement du signal activable, statut : établi

Pour $`p>2/3`$, l'unique solution dans $`(0,1)`$ de

```math
M_{\rm future}(\beta)=M_{\rm false}
```

est

```math
\boxed{
\beta_{\rm act}(p)
=
1-\frac{\log2}{u_p}.
}
\qquad\text{(3.2)}
```

Pour $`\beta>\beta_{\rm act}(p)`$, les fausses arêtes sont plus nombreuses
que les vraies arêtes qui sonneront encore avant la censure.

#### Preuve

L'égalité équivaut à

```math
pe^{-u_p\beta}=2(1-p)=2pe^{-u_p}.
```

Après logarithme,

```math
-u_p\beta=\log2-u_p,
```

d'où (3.2). La solution est positive exactement lorsque
$`u_p>\log2`$, soit $`p>2/3`$.

### Contre-audit

Ce croisement n'est **pas** un seuil d'information. Les vraies censurées et
les fausses ont exactement la même masse $`1-p`$ et se compensent dans le
biais signé. Au temps $`\beta_{\rm act}`$, les trois catégories résiduelles

```math
(\text{vraie future},\text{vraie censurée},\text{fausse})
```

ont chacune probabilité conditionnelle $`1/3`$. La probabilité totale d'une
arête vraie est donc encore $`2/3`$.

Le fait que les fausses dépassent les seules vraies **futures** ne signifie
donc ni « aucune information », ni « flip équitable ».

## 4. Loi exacte d'une arête de frontière

Conditionnons une arête par le fait qu'elle n'est pas ouverte avant
$`\beta`$. Le normalisateur est

```math
1-q_p(\beta)
=
1-p+pe^{-u_p\beta}.
\qquad\text{(4.1)}
```

Sa probabilité d'être vraie vaut

```math
\boxed{
s_p(\beta)
=
\frac{pe^{-u_p\beta}}
{1-p+pe^{-u_p\beta}}
=
\frac1{1+e^{-u_p(1-\beta)}}.
}
\qquad\text{(4.2)}
```

Son biais signé est

```math
\boxed{
h_p(\beta)
=
2s_p(\beta)-1
=
\tanh\left(\frac{u_p(1-\beta)}2\right).
}
\qquad\text{(4.3)}
```

### Proposition 4.1 — dégradation continue, statut : établi

Sur $`[0,1]`$,

```math
s_p(0)=p,
\qquad
s_p(1)=\frac12,
\qquad
s_p'(\beta)<0.
\qquad\text{(4.4)}
```

En particulier, il n'existe pas de transition abrupte locale : le canal se
dégrade continûment de $`p`$ vers $`1/2`$.

#### Preuve

La forme logistique (4.2) donne les valeurs aux bords. Sa dérivée est

```math
s_p'(\beta)
=
-u_p s_p(\beta)(1-s_p(\beta))<0.
```

Le biais (4.3) est précisément la masse conditionnelle de la catégorie
« vraie future ». Les catégories vraie censurée et fausse ont les mêmes
masses et s'annulent dans le contraste vrai/faux.

## 5. Spécialisation au seuil de percolation

Sur la grille triangulaire,

```math
q_c=2\sin(\pi/18)=0.347296355334\ldots.
\qquad\text{(5.1)}
```

Le temps critique satisfait

```math
q_p(\beta_c)=q_c,
\qquad
\beta_c(p)
=
-\frac1{u_p}\log\left(1-\frac{q_c}{p}\right).
\qquad\text{(5.2)}
```

Il appartient à $`[0,1]`$ exactement lorsque

```math
p\ge p_{\rm SW}:=\frac{1+q_c}{2}=0.673648177667\ldots.
\qquad\text{(5.3)}
```

Les masses critiques non conditionnelles sont

```math
\begin{array}{c|c}
\text{catégorie}&\text{masse}\cr
\hline
\text{vraie précoce}&q_c\cr
\text{vraie future}&2p-1-q_c\cr
\text{vraie censurée}&1-p\cr
\text{fausse}&1-p.
\end{array}
\qquad\text{(5.4)}
```

Sur une arête de frontière critique,

```math
s_c(p)
=
\frac{p-q_c}{1-q_c},
\qquad
h_c(p)
=
\frac{2p-1-q_c}{1-q_c}.
\qquad\text{(5.5)}
```

Comme $`h_p(\beta)`$ décroît, $`\beta_c`$ maximise la qualité par arête parmi
les temps $`\beta\ge\beta_c`$, à coupe fixée. C'est le sens rigoureux le plus
simple de « la paire fusionne dans le cas le plus favorable ».

### Diagnostic « future contre fausse » au seuil

L'égalité

```math
2p-1-q_c=1-p
```

donne

```math
\boxed{
p_{\partial,\rm act}
=
\frac{2+q_c}{3}
=
0.782432118445\ldots.
}
\qquad\text{(5.6)}
```

C'est aussi l'unique valeur telle que

```math
\beta_c(p)=\beta_{\rm act}(p).
```

Mais

```math
p_{\partial,\rm act}
<
p_{\rm info}
=
\frac{1+\sqrt{q_c}}2
=
0.794659275831\ldots,
\qquad\text{(5.7)}
```

où $`p_{\rm info}`$ est la borne d'information-percolation déjà connue.
Même si le diagnostic (5.6) était converti en obstruction globale, il ne
battrait donc pas la meilleure borne actuelle.

### Valeurs à $`p=0.8`$

On a

```math
\beta_c=0.410716539196\ldots,
\qquad
\beta_{\rm act}=0.5.
\qquad\text{(5.8)}
```

À $`\beta_c`$, les masses non conditionnelles sont

```math
(q_c, 0.252703644666\ldots, 0.2, 0.2).
```

Après conditionnement par la frontière, les trois masses résiduelles sont

```math
(0.387164445505\ldots,
  0.306417777248\ldots,
  0.306417777248\ldots),
```

et

```math
s_c(0.8)=0.693582222752\ldots.
\qquad\text{(5.9)}
```

Le bucket critique reste donc nettement biaisé vers la vérité. Une preuve à
$`p=0.8`$ ne peut pas reposer sur une absence locale de majorité au LCA.

## 6. Pourquoi les premiers clusters sont informatifs

Dans la jauge plantée, seules les vraies arêtes possèdent une horloge finie.
Le graphe ouvert à $`\beta`$ a, sous la loi annealed, la loi d'une percolation
de paramètre $`q_p(\beta)`$. Toute connexion ouverte fournit un chemin dont
chaque arête est conforme à la réplique de référence.

Pour $`q_p(\beta)<q_c`$ à distance fixe du seuil, les composantes sont petites
et une paire macroscopiquement éloignée n'est connectée qu'avec probabilité
tendant vers zéro. Le premier temps non sous-critique auquel une telle paire
peut raisonnablement apparaître dans une même composante est donc
$`\beta_c`$.

Cela justifie le scénario favorable :

1. imposer que $`i,j`$ soient connectés dès la fenêtre critique ;
2. leur donner ainsi les canaux de frontière postcritiques les plus
   informatifs ;
3. montrer que même cet oracle favorable perd la corrélation.

La troisième ligne reste une réduction à prouver sur la géométrie : le
conditionnement critique est rare pour une paire uniforme et ne peut pas
être substitué à la moyenne globale sans domination favorable.

## 7. Canal exact d'un bucket de fusion

Considérons une coupe homogène de taille $`m`$ fusionnant au niveau
$`\beta`$. Le dendrogramme non marqué révèle qu'une arête satisfaite a gagné la
course, sans révéler son identité. Pour la parité latente
$`X\in\{-1,+1\}`$, le compte $`K`$ suit

```math
K\mid X=+1
\sim
1+\mathrm{Bin}(m-1,s_p(\beta)),
\qquad\text{(7.1)}
```

```math
K\mid X=-1
\sim
\mathrm{Bin}(m-1,1-s_p(\beta)).
\qquad\text{(7.2)}
```

Pour $`1\le k\le m-1`$, le log-rapport local vaut

```math
\ell_{m,k}(\beta;p)
=
\log\frac{k}{m-k}
+u_p(1-\beta)(2k-m).
\qquad\text{(7.3)}
```

Avec le message ancestral $`B_u`$, le log-rapport exact devient

```math
L_u=B_u+\ell_{m,K}(\beta;p).
\qquad\text{(7.4)}
```

Une majorité $`K>m/2`$ rend $`\ell_{m,K}>0`$, mais ne contrôle pas le signe
de $`L_u`$ si $`B_u`$ est négatif. « Plus d'arêtes vraies » est donc un
certificat local, pas encore une indépendance vis-à-vis de l'extérieur.

## 8. Deux vrais régimes de perte locale

Définissons la fiabilité oracle sans message extérieur

```math
\Gamma_m(\beta;p)
=
\mathbb Eleft[
\tanh^2\left(\frac{\ell_{m,K}(\beta;p)}2\right)
\right].
\qquad\text{(8.1)}
```

La probabilité moyenne de conserver la parité vraie lors du heat bath est

```math
\overline P_m(\beta;p)
=
\frac{1+\Gamma_m(\beta;p)}2.
\qquad\text{(8.2)}
```

### Théorème 8.1 — une grande coupe précoce est très informative, établi

Pour tout $`p>1/2`$ et tout $`\beta<1`$ fixé,

```math
\Gamma_m(\beta;p)\longrightarrow1
\qquad(m\to\infty).
\qquad\text{(8.3)}
```

Plus précisément, l'erreur décroît exponentiellement à l'échelle

```math
mD\left(\frac12\middle\|s_p(\beta)\right).
```

#### Idée de preuve

Comme $`s_p(\beta)>1/2`$ est fixé,

```math
\frac Km\longrightarrow s_p(\beta)>\frac12
```

en probabilité. Le terme logarithmique de (7.3) converge vers une constante
positive, tandis que le second terme est positif et de taille linéaire en
$`m`$. La probabilité d'une majorité erronée est exponentiellement petite par
Chernoff.

Ce résultat contre-audite une autre version trop forte de l'intuition : une
grande coupe critique n'est pas noyée par les fausses arêtes ; à $`p=0.8`$,
elle devient au contraire presque parfaitement informative.

### Théorème 8.2 — identité terminale, statut : établi

Au bord $`\beta=1`$,

```math
\boxed{
\Gamma_m(1;p)=\frac1m,
\qquad
\overline P_m(1;p)=\frac12+\frac1{2m}.
}
\qquad\text{(8.4)}
```

#### Preuve

À $`\beta=1`$, $`s_p(1)=1/2`$ et

```math
K=1+\mathrm{Bin}(m-1,1/2).
```

De plus,

```math
\tanh\left(\frac12\log\frac K{m-K}\right)
=
\frac{2K-m}{m}.
```

Or

```math
\mathbb E[(2K-m)^2]
=
\mathrm{Var}(2K-m)
+\mathbb E[2K-m]^2
=(m-1)+1=m.
```

Après division par $`m^2`$, on obtient (8.4).

Une seule grande coupe terminale, avec message extérieur neutre, rend donc
la parité asymptotiquement équitable.

### Théorème 8.3 — fenêtre terminale, statut : établi dans l'oracle local

Le régime non trivial est

```math
1-\beta\asymp m^{-1/2}.
```

Si

```math
u_p(1-\beta_m)=\frac a{\sqrt m},
```

alors, pour $`Z\sim\mathcal N(0,1)`$,

```math
\ell_{m,K}(\beta_m;p)
\Longrightarrow
aZ+\frac{a^2}{2},
\qquad\text{(8.5)}
```

et

```math
\Gamma_m(\beta_m;p)
\longrightarrow
\mathbb E\left[
\tanh^2\left(\frac{aZ+a^2/2}{2}\right)
\right].
\qquad\text{(8.6)}
```

Le paramètre simple qui distingue connexion et information est donc

```math
m h_p(\beta)^2,
```

et non le seul bilan du nombre moyen d'arêtes.

## 9. Le bucket $`m=2`$ : mécanisme simple d'obstruction

Prenons $`m=2`$ et notons $`s=s_p(\beta)`$. Les deux lois du compte sont

```math
\begin{array}{c|ccc}
&K=0&K=1&K=2\cr
\hline
X=+1&0&1-s&s\cr
X=-1&s&1-s&0.
\end{array}
\qquad\text{(9.1)}
```

Avec probabilité $`s`$, le compte révèle parfaitement la parité. Avec
probabilité $`1-s`$, l'événement $`K=1`$ ne donne aucune information locale.
Il s'agit exactement d'un canal d'effacement.

### Lemme 9.1 — fiabilité exacte, statut : établi

Sans message extérieur,

```math
\boxed{
\Gamma_2(\beta;p)=s_p(\beta).
}
\qquad\text{(9.2)}
```

Avec un message ancestral $`B`$, l'issue effacée conserve seulement ce
message. La fiabilité devient

```math
\boxed{
\kappa_2(B;\beta,p)
=
s_p(\beta)
+(1-s_p(\beta))\tanh^2(B/2).
}
\qquad\text{(9.3)}
```

Comme $`s_p(\beta)\le p`$ sur $`[0,1]`$, si $`|B|\le B_0<\infty`$,

```math
\kappa_2(B;\beta,p)
\le
\overline\kappa_2(p,B_0)
:=
p+(1-p)\tanh^2(B_0/2)
<1.
\qquad\text{(9.4)}
```

Cette borne est uniforme sur **tous** les niveaux, y compris ceux bien avant
$`\beta_c`$.

### Théorème 9.2 — corridor factorisé, statut : établi

Supposons qu'un corridor comporte $`N`$ parités indépendantes uniformes,
chacune observée par un bucket $`m=2`$ au niveau $`\beta_r`$. Pour la parité
produit entre les endpoints,

```math
\boxed{
A_N
=
\prod_{r=1}^N s_p(\beta_r)
\le
p^N.
}
\qquad\text{(9.5)}
```

En particulier, pour $`p<1`$ fixé,

```math
N\to\infty
\quad\Longrightarrow\quad
A_N\to0.
```

À $`p=0.8`$ et si tous les blocs sont critiques et neutres,

```math
A_N=s_c(0.8)^N.
```

Dix blocs donnent

```math
A_{10}=0.025761997386\ldots,
```

et quarante blocs

```math
A_{40}=4.4047181845\,10^{-7}.
```

## 10. Théorème simple conditionnel sur la grille

Le produit (9.5) n'est pas directement valide sur la grille : les parités et
les messages sont corrélés par les cycles et les ancêtres. Il suggère
cependant un lemme géométrique beaucoup plus simple que la description
complète de tous les buckets.

### Hypothèse SB — buckets bornés screenés

Sous la loi Palm d'une paire lointaine dont le LCA est dans la fenêtre
critique, il existe sur les deux bras descendants un nombre $`N_L`$ de blocs
disjoints tels que :

1. chaque bloc contient un bucket séparant la parité de l'endpoint, avec
   $`2\le m\le M`$ pour un $`M`$ fixe ;
2. les routes latérales sont screenées par l'état de bord du bloc ;
3. le message extérieur satisfait $`|B_r|\le B_0`$ ;
4. le transfert répliqué conditionnel de chaque bloc a un coefficient
   $`\chi^2`$ au plus $`\kappa(p,M,B_0)<1`$ ;
5. $`N_L\to\infty`$ en probabilité.

Pour le sous-cas $`m=2`$, le coefficient de la ligne 4 est explicitement
majoré par (9.4).

### Théorème 10.1 — obstruction sous SB, statut : conditionnel

Sous SB et sous une réduction favorable contrôlant les paires postcritiques,

```math
\mathbb E[A_{I_LJ_L}]longrightarrow0,
```

donc la weak recovery est impossible.

#### Preuve

La composition conditionnelle des transferts et la sous-multiplicativité
donnent

```math
A_{I_LJ_L}
\le
\kappa(p,M,B_0)^{N_L}
+\varepsilon_L,
\qquad\text{(10.1)}
```

avec $`\mathbb E\varepsilon_L=o(1)`$. Comme $`0<\kappa<1`$ et
$`N_L\to\infty`$ en probabilité,

```math
\kappa^{N_L}\longrightarrow0
```

en probabilité et dans $`L^1`$, car cette variable est bornée par un. Les
paires précoces ont masse $`o(1)`$ et les paires dans des racines distinctes
ont persistance nulle. La réduction pairwise du fichier 20 conclut.

### Forme quantitative

S'il existe $`c>0`$ tel que

```math
\mathbb P(N_L<c\log L)\longrightarrow0,
```

alors

```math
\mathbb E[A_{I_LJ_L}]
\le
o(1)+L^{-c|\log\kappa|}.
\qquad\text{(10.2)}
```

## 11. Ce que cette voie simple change dans le programme

La cible géométrique prioritaire peut être affaiblie. Le fichier 25 montre
que le bon premier tri porte sur la charge
$`\mathcal J_r=m_rh_p(\beta_r)^2`$. Avant de construire le transfert complet
d'une bande, il faut essayer de prouver :

```math
\boxed{
\text{sous Palm critique, le corridor contient }
N_L\to\infty
\text{ coupes de charge bornée et screenées.}
}
\qquad\text{(11.1)}
```

Une version minimale chercherait seulement des motifs triangulaires donnant
$`m=2`$. Une version robuste autoriserait des tailles variables sous la
contrainte $`m_rh_p(\beta_r)^2\le M`$ et certifierait leur coefficient exact,
avec l'arête gagnante et l'état de bord. Les buckets bornés restent un
sous-cas commode, pas la définition géométrique de la perte.

Si (11.1) échoue parce que les interfaces ou messages de bord ne peuvent pas
être isolés, alors le transfert de bande de largeur deux redevient nécessaire.
La bande est donc maintenant la seconde ligne de défense, pas un préalable
logique au test du mécanisme simple.

## 12. Audit et contre-audit

| affirmation | verdict | raison |
|---|---|---|
| Au début, seules les vraies arêtes sonnent dans la jauge plantée | Établi annealed | les fausses ont horloge $`+\infty`$ |
| Les petits clusters précoces sont informatifs | Établi dans le couplage | tout chemin ouvert est conforme à la réplique de référence |
| Avant $`1`$, les fausses dépassent toutes les vraies non activées | Faux | différence strictement positive dans (2.3) |
| Les fausses dépassent les vraies qui sonneront encore avant $`1`$ | Vrai après $`\beta_{\rm act}`$ | mais les vraies censurées restent présentes |
| $`\beta_{\rm act}`$ est un seuil de perte d'information | Faux | à ce temps, une arête fermée est encore vraie avec probabilité $`2/3`$ |
| Le diagnostic critique correspondant bat la borne connue | Faux | $`0.782432<0.794659`$ |
| Une grande coupe critique à $`p=0.8`$ est peu informative | Faux | à temps fixé $`<1`$, sa fiabilité tend vers un |
| Une grande coupe terminale est peu informative | Établi si $`B=0`$ | fiabilité exacte $`1/m`$ |
| Un bucket $`m=1`$ contracte | Faux | l'arête gagnante révèle parfaitement la parité |
| Un bucket $`m=2`$ contracte | Établi | canal d'effacement de coefficient $`s_p(\beta)`$ |
| Beaucoup de buckets $`m=2`$ suffisent toujours | Faux sans screening | un message ancestral divergent ou une route latérale peut conserver la parité |
| L'absence d'information locale implique l'absence de corrélation globale | Faux | l'extérieur peut déjà connaître la parité via $`B_u`$ |
| Les arêtes fausses internes votent au nœud courant | Faux | seules les arêtes de la coupe $`E_u`$ entrent dans $`\Lambda_u`$ |
| Le LCA critique est géométriquement le plus favorable | Établi sur cactus, ouvert sur grille | la qualité à taille fixée est ordonnée, pas les tailles des coupes |

## Conclusion

L'intuition du « signal qui s'épuise » est correcte si elle est formulée en
termes de biais continu $`h_p(\beta)`$ et de rapport signal sur bruit
$`m h_p(\beta)^2`$. Elle est fausse sous la forme d'une majorité soudaine de
fausses arêtes avant $`\beta=1`$.

Deux mécanismes simples et rigoureux restent disponibles :

1. une grande coupe fusionnant dans la fenêtre terminale
   $`1-\beta\asymp m^{-1/2}`$ ;
2. l'accumulation de buckets bornés ambigus, spécialement $`m=2`$, le long du
   corridor descendant.

Pour attaquer $`p=0.8`$, le second est le plus prometteur. Le scénario où
$`i,j`$ fusionnent dès $`\beta_c`$ reste le test favorable le plus sévère : le
bucket du LCA y est fortement informatif, mais la distance peut encore être
exploitée si le corridor contient un nombre divergent de motifs ambigus
screenés.
