# Oracle de fusion critique et seuil local

Ce fichier pousse jusqu'au bout l'hypothèse favorable suivante : deux sommets
lointains sont étudiés conditionnellement au fait que leur fusion de Kruskal a
lieu au seuil géométrique de percolation. Le calcul local se ferme exactement
dans le GSBM homogène sur la grille triangulaire.

Le verdict comporte trois parties distinctes.

> **Périmètre.** Ce fichier calibre le bucket $u$. Le problème principal des
> slides 31--33 — estimer tous les $`\Lambda_v`$ pour $`v\succ u`$ — est
> désormais formulé et poussé dans le
> [fichier 10](10_ANCESTRAL_LAMBDA_ESTIMATION.md). Aucun résultat local de ce
> fichier ne permet de poser $`B_u=0`$ dans le heat bath hiérarchique complet.

1. **Résultat exact local.** Dans le modèle sans message ancestral
   $`B_u=0`$, le seuil de la grande coupe critique est exactement la borne
   Swendsen--Wang
   $`p_{\mathrm{SW}}=(1+q_c)/2`$. Au seuil,
   $`\Gamma_m^c=1/m`$ ; pour tout $`p>p_{\mathrm{SW}}`$ fixé,
   $`\Gamma_m^c\to1`$ exponentiellement vite. La fenêtre de transition est
   $`p-p_{\mathrm{SW}}\asymp m^{-1/2}`$ et sa limite est explicite.
2. **Contre-audit global.** Cette fiabilité est conditionnelle et oracle. Une
   fenêtre critique de masse $`o(1)`$ parmi les paires contribue au plus
   $`o(1)`$ à l'overlap, même si sa fiabilité conditionnelle tend vers $1$.
   Elle ne mesure donc pas, à elle seule, le seuil de weak recovery.
3. **Usage correct du cas favorable.** Le
   [fichier 12](12_FAVORABLE_HIERARCHICAL_REDUCTION.md) formule le lemme de
   domination qui remplace toutes les expériences postcritiques par l'oracle
   où la paire lointaine du même arbre se sépare au seuil. Sous ce lemme,
   l'échec de l'oracle favorable donne bien une impossibilité globale, sans
   demander que la masse réelle des temps LCA se concentre au seuil.

Autrement dit, se placer au temps critique est le bon oracle extrémal pour la
voie hiérarchique, mais son caractère « le plus favorable » doit être prouvé
sur le log-rapport complet. La quantité à contrôler n'est pas seulement le
bucket local : elle contient tous les $`\Lambda_v`$ pour $`v\succ u`$.

## 1. Cadre fini et convention de conditionnement

Soit $`G_n=(V_n,E_n)`$, avec $`|V_n|=n`$, et soient $`I_n,J_n`$ uniformes et
indépendants dans $`V_n`$. Pour le dendrogramme non marqué de Kruskal, posons

```math
\beta_{ij}
=
\inf\{t:i\leftrightarrow j\text{ dans }G_t\},
```

avec la valeur $`+\infty`$ lorsque les deux sommets appartiennent à des racines
distinctes à la coupe $`1`$. Le score LCA étendu est

```math
\eta_{ij}^{\mathrm{LCA}}
=
\mathbf1_{\{\beta_{ij}\le1\}}
\tanh^2\left(\frac{L_{u_{ij}}}{2}\right).
```

Les horloges ayant des lois continues, pour tout $`n`$ fini et tout temps
déterministe $`t`$,

```math
\mathbb P(\beta_{I_nJ_n}=t,\ I_n\ne J_n)=0.
```

L'expression « fusion exactement à $`\beta_c`$ » doit donc désigner l'un des
trois objets suivants :

- une fenêtre $`\beta_c<\beta_{ij}\le\beta_c+\delta`$ ;
- une désintégration régulière par rapport à la densité du temps de fusion ;
- une mesure de Palm pondérée par le flux des fusions de Kruskal.

Dans ce fichier, les énoncés globaux utilisent une fenêtre. Le calcul local est
la valeur de la désintégration au temps $`\beta_c`$.

## 2. Identité de masse : premier contre-audit

Rappelons

```math
S_n(t)
=
\mathbb P(I_n\leftrightarrow J_n\text{ dans }G_t)
=
\frac1{n^2}
\mathbb E\sum_{C\in\Pi_t}|C|^2.
```

Pour $`0<\delta\le1-\beta_c`$, définissons

```math
\mathcal A_{n,\delta}^c
=
\{\beta_c<\beta_{I_nJ_n}\le\beta_c+\delta\},
```

```math
\rho_{n,\delta}^c
=
\mathbb P(\mathcal A_{n,\delta}^c),
\qquad
\mathcal C_{n,\delta}^c
=
\mathbb E\left[
\eta_{I_nJ_n}^{\mathrm{LCA}}
\mathbf1_{\mathcal A_{n,\delta}^c}
\right].
```

### Proposition 2.1 — statut : établi, volume fini

On a exactement

```math
\rho_{n,\delta}^c
=
S_n(\beta_c+\delta)-S_n(\beta_c)
```

et

```math
\boxed{
0
\le
\mathcal C_{n,\delta}^c
\le
\rho_{n,\delta}^c.
}
```

La première identité vient de
$`\{\beta_{ij}\le t\}=\{i\leftrightarrow j\text{ dans }G_t\}`$. La seconde
vient de $`0\le\eta_{ij}^{\mathrm{LCA}}\le1`$.

Lorsque $`\rho_{n,\delta}^c>0`$, la fiabilité conditionnelle

```math
\overline\Gamma_{n,\delta}^c
=
\mathbb E\left[
\eta_{I_nJ_n}^{\mathrm{LCA}}
\mid
\mathcal A_{n,\delta}^c
\right]
```

vérifie donc

```math
\mathcal C_{n,\delta}^c
=
\rho_{n,\delta}^c\,
\overline\Gamma_{n,\delta}^c.
```

### Corollaire 2.2 — une fenêtre rare ne porte pas la weak recovery

Pour toute suite $`\delta_n`$ telle que

```math
S_n(\beta_c+\delta_n)\longrightarrow0,
```

on a

```math
\mathcal C_{n,\delta_n}^c\longrightarrow0,
```

même si $`\overline\Gamma_{n,\delta_n}^c\to1`$. Ce corollaire ne suppose
aucun exposant critique. Pour l'appliquer à une fenêtre proche-critique
particulière, il reste à prouver l'hypothèse sur $`S_n`$ dans l'exhaustion et
avec les conditions de bord choisies.

Sous A1, la borne de bande se décompose plus précisément en

```math
Q_n
\le
S_n(\beta_c)
+
\mathcal C_{n,\delta}^c
+
\mathcal M_n((\beta_c+\delta,1]).
```

Ainsi, concentrer toute l'analyse sur une fenêtre critique rétrécissante exige
un théorème supplémentaire montrant que le dernier terme contient encore une
masse macroscopique, ou qu'une dynamique ultérieure transforme les pivots
critiques en information cohérente sur une masse macroscopique de sommets.

Il existe une seconde voie, prioritaire ici : ne pas supposer la concentration
des temps réels, mais démontrer que leur expérience hiérarchique est dominée
par celle d'une paire lointaine qui se sépare au seuil. Le théorème 3.1 du
[fichier 12](12_FAVORABLE_HIERARCHICAL_REDUCTION.md) absorbe alors en une seule
borne la fenêtre et tout le reste tardif. La preuve de cette domination porte
nécessairement sur le message $`B_u`$ construit avec tous les taux ancestraux.

## 3. Les paramètres critiques se simplifient exactement

Dans le GSBM homogène, posons

```math
u_p=\log\frac p{1-p},
\qquad
q_p(t)=p(1-e^{-u_pt}),
\qquad
q_c=2\sin(\pi/18).
```

La coupe critique appartient à $`[0,1]`$ si et seulement si

```math
p\ge p_{\mathrm{SW}}
:=
\frac{1+q_c}{2}
=
0.6736481777\ldots,
```

et vaut

```math
\beta_c(p)
=
q_p^{-1}(q_c)
=
-\frac1{u_p}\log\left(1-\frac{q_c}{p}\right).
```

Pour une arête encore fermée au temps $`t`$, introduisons

```math
s_p(t)
=
\frac{pe^{-u_pt}}{1-p+pe^{-u_pt}},
\qquad
h_p(t)=2s_p(t)-1,
\qquad
a_p(t)=u_p(1-t).
```

### Proposition 3.1 — statut : établi, identités algébriques

Au temps critique, pour tout $`p\ge p_{\mathrm{SW}}`$,

```math
\boxed{
s_c(p)
:=
s_p(\beta_c)
=
\frac{p-q_c}{1-q_c},
}
```

```math
\boxed{
h_c(p)
:=
h_p(\beta_c)
=
\frac{2p-1-q_c}{1-q_c}
=
\frac{2(p-p_{\mathrm{SW}})}{1-q_c},
}
```

et

```math
\boxed{
a_c(p)
:=
a_p(\beta_c)
=
\log\frac{p-q_c}{1-p}
=
2\,\mathrm{artanh}(h_c(p)).
}
```

### Preuve

L'identité $`q_p(\beta_c)=q_c`$ donne

```math
pe^{-u_p\beta_c}=p-q_c.
```

La substitution dans $`s_p`$ donne la première formule, puis
$`h_c=2s_c-1`$ donne la deuxième. Enfin,

```math
u_p(1-\beta_c)
=
\log\frac p{1-p}
+
\log\left(1-\frac{q_c}{p}\right)
=
\log\frac{p-q_c}{1-p}.
```

Comme $`s_c=(1+h_c)/2`$, ce dernier terme est aussi
$`\log((1+h_c)/(1-h_c))=2\,\mathrm{artanh}(h_c)`$.

Au point $`p=p_{\mathrm{SW}}`$, ces identités donnent exactement

```math
\beta_c=1,
\qquad
s_c=\frac12,
\qquad
h_c=a_c=0.
```

## 4. Canal local critique exact

Fixons un bucket de $`m\ge1`$ arêtes qui fusionne à $`\beta_c`$, et
conditionnons par son squelette non marqué. L'arête gagnante latente est
satisfaite. Le nombre $`K`$ d'arêtes satisfaites suit donc

```math
K
\stackrel d=
1+\mathrm{Bin}(m-1,s_c).
```

Dans le modèle local $`B_u=0`$, le log-rapport des poids de parité vaut, pour
$`1\le k\le m-1`$,

```math
\ell_{m,k}^c(p)
=
\log\frac{k}{m-k}
+
a_c(p)(2k-m),
```

avec les conventions $`\ell_{m,0}^c=-\infty`$ et
$`\ell_{m,m}^c=+\infty`$.

### Contre-audit indépendant par les deux expériences de parité

Sous la parité vraie $`+`$,

```math
P_+(K=k)
=
\binom{m-1}{k-1}
s_c^{k-1}(1-s_c)^{m-k},
\qquad
1\le k\le m.
```

Sous la parité opposée, toutes les satisfactions sont complémentées, donc

```math
P_-(K=k)=P_+(K=m-k).
```

Pour $`1\le k\le m-1`$, un calcul direct donne

```math
\log\frac{P_+(K=k)}{P_-(K=k)}
=
\log\frac{k}{m-k}
+
(2k-m)\log\frac{s_c}{1-s_c}
=
\ell_{m,k}^c(p).
```

La formule issue du heat bath hiérarchique est donc exactement le LLR de
l'expérience statistique conditionnelle ; le facteur
$`\log(k/(m-k))`$ est la correction due à l'arête gagnante non marquée.

### Définition 4.1 — fiabilité de l'oracle local critique

Posons

```math
\boxed{
\Gamma_m^c(p)
=
\sum_{k=1}^m
\binom{m-1}{k-1}
s_c^{k-1}(1-s_c)^{m-k}
\tanh^2\left(\frac{\ell_{m,k}^c(p)}2\right).
}
```

Cette quantité est exacte, finie et comprise dans $`[0,1]`$. Elle est la
contraction $`L^2`$ de l'expérience binaire locale conditionnelle au bucket et
au temps critiques. Elle n'est pas encore la contraction non oracle du modèle
après marginalisation de tout le dendrogramme.

La probabilité moyenne de choisir les deux états de parité paire est reliée
exactement à cette fiabilité par

```math
\overline P_m^c(p)
=
\frac{1+\Gamma_m^c(p)}2.
```

La preuve, l'équivalent précis de son déficit et les hypothèses CUT/ANC
nécessaires pour passer d'une grande coupe à une paire lointaine de la grille
sont donnés dans le
[fichier 15](15_CRITICAL_GIANT_PAIR_FLIP.md).

### Proposition 4.2 — valeur exacte au seuil géométrique

Pour tout $`m\ge1`$,

```math
\boxed{
\Gamma_m^c(p_{\mathrm{SW}})=\frac1m.
}
```

### Preuve

Au seuil, $`s_c=1/2`$ et $`a_c=0`$. Si $`X=2K-m`$, alors

```math
\tanh\left(\frac12\log\frac K{m-K}\right)
=
\frac{2K-m}{m}
=
\frac Xm.
```

De plus,

```math
X
=
1+\sum_{r=1}^{m-1}\varepsilon_r,
```

où les $`\varepsilon_r`$ sont i.i.d. uniformes sur $`\{-1,+1\}`$. Ainsi
$`\mathbb E[X^2]=m`$, d'où $`\Gamma_m^c=\mathbb E[X^2]/m^2=1/m`$.
La preuve couvre aussi $`m=1`$ avec la convention de LLR infini.

## 5. Dichotomie de grande coupe et borne quantitative

### Théorème 5.1 — statut : établi

Pour tout $`p>p_{\mathrm{SW}}`$ fixé,

```math
\Gamma_m^c(p)\longrightarrow1
\qquad(m\to\infty).
```

Plus précisément, pour $`m\ge2`$, en posant
$`h=h_c(p)>0`$ et $`a=a_c(p)>0`$,

```math
\boxed{
1-\Gamma_m^c(p)
\le
\exp\left(-\frac{(m-1)h^2}{8}\right)
+
4\exp\left(-\frac{(m-1)ah}{2}\right).
}
```

### Preuve

Écrivons $`X=2K-m=1+\sum_{r=1}^{m-1}\varepsilon_r`$, où
$`\mathbb E\varepsilon_r=h`$. Sur l'événement

```math
\mathcal G_m
=
\left\{X\ge\frac{(m-1)h}{2}\right\},
```

on a $`K>m/2`$, donc $`\log(K/(m-K))\ge0`$, et par conséquent

```math
\ell_{m,K}^c
\ge
\frac{a(m-1)h}{2}.
```

Pour $`x\ge0`$,

```math
1-\tanh^2(x/2)
=
\frac1{\cosh^2(x/2)}
\le
4e^{-x}.
```

Une inégalité de Hoeffding donne par ailleurs

```math
\mathbb P(\mathcal G_m^c)
\le
\exp\left(-\frac{(m-1)h^2}{8}\right).
```

On borne la perte par $`1`$ sur $`\mathcal G_m^c`$ et par
$`4e^{-a(m-1)h/2}`$ sur $`\mathcal G_m`$.

### Corollaire 5.2 — le seuil local favorable retombe sur SW

Si l'on définit le bord du régime informatif de grande coupe par

```math
p_{\mathrm{loc}}^c
:=
\inf\left\{
p\ge p_{\mathrm{SW}}:
\liminf_{m\to\infty}\Gamma_m^c(p)>0
\right\},
```

alors

```math
\boxed{
p_{\mathrm{loc}}^c=p_{\mathrm{SW}}.
}
```

Il faut lire correctement cette égalité : au point frontière lui-même,
$`\Gamma_m^c=1/m\to0`$ ; pour tout $`p>p_{\mathrm{SW}}`$, la limite vaut
$`1`$. Pour $`m`$ fixé, la fiabilité est déjà positive au seuil, et le cas
$`m=1`$ donne même $`\Gamma_1^c=1`$ pour tout $`p`$. Ces faits montrent que
l'oracle **local avec $`B_u=0`$** est trop favorable pour localiser le seuil
réel de weak recovery. Ils ne dispensent pas d'estimer le message ancestral
du véritable oracle hiérarchique.

## 6. Fenêtre critique en $p$ : limite exacte

La borne précédente indique l'échelle $`m h_c(p)^2`$. Comme $`h_c`$ est
linéaire en $`p`$, la fenêtre non triviale est exactement de largeur
$`m^{-1/2}`$ autour de $`p_{\mathrm{SW}}`$.

### Théorème 6.1 — statut : établi

Fixons $`\alpha\ge0`$ et posons, pour $`m`$ assez grand,

```math
p_m
=
p_{\mathrm{SW}}
+
\frac{(1-q_c)\alpha}{2\sqrt m}.
```

Si $`Z\sim\mathcal N(0,1)`$, alors

```math
\boxed{
\Gamma_m^c(p_m)
\longrightarrow
\Psi(\alpha)
:=
\mathbb E\left[
\tanh^2(\alpha Z+\alpha^2)
\right].
}
```

### Preuve

Les identités critiques donnent exactement

```math
h_c(p_m)=\frac\alpha{\sqrt m},
\qquad
a_c(p_m)=2\,\mathrm{artanh}\left(\frac\alpha{\sqrt m}\right).
```

Avec $`X_m=2K-m`$, on a

```math
\frac{\mathbb E X_m}{\sqrt m}
=
\frac1{\sqrt m}+\frac{(m-1)\alpha}{m}
\longrightarrow\alpha,
```

et

```math
\mathrm{Var}\left(\frac{X_m}{\sqrt m}\right)
=
\frac{m-1}{m}\left(1-\frac{\alpha^2}{m}\right)
\longrightarrow1.
```

Les incréments sont uniformément bornés, donc la condition de Lindeberg est
automatique. Le théorème central limite pour tableaux triangulaires donne

```math
\frac{X_m}{\sqrt m}\Longrightarrow Z+\alpha.
```

Or

```math
\frac{\ell_{m,K}^c(p_m)}2
=
\mathrm{artanh}\left(\frac{X_m}{m}\right)
+
\mathrm{artanh}\left(\frac\alpha{\sqrt m}\right)X_m.
```

Le premier terme converge vers zéro en probabilité ; le second converge en
loi vers $`\alpha(Z+\alpha)`$. La fonction $`\tanh^2`$ étant continue et
bornée, la convergence des espérances suit.

Ce résultat est cohérent avec la fenêtre terminale du fichier 07 : son
paramètre $`a=u_p(1-\beta)`$ vaut ici $`2\alpha/\sqrt m+o(m^{-1/2})`$.

## 7. Tous les $`\Lambda_v`$ au-dessus d'une fusion critique

Supposons désormais $`\beta_u=\beta_c`$ et considérons un ancêtre strict
$`v\succ u`$. Alors

```math
\beta_c<\beta_v\le1,
```

donc, par décroissance de $`s_p`$,

```math
\frac12
\le
s_v:=s_p(\beta_v)
<
s_c.
```

Pour les trois groupes du bucket ancestral, le fichier 08 donne, conditionnellement
au squelette non marqué et à la catégorie gagnante $`G_v`$,

```math
K_{v,r}
=
\delta_{v,r}
+
\mathrm{Bin}(n_{v,r},s_v),
```

où

```math
\delta_{v,r}=\mathbf1_{\{G_v=r\}},
\qquad
n_{v,r}=m_{v,r}-\delta_{v,r}.
```

### Proposition 7.1 — sandwich critique des comptes

Sur un même espace muni de variables uniformes communes, on peut construire

```math
K_{v,r}^-
=
\delta_{v,r}+\mathrm{Bin}(n_{v,r},1/2),
```

```math
K_{v,r}^+
=
\delta_{v,r}+\mathrm{Bin}(n_{v,r},s_c),
```

de sorte que, presque sûrement,

```math
\boxed{
K_{v,r}^-
\le
K_{v,r}
\le
K_{v,r}^+.
}
```

La construction consiste à écrire chaque Bernoulli comme
$`\mathbf1_{\{U_e\le s\}}`$. Elle est simultanée pour tous les groupes et,
conditionnellement au squelette et aux catégories gagnantes, pour tous les
ancêtres.

Les quatre taux sont ensuite exactement

```math
\boxed{
\Lambda_v^{ab}
=
u_p\left[
K_{v,0}
+
c_a(K_{v,1};m_{v,1})
+
c_b(K_{v,2};m_{v,2})
\right],
}
```

avec

```math
c_0(k;m)=k,
\qquad
c_1(k;m)=m-k.
```

Le sandwich doit donc être complémenté pour les groupes retournés. Si
$`\underline K_{v,r}\le K_{v,r}\le\overline K_{v,r}`$, alors

```math
c_0(K_{v,r};m_{v,r})
\in
[\underline K_{v,r},\overline K_{v,r}],
```

mais

```math
c_1(K_{v,r};m_{v,r})
\in
[m_{v,r}-\overline K_{v,r},
  m_{v,r}-\underline K_{v,r}].
```

Cette règle donne des intervalles certifiés pour chaque
$`\Lambda_v^{00},\Lambda_v^{01},\Lambda_v^{10},\Lambda_v^{11}`$, puis pour
les quatre log-poids ancestraux.

### Contre-audit de monotonie

Le fait que $`s_v\le s_c`$ ne permet pas d'ordonner composante par composante
les quatre taux : retourner un groupe renverse l'ordre de son compte. Il ne
permet pas non plus d'ordonner directement la fiabilité complète. En effet,

```math
L_u
=
B_u+\ell_{m,K}^c,
```

et le message ancestral peut renforcer ou annuler le message local. Si un
calcul certifie $`B_u\in[\underline B,\overline B]`$ et si
$`\ell=\ell_{m,K}^c`$, posons

```math
d
=
\inf\{|x|:x\in[\ell+\underline B,\ell+\overline B]\},
```

```math
M
=
\max\{|\ell+\underline B|,|\ell+\overline B|\}.
```

Cette annulation est un énoncé point par point. Le modèle $`B_u=0`$ reste une
calibration, pas un préfixe de la filtration réelle. La
[calcul des frontières critiques](14_CRITICAL_COMPONENT_BOUNDARY.md)
donne le contrôle qui subsiste sans annulation : chaque ancêtre possède un
biais $`h_p(\beta_v)<h_p(\beta_c)`$, mais la préférence finale est décidée
par la somme des deux poids pairs contre celle des deux poids impairs.

Alors seulement

```math
\boxed{
\tanh^2(d/2)
\le
\eta_u
\le
\tanh^2(M/2).
}
```

En particulier, la borne inférieure vaut zéro dès que l'intervalle ancestral
peut contenir $`-\ell`$. Toute affirmation selon laquelle la fusion à
$`\beta_c`$ maximise le canal hiérarchique complet demande donc un ordre de
Blackwell, une SDPI ou une preuve de monotonie supplémentaire. Ce qui est
déjà rigoureux est plus limité : $`h_p(t)`$, $`a_p(t)`$ et l'échelle locale
$`m h_p(t)^2`$ sont maximaux à $`t=\beta_c`$ parmi les ancêtres de cette
fusion.

## 8. Ce que l'oracle critique permet — et ne permet pas

Trois opérations différentes ne doivent pas être confondues.

| Opération | Relation avec le problème original |
|---|---|
| Révéler $`D`$ sans changer sa loi | véritable oracle, donc information supplémentaire |
| Conditionner par $`\beta_{ij}\simeq\beta_c`$ | change la population de paires ; aucun ordre global automatique |
| Poser $`B_u=0`$ | sous-modèle local soluble ; ni majorant ni minorant universel du message complet |

Par conséquent :

- l'impossibilité pour un véritable oracle dominant impliquerait
  l'impossibilité originale ;
- le succès d'un oracle ne donne aucune condition suffisante pour le modèle
  non oracle ;
- le succès ou l'échec après conditionnement sur une classe rare de paires ne
  conclut rien sans remettre son facteur de masse $`\rho_{n,\delta}^c`$ ;
- la valeur $`p_{\mathrm{loc}}^c=p_{\mathrm{SW}}`$ est un résultat de
  calibration, pas le seuil de weak recovery de la grille triangulaire.

Le contraste numérique est instructif :

```math
p_{\mathrm{loc}}^c
=
p_{\mathrm{SW}}
=
0.673648\ldots
<
p_{\mathrm{info}}
=
0.794659\ldots,
```

tandis que $`0.835806\ldots`$ reste seulement le repère numérique du point
multicritique de Nishimori dans ce dossier. La borne
d'information-percolation interdit déjà la weak recovery pour
$`p<p_{\mathrm{info}}`$. Le critère local qui déclare toute grande coupe
critique informative dès $`p>p_{\mathrm{SW}}`$ ne peut donc pas, à lui seul,
améliorer cette borne d'impossibilité ; il est trop informatif.

## 9. Formalisation correcte de la prochaine cible

La quantité critique utile n'est pas $`\overline\Gamma_{n,\delta}^c`$ seule,
mais la mesure non normalisée

```math
\mathfrak C_{n,p}(dt)
=
\frac2{n^2}
\mathbb E_{\mathscr D}
\sum_u
|C_{u,1}||C_{u,2}|\,
\Gamma_u(\mathscr D)\,
\delta_{\beta_u}(dt).
```

Elle doit ensuite être remplacée, pour une condition suffisante, par une
contraction où $`D`$ est marginalisé et par un score dont les signes sont
cohérents entre les paires. Un programme rigoureux comporte donc les étapes
suivantes.

1. **Géométrie.** Déterminer la loi pondérée par
   $`|C_{u,1}||C_{u,2}|`$ de
   $`(\beta_u,m_{v,r},\beta_v)_{v\succeq u}`$, sans conditionner gratuitement
   sur une paire rare.
2. **Marques.** Utiliser la loi groupée exacte et le sandwich de la section 7
   pour tous les $`\Lambda_v^{ab}`$.
3. **Marginalisation.** Calculer une contraction de bundle après oubli du
   dendrogramme ; elle doit redonner $`(2p-1)^2`$ pour une arête.
4. **Agrégation.** Contrôler la cohérence signée ou un rayon spectral global,
   et non la seule moyenne de $`\eta_u`$.
5. **Fenêtre.** Vérifier séparément la masse
   $`S_n(\beta_c+\delta_n)-S_n(\beta_c)`$ et la fiabilité conditionnelle.

Sur cactus, les étapes 1 à 3 sont accessibles par récursion exacte. Sur bandes
de largeur fixée, elles se prêtent à une matrice de transfert par intervalles.
Sur la grille entière, la loi du squelette proche-critique et la cohérence
signée restent les deux verrous.

## 10. Statut et contre-audits reproductibles

| Énoncé | Statut | Contre-audit |
|---|---|---|
| Identités $`s_c,h_c,a_c`$ | Établi | substitution dans $`q_p(\beta_c)=q_c`$ |
| Loi $`K=1+\mathrm{Bin}(m-1,s_c)`$ | Établi conditionnellement au squelette non marqué | course exponentielle du fichier 08 |
| LLR local $`\ell_{m,k}^c`$ | Établi | rapport indépendant $`P_+(k)/P_-(k)`$ |
| $`\Gamma_m^c(p_{\mathrm{SW}})=1/m`$ | Établi | identité de second moment de la marche symétrique |
| Convergence exponentielle pour $`p>p_{\mathrm{SW}}`$ | Établi | Hoeffding plus borne exacte de $`\mathrm{sech}^2`$ |
| Limite $`m^{-1/2}`$ | Établi | TCL triangulaire et calcul numérique indépendant |
| Sandwich des $`\Lambda_v^{ab}`$ | Établi conditionnellement au squelette | couplage par uniformes et complément des groupes retournés |
| Fusion critique maximisant tout le canal ancestral | Non établi | faux raisonnement monotone identifié ; annulation possible par $`B_u`$ |
| Seuil réel de weak recovery égal à $`p_{\mathrm{SW}}`$ | Exclu par la borne $`p_{\mathrm{info}}>p_{\mathrm{SW}}`$ | l'oracle oublie la contraction non oracle et la masse |

Le script
[computations/critical_merger_oracle.py](computations/critical_merger_oracle.py)
recalcule les identités, la somme finie $`\Gamma_m^c`$, la borne exponentielle
et la limite de fenêtre. Les tests unitaires vérifient les mêmes résultats par
des représentations algébriques indépendantes.
