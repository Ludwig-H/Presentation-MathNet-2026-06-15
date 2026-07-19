# Certificat rationnel A0 à $`p=0.805`$

> [!NOTE]
> **Jalon rigoureux archivé.** Ce certificat est subsumé par le
> [certificat P809439](../../results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).

> [!NOTE]
> Ce jalon exact est conservé pour sa marge simple $`1/200`$. La meilleure
> borne quantitative du dossier est le [certificat
> P809439](../../results/non_hierarchical/34_CERTIFICAT_RATIONNEL_P809439.md).

## 1. Résultat

Fixons

```math
p_0=\frac{161}{200},
\qquad
q_0=2p_0-1=\frac{61}{100},
```

et le canal d'effacement multi-état

```math
(a_0,s_0,e_0)=\frac1{125}(41,14,42).
\qquad\text{(1.1)}
```

Le calcul rationnel de ce fichier établit le lemme suivant.

**Théorème A0.** Pour tout a priori $`\mu`$ sur les quatre états relatifs et
toute fonction $`f`$,

```math
Q_{E_{a_0,s_0,e_0}}(\mu,f)
-Q_{Y_{q_0}}(\mu,f)
\ge
\frac1{200}\,\mathrm{Var}_\mu(f).
\qquad\text{(1.2)}
```

Il s'agit d'un certificat exhaustif : il n'utilise ni maillage du simplexe,
ni optimisation flottante, ni hypothèse sur la dynamique hiérarchique. Les
calculs sont effectués avec des fractions exactes et quatre suites de Sturm.

La preuve avait d'abord été formulée comme un problème de positivité de
mineurs rationnels à trois variables. Une annulation supplémentaire ramène
en fait tout le secteur polarisé à trois inégalités univariées et à la
dominance diagonale.

## 2. Profil entier du canal physique

Posons

```math
H(t)=(1521+24400t)(25921-24400t).
\qquad\text{(2.1)}
```

On a $`H(t)>0`$ sur $`[0,1]`$ et le profil exact du triangle vaut

```math
c_{q_0}(t)
=
\frac{14884t(21163-14884t)}{H(t)}.
\qquad\text{(2.2)}
```

Après soustraction du membre droit de (1.2), le coefficient du terme de
variance est

```math
\bar a
:=
a_0-\frac1{200}
=
\frac{323}{1000}.
\qquad\text{(2.3)}
```

Définissons enfin

```math
d(t)
:=
\frac{\bar a-c_{q_0}(t)}{t}
=
\frac{
12734546643-122688812000t+29232176000t^2
}{1000tH(t)}.
\qquad\text{(2.4)}
```

## 3. Secteur non polarisé

Le lemme des trois projections du
[fichier 11](../../results/non_hierarchical/11_TRIANGLE_BLOCK_SDPI.md) minore leur somme par
$`4\sum_x\mu_x^2f_x^2`$. Il suffit donc de vérifier

```math
\bar a+4s_0t-c_{q_0}(t)\ge0
\qquad
(0\le t\le1/2).
\qquad\text{(3.1)}
```

Le membre gauche est $`P_{\mathrm{np}}(t)/(1000H(t))`$, où

```math
P_{\mathrm{np}}(t)
=
12734546643
-105026035232t
+295953456000t^2
-266721280000t^3.
\qquad\text{(3.2)}
```

La suite de Sturm de $`P_{\mathrm{np}}`$ compte zéro racine sur
$`[0,1/2]`$. Comme $`P_{\mathrm{np}}(0)>0`$, l'inégalité (3.1) est stricte.
Cela démontre (1.2) dès que $`\max_x\mu_x\le1/2`$.

## 4. Séparation des masses faibles et dominantes

Supposons désormais $`\mu_0>1/2`$. Deux certificats de Sturm donnent

```math
d(t)\ge d(1/2)
\qquad
(0<t\le1/2),
\qquad\text{(4.1)}
```

et

```math
d'(t)<0
\qquad
(1/2\le t\le1).
\qquad\text{(4.2)}
```

Pour (4.1), après extraction du facteur positif $`1-2t`$, le polynôme à
certifier est

```math
P_{\mathrm{tail}}(t)
=
174730714488603
-1096598972008000t
+1792103257120000t^2.
\qquad\text{(4.3)}
```

Il n'a aucune racine sur $`[0,1/2]`$ et sa valeur en zéro est positive. Pour
(4.2), le numérateur de $`-d'(t)`$ est

```math
P_{\mathrm{dec}}(t)
=
502070211154001763
+15163279378752960000t
-96941433303509456000t^2
+146088022224640000000t^3
-17403668303360000000t^4.
\qquad\text{(4.4)}
```

La suite de Sturm compte zéro racine sur $`[1/2,1]`$ et
$`P_{\mathrm{dec}}(1/2)>0`$. Par conséquent, pour chacune des trois masses de
queue,

```math
d(\mu_i)\ge d(1/2)\ge d(\mu_0).
\qquad\text{(4.5)}
```

## 5. Positivité des termes hors diagonale

Le dernier certificat univarié est

```math
d(t)+\frac{s_0}{t(1-t)}>0
\qquad
(1/2\le t<1).
\qquad\text{(5.1)}
```

Après multiplication par le dénominateur positif
$`1000tH(t)(1-t)`$, son numérateur est

```math
P_{\mathrm{off}}(t)
=
17150240835
-68743038643t
+85240668000t^2
-29232176000t^3.
\qquad\text{(5.2)}
```

Il n'a aucune racine sur $`[1/2,1]`$ et sa valeur en $`1/2`$ est positive.

Pour $`\{i,j,k\}=\{1,2,3\}`$, écrivons

```math
D_k
=(\mu_0+\mu_k)(1-\mu_0-\mu_k)
=w(1-w),
\qquad
w=\mu_i+\mu_j.
\qquad\text{(5.3)}
```

Comme $`0<w\le1-\mu_0<1/2`$ et que $`w(1-w)`$ est croissante sur
$`[0,1/2]`$,

```math
D_k\le\mu_0(1-\mu_0).
\qquad\text{(5.4)}
```

Ainsi chaque terme hors diagonale de la matrice polarisée vérifie

```math
M_{ij}
=d(\mu_0)+\frac{s_0}{D_k}
\ge
d(\mu_0)+\frac{s_0}{\mu_0(1-\mu_0)}
>0.
\qquad\text{(5.5)}
```

## 6. Dominance diagonale

L'expression exacte des termes diagonaux est

```math
M_{ii}
=
d(\mu_0)+d(\mu_i)
+s_0\sum_{\ell\ne i}\frac1{D_\ell}.
\qquad\text{(6.1)}
```

Les deux termes hors diagonale de la même ligne contiennent exactement les
deux fractions de (6.1). Elles s'annulent, et il reste l'identité

```math
M_{ii}-M_{ij}-M_{ik}
=
d(\mu_i)-d(\mu_0)
\ge0.
\qquad\text{(6.2)}
```

Les termes hors diagonale étant positifs par (5.5), le membre gauche de
(6.2) est le résidu exact de dominance diagonale
$`M_{ii}-\sum_{j\ne i}|M_{ij}|`$. La matrice symétrique $`M`$ est donc à
diagonale dominante avec diagonale positive. Le théorème de Gershgorin donne
$`M\succeq0`$.

Cette preuve traite les a priori strictement positifs. Une masse nulle retire
simplement l'état correspondant du support ; de façon équivalente, les formes
quadratiques des deux canaux finis sont continues jusqu'aux faces du
simplexe. Les faces ternaires, binaires et les sommets sont donc incluses.
Sur la face binaire uniforme, la marge renforcée vaut exactement

```math
\frac{547}{1000}-\frac{7442}{13721}
=
\frac{63387}{13721000}
>0.
\qquad\text{(6.3)}
```

Cela achève la preuve de (1.2).

## 7. Vérification autonome

Le module
[rational_a0_less_noisy_certificate.py](../../computations/rational_a0_less_noisy_certificate.py)
construit les quatre polynômes depuis les paramètres rationnels, calcule les
suites de Sturm et vérifie les identités de dominance diagonale. Les
[tests associés](../../computations/test_rational_a0_less_noisy_certificate.py)
recalculent aussi directement les deux formes quadratiques sur des a priori
intérieurs et sur les faces.

Depuis la racine du dépôt :

```bash
python3 \
  research/hierarchical-swendsen-wang/computations/rational_a0_less_noisy_certificate.py
python3 -m unittest \
  research/hierarchical-swendsen-wang/computations/test_rational_a0_less_noisy_certificate.py \
  -v
```

La sortie du certificateur commence par

```text
status: CERTIFIED_PSD
scope: exhaustive
unresolved_regions: 0
method: exact Sturm certificates plus diagonal dominance
```

Ce statut certifie le lemme less-noisy local A0. Le passage à une borne de
weak recovery utilise ensuite la tensorisation facteur par facteur et le
régime strictement sous-critique de Chayes--Lei ; il est indépendant de la
dynamique hiérarchique.
