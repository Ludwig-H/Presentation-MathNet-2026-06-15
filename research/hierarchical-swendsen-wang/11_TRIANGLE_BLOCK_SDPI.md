# Canal de triangle, information latérale et candidat multi-état

> **Statut dans le programme.** Ce fichier est un audit auxiliaire d'un canal
> physique. Il ne remplace pas la dynamique hiérarchique des slides 31--33.
> La voie prioritaire est maintenant le
> [corridor collapsed critique](20_COLLAPSED_CORRIDOR_BLACKWELL.md), qui
> conserve tous les $`\Lambda_v^{ab}`$ au-dessus du LCA d'une paire lointaine.

Ce fichier répond à la question quantitative suivante : peut-on dépasser la
borne d'impossibilité

```math
p_{\mathrm{info}}
=\frac{1+\sqrt{2\sin(\pi/18)}}2
=0.7946592758\ldots
```

sur le GSBM triangulaire en exploitant simultanément les trois observations
d'un triangle, puis les quatre poids construits à partir des
$`\Lambda_v`$ ancestraux ?

La réponse actuelle est précise mais en deux parties.

1. **Établi.** Regrouper un triangle en un seul facteur scalaire ne peut pas
   améliorer information-percolation. La contraction uniforme séduisante
   $`\eta_\triangle`$ n'est pas stable sous information latérale. La vraie
   contraction SDPI globale du facteur vaut
   $`2q^2/(1+q^2)`$ et conduit seulement à
   $`p<(1+1/\sqrt3)/2=0.788675\ldots`$.
2. **À prouver.** Un canal d'effacement multi-état — vide, une relation
   révélée, ou triangle entièrement révélé — conserve la géométrie des trois
   modes. Son enveloppe exacte produit le candidat conditionnel
   $`p_\star=0.8099092892\ldots`$. L'inégalité quadratique requise par le
   critère less-noisy est démontrée pour tous les a priori quatre états non
   polarisés. Il reste à la démontrer lorsqu'un état a masse strictement
   supérieure à $`1/2`$ ; sans ce lemme, $`p_\star`$ n'est **pas** une
   nouvelle borne rigoureuse.

Cette séparation est importante pour la chaîne ancestrale. Les slides 31--33
imposent précisément de calculer les quatre poids du heat bath, donc aussi
leur polarisation ; remplacer ces poids par l'a priori uniforme supprime la
difficulté au lieu de la résoudre.

## 1. Partition canonique de la grille triangulaire

Écrivons les sommets de la grille comme $`(i,j)\in\mathbb Z^2`$ avec les
trois directions $`(1,0),(0,1),(-1,1)`$, et posons

```math
T_{i,j}=\{(i,j),(i+1,j),(i,j+1)\}.
```

Les triangles montants $`T_{i,j}`$ sont disjoints par arêtes et couvrent
chaque arête exactement une fois. Les trois observations contenues dans deux
triangles distincts sont donc conditionnellement indépendantes sachant les
spins. Deux triangles montants sont adjacents lorsqu'ils partagent un sommet ;
leur graphe d'intersection est encore une grille triangulaire.

On peut ainsi utiliser les triangles montants comme nœuds facteurs de degré
trois. Le théorème multi-terminal de Polyanskiy--Wu s'applique bien à de tels
facteurs ; en revanche, son paramètre est la contraction SDPI **globale** du
canal facteur, pas sa contraction sous l'a priori uniforme.

## 2. Canal relatif exact d'un triangle

Posons

```math
q=2p-1\in[0,1),
```

et, modulo le flip global,

```math
A=\sigma_1\sigma_2,
\qquad
B=\sigma_1\sigma_3.
```

Le canal d'entrée $`X=(A,B)\in G:=\{\pm1\}^2`$ est

```math
Y_1=AZ_1,
\qquad
Y_2=BZ_2,
\qquad
Y_3=ABZ_3,
```

où les $`Z_r`$ sont indépendants et
$`\mathbb E Z_r=q`$. Le syndrome

```math
S=Y_1Y_2Y_3=Z_1Z_2Z_3
```

est indépendant de $X$. Ses probabilités sont

```math
w_+=\frac{1+q^3}{2},
\qquad
w_-=\frac{1-q^3}{2}.
```

Conditionnellement à $`S=\pm1`$, et après une translation connue de la
sortie, le canal $`X\mapsto(Y_1,Y_2)`$ est le canal symétrique à quatre
symboles

```math
T_\lambda(y\mid x)
=
\lambda\mathbf1_{\{y=x\}}+\frac{1-\lambda}{4}
```

avec

```math
\lambda_+=\frac q{1-q+q^2},
\qquad
\lambda_-=-\frac q{1+q+q^2}.
```

Soit $`\mu=(\mu_x)_{x\in G}`$ un a priori quelconque et soit
$`f:G\to\mathbb R`$ tel que $`\sum_x\mu_xf_x=0`$. Un calcul direct dans les
deux canaux $`T_{\lambda_+},T_{\lambda_-}`$ donne l'identité quadratique
exacte

```math
\boxed{
\mathrm{Var}\bigl(\mathbb E[f(X)\mid Y]\bigr)
=
\sum_{x\in G}\mu_x c_q(\mu_x)f_x^2,
}
```

où

```math
\boxed{
c_q(t)
=
2q^2t\left[
\frac{1+q}{(1-q)^2+4qt}
+
\frac{1-q}{(1+q)^2-4qt}
\right].
}
```

Cette formule est valable pour $`0\le t\le1`$. Lorsque $`t>1/2`$, la valeur
$`c_q(t)`$ peut dépasser $1$ ; ce n'est pas une contraction utilisable seule,
car la contrainte de moyenne nulle couple nécessairement ce coin aux trois
autres.

### Preuve courte de l'identité

Pour $`T_\lambda`$, la masse de sortie en $x$ vaut

```math
\nu_x=\frac{1-\lambda}{4}+\lambda\mu_x.
```

La moyenne de $f$ étant nulle,

```math
\mathbb E[f(X)\mid Y=x]
=\frac{\lambda\mu_xf_x}{\nu_x}.
```

Il suffit de sommer $`\lambda^2\mu_x^2f_x^2/\nu_x`$ avec les poids
$`w_+,w_-`$, puis d'utiliser les quatre identités

```math
w_+=\frac{(1+q)(1-q+q^2)}2,
\quad
w_- =\frac{(1-q)(1+q+q^2)}2,
```

```math
1-\lambda_+=\frac{(1-q)^2}{1-q+q^2},
\quad
1-\lambda_-=\frac{(1+q)^2}{1+q+q^2}.
```

## 3. Deux calibrations et la contraction globale

Sous l'a priori uniforme $`\mu_x=1/4`$,

```math
\boxed{
c_q(1/4)
=\eta_\triangle(q)
=\frac{q^2(1+2q^2)}{1+q^2+q^4}.
}
```

Si une relation indépendante est déjà révélée par l'extérieur, l'a priori
est uniforme sur deux états. Les deux arêtes variables sont alors deux
observations BSC indépendantes du même bit et

```math
\boxed{
c_q(1/2)
=\gamma_2(q)
=\frac{2q^2}{1+q^2}.
}
```

Pour $`0<q<1`$,

```math
\gamma_2(q)-\eta_\triangle(q)
=
\frac{q^2(1-q^2)}
{(1+q^2)(1+q^2+q^4)}
>0.
```

L'information latérale augmente donc strictement la contraction.

### Proposition 3.1 — coefficient SDPI global exact

Pour le canal de triangle,

```math
\boxed{
\eta_{\chi^2}^{\mathrm{glob}}(q)
=
\eta_{\mathrm{KL}}^{\mathrm{glob}}(q)
=
\frac{2q^2}{1+q^2}.
}
```

La seconde égalité utilise l'identité générale entre les coefficients globaux
KL et $`\chi^2`$ d'un canal fini. Voici une preuve autonome de la première.

La fonction $`c_q`$ est strictement croissante : chacune des deux fonctions
$`t/[(1\mp q)^2\pm4qt]`$ apparaissant dans sa définition a une dérivée
strictement positive. Posons $`g=c_q(1/2)`$. Si toutes les masses
$`\mu_x`$ sont au plus $`1/2`$, l'identité de la section 2 donne immédiatement
une contraction au plus $g$.

Supposons donc $`\mu_1=t>1/2`$, posons $`s_i=\mu_i`$ pour $`i>1`$ et

```math
d_+(t)=c_q(t)-g,
\qquad
d_-(s)=g-c_q(s).
```

Avec $`a=1+q^2`$ et $`w\in(0,1)`$, un calcul direct donne

```math
g-c_q\left(\frac{1-w}{2}\right)
=
\frac{2q^2(1-q^2)w(a-2q^2w)}
{a(a^2-4q^2w^2)},
```

```math
c_q\left(\frac{1+w}{2}\right)-g
=
\frac{2q^2(1-q^2)w(a+2q^2w)}
{a(a^2-4q^2w^2)}.
```

Définissons $`h(s)=s/d_-(s)`$ sur $`(0,1/2)`$. Comme $`c_q`$ est
croissante, $`h(s)/s`$ est croissante ; $h$ est donc superadditive. Les deux
formules précédentes montrent en outre

```math
h(1-t)\le\frac{t}{d_+(t)},
```

car cette inégalité se réduit à
$`(2t-1)(1-q^2)\ge0`$. Enfin,
$`tf_1=-\sum_{i>1}s_if_i`$ et Cauchy--Schwarz donnent

```math
t^2f_1^2
\le
\left(\sum_{i>1}\frac{s_i}{d_-(s_i)}\right)
\left(\sum_{i>1}s_id_-(s_i)f_i^2\right)
\le
\frac{t}{d_+(t)}
\sum_{i>1}s_id_-(s_i)f_i^2.
```

Ainsi la contribution positive du coin dominant est compensée par les trois
autres. L'égalité est atteinte par l'a priori uniforme sur deux états, ce qui
achève la preuve.

## 4. Contre-audit de la percolation scalaire

Si chaque triangle facteur était ouvert avec la contraction uniforme
$`\eta_\triangle(q)`$, la percolation de ses centres aurait seuil $`1/2`$ et
on obtiendrait formellement

```math
\eta_\triangle(q)=\frac12
\quad\Longleftrightarrow\quad
3q^4+q^2-1=0,
```

donc

```math
p=\frac12\left[
1+\sqrt{\frac{\sqrt{13}-1}{6}}
\right]
=0.8294914816\ldots
```

Cette valeur est **fausse comme borne démontrée** : le théorème facteur de
Polyanskiy--Wu demande la contraction globale de la proposition 3.1. Le vrai
test scalaire est

```math
\frac{2q^2}{1+q^2}<\frac12
\quad\Longleftrightarrow\quad
q^2<\frac13,
```

soit seulement

```math
p<\frac12\left(1+\frac1{\sqrt3}\right)
=0.7886751346\ldots
```

Cette borne est strictement plus faible que $`p_{\mathrm{info}}`$. Le
regroupement scalaire des triangles est donc une impasse quantitative, même si
le calcul du canal de triangle est exact.

## 5. Canal d'effacement multi-état

Pour conserver les trois modes, introduisons un canal $`E_{a,s,e}`$ qui,
indépendamment de l'entrée :

- révèle les deux bits relatifs avec probabilité $a$ ;
- choisit chacune des trois relations non triviales et la révèle exactement
  avec probabilité $s$ ;
- ne révèle rien avec probabilité $e$ ;
- satisfait $`a+3s+e=1`$.

Géométriquement, ce sont respectivement un triangle plein, chacun des trois
états à une arête et le triangle vide du modèle de Chayes--Lei.

Pour un a priori $\mu$ et une fonction centrée $f$, sa forme quadratique est

```math
Q_E(\mu,f)
=
a\,\mathrm{Var}_\mu(f)
+s\sum_{\chi\ne1}
\mathrm{Var}_\mu\bigl(\mathbb E[f\mid\chi(X)]\bigr).
```

### Lemme 5.1 — borne diagonale des trois projections

On a exactement

```math
\boxed{
\sum_{\chi\ne1}
\mathrm{Var}_\mu\bigl(\mathbb E[f\mid\chi(X)]\bigr)
\ge
4\sum_{x\in G}\mu_x^2f_x^2.
}
```

Pour le voir, écrivons $`g_x=\mu_xf_x`$, donc $`\sum_xg_x=0`$. Les trois
bipartitions de $G$ donnent des termes

```math
\frac{(g_i+g_j)^2}{m_{ij}(1-m_{ij})}
\ge4(g_i+g_j)^2.
```

En choisissant les trois paires contenant un état fixé,

```math
\sum_{j\ne i}(g_i+g_j)^2=\sum_xg_x^2,
```

ce qui prouve le lemme.

### Corollaire 5.2 — secteur non polarisé

Si

```math
\max_x\mu_x\le\frac12
```

et si

```math
\boxed{
a+4st\ge c_q(t)
\qquad(0\le t\le1/2),
}
```

alors $`Q_E(\mu,f)\ge Q_Y(\mu,f)`$ pour toute fonction centrée $f$. C'est une
comparaison $`\chi^2`$ exacte, uniforme sur tous les a priori non polarisés.

Les seules contraintes aux extrémités,

```math
a+s\ge\eta_\triangle(q),
\qquad
a+2s\ge\gamma_2(q),
```

ne suffisent pas. À $`p=p_{\mathrm{info}}`$, le choix qui met ces deux
inégalités à égalité donne

```math
a=0.2863093150\ldots,
\qquad
s=0.1146181435\ldots,
```

mais, à $`t=0.3136932605\ldots`$,

```math
a+4st-c_q(t)=-0.0012919051\ldots<0.
```

C'est le second contre-audit : même après avoir traité l'information latérale
binaire, il faut contrôler tout le profil quatre états.

## 6. Enveloppe optimale et constante conditionnelle

Ancrons la droite au cas binaire exact $`a+2s=\gamma_2(q)`$ et définissons

```math
s(q)
:=
\inf_{0\le t<1/2}
\frac{\gamma_2(q)-c_q(t)}{2-4t},
```

```math
a(q):=\gamma_2(q)-2s(q),
\qquad
e(q):=1-a(q)-3s(q).
```

Alors $`a(q)+4s(q)t\ge c_q(t)`$ sur $`[0,1/2]`$ par construction, et le
corollaire 5.2 s'applique à tout a priori non polarisé.

Cette enveloppe ne nécessite pas d'optimisation numérique. Avec $`r=q^2`$ et
$`x=2t-1`$,

```math
\frac{\gamma_2(q)-c_q(t)}{2-4t}
=
\frac{r(1-r)(1+r+2rx)}
{(1+r)((1+r)^2-4rx^2)}.
```

Sa dérivée s'annule en l'unique point intérieur

```math
x(q)=-\frac{1+r}{2(1+\sqrt{1-r})},
\qquad
t(q)=\frac12-\frac{1+r}{4(1+\sqrt{1-r})}.
```

Les valeurs $`s(q),a(q),e(q)`$ sont donc des expressions fermées.

Le seuil autodual du modèle triangulaire multi-état est $`a=e`$, soit

```math
2a(q)+3s(q)=1.
```

Dans la branche où l'infimum est atteint à l'intérieur, cette équation et la
condition de tangence s'écrivent

```math
4s=c_q'(t),
\qquad
c_q(t)=\gamma_2(q)-2s+4st,
\qquad
s=2\gamma_2(q)-1.
```

Elles s'éliminent exactement. La valeur pertinente
$`q_\star\in(0.61,0.63)`$ se voit sans calcul formel opaque. Sur la surface
autoduale,

```math
s=\frac{3r-1}{1+r},
\qquad
a=e=\frac{2(1-2r)}{1+r}.
```

Après multiplication par un dénominateur strictement positif, le numérateur
de $`a+4st-c_q(t)`$ est

```math
x\left[
(1+r)(4r^2+r-1)
-2r^2(1-r)x
-4r(3r-1)x^2
\right].
```

La tangence intérieure impose que le facteur quadratique ait une racine
double. Son discriminant nul donne exactement

```math
q^{10}+46q^8+45q^6-20q^4-12q^2+4=0.
```

Cette équation a une unique racine dans $`(0.61,0.63)`$ : le polynôme change
de signe aux extrémités et sa dérivée en $r$ y est strictement positive. La
masse de tangence s'écrit exactement

```math
t_\star
=\frac12-\frac{r_\star(1-r_\star)}{8(3r_\star-1)}.
```

Numériquement,

```math
q_\star=0.619818578503838\ldots,
\qquad
p_\star^{\mathrm{cond}}
=\frac{1+q_\star}{2}
=0.809909289251919\ldots,
```

et

```math
a_\star=e_\star=0.334711792547597\ldots,
```

```math
s_\star=0.110192138301602\ldots,
\qquad
t_\star=0.306110261668054\ldots.
```

Les deux conditions auxiliaires utilisées dans le régime de Chayes--Lei sont
largement satisfaites :

```math
a_\star e_\star>2s_\star^2,
\qquad
a_\star+e_\star>
\frac{2\sqrt2}{3+2\sqrt2}.
```

### Lemme manquant $`P_\star`$

Il suffit d'établir, pour tout a priori $\mu$ sur $G$ et toute fonction
centrée $f$,

```math
\boxed{
Q_{Y_{q_\star}}(\mu,f)
\le
Q_{E_{a_\star,s_\star,e_\star}}(\mu,f).
}
```

Le critère matriciel de Makur--Polyanskiy montre que $`P_\star`$ pour tous
$`\mu,f`$ est exactement la relation less-noisy
$`E_{a_\star,s_\star,e_\star}\succeq_{\mathrm{ln}}Y_{q_\star}`$. Il n'y a
donc pas une seconde étape informationnelle cachée après ce lemme.

Le corollaire 5.2 prouve $`P_\star`$ lorsque
$`\max_x\mu_x\le1/2`$. Dans le secteur restant, réordonnons les états de sorte
que $`\mu_0>1/2`$, posons $`g_x=\mu_xf_x`$ et éliminons
$`g_0=-g_1-g_2-g_3`$. Avec

```math
d_i:=\frac{a_\star-c_{q_\star}(\mu_i)}{\mu_i},
\qquad
D_k:=(\mu_0+\mu_k)(1-\mu_0-\mu_k),
```

la différence des deux formes quadratiques est exactement $`g^\mathsf TMg`$,
où

```math
\boxed{
M
=d_0\mathbf1\mathbf1^\mathsf T
+\mathrm{diag}(d_1,d_2,d_3)
+s_\star\sum_{k=1}^3\frac{v_kv_k^\mathsf T}{D_k},
}
```

```math
v_1=(0,1,1)^\mathsf T,
\qquad
v_2=(1,0,1)^\mathsf T,
\qquad
v_3=(1,1,0)^\mathsf T.
```

Ainsi $`P_\star`$ est réduit à la positivité d'une matrice rationnelle
$`3\times3`$ sur le simplexe $`\mu_0>1/2`$. Les calibrations binaire et
uniforme sont exactes. Le contre-audit à graine fixe vérifie tous les mineurs
principaux de $M$ pour les a priori polarisés échantillonnés et n'a pas trouvé
de contre-exemple ; cet audit fini n'est pas une preuve uniforme.

Il n'est pas nécessaire de prouver séparément le lemme pour chaque
$`q<q_\star`$. En effet, $`Y_q`$ s'obtient à partir de $`Y_{q_\star}`$ en
multipliant chacune des trois sorties par un bruit BSC indépendant de moyenne
$`q/q_\star`$. C'est une dégradation explicite.

**Théorème conditionnel.** Si $`P_\star`$ est vrai, la transitivité de l'ordre
less-noisy et le théorème de comparaison de Polyanskiy--Wu permettent de
remplacer chaque triangle observé, pour tout $`q\le q_\star`$, par le même
canal critique $`E_{a_\star,s_\star,e_\star}`$. Les hypothèses de
Chayes--Lei vérifiées ci-dessus et $`a_\star=e_\star`$ impliquent l'absence de
composante infinie ; la probabilité de relier deux points dont la distance
tend vers l'infini s'annule. On obtiendrait donc

```math
p<p_\star^{\mathrm{cond}}
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
```

Ni $`P_\star`$ ni cette conclusion ne sont actuellement classés « établis ».

## 7. Où interviennent les $`\Lambda_v`$ ancestraux ?

Pour le heat bath au LCA $u$, les slides 31--33 donnent les quatre log-poids

```math
\ell_u^{ab}
=
\log\mu_0(\sigma^{ab})
+
\sum_{v\succeq u}
\left[
\log\Lambda_v(\sigma^{ab})
+(1-\beta_v)\Lambda_v(\sigma^{ab})
\right],
```

et

```math
\pi_u^{ab}
=
\frac{e^{\ell_u^{ab}}}
{\sum_{c,d}e^{\ell_u^{cd}}}.
```

Le fichier 10 estime précisément tous les
$`\Lambda_v(\sigma^{ab})`$. Le nouvel objectif utile est donc de contrôler la
polarisation

```math
M_u:=\max_{a,b}\pi_u^{ab},
```

et non seulement un odds binaire ou la moyenne de $`\eta_u`$.

Le critère exact $`M_u\le1/2`$ est

```math
e^{\ell_u^{\widehat a\widehat b}}
\le
\sum_{(a,b)\ne(\widehat a,\widehat b)}e^{\ell_u^{ab}},
```

où $`(\widehat a,\widehat b)`$ maximise $`\ell_u^{ab}`$. Un certificat simple
est

```math
\max_{a,b}\ell_u^{ab}-\min_{a,b}\ell_u^{ab}\le\log3.
```

Dans les coordonnées de Walsh du fichier 10, une condition suffisante est

```math
|h_{u,1}^{\mathrm{tot}}|
+|h_{u,2}^{\mathrm{tot}}|
+|J_u^{\mathrm{tot}}|
\le\frac{\log3}{2}.
```

Les moments exacts, les bornes de Bernstein et le certificat de queue
$`\mathcal R_u`$ peuvent maintenant servir à borner cet événement. Mais deux
garde-fous subsistent.

1. Une fusion critique favorable peut être fortement polarisée ; on ne peut
   pas supposer $`M_u\le1/2`$ gratuitement.
2. L'a priori $\mu$ du canal facteur physique et la loi de flip
   $`\pi_u^{ab}`$ du heat bath LCA ne sont pas automatiquement le même objet.
   Un argument de marginalisation de $D$ ou d'entrelacement des deux canaux est
   nécessaire avant d'utiliser une borne sur les $`\Lambda_v`$ dans le
   théorème facteur.

La voie hiérarchique peut donc produire un gain, mais seulement si elle ferme
l'un des deux énoncés suivants :

- le lemme $`P_\star`$ sur le secteur polarisé ;
- une domination adaptée qui paie explicitement la probabilité et la
  géométrie des nœuds où $`M_u>1/2`$.

## 8. Statut et ordre de travail

| Objet | Statut | Conséquence |
|---|---|---|
| Profil $`c_q(t)`$ | Établi, énumération indépendante | contraction exacte sous tout a priori fixé |
| SDPI globale $`2q^2/(1+q^2)`$ | Établi | le triangle scalaire ne bat pas la baseline |
| Borne des trois projections | Établi | inégalité $`\chi^2`$ pour $`\max\mu\le1/2`$ |
| Enveloppe $`a(q),s(q),e(q)`$ | Établie comme optimisation locale | candidat algébrique $`0.809909\ldots`$ |
| Lemme $`P_\star`$, masse dominante | À prouver | positivité rationnelle $`3\times3`$ |
| Passage des $`\pi_u^{ab}`$ LCA au canal facteur | À prouver | verrou propre aux $`\Lambda_v`$ |
| Nouvelle borne weak recovery | Conditionnelle | ne pas citer $`0.809909\ldots`$ comme théorème |

Ordre recommandé : prouver ou réfuter $`P_\star`$ en certifiant les mineurs de
la matrice $M$ sur le simplexe polarisé. Le critère de Makur--Polyanskiy donne
alors directement la comparaison less-noisy. Seulement ensuite appliquer
Chayes--Lei et revenir à la loi des $`\Lambda_v`$ le long d'une paire
critique.

## 9. Vérification reproductible

Le module
[computations/triangle_block_sdpi.py](computations/triangle_block_sdpi.py)
calcule les deux contractions, le profil $`c_q`$, l'enveloppe affine et la
racine algébrique. Les tests
[computations/test_triangle_block_sdpi.py](computations/test_triangle_block_sdpi.py)
contre-auditent la formule par énumération directe des $`4\times8`$
probabilités, vérifient le défaut de la droite naïve, la borne des trois
projections et les conditions de Chayes--Lei. Ils cherchent aussi des
contre-exemples à $`P_\star`$ dans le secteur polarisé en testant les mineurs
principaux de la matrice exacte $`3\times3`$ ; ce dernier test est uniquement
diagnostique.

Sources primaires :

- [Polyanskiy--Wu, information-percolation multi-terminale](https://arxiv.org/abs/1806.04195) ;
- [Abbe--Boix, borne $`\chi^2`$ sur les graphes](https://arxiv.org/abs/1806.03227) ;
- [Makur--Polyanskiy, critère less-noisy par $`\chi^2`$](https://arxiv.org/abs/1609.06877) ;
- [Chayes--Lei, modèle random-cluster triangulaire](https://arxiv.org/abs/cond-mat/0508254).
