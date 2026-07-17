# Réduction hiérarchique au cas favorable critique

Ce fichier formalise la voie prioritaire pour la weak recovery. On suit deux
sommets lointains $`i,j`$ dans le dendrogramme d'horloges exponentielles. Le
cas géométriquement le plus favorable est celui où ils appartiennent au même
arbre à la coupe $`1`$ et où, en descendant cet arbre, leurs branches se
séparent au seuil de percolation. De façon équivalente, leur LCA $`u`$ est une
fusion de Kruskal au temps $`\beta_u\simeq\beta_c`$.

La conclusion rigoureuse comporte trois niveaux.

1. **Établi, volume fini.** Les paires proches sont négligeables et les paires
   dans deux racines distinctes ont un score LCA nul. La borne de weak
   recovery se réduit donc aux paires lointaines du même arbre.
2. **Établi sous un lemme de domination explicite.** Si l'expérience LCA
   postcritique est moins informative que l'oracle où la séparation a lieu au
   seuil, alors toute la contribution postcritique est majorée par cet oracle
   favorable. Il n'est pas nécessaire que les temps LCA réels se concentrent
   en $`\beta_c`$.
3. **Problème central.** Dans cet oracle, la fiabilité dépend du log-rapport
   complet

   ```math
   L_u
   =
   \ell_u^{\mathrm{crit}}+B_u,
   ```

   où $`B_u`$ utilise les quatre taux
   $`\Lambda_v(\sigma^{ab})`$ de **chaque** ancêtre $`v\succ u`$. Poser
   $`B_u=0`$ ne traite donc pas la dynamique hiérarchique des slides 31--33.

Le résultat recherché est une impossibilité : si même cette expérience
oracle favorable contracte la parité, et si elle domine les expériences
postcritiques réelles, alors la weak recovery est impossible. Le succès de
l'oracle ne suffit en revanche pas à construire un estimateur non oracle.

## 1. Cadre fini et expérience favorable

Soit $`G_L=(V_L,E_L)`$ une exhaustion finie, $`n_L=|V_L|`$, et supposons
l'a priori binaire i.i.d. uniforme. Le dendrogramme $`D_L`$ est construit avec
les horloges exponentielles censurées à la coupe $`1`$. Pour deux sommets,
posons

```math
\beta_{ij}
=
\inf\{t:i\leftrightarrow j\text{ dans la coupe }t\},
```

avec $`\beta_{ij}=+\infty`$ lorsque $`i,j`$ appartiennent à deux racines
distinctes et la convention $`\beta_{ii}=0`$. Lorsque $`i\ne j`$ et
$`\beta_{ij}\le1`$, leur LCA est noté $`u_{ij}`$.

Soient $`I_L,J_L`$ indépendants et uniformes dans $`V_L`$. Le score étendu
est

```math
\eta_{I_LJ_L}^{\mathrm{LCA}}
=
\begin{cases}
1,&I_L=J_L,\\
\tanh^2(L_{u_{I_LJ_L}}/2),&I_L\ne J_L,\ \beta_{I_LJ_L}\le1,\\
0,&\beta_{I_LJ_L}=+\infty.
\end{cases}
```

La projection au LCA donne

```math
Q_L
\le
\mathbb E\left[\eta_{I_LJ_L}^{\mathrm{LCA}}\right].
```

Choisissons $`r_L\to\infty`$ et posons

```math
b_L
:=
\mathbb P\bigl(d(I_L,J_L)<r_L\bigr).
```

Sur toute exhaustion de dimension fixe pour laquelle $`r_L`$ croît plus
lentement que le diamètre, on peut choisir $`r_L`$ de sorte que $`b_L\to0`$.
Cette propriété géométrique doit être vérifiée pour les conditions de bord
retenues.

Pour $`\varepsilon>0`$, notons

```math
a_\varepsilon
=
\max(0,\beta_c-\varepsilon),
\qquad
c_\varepsilon
=
\min(1,\beta_c+\varepsilon)
```

On suppose ici $`\beta_c\in[0,1]`$. Dans le GSBM triangulaire homogène, cela
équivaut à $`p\ge p_{\mathrm{SW}}`$ ; en dessous, le seuil géométrique n'est
pas atteint avant la censure et l'expérience favorable critique est vide.

Définissons alors l'événement favorable

```math
\mathcal F_{L,\varepsilon}
=
\left\{
d(I_L,J_L)\ge r_L,
\ a_\varepsilon\le\beta_{I_LJ_L}\le c_\varepsilon
\right\}.
```

La condition $`\beta_{I_LJ_L}\le c_\varepsilon\le1`$ implique déjà que les
deux sommets sont dans le même arbre. Comme les temps sont continus, une
fusion exactement à $`\beta_c`$ a probabilité nulle en volume fini. Le cas
« au seuil » désigne donc cette fenêtre suivie de la limite
$`L\to\infty`$, puis $`\varepsilon\downarrow0`$, ou une désintégration de Palm
équivalente.

## 2. Réduction exacte : loin, même arbre, puis temps critique

Rappelons le second moment géométrique

```math
S_L(t)
=
\mathbb P(\beta_{I_LJ_L}\le t)
=
\frac1{n_L^2}
\mathbb E\sum_{C\in\Pi_t}|C|^2.
```

Posons

```math
\rho_{L,\varepsilon}^{\mathrm{fav}}
=
\mathbb P(\mathcal F_{L,\varepsilon}),
```

```math
\Gamma_{L,\varepsilon}^{\mathrm{fav}}
=
\mathbb E\left[
\eta_{I_LJ_L}^{\mathrm{LCA}}
\mid\mathcal F_{L,\varepsilon}
\right]
```

lorsque la fenêtre a une masse non nulle, et

```math
\Delta_{L,\varepsilon}^{\mathrm{late}}
=
\mathbb E\left[
\eta_{I_LJ_L}^{\mathrm{LCA}}
\mathbf1_{\{
d(I_L,J_L)\ge r_L,
\ c_\varepsilon<\beta_{I_LJ_L}\le1
\}}
\right].
```

### Proposition 2.1 — décomposition favorable, statut : établi

Sous les hypothèses de la projection LCA,

```math
\boxed{
Q_L
\le
b_L
+S_L(a_\varepsilon)
+\rho_{L,\varepsilon}^{\mathrm{fav}}
 \Gamma_{L,\varepsilon}^{\mathrm{fav}}
+\Delta_{L,\varepsilon}^{\mathrm{late}}.
}
```

### Preuve

On partitionne les paires en quatre classes : distance inférieure à $`r_L`$,
distance supérieure avec fusion avant $`a_\varepsilon`$, fusion dans la
fenêtre, puis fusion après $`c_\varepsilon`$. Les deux premières contributions
sont majorées respectivement par $`b_L`$ et $`S_L(a_\varepsilon)`$, car
$`0\le\eta^{\mathrm{LCA}}\le1`$. La contribution de la fenêtre est exactement
$`\rho^{\mathrm{fav}}\Gamma^{\mathrm{fav}}`$. Les racines distinctes ont
$`\eta^{\mathrm{LCA}}=0`$ par définition. Cela donne l'inégalité annoncée.

Cette proposition justifie sans hypothèse supplémentaire les deux premières
restrictions demandées : les sommets peuvent être pris lointains, et seuls
ceux du même arbre contribuent. Elle montre aussi exactement ce qu'il faut
ajouter pour remplacer toutes les fusions postcritiques par une séparation au
seuil.

## 3. Le lemme « cas le plus favorable »

Écrivons

```math
f(x)=\tanh^2(x/2).
```

Sous la loi de la paire lointaine conditionnée par
$`a_\varepsilon\le\beta_{I_LJ_L}\le1`$, notons $`L_L^{\mathrm{post}}`$ le
log-rapport hiérarchique complet. Sous la loi favorable
$`\mathcal F_{L,\varepsilon}`$, notons $`L_L^{\mathrm c}`$ le même log-rapport.
Ces deux variables comprennent le bucket du LCA et tous ses ancêtres.

### Hypothèse HF — domination hiérarchique favorable

Il existe un couplage de $`L_L^{\mathrm{post}}`$ et $`L_L^{\mathrm c}`$ tel
que, pour des erreurs $`\zeta_{L,\varepsilon},\delta_{L,\varepsilon}\ge0`$,

```math
\mathbb P\left(
|L_L^{\mathrm{post}}|
>
|L_L^{\mathrm c}|+\zeta_{L,\varepsilon}
\right)
\le
\delta_{L,\varepsilon}.
```

La version exacte du « cas le plus favorable » correspond à
$`\zeta_{L,\varepsilon}=\delta_{L,\varepsilon}=0`$. Cette formulation est plus
forte qu'une comparaison des moyennes de $`\Lambda_v`$ et plus faible qu'une
domination point par point de tous les squelettes. Elle vise directement la
quantité de weak recovery.

La dérivée de $`f`$ vérifie

```math
\sup_{x\in\mathbb R}|f'(x)|
=
\frac{2}{3\sqrt3}
=:
c_{\mathrm{rel}}.
```

### Théorème 3.1 — réduction à l'oracle critique, statut : établi sous HF

Sous HF,

```math
\boxed{
Q_L
\le
b_L
+S_L(a_\varepsilon)
+\Gamma_{L,\varepsilon}^{\mathrm{fav}}
+c_{\mathrm{rel}}\zeta_{L,\varepsilon}
+\delta_{L,\varepsilon}.
}
```

En particulier, si l'on peut choisir les limites de sorte que

```math
b_L+S_L(a_\varepsilon)
+\Gamma_{L,\varepsilon}^{\mathrm{fav}}
+\zeta_{L,\varepsilon}
+\delta_{L,\varepsilon}
\longrightarrow0,
```

alors $`Q_L\to0`$ et la weak recovery est impossible.

### Preuve

Sur l'événement de bon couplage, $`f`$ est croissante en $`|x|`$ et
$`c_{\mathrm{rel}}`$-Lipschitz, donc

```math
f(L_L^{\mathrm{post}})
\le
f(L_L^{\mathrm c})
+c_{\mathrm{rel}}\zeta_{L,\varepsilon}.
```

Sur l'événement exceptionnel, on utilise $`f\le1`$. La fiabilité moyenne de
toutes les paires lointaines postcritiques est donc au plus
$`\Gamma_{L,\varepsilon}^{\mathrm{fav}}
+c_{\mathrm{rel}}\zeta_{L,\varepsilon}
+\delta_{L,\varepsilon}`$. Sa masse est au plus $`1`$. On ajoute les paires
proches et les fusions sous-critiques comme dans la proposition 2.1.

Ce théorème donne le sens précis de l'expérience favorable : on ne prétend
pas que tous les LCA réels sont critiques ; on démontre qu'on ne perd rien
pour une borne d'impossibilité en les remplaçant par l'expérience critique la
plus informative.

## 4. Le log-rapport exact de la dynamique hiérarchique

Pour $`u=\mathrm{LCA}(i,j)`$, retournons indépendamment les deux fils par
$`a,b\in\{0,1\}`$. Les slides 31--33 donnent les quatre poids

```math
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.
```

Avec

```math
\phi_v(x)=\log x+(1-\beta_v)x
```

pour $`x>0`$, posons

```math
\Phi_u^{ab}
=
\log\mu_0(\sigma^{ab})
+\sum_{v\succ u}\phi_v(\Lambda_v^{ab}).
```

Alors

```math
B_u
=
\mathrm{LSE}(\Phi_u^{00},\Phi_u^{11})
-\mathrm{LSE}(\Phi_u^{01},\Phi_u^{10})
```

et

```math
\boxed{
L_u
=
B_u
+\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u).
}
```

Cette écriture logarithmique vaut directement lorsque les deux taux locaux
sont strictement positifs. Si l'un des quatre poids s'annule, on conserve les
$`q_u^{ab}`$ et les deux `log-sum-exp` dans les réels étendus ; la fiabilité
reste bien définie dans $`[0,1]`$.

La fiabilité favorable est donc exactement

```math
\boxed{
\Gamma_{L,\varepsilon}^{\mathrm{fav}}
=
\mathbb E_{L,\varepsilon}^{\star}
\left[
\tanh^2\left(
\frac{\ell_u^{\mathrm{crit}}+B_u}{2}
\right)
\right],
}
```

où

```math
\ell_u^{\mathrm{crit}}
=
\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u)
```

et $`\mathbb E_{L,\varepsilon}^{\star}`$ est la loi biaisée par le nombre de
paires lointaines séparées par $`u`$. Cette identité est le cœur de la voie
hiérarchique. Le calcul local $`B_u=0`$ n'est qu'une calibration.

## 5. Ce qu'il faut estimer pour chaque ancêtre

Fixons un ancêtre $`v\succ u`$. Son bucket se partage en trois groupes : les
arêtes incidentes à la partie du fils ancestral hors de $`C_u`$, celles
incidentes au premier fils $`C_1`$ de $`u`$, et celles incidentes au second
fils $`C_2`$. Avec

```math
T_{v,r}
=
\sum_{e\in E_v^{(r)}}|W_e|,
\qquad
\lambda_{v,r}
=
\sum_{e\in E_v^{(r)}}
|W_e|\mathbf1_{\{e\text{ satisfaite}\}},
```

et $`X_{v,r}=2\lambda_{v,r}-T_{v,r}`$, on a l'identité déterministe

```math
\boxed{
\Lambda_v^{ab}
=
\frac12\left[
T_{v,0}+T_{v,1}+T_{v,2}
+X_{v,0}
+(-1)^aX_{v,1}
+(-1)^bX_{v,2}
\right].
}
```

Il faut donc estimer les trois déséquilibres, et pas seulement
$`\Lambda_v^{00}`$.

Dans le GSBM triangulaire homogène, conditionnellement au squelette non marqué,
notons $`m_{v,r}=|E_v^{(r)}|`$, $`m_v=\sum_rm_{v,r}`$ et

```math
s_v
=
\frac1{1+e^{-u_p(1-\beta_v)}},
\qquad
u_p=\log\frac p{1-p}.
```

Si $`G_v`$ est le groupe de l'arête gagnante, alors

```math
\mathbb P(G_v=r\mid\mathscr S_u)
=
\frac{m_{v,r}}{m_v}
```

et, conditionnellement à $`G_v`$, les trois comptes satisfaits sont
indépendants avec

```math
K_{v,r}
\stackrel d=
\mathbf1_{\{G_v=r\}}
+\mathrm{Bin}\left(
m_{v,r}-\mathbf1_{\{G_v=r\}},s_v
\right).
```

On obtient $`\lambda_{v,r}=u_pK_{v,r}`$ puis les quatre
$`\Lambda_v^{ab}`$ par la formule affine. Pour un nombre fini d'ancêtres,
c'est un noyau exact et calculable par somme finie. Le verrou asymptotique est
la loi jointe

```math
\left(
\beta_v,m_{v,0},m_{v,1},m_{v,2}
\right)_{v\succ u}
```

sous le biais d'une paire lointaine critique.

## 6. Deux lemmes de fermeture pour la chaîne des $`\Lambda_v`$

Ordonnons les ancêtres $`u=v_0\prec v_1\prec\cdots`$. Soit $`B_u^{(K)}`$ le
message obtenu en conservant les $`K`$ premiers ancêtres et

```math
L_u^{(K)}
=
\ell_u^{\mathrm{crit}}+B_u^{(K)},
\qquad
\eta_u^{(K)}
=
\tanh^2(L_u^{(K)}/2).
```

Le [certificat de queue](10_ANCESTRAL_LAMBDA_ESTIMATION.md) fournit une
fonctionnelle $`\mathcal R_u^{(>K)}`$ des quatre coins ancestraux telle que

```math
|B_u-B_u^{(K)}|
\le
\mathcal R_u^{(>K)}.
```

### Lemme 6.1 — transport exact de la troncature, statut : établi

On a point par point

```math
\boxed{
|\eta_u-\eta_u^{(K)}|
\le
\min\left(
1,
\frac{2}{3\sqrt3}\mathcal R_u^{(>K)}
\right).
}
```

Par conséquent,

```math
\left|
\Gamma_{L,\varepsilon}^{\mathrm{fav}}
-
\mathbb E_{L,\varepsilon}^{\star}[\eta_u^{(K)}]
\right|
\le
\mathbb E_{L,\varepsilon}^{\star}
\left[
\min\left(
1,
\frac{2\mathcal R_u^{(>K)}}{3\sqrt3}
\right)
\right].
```

Ce lemme transforme la sommabilité de tous les $`\Lambda_v`$ éloignés en une
erreur certifiée sur la quantité de weak recovery.

### Lemme 6.2 — critère quadratique du message, statut : établi

Lorsque $`L_u`$ est fini presque sûrement, $`|\tanh x|\le|x|`$ donne

```math
\boxed{
\Gamma_{L,\varepsilon}^{\mathrm{fav}}
\le
\frac14
\mathbb E_{L,\varepsilon}^{\star}[L_u^2].
}
```

Il suffit donc, pour une borne d'impossibilité dans l'oracle favorable, de
montrer que le message hiérarchique total tend vers zéro dans $`L^2`$. Une
version en probabilité suffit également par convergence dominée, puisque
$`0\le\eta_u\le1`$.

Les lemmes 6.1 et 6.2 donnent un programme fermé : calculer exactement les
$`K`$ premiers ancêtres, contrôler leur message total, puis envoyer
$`K\to\infty`$ avec le certificat de queue.

### Lemme 6.3 — test de non-contraction, statut : établi

Point par point, l'inégalité triangulaire inverse donne

```math
\boxed{
\eta_u
\ge
\tanh^2\left(
\frac{\bigl||\ell_u^{\mathrm{crit}}|-|B_u|\bigr|}{2}
\right).
}
```

En outre, si $`\Gamma_{L,\varepsilon}^{\mathrm{fav}}\to0`$, alors

```math
B_u+\ell_u^{\mathrm{crit}}
\longrightarrow0
```

en probabilité sous la loi favorable. En effet, pour tout $`x>0`$,

```math
\mathbb P_{L,\varepsilon}^{\star}(|L_u|>x)
\le
\frac{\Gamma_{L,\varepsilon}^{\mathrm{fav}}}
{\tanh^2(x/2)}.
```

Ce lemme est un contre-audit décisif. Si le message local critique diverge,
une contraction favorable ne peut avoir lieu que si le message construit par
les $`\Lambda_v`$ ancestraux le compense avec la même échelle et le signe
opposé, à une erreur $`o_{\mathbb P}(1)`$. Une simple borne $`B_u=O(1)`$ ne
suffirait pas.

## 7. Conséquence recherchée sur le GSBM triangulaire

Fixons $`p_\star`$. Pour démontrer l'absence de weak recovery pour tout
$`p<p_\star`$ par cette voie, il suffit d'établir uniformément dans cet
intervalle les quatre entrées suivantes.

1. **Sous-criticité géométrique.** Pour tout $`\varepsilon>0`$ fixé,
   $`S_L(\beta_c(p)-\varepsilon)\to0`$ avec l'exhaustion et les conditions de
   bord choisies.
2. **Domination favorable HF.** Toute paire lointaine du même arbre qui se
   sépare après le seuil est, pour la parité, moins informative que la paire
   oracle séparée au seuil, à une erreur qui tend vers zéro.
3. **Limite ancestrale.** La loi des premiers groupes
   $`(m_{v,0},m_{v,1},m_{v,2},\beta_v)`$ converge sous le biais de la paire
   critique, les quatre coins proches de zéro sont traités exactement et la
   queue $`\mathcal R_u^{(>K)}`$ est sommable.
4. **Contraction critique.** Le calcul fini obtenu satisfait

   ```math
   \lim_{\varepsilon\downarrow0}
   \limsup_{L\to\infty}
   \Gamma_{L,\varepsilon}^{\mathrm{fav}}(p)
   =0.
   ```

Le théorème 3.1 donnerait alors l'impossibilité pour $`p<p_\star`$. Une borne
strictement meilleure que l'information-percolation correspondrait à

```math
p_\star
>
\frac{1+\sqrt{2\sin(\pi/18)}}2
=
0.794659\ldots.
```

Cette possibilité n'est pas encore démontrée. Le calcul local du fichier 09
constitue un contre-audit important : si $`B_u=0`$ et si le bucket critique a
$`m\to\infty`$ arêtes, alors sa fiabilité tend vers $`1`$ pour tout
$`p>p_{\mathrm{SW}}`$ fixé. Une amélioration ne peut donc provenir du seul
bucket $`u`$. Elle exige le calcul de la loi réelle des tailles critiques et du
message formé par tous les $`\Lambda_v`$ ancestraux. Par le lemme 6.3, si la
coupe critique est grande, ce message doit compenser le LLR local et pas
seulement rester borné. S'il renforce le signal local au lieu de le contracter,
cette mise à jour LCA à un pas ne pourra pas améliorer la borne ; ce résultat
négatif serait lui-même exact.

## 8. Comment prouver HF avec les $`\Lambda_v`$

HF ne doit pas être attaquée par une monotonie informelle en $`\beta_u`$.
L'objet à comparer est le vecteur quatre états complet. Une stratégie
certifiable est la suivante.

1. Tronquer les deux chaînes à $`K`$ ancêtres et conserver exactement les
   configurations où un taux vaut zéro.
2. Coupler les squelettes postcritique et critique au niveau des vecteurs
   $`(\beta_v,m_{v,0},m_{v,1},m_{v,2})_{v\le K}`$.
3. Utiliser le même aléa de gagnante et les mêmes uniformes pour les marques
   résiduelles, puis transporter ce couplage vers les quatre
   $`\Lambda_v^{ab}`$.
4. Calculer les deux log-rapports par `log-sum-exp` et certifier
   $`|L^{\mathrm{post}}|\le|L^{\mathrm c}|+\zeta_K`$ par arithmétique
   d'intervalles.
5. Ajouter les deux queues avec le lemme 6.1 et faire tendre $`K`$ vers
   l'infini.

Sur un cactus de triangles, l'étape 2 relève d'une récurrence finie. Sur une
bande de largeur fixée, elle relève d'une matrice de transfert. Sur la grille
entière, elle demande des estimations multi-bras pour les trois interfaces
ancestrales. Ce sont ces objets géométriques, et non un canal de triangle
isolé, qui déterminent la voie hiérarchique demandée.

## 9. Contre-audits obligatoires

1. **Conditionnement contre domination.** Conditionner sur
   $`\mathcal F_{L,\varepsilon}`$ ne donne pas automatiquement une borne
   supérieure sur les autres paires. C'est exactement le contenu de HF.
2. **Même arbre.** Cette restriction est rigoureuse pour la borne LCA parce
   que deux racines distinctes ont un score nul ; elle ne signifie pas qu'une
   composante géante suffit à la weak recovery.
3. **Temps exact.** En volume fini, l'événement
   $`\beta_u=\beta_c`$ a probabilité nulle. Toute preuve doit annoncer la
   fenêtre ou la mesure de Palm et l'ordre des limites.
4. **Tous les ancêtres.** Le facteur $`B_u`$ ne peut être supprimé. La
   criticité de $`u`$ n'implique pas celle de $`v\succ u`$.
5. **Quatre taux.** Une estimation de $`\Lambda_v^{00}`$ seule ne contrôle ni
   les flips impairs, ni le couplage de Walsh, ni les coins nuls.
6. **Non-linéarité.** Remplacer un taux par sa moyenne à l'intérieur de
   $`\log\Lambda_v`$ est invalide.
7. **Biais de paire.** La loi pertinente pondère un nœud par le nombre de
   paires lointaines séparées par ses deux fils. Un nœud uniforme du
   dendrogramme a une autre loi.
8. **Portée logique.** L'échec de l'oracle favorable, combiné à HF, prouve une
   impossibilité globale. Son succès ne prouve pas la weak recovery et ne
   fournit pas d'algorithme.

## 10. Statut des lemmes

| Élément | Statut | Rôle |
|---|---|---|
| Réduction aux paires lointaines du même arbre | Établi, volume fini | élimine les cas sans contribution LCA |
| Décomposition critique de la proposition 2.1 | Établi, volume fini | sépare géométrie, fenêtre et reste tardif |
| Formule de $`L_u`$ avec tous les $`\Lambda_v^{ab}`$ | Établi, volume fini | définit l'oracle hiérarchique exact |
| Noyau conditionnel des marques ancestrales | Établi, volume fini | calcule tout préfixe fini de la chaîne |
| Transport de la queue vers $`\eta_u`$ | Établi | donne une erreur certifiée sur la fiabilité |
| Nécessité d'une compensation $`B_u\simeq-\ell_u^{\mathrm{crit}}`$ si l'oracle contracte | Établi | contre-audit des grandes coupes critiques |
| Domination favorable HF | À prouver | autorise le remplacement de toutes les paires par le cas critique |
| Convergence du squelette critique et sommabilité | À prouver | ferme la limite ancestrale |
| Nouvelle borne triangulaire $`p_\star>0.794659\ldots`$ | Ouvert | conséquence seulement après fermeture des deux lignes précédentes |

## 11. Calibration autoduale de Nishimori

Le [fichier 13](13_NISHIMORI_HIERARCHICAL_CLOCKS.md) montre que l'équation de
face de Nishimori--Ohzeki est exactement

```math
H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3)=1\ \text{bit},
```

et qu'elle se réalise par une course conditionnelle de quatre horloges. Son
niveau zéro redonne $`0.835805792367\ldots`$. Le pont avec le présent théorème
favorable n'est toutefois pas automatique : il faut construire le défaut
autodual sur les mêmes buckets ancestraux que
$`\Gamma_{L,\varepsilon}^{\mathrm{fav}}`$, prouver que son signe contracte la
fiabilité collapsed après oubli de $D$, puis appliquer HF. Cette exigence est
le lemme NH3 ; elle empêche d'annoncer la racine de face comme une borne déjà
démontrée.
