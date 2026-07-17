# Cas test : GSBM homogène sur la grille triangulaire

Ce cas est le premier banc d'essai parce que plusieurs quantités de percolation sont explicites. Il faut toutefois fixer une exhaustion finie — idéalement un tore triangulaire — avant tout énoncé asymptotique rigoureux.

## 1. Modèle

On prend
```math
f_{\mathrm{in}}\equiv p,
\qquad
f_{\mathrm{out}}\equiv1-p,
\qquad
p\in[1/2,1).
```
Toutes les interactions observées ont le même module
```math
u_p=\log\frac p{1-p}.
```
Les calculs conditionnels à une fusion supposent $p>1/2$, donc $`u_p>0`$. Au point dégénéré $p=1/2$, toutes les horloges ont taux nul, le dendrogramme ne contient aucune fusion et les formules conditionnelles de type $0/0$ ne sont pas utilisées.

Sous la loi générative, conditionnellement à la vérité $\Sigma$ mais après moyenne sur l'observation $O$, le signe observé d'une arête est satisfait par $\Sigma$ avec probabilité $p$, indépendamment des autres arêtes.

## 2. Loi exacte de chaque coupe du dendrogramme

Une arête est présente à la coupe $t$ si elle est satisfaite et si son horloge est inférieure à $t$. Par conséquent,
```math
\boxed{
q_p(t)
=
\mathbb P(\xi_e\le t\mid\Sigma)
=
p(1-e^{-u_pt})
=
p\left[
1-\left(\frac{1-p}{p}\right)^t
\right].
}
```
Sous cette même loi annealed, conditionnellement à $\Sigma$, les indicatrices sont indépendantes. Ainsi $`\Pi_t`$ a exactement la loi d'une percolation indépendante par arêtes de paramètre $`q_p(t)`$, et les différentes coupes sont couplées de façon monotone par les mêmes horloges.

Cette identité ne vaut pas conditionnellement à une observation $O$ fixée. Elle vaut aussi lorsqu'on remplace la vérité par une réplique postérieure à l'équilibre, après moyenne jointe sur $O$, par l'identité de Nishimori.

En particulier,
```math
q_p(1)=2p-1,
```
ce qui redonne la percolation de Swendsen–Wang du chapitre 11.

Soit
```math
q_c^{\mathrm{bond}}(\mathbb T)=2\sin(\pi/18).
```

Le temps auquel la filtration brute atteint le seuil de percolation est

```math
\beta_c(p):=t_c(p)
=
\frac{
\log\left(1-q_c^{\mathrm{bond}}/p\right)
}{
\log\left((1-p)/p\right)
},
```

lorsque le membre de droite est défini. On a $`\beta_c(p)\le1`$ si et seulement si

```math
p\ge
\frac12+\sin(\pi/18)
=0.673648\ldots
```
La distribution des temps de coalescence est donc reliée aux fonctions de connexion de la percolation :
```math
\mathbb P_{\mathrm{ann}}(\beta_{ij}\le t\mid\Sigma)
=
\mathbb P_{q_p(t)}(i\leftrightarrow j).
```
Cette identité rend toute la **géométrie brute** du dendrogramme calculable. La nouvelle information doit venir des probabilités de flips, pas seulement de $`\Pi_t`$.

## 3. Trois bornes rigoureuses déjà disponibles

### Swendsen–Wang par arêtes
```math
p_c^{\mathrm{edge}}
=
\frac12+\sin(\pi/18)
=0.673648\ldots
```
En dessous, le graphe gelé est sous-critique et la borne du chapitre 11 interdit la weak recovery, après formalisation du passage en volume fini.

### Dynamique triangulaire d'ordre supérieur

Pour les triangles disjoints d'une couleur, posons
```math
\alpha_p=1-e^{-2u_p}
=
\frac{2p-1}{p^2}.
```
Après moyenne sur les observations, les états locaux de percolation ont les probabilités
```math
a=p(2p-1)
```
pour le triangle plein,
```math
s=(1-p)(2p-1)
```
pour chacun des trois états à une arête, et
```math
e=4(1-p)^2
```
pour le triangle vide. On vérifie $a+3s+e=1$.

La condition autoduale $a=e$ donne
```math
p_c^\triangle
=
\frac{7-\sqrt{17}}4
=0.719224\ldots
```
Le théorème de Chayes–Lei demande aussi
```math
ae\ge2s^2,
\qquad
a+e>\frac{2\sqrt2}{3+2\sqrt2},
```
ainsi que l'isotropie et l'indépendance entre les triangles choisis. Ces hypothèses doivent accompagner l'utilisation du seuil.

L'intervalle de stricte amélioration sur la dynamique par arêtes est
```math
\left(
\frac12+\sin(\pi/18),
\frac{7-\sqrt{17}}4
\right).
```
La borne gauche est ouverte : au point critique bidimensionnel, il n'y a pas de composante de densité positive sous les conditions de bord usuelles.

### Information-percolation

La contraction $\chi^2$ du canal d'une arête vaut
```math
\eta=(2p-1)^2.
```
La sous-criticité de la percolation d'information donne
```math
p<
p_c^{\mathrm{info}}
:=
\frac{1+\sqrt{2\sin(\pi/18)}}2
=0.794659\ldots
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
```
Cette borne est la baseline rigoureuse à dépasser.

## 4. Seuil informationnel attendu

Après jauge par la vérité, le modèle est relié à l'Ising $\pm J$ sur la ligne de Nishimori. Les calculs de dualité/répliques et les simulations de la littérature situent le point multicritique triangulaire vers
```math
p_{\mathrm N}=0.8358058\ldots
```
Cette valeur est **conjecturale**, pas un seuil rigoureux à utiliser dans une preuve.

Le paysage de travail est donc :

| Méthode | Seuil $p$ | Statut |
|---|---:|---|
| FK / Swendsen–Wang par arêtes | $0.673648\ldots$ | Rigoureux |
| Dynamique triangulaire | $0.719224\ldots$ | Rigoureux sous les hypothèses de Chayes–Lei |
| Information-percolation $\chi^2$ | $0.794659\ldots$ | Rigoureux |
| Point de Nishimori triangulaire | $0.8358058\ldots$ | Conjecture + numérique |

## 5. Poids exact d'une fusion locale

Considérons une coupe $`E_u`$ contenant $m$ arêtes, dont $k$ sont satisfaites par la configuration courante. Alors
```math
T_u=mu_p,
\qquad
\Lambda_u=ku_p.
```
Si l'a priori et les facteurs ancêtres s'annulent dans le rapport de parité, la log-vraisemblance locale vaut
```math
\boxed{
L_{m,k,\beta}^{\mathrm{loc}}
=
\log\frac{k}{m-k}
+(1-\beta)u_p(2k-m).
}
```
La fiabilité associée est
```math
\rho_{m,k,\beta}
=
\left|\tanh\frac{L_{m,k,\beta}^{\mathrm{loc}}}{2}\right|.
```
Pour une coupe fixée ayant $k\ge1$ arêtes satisfaites,
```math
\beta=\min_{1\le j\le k}\xi_j
\sim\mathrm{Exp}(ku_p).
```
Conditionnellement à une fusion avant $1$,
```math
f(\beta\mid\beta\le1,k)
=
\frac{ku_pe^{-ku_p\beta}}{1-e^{-ku_p}},
\qquad 0\le\beta\le1.
```
Cette loi ne peut pas être appliquée naïvement à une coupe **sélectionnée par Kruskal** : la sélection introduit un biais. Les premiers calculs doivent donc distinguer :

1. une coupe déterministe ;
2. une fusion du dendrogramme conditionnée par son passé ;
3. une fusion typique sous la loi stationnaire.

## 6. Identité à deux répliques

Dans ce modèle,
```math
Q_n
=
\frac1{n^2}
\sum_{i,j}
\mathbb E\left[
\left\langle\sigma_i\sigma_j\right\rangle_O^2
\right].
```
Sur la ligne de Nishimori, des identités supplémentaires relient les corrélations plantées et les corrélations de répliques. La stratégie correcte est de contrôler cette somme, directement ou via $`H_S`$, plutôt que de chercher seulement une composante géante dans $`\Pi_1`$.

## 7. Programme LCA explicite

Pour $`u=C_1\mathbin{\dot\cup}C_2`$, le poids informationnel de la fusion est
```math
\eta_u
=
\tanh^2\left[
\frac12\left\{
B_u
+\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u)
\right\}
\right],
```
où $`T_u=|E_u|u_p`$, $`\Lambda_u`$ est le poids total des liens satisfaits entre les deux fils et $`B_u`$ est le message exact de l'a priori et des ancêtres. La borne à calculer est
```math
H_n^{\mathrm{LCA}}
=
\frac1{n^2}\mathbb E\left[
n+2\sum_{u\in D}|C_{u,1}||C_{u,2}|\eta_u
\right].
```
Elle vérifie
```math
Q_n\le H_n^{\mathrm{LCA}}
\le
\frac1{n^2}\mathbb E\sum_{R\text{ racine}}|R|^2.
```
### Géométrie du temps de fusion

Pour une paire $(i,j)$,
```math
\beta_{ij}
=
\min_{\gamma:i\leadsto j}\max_{e\in\gamma}\xi_e,
\qquad
\mathbb P(\beta_{ij}\le t)
=
\tau_{ij}(q_p(t)),
```
où $`\tau_{ij}(q)=\mathbb P_q(i\leftrightarrow j)`$. Sur un graphe fini, la densité de $`\beta_{ij}`$ peut donc s'écrire par dérivation de la fonction de connexion, ou par une somme d'événements pivotaux. C'est le point d'entrée vers les outils de percolation proche-critique et de minimum spanning forest.

### Loi des liens encore fermés au temps $t$

Conditionnellement au fait qu'un lien ne soit pas encore ouvert juste avant $t$,
```math
s_p(t)
:=
\mathbb P(\text{lien satisfait}\mid\text{fermé à }t^-)
=
\frac{pe^{-u_pt}}{1-p+pe^{-u_pt}}
=
\mathrm{logistic}(u_p(1-t)).
```
Pour une coupe conditionnée de $m$ liens fusionnant à $t$, l'arête ouvrante est satisfaite et
```math
k\stackrel d=1+\mathrm{Bin}(m-1,s_p(t)).
```
Cette identité est annealed. Plus précisément, elle reste exacte conditionnellement au squelette de Kruskal non marqué, au bucket $`E_u`$ et à son temps : la géométrie aléatoire de $`E_u`$ reste biaisée par le passé, mais les marques des $m-1$ arêtes non gagnantes sont alors indépendantes de paramètre $`s_p(t)`$. La version groupée simultanée sur tous les ancêtres est démontrée dans [08_ANCESTRAL_LAMBDA_CHAIN.md](08_ANCESTRAL_LAMBDA_CHAIN.md).

Au temps $t=1$, dans le modèle local $`B_u=0`$,

```math
\mathbb E[\eta_u\mid m,\beta_u=1]=\frac1m.
```

Pour $p>1/2$ et $t<1$ fixés, la même fiabilité locale tend au contraire vers $1$ lorsque $m\to\infty$. Les grandes coupes ne randomisent donc fortement la parité que dans une fenêtre proche de la coupure $1$. Le [calcul critique exact](09_CRITICAL_MERGER_ORACLE.md) donne désormais une borne exponentielle et identifie la fenêtre $`p-p_{\mathrm{SW}}\asymp m^{-1/2}`$. Comme les coalescences longue portée apparaissent autour de $`\beta_c(p)<1`$ dès que la percolation est supercritique, un seul pas LCA reste trop proche de la borne percolative. Il faut alors étudier la masse des paires, les itérations avec rafraîchissement de $D$ et la loi complète du message $`B_u`$.

### Flux de susceptibilité informationnel

Chaque ouverture qui joint deux composantes $A$ et $B$ augmente le second moment géométrique de $2|A||B|$. Le score LCA remplace cet accroissement par
```math
2|A||B|\,\eta_u.
```

Le calcul recherché sur la grille est donc une mesure de fusions, indexée par $(t,m)$, multipliée par la fiabilité conditionnelle exacte. Cette décomposition sépare :

1. la fréquence et la taille des fusions, problème de percolation/Kruskal ;
2. le canal de parité à la fusion, problème de LLR et de messages d'ancêtres.

### Coupe critique et temps informationnel

Abrégeons

```math
q_c=2\sin(\pi/18).
```

Le temps critique géométrique est

```math
\beta_c(p)
=
q_p^{-1}(q_c)
=
-\frac1{u_p}
\log\left(1-\frac{q_c}{p}\right).
```

La contraction $`\chi^2`$ d'une arête observée vaut $`\gamma_p=(2p-1)^2`$. Pour placer la borne d'information-percolation sur la même échelle que le dendrogramme, définissons

```math
t_\chi(p)
=
q_p^{-1}(\gamma_p)
=
-\frac1{u_p}
\log\left(
1-\frac{(2p-1)^2}{p}
\right).
```

La convention est $`q_p(t_\chi)=(2p-1)^2`$, et non $`q_p(t_\chi)^2=q_c`$.

Alors

```math
t_\chi(p)>\beta_c(p)
\quad\Longleftrightarrow\quad
(2p-1)^2>q_c.
```

On retrouve donc exactement

```math
p>p_c^{\mathrm{info}}
=
\frac{1+\sqrt{q_c}}2
=0.7946592758\ldots
```

comme condition nécessaire. Au seuil, $`t_\chi=\beta_c=0.4245677743\ldots`$. Pour $`p\ge p_c^{\mathrm{info}}`$, après contraction des composantes critiques, toute connexion nouvellement créée dans le graphe d'information-percolation peut être représentée en loi par la sous-bande $`\beta_c<\xi_e\le t_\chi`$. Cette représentation est annealed et marginale ; elle ne transforme pas $D$ révélé en information observée supplémentaire.

Trois tests doivent rester distincts :

```math
\begin{array}{c|c|c}
\text{objet}&\text{condition}&\text{seuil en }p\\
\hline
\text{quotient jusqu'à }1&q_p(1)>q_c&0.673648\ldots\\
\text{quotient informationnel}&(2p-1)^2>q_c&0.794659\ldots\\
\text{bande pure }(\beta_c,1]&q_p(1)-q_c>q_c&0.847296\ldots
\end{array}
```

Le premier redonne Swendsen--Wang. Le dernier est trop exigeant : il interdit d'utiliser les composantes déjà construites sous $`\beta_c`$.

### Décomposition exacte du vote de bande conditionnellement au squelette

Pour un lien encore fermé au temps $\beta$, posons

```math
h_p(\beta)
=
\mathbb P(\beta<\xi_e\le1\mid\xi_e>\beta)
=
\tanh\left(\frac{u_p(1-\beta)}2\right).
```

Conditionnellement au squelette de Kruskal non marqué, à un bucket de $m$ liens fusionnant à $\beta$, une arête gagnante latente est uniforme et satisfaite. Pour les $m-1$ autres liens, notons $R$ le nombre de liens satisfaits dans la bande, $S$ le nombre de liens satisfaits après la censure $1$, et $U$ le nombre de liens insatisfaits. On a exactement

```math
(R,S,U)
\sim
\mathrm{Mult}\left(
m-1;
h_p(\beta),
\frac{1-h_p(\beta)}2,
\frac{1-h_p(\beta)}2
\right),
```

et

```math
k=1+R+S.
```

Ainsi, $R$ est exactement le signal biaisé produit par la bande, tandis que $S-U$ est un bruit symétrique. L'existence d'une arête de bande est gouvernée par $`m h_p(\beta)`$ ; la transmission de la parité est gouvernée par le rapport signal sur bruit quadratique

```math
m h_p(\beta)^2.
```

Pour une grande coupe, la fenêtre locale non triviale est $`1-\beta\asymp m^{-1/2}`$. Le biais de Kruskal porte donc sur la loi géométrique de $`(E_u,m,\beta)`$, pas sur le noyau conditionnel des marques. Le prochain objet à estimer est la loi jointe de ce squelette et du message ancestral $`B_u`$, et non une probabilité de chemin tardif non pondérée.

### Calibrations finies annealed

Dans les quatre exemples suivants, $`q=q_p(1)=2p-1`$.

- **Une arête.** La valeur propre centrée du noyau pair-spécifique vaut $\lambda=q/(1+q)$, donc
```math
    A_{ij}^{(m)}
    =
    q^2+(1-q^2)\left(\frac{q}{1+q}\right)^m.
```
    En particulier $`A_{ij}^{(1)}=q`$ et $`A_{ij}^{(m)}\downarrow q^2=c_{ij}^2`$.
- **Chemin de longueur $\ell$.** Chaque coupe a $m=1$, donc la borne à un pas vaut $`A_{ij}^{(1)}=q^\ell`$, tandis que la valeur exacte est $`c_{ij}^2=q^{2\ell}`$. Ce gap montre pourquoi les itérations sont nécessaires même sans cycles.
- **Triangle isolé.** La probabilité de connexion FK vaut $q+q^2-q^3$, mais la borne LCA à un pas se simplifie en $`A_{ij}^{(1)}=q`$. La corrélation exacte à deux répliques est
```math
    \eta_\triangle
    =
    \frac12\left[
    \frac{(q+q^2)^2}{1+q^3}
    +
    \frac{(q-q^2)^2}{1-q^3}
    \right].
```
    Elle se simplifie en
```math
    \eta_\triangle
    =
    \frac{q^2(1+2q^2)}{1+q^2+q^4}.
```
- **Cactus de $L$ triangles.** La fonction de connexion brute se factorise en $`[q+q^2-q^3]^L`$, tandis que la corrélation exacte vaut $`\eta_\triangle^L`$. C'est le premier banc d'essai pour vérifier numériquement et analytiquement $`A_{ij}^{(m)}\downarrow\eta_\triangle^L`$.

Les démonstrations générales et la distinction entre couplage depuis la vérité et réplique indépendante sont dans [06_LCA_SPIN_CORRELATION.md](06_LCA_SPIN_CORRELATION.md).

## 8. Cas calculables à traiter dans cet ordre

1. **Arbre fini à degrés bornés.** Kruskal n'introduit aucun cycle ; comparer exactement capacité hiérarchique et reconstruction broadcast.
2. **Cactus de triangles.** Premier modèle avec interactions triangulaires mais dépendances contrôlables.
3. **Bandes triangulaires de largeur fixe.** Matrice de transfert et seuils numériques certifiables.
4. **Tore triangulaire fini.** Énumération exacte pour petites tailles, puis bornes avec conditions de bord périodiques.
5. **Limite planaire.** Seulement après une domination uniforme et un contrôle des effets de bord.

## 9. Objectif quantitatif

Le premier succès non trivial est une constante $`p_\star`$ telle que
```math
p_\star>0.794659\ldots
```
et une preuve que $`p<p_\star`$ interdit la weak recovery. Une condition nécessaire et suffisante explicite sur toute la grille triangulaire toucherait au point multicritique de Nishimori et constitue un objectif de long terme.

## 10. Audit du canal de triangle

Les triangles montants couvrent chaque arête exactement une fois. Il est donc
légitime de regrouper les trois observations d'un triangle en un facteur
indépendant conditionnellement aux spins. Ce regroupement ne donne toutefois
pas la contraction uniforme comme paramètre SDPI global.

Avec $`q=2p-1`$, les deux valeurs exactes sont

```math
\eta_\triangle(q)
=\frac{q^2(1+2q^2)}{1+q^2+q^4}
```

sous l'a priori uniforme, et

```math
\gamma_2(q)=\frac{2q^2}{1+q^2}
```

pour la SDPI globale. La seconde valeur est atteinte lorsque l'extérieur a
déjà révélé une relation et que les deux arêtes variables deviennent deux
observations indépendantes du même bit. Une percolation scalaire de triangles
aurait donc le seuil

```math
\gamma_2(q)=\frac12
\quad\Longleftrightarrow\quad
p=\frac12\left(1+\frac1{\sqrt3}\right)
=0.788675\ldots,
```

moins bon que la baseline par arêtes.

Le [calcul multi-état complet](11_TRIANGLE_BLOCK_SDPI.md) conserve séparément
le triangle plein, chacune des trois relations révélées et l'état vide. Il
produit le candidat algébrique

```math
p_\star^{\mathrm{cond}}=0.8099092892\ldots,
```

mais seulement sous un lemme less-noisy encore ouvert pour les a priori ayant
un atome de masse $`>1/2`$. Ce lemme est réduit à la positivité d'une matrice
rationnelle $`3\times3`$ au seul point algébrique critique. Le gain numérique
n'est donc pas encore une borne rigoureuse. La fonction exacte $`c_q(t)`$ de
ce fichier isole un problème de polarisation pour un facteur physique, mais ne
transforme pas la difficulté des $`\Lambda_v`$ ancestraux. La voie
hiérarchique prioritaire est le
[cas favorable critique](12_FAVORABLE_HIERARCHICAL_REDUCTION.md) : paire
lointaine du même arbre, séparation en $`\beta_c`$, puis calcul du message
formé par tous les ancêtres et preuve de la domination HF.
