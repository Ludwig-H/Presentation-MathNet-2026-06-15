# Cas test : GSBM homogène sur la grille triangulaire

Ce cas est le premier banc d'essai parce que plusieurs quantités de percolation sont explicites. Il faut toutefois fixer une exhaustion finie — idéalement un tore triangulaire — avant tout énoncé asymptotique rigoureux.

## 1. Modèle

On prend

\[
f_{\mathrm{in}}\equiv p,
\qquad
f_{\mathrm{out}}\equiv1-p,
\qquad
p\in[1/2,1).
\]

Toutes les interactions observées ont le même module

\[
u_p=\log\frac p{1-p}.
\]

Sous la loi générative, conditionnellement à la vérité \(\Sigma\) mais après moyenne sur l'observation \(O\), le signe observé d'une arête est satisfait par \(\Sigma\) avec probabilité \(p\), indépendamment des autres arêtes.

## 2. Loi exacte de chaque coupe du dendrogramme

Une arête est présente à la coupe \(t\) si elle est satisfaite et si son horloge est inférieure à \(t\). Par conséquent,

\[
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
\]

Sous cette même loi annealed, conditionnellement à \(\Sigma\), les indicatrices sont indépendantes. Ainsi \(\Pi_t\) a exactement la loi d'une percolation indépendante par arêtes de paramètre \(q_p(t)\), et les différentes coupes sont couplées de façon monotone par les mêmes horloges.

Cette identité ne vaut pas conditionnellement à une observation \(O\) fixée. Elle vaut aussi lorsqu'on remplace la vérité par une réplique postérieure à l'équilibre, après moyenne jointe sur \(O\), par l'identité de Nishimori.

En particulier,

\[
q_p(1)=2p-1,
\]

ce qui redonne la percolation de Swendsen–Wang du chapitre 11.

Soit

\[
q_c^{\mathrm{bond}}(\mathbb T)=2\sin(\pi/18).
\]

Le temps auquel la filtration brute atteint le seuil de percolation est

\[
t_c(p)
=
\frac{
\log\!\left(1-q_c^{\mathrm{bond}}/p\right)
}{
\log\!\left((1-p)/p\right)
},
\]

lorsque le membre de droite est défini. On a \(t_c(p)\le1\) si et seulement si

\[
p\ge
\frac12+\sin(\pi/18)
=0.673648\ldots
\]

La distribution des temps de coalescence est donc reliée aux fonctions de connexion de la percolation :

\[
\mathbb P_{\mathrm{ann}}(\beta_{ij}\le t\mid\Sigma)
=
\mathbb P_{q_p(t)}(i\leftrightarrow j).
\]

Cette identité rend toute la **géométrie brute** du dendrogramme calculable. La nouvelle information doit venir des probabilités de flips, pas seulement de \(\Pi_t\).

## 3. Trois bornes rigoureuses déjà disponibles

### Swendsen–Wang par arêtes

\[
p_c^{\mathrm{edge}}
=
\frac12+\sin(\pi/18)
=0.673648\ldots
\]

En dessous, le graphe gelé est sous-critique et la borne du chapitre 11 interdit la weak recovery, après formalisation du passage en volume fini.

### Dynamique triangulaire d'ordre supérieur

Pour les triangles disjoints d'une couleur, posons

\[
\alpha_p=1-e^{-2u_p}
=
\frac{2p-1}{p^2}.
\]

Après moyenne sur les observations, les états locaux de percolation ont les probabilités

\[
a=p(2p-1)
\]

pour le triangle plein,

\[
s=(1-p)(2p-1)
\]

pour chacun des trois états à une arête, et

\[
e=4(1-p)^2
\]

pour le triangle vide. On vérifie \(a+3s+e=1\).

La condition autoduale \(a=e\) donne

\[
p_c^\triangle
=
\frac{7-\sqrt{17}}4
=0.719224\ldots
\]

Le théorème de Chayes–Lei demande aussi

\[
ae\ge2s^2,
\qquad
a+e>\frac{2\sqrt2}{3+2\sqrt2},
\]

ainsi que l'isotropie et l'indépendance entre les triangles choisis. Ces hypothèses doivent accompagner l'utilisation du seuil.

L'intervalle de stricte amélioration sur la dynamique par arêtes est

\[
\left(
\frac12+\sin(\pi/18),
\frac{7-\sqrt{17}}4
\right).
\]

La borne gauche est ouverte : au point critique bidimensionnel, il n'y a pas de composante de densité positive sous les conditions de bord usuelles.

### Information-percolation

La contraction \(\chi^2\) du canal d'une arête vaut

\[
\eta=(2p-1)^2.
\]

La sous-criticité de la percolation d'information donne

\[
p<
p_c^{\mathrm{info}}
:=
\frac{1+\sqrt{2\sin(\pi/18)}}2
=0.794659\ldots
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
\]

Cette borne est la baseline rigoureuse à dépasser.

## 4. Seuil informationnel attendu

Après jauge par la vérité, le modèle est relié à l'Ising \(\pm J\) sur la ligne de Nishimori. Les calculs de dualité/répliques et les simulations de la littérature situent le point multicritique triangulaire vers

\[
p_{\mathrm N}=0.8358058\ldots
\]

Cette valeur est **conjecturale**, pas un seuil rigoureux à utiliser dans une preuve.

Le paysage de travail est donc :

| Méthode | Seuil \(p\) | Statut |
|---|---:|---|
| FK / Swendsen–Wang par arêtes | \(0.673648\ldots\) | Rigoureux |
| Dynamique triangulaire | \(0.719224\ldots\) | Rigoureux sous les hypothèses de Chayes–Lei |
| Information-percolation \(\chi^2\) | \(0.794659\ldots\) | Rigoureux |
| Point de Nishimori triangulaire | \(0.8358058\ldots\) | Conjecture + numérique |

## 5. Poids exact d'une fusion locale

Considérons une coupe \(E_u\) contenant \(m\) arêtes, dont \(k\) sont satisfaites par la configuration courante. Alors

\[
T_u=mu_p,
\qquad
\Lambda_u=ku_p.
\]

Si l'a priori et les facteurs ancêtres s'annulent dans le rapport de parité, la log-vraisemblance locale vaut

\[
\boxed{
L_{m,k,\beta}^{\mathrm{loc}}
=
\log\frac{k}{m-k}
+(1-\beta)u_p(2k-m).
}
\]

La fiabilité associée est

\[
\rho_{m,k,\beta}
=
\left|\tanh\frac{L_{m,k,\beta}^{\mathrm{loc}}}{2}\right|.
\]

Pour une coupe fixée ayant \(k\ge1\) arêtes satisfaites,

\[
\beta=\min_{1\le j\le k}\xi_j
\sim\operatorname{Exp}(ku_p).
\]

Conditionnellement à une fusion avant \(1\),

\[
f(\beta\mid\beta\le1,k)
=
\frac{ku_pe^{-ku_p\beta}}{1-e^{-ku_p}},
\qquad 0\le\beta\le1.
\]

Cette loi ne peut pas être appliquée naïvement à une coupe **sélectionnée par Kruskal** : la sélection introduit un biais. Les premiers calculs doivent donc distinguer :

1. une coupe déterministe ;
2. une fusion du dendrogramme conditionnée par son passé ;
3. une fusion typique sous la loi stationnaire.

## 6. Identité à deux répliques

Dans ce modèle,

\[
Q_n
=
\frac1{n^2}
\sum_{i,j}
\mathbb E\left[
\left\langle\sigma_i\sigma_j\right\rangle_O^2
\right].
\]

Sur la ligne de Nishimori, des identités supplémentaires relient les corrélations plantées et les corrélations de répliques. La stratégie correcte est de contrôler cette somme, directement ou via \(H_S\), plutôt que de chercher seulement une composante géante dans \(\Pi_1\).

## 7. Cas calculables à traiter dans cet ordre

1. **Arbre fini à degrés bornés.** Kruskal n'introduit aucun cycle ; comparer exactement capacité hiérarchique et reconstruction broadcast.
2. **Cactus de triangles.** Premier modèle avec interactions triangulaires mais dépendances contrôlables.
3. **Bandes triangulaires de largeur fixe.** Matrice de transfert et seuils numériques certifiables.
4. **Tore triangulaire fini.** Énumération exacte pour petites tailles, puis bornes avec conditions de bord périodiques.
5. **Limite planaire.** Seulement après une domination uniforme et un contrôle des effets de bord.

## 8. Objectif quantitatif

Le premier succès non trivial est une constante \(p_\star\) telle que

\[
p_\star>0.794659\ldots
\]

et une preuve que \(p<p_\star\) interdit la weak recovery. Une condition nécessaire et suffisante explicite sur toute la grille triangulaire toucherait au point multicritique de Nishimori et constitue un objectif de long terme.
