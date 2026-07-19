# Projection collapsed et tensorisation scalaire restreinte

> [!CAUTION]
> Les théorèmes de projection de la section 2 restent exacts. Le théorème 4.1
> reste exact pour le **surrogate produit mono-bit** défini en section 3, mais
> son ancienne identification au corridor collapsed réel est fausse. Une
> fusion multiport à gagnante marginalisée fournit un contre-exemple de même
> squelette et de même taille ; voir le
> [fichier 29](../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md).

Cette note affine la voie favorable des fichiers 18--19. Elle répond à deux
questions distinctes.

1. Quelle dynamique hiérarchique donne le certificat $`L^2`$ le plus fort
   pour une paire fixée ?
2. Comment composer les transferts multiports aux rangs réalisés sans
   réutiliser séquentiellement les mêmes buckets dans les messages
   descendants ?

La mise à jour conjointe des orientations du **corridor** de la paire reste
un heat bath hiérarchique exact et fournit l'enveloppe $`L^2`$ de la section
2. En revanche, l'ordre de Blackwell du fichier 19 ne se tensorise sur tout
le bloc que lorsque chaque bucket dépend réellement d'un unique bit latent.

Le résultat ne constitue pas encore une preuve d'impossibilité à $`p=0.8`$ :
la géométrie Palm du corridor triangulaire et le contrôle de son état de bord
restent ouverts.

## 1. Le corridor hiérarchique d'une paire

Fixons $`O,D`$ et

```math
\pi_D(d\sigma)=\nu_O(d\sigma\mid D).
```

Pour une paire $`i,j`$ appartenant à une même racine, notons
$`\mathcal C_{ij}`$ l'union des deux chaînes du dendrogramme qui joignent
$`i,j`$ à leur LCA, prolongée vers la racine lorsque les facteurs ancestraux
sont conservés. Les degrés de liberté sont les orientations relatives des
enfants aux nœuds de ce corridor.

Soit $`\mathcal A_{ij}`$ la tribu qui conserve :

- les spins et orientations hors du corridor ;
- le dendrogramme non marqué $`D`$ ;
- toutes les données qui ne sont pas rééchantillonnées par le bloc.

De façon intrinsèque, si $`\mathcal G_u`$ est la tribu conservée par le heat
bath local en $`u`$, on prend

```math
\mathcal A_{ij}
=
\bigcap_{u\in\mathcal C_{ij}}\mathcal G_u.
\qquad\text{(1.1)}
```

La dynamique **collapsed corridor** est le heat bath de tous les degrés de
liberté du corridor en un seul bloc. Son opérateur est

```math
P_{\mathcal C}g
=
\mathbb E_{\pi_D}[g\mid\mathcal A_{ij}].
\qquad\text{(1.2)}
```

Il s'agit d'une mise à jour hiérarchique valide, mais pair-spécifique. Cette
dépendance en $`(i,j)`$ interdit d'appliquer directement la matrice de Gram
d'un sweep commun du fichier 18. Le théorème 2.2 ci-dessous donne la
globalisation pairwise correcte.

## 2. Optimalité $`L^2`$ du bloc collapsed

Pour chaque nœud $`u\in\mathcal C_{ij}`$, soit

```math
P_u=\mathbb E_{\pi_D}[\,\cdot\mid\mathcal G_u]
```

le heat bath local. On a $`\mathcal A_{ij}\subseteq\mathcal G_u`$ : un
update local conserve davantage d'information que l'update conjoint de tout
le corridor.

### Théorème 2.1 — enveloppe collapsed, statut : établi

Pour tout programme systématique

```math
K=P_{u_1}\cdots P_{u_M},
\qquad u_r\in\mathcal C_{ij},
```

et tout $`g\in L^2(\pi_D)`$,

```math
\boxed{
\|P_{\mathcal C}g\|_2^2
\le
\|Kg\|_2^2.
}
\qquad\text{(2.1)}
```

Plus précisément,

```math
\boxed{
\|Kg\|_2^2
=
\|P_{\mathcal C}g\|_2^2
+\|K(I-P_{\mathcal C})g\|_2^2.
}
\qquad\text{(2.2)}
```

#### Preuve

Les tribus sont emboîtées, donc

```math
P_{\mathcal C}P_u=P_{\mathcal C},
\qquad
P_uP_{\mathcal C}=P_{\mathcal C}.
```

Par composition,

```math
P_{\mathcal C}K=P_{\mathcal C},
\qquad
KP_{\mathcal C}=P_{\mathcal C}.
```

Écrivons $`g=P_{\mathcal C}g+h`$ avec
$`h=(I-P_{\mathcal C})g`$. Alors
$`Kh\perp\mathrm{Ran}(P_{\mathcal C})`$, puisque
$`P_{\mathcal C}Kh=0`$. Ainsi

```math
Kg=P_{\mathcal C}g+Kh
```

est une somme orthogonale, ce qui donne (2.2).

### Conséquence

Pour $`f_{ij}(\sigma)=\sigma_i\sigma_j`$,

```math
\mathbb E[H_{\mathcal C}(i,j)^2\mid O,D]
=
\|P_{\mathcal C}f_{ij}\|_2^2.
\qquad\text{(2.3)}
```

Le bloc collapsed est donc le certificat le plus contractant parmi les
sweeps qui n'utilisent que les mêmes coordonnées du corridor. Cela ne dit pas
qu'il se mélange plus vite comme algorithme global ; (2.1) est un ordre sur
la persistance $`L^2`$ d'une paire après un bloc fixé.

Le [fichier 22](../results/hierarchical/22_LCA_VS_FULL_HIERARCHY.md) spécialise cette enveloppe à la
question de profondeur. Si $`u`$ est le LCA, le corridor descendant collapsed
vérifie l'identité plus précise

```math
\|P_uf_{ij}\|_2^2
=
\|P_{\downarrow}f_{ij}\|_2^2
+\|(P_u-P_{\downarrow})f_{ij}\|_2^2.
```

Il distingue aussi le bloc collapsed d'un unique passage séquentiel :
bottom-up est comparable au LCA seul, tandis que top-down exige un contrôle
supplémentaire de la non-commutation.

### Théorème 2.2 — globalisation pair-spécifique, statut : établi

Pour chaque paire $`i,j`$, choisissons son propre corridor collapsed et
posons

```math
A_{ij}(O,D)
=
\|P_{\mathcal C_{ij}}f_{ij}\|_{L^2(\pi_D)}^2.
```

Pour deux racines distinctes, $`P_{\mathcal C_{ij}}`$ inclut les deux
recolorations globales et $`A_{ij}=0`$. Sur la diagonale, on pose
$`A_{ii}=1`$.

Si $`I_n,J_n`$ sont uniformes et indépendants, alors

```math
\boxed{
Q_n
\le
\mathbb E[A_{I_nJ_n}(O,D)].
}
\qquad\text{(2.4)}
```

Par conséquent,

```math
\mathbb E[A_{I_nJ_n}(O,D)]\longrightarrow0
\quad\Longrightarrow\quad
\text{pas de weak recovery}.
\qquad\text{(2.5)}
```

#### Preuve

Écrivons

```math
c_{ij}(O)=\mathbb E_{\mu_O}[f_{ij}].
```

La désintégration de la loi jointe puis Jensen donnent

```math
c_{ij}(O)^2
\le
\mathbb E_{D\mid O}
\left[
\mathbb E_{\pi_D}[f_{ij}]^2
\right].
```

Comme $`P_{\mathcal C_{ij}}`$ est une espérance conditionnelle,

```math
\mathbb E_{\pi_D}[f_{ij}]
=
\mathbb E_{\pi_D}[P_{\mathcal C_{ij}}f_{ij}],
```

et une seconde application de Jensen donne

```math
\mathbb E_{\pi_D}[f_{ij}]^2
\le
\|P_{\mathcal C_{ij}}f_{ij}\|_2^2.
```

On moyenne en $`O,D`$, puis on somme sur les paires. L'identité à deux
répliques écrit $`Q_n=n^{-2}\sum_{i,j}\mathbb E[c_{ij}(O)^2]`$, d'où (2.4).
Cette preuve n'exige pas que les différentes quantités $`A_{ij}`$ forment une
matrice positive semi-définie commune.

## 3. Surrogate de comptes mono-bit sur un squelette fixé

Fixons maintenant une expérience factorisée de $`h`$ buckets. Pour
$`1\le r\le h`$, notons :

- $`X_r\in\{-1,+1\}`$ sa parité latente ;
- $`m_r`$ la taille du bucket ;
- $`t_r`$ son niveau ;
- $`K_r`$ son nombre de liens satisfaits dans l'orientation de référence.

Conditionnellement au squelette non marqué et au vecteur
$`X=(X_1,\ldots,X_h)`$, les marques de buckets distincts sont indépendantes.
Avec $`s_r=s_p(t_r)`$,

```math
\begin{aligned}
K_r\mid X_r=+1
&\ \overset{\mathrm d}=\ 1+\mathrm{Bin}(m_r-1,s_r),\\
K_r\mid X_r=-1
&\ \overset{\mathrm d}=\ \mathrm{Bin}(m_r-1,1-s_r).
\end{aligned}
\qquad\text{(3.1)}
```

L'hypothèse décisive est que toutes les arêtes du bucket $`r`$ codent le même
bit $`X_r`$. Elle correspond notamment à des relations descendantes figées
ou révélées. Elle ne décrit pas un bucket ancestral multiport lorsque des
flips descendants font varier séparément plusieurs groupes d'incidence.

La loi a priori $`\rho(dx)`$ de $`X`$ peut être arbitraire : les parités
peuvent être corrélées par un état de bord, des branches latérales ou des
contraintes de cycles. Seule l'indépendance des **observations de buckets
sachant $`X`$** est utilisée ci-dessous.

Pour une cible $`F(X)\in[-1,1]`$, la persistance collapsed vaut

```math
\mathscr R(\rho,F;(m_r,t_r)_{r\le h})
=
\mathbb E\left[
\mathbb E[F(X)\mid K_1,\ldots,K_h]^2
\right].
\qquad\text{(3.2)}
```

Pour la paire, $`F`$ est le produit des parités qui séparent $`i`$ de $`j`$.

## 4. Tensorisation exacte de Blackwell

Comparons deux expériences sur le **même** squelette et les mêmes tailles.
Dans l'expérience critique, le bucket $`r`$ a un paramètre
$`s_r^{\mathrm c}`$ ; dans l'expérience tardive,

```math
\frac12\le s_r^{\mathrm l}\le s_r^{\mathrm c}\le1.
\qquad\text{(4.1)}
```

### Théorème 4.1 — surrogate produit mono-bit, statut : établi

Pour toute loi $`\rho`$ sur les parités, toute cible bornée $`F`$ et toute
famille satisfaisant (4.1),

```math
\boxed{
\mathscr R_{\mathrm l}(\rho,F)
\le
\mathscr R_{\mathrm c}(\rho,F).
}
\qquad\text{(4.2)}
```

Cette domination ne suppose ni l'indépendance des $`X_r`$, ni la neutralité
du message de bord.

#### Preuve

Le théorème 4.2 du fichier 19 fournit, pour chaque bucket, un noyau
stochastique $`G_r`$ indépendant de $`X_r`$ tel que

```math
P_{s_r^{\mathrm l}}^x
=
P_{s_r^{\mathrm c}}^xG_r,
\qquad x\in\{-1,+1\}.
```

Conditionnellement à $`X`$, les comptes sont indépendants. Le produit

```math
G=\bigotimes_{r=1}^hG_r
```

transforme donc le vecteur de comptes critique en le vecteur tardif, pour
chaque valeur de $`X`$. On a la chaîne de Markov

```math
X\longrightarrow K^{\mathrm c}\longrightarrow K^{\mathrm l}.
```

Par la propriété de tour,

```math
\mathbb E[F(X)\mid K^{\mathrm l}]
=
\mathbb E[
\mathbb E[F(X)\mid K^{\mathrm c}]
\mid K^{\mathrm l}],
```

et Jensen conditionnelle donne (4.2).

### Portée exacte et limite multiport

Le théorème 4.1 résout la comparaison critique/postcritique dans
l'**expérience produit mono-bit** de la section 3. Il supprime, dans ce
surrogate :

- les contre-exemples pointwise d'anti-alignement ;
- le feedback séquentiel où un descendant réutilise un compte ancestral.

Il ne décrit pas le corridor collapsed réel dès qu'un même facteur ancestral
dépend de plusieurs caractères du vecteur de flips. Avec deux incidences et
une gagnante marginalisée, les canaux critique et tardif ne sont pas ordonnés
par Blackwell ; sous un bord polarisé, le second moment de la parité cible
peut même être plus grand au temps tardif. Le certificat est dans le fichier
29. Cette obstruction apparaît sans changer ni squelette ni taille.

### Contre-lemme 4.2 — changement de taille, statut : établi

L'hypothèse « mêmes tailles » ne peut pas être retirée. Pour une expérience
binaire $`\mathcal E`$ sous prior uniforme, soit

```math
\Pi_{\mathcal E}=\mathbb P(X=+1\mid K),
\qquad
C_{\mathcal E}(z)=\mathbb E[(\Pi_{\mathcal E}-z)_+].
\qquad\text{(4.3)}
```

Deux expériences binaires ont une moyenne postérieure égale à $`1/2`$.
Ainsi, $`\mathcal E_A`$ Blackwell-domine $`\mathcal E_B`$ si et seulement si

```math
C_{\mathcal E_A}(z)\ge C_{\mathcal E_B}(z)
\quad\text{pour tout }z\in[0,1],
\qquad\text{(4.4)}
```

c'est-à-dire si la première postérieure domine la seconde dans l'ordre
convexe.

Prenons $`p=4/5`$ et comparons :

```math
\mathcal E_A=\mathcal E_{4,s_c}
\quad\text{au niveau critique},
\qquad
\mathcal E_B=\mathcal E_{2,s_p(4/5)}
\quad\text{au niveau }t=4/5.
```

Ces deux expériences sont **incomparables** au sens de Blackwell. Un
certificat par arithmétique rationnelle donne

```math
C_A(z_A)-C_B(z_A)
\in[-0.007184305272,-0.007184305271]
\qquad\text{(4.5)}
```

au point

```math
z_A
=
\frac{(1-s_c)^2}{(1-s_c)^2+3s_c^2}
\in[0.061085324056076,0.061085324056078],
```

tandis que

```math
C_B(1/2)-C_A(1/2)
\in[-0.044555124600,-0.044555124599].
\qquad\text{(4.6)}
```

La première inégalité interdit $`A\succeq B`$ et la seconde interdit
$`B\succeq A`$.

#### Certification

On utilise seulement les inclusions rationnelles

```math
q_\triangle
\in
[0.347296355333860,0.347296355333861],
```

certifiées par $`q^3-3q+1=0`$ et le changement de signe aux deux bornes, et

```math
4^{-1/5}
\in
[0.757858283255199,0.757858283255200],
```

certifiée en élevant les deux bornes à la puissance cinq. On propage ensuite
des intervalles de fractions dans

```math
s_c=\frac{4/5-q_\triangle}{1-q_\triangle},
\qquad
s_p(4/5)=\frac{4y^4}{1+4y^4},
\quad y=4^{-1/5}.
```

Le calcul reproductible est
[`p_eight_cross_size_incomparability_certificate`](../computations/favorable_time_comparison.py).
Aucun arrondi flottant n'intervient dans le signe de (4.5)--(4.6).

Ce contre-lemme montre déjà que la taille doit être conservée dans le
surrogate mono-bit. Le contre-exemple multiport du fichier 29 va plus loin :
même conserver le squelette, la taille et les incidences ne suffit pas si
plusieurs caractères latents varient dans un bucket.

### Lemme 4.3 — séparation géométrie/canal au seuil, statut : établi

Sous la loi jointe annealed de Nishimori, jaugeons les observations par la
réplique postérieure utilisée pour construire la hiérarchie. Pour chaque
arête, l'indicateur de satisfaction est alors Bernoulli $`p`$, indépendamment
des autres arêtes. Posons

```math
q_p(t)=p(1-e^{-u_pt}),
\qquad
u_p=\log\frac p{1-p}.
\qquad\text{(4.7)}
```

Si $`T_e`$ est l'horloge censurée de l'arête et
$`R_e=q_p(T_e)`$ pour $`T_e\le1`$, alors, pour
$`0\le r\le2p-1`$,

```math
\mathbb P(R_e\le r)=r.
\qquad\text{(4.8)}
```

Les rangs $`R_e`$ sont indépendants. Pour tout
$`p>(1+q_\triangle)/2`$, le graphe ouvert à
$`r=q_\triangle`$ est donc exactement une percolation de Bernoulli de
paramètre $`q_\triangle`$. La forme non marquée de la forêt de Kruskal
jusqu'au seuil critique, ses tailles de buckets descendantes et tout
conditionnement de paire mesurable par rapport à cette forêt ont une loi qui
ne dépend pas de $`p`$. Le paramètre statistique reste, lui,

```math
s_c(p)=\frac{p-q_\triangle}{1-q_\triangle}.
\qquad\text{(4.9)}
```

#### Preuve

Une arête satisfaite a une horloge exponentielle de taux $`u_p`$ ; une arête
insatisfaite n'ouvre jamais. Par conséquent

```math
\mathbb P(T_e\le t)=p(1-e^{-u_pt})=q_p(t).
```

Le changement de variable monotone $`r=q_p(t)`$ donne (4.8). À
$`r=q_\triangle`$, les indicatrices
$`\mathbf1_{\{R_e\le q_\triangle\}}`$ sont donc i.i.d. Bernoulli
$`q_\triangle`$. Enfin, conditionnellement à $`R_e>q_\triangle`$, la
probabilité que l'arête soit satisfaite vaut (4.9).

Pour une coupe **fixée** de taille $`m`$, la densité de son premier rang est

```math
m(1-r)^{m-1}\,dr.
\qquad\text{(4.10)}
```

Au seuil, les horloges pénalisent donc une grande coupe par le facteur
$`(1-q_\triangle)^{m-1}`$. Cela ne prouve pas la tension des tailles sous la
loi Palm : le nombre de coupes candidates, le biais
$`|C_1||C_2|`$ de la paire et les contraintes imposées par les buckets
précédents peuvent compenser ce facteur. La conséquence rigoureuse de ce
lemme est une réduction : la géométrie critique descendante peut être
étudiée une seule fois en coordonnée $`q`$, puis le canal (4.9) évalué au
$`p`$ désiré.

## 5. Corridor factorisé et critère de perte

Supposons temporairement les $`X_r`$ indépendants et uniformes, et prenons

```math
F(X)=\prod_{r=1}^hX_r.
```

Écrivons $`\Gamma_{m_r}(t_r;p)`$ pour la fiabilité de l'expérience binaire
du bucket $`r`$.

### Corollaire 5.1 — produit exact, statut : établi sous factorisation

```math
\boxed{
\mathscr R
=
\prod_{r=1}^h\Gamma_{m_r}(t_r;p).
}
\qquad\text{(5.1)}
```

Par conséquent,

```math
\mathscr R\longrightarrow0
\quad\Longleftrightarrow\quad
\sum_{r=1}^h-log\Gamma_{m_r}(t_r;p)
\longrightarrow+\infty.
\qquad\text{(5.2)}
```

#### Preuve

La loi jointe et la postérieure se factorisent par bucket. Ainsi

```math
\mathbb E[F(X)\mid K_1,\ldots,K_h]
=
\prod_{r=1}^h\mathbb E[X_r\mid K_r].
```

Le carré et l'espérance donnent (5.1).

Cette formule redonne PATH-FAC, mais avec une interprétation désormais
précise : elle est exacte pour le corridor collapsed factorisé, pas pour un
sweep séquentiel quelconque sur la grille frustrée.

## 6. Certificat numérique exact à $`p=0.8`$

Au niveau critique,

```math
s_c(0.8)=0.693582222752\ldots.
```

Pour $`m=2`$, le canal de compte est un effacement et

```math
\Gamma_2(\beta_c;0.8)=s_c(0.8).
```

Un corridor neutre contenant $`N`$ blocs indépendants de taille deux donne
donc exactement

```math
\mathscr R_N=s_c^N
=\exp(-0.365885484247\ldots N).
\qquad\text{(6.1)}
```

| $`N`$ | $`s_c^N`$ |
|---:|---:|
| 5 | $`0.160505443478`$ |
| 10 | $`0.025761997386`$ |
| 20 | $`6.63680509318\,10^{-4}`$ |
| 40 | $`4.40471818450\,10^{-7}`$ |

Avec un message de bord borné par $`|B|\le b`$, le coefficient local du
fichier 19 est

```math
\kappa_2(b)
=
s_c+(1-s_c)\tanh^2(b/2)<1
\qquad\text{(6.2)}
```

pour tout $`b<\infty`$. Multiplier ces coefficients exige encore un
découplage ou un transfert de blocs ; la seule borne $`|B|\le b`$ n'introduit
pas cette indépendance.

Le script
[`collapsed_corridor_transfer.py`](../computations/collapsed_corridor_transfer.py)
énumère (3.2). À $`p=0.8`$, sur les tailles $`(2,3,2,4)`$, il donne :

| prior des parités | tout critique | niveaux $`(0.55,0.70,0.85,1)`$ |
|---|---:|---:|
| uniforme | $`0.232015050844`$ | $`0.047131567858`$ |
| chaîne d'Ising, interaction $`0.6`$ | $`0.426226710965`$ | $`0.221677424071`$ |

Ces nombres contre-auditent la tensorisation du **surrogate factorisé** avec
un prior indépendant puis corrélé. Ils ne sont ni des estimations du tore
triangulaire, ni des calculs du corridor collapsed multiport.

## 7. Théorème conditionnel ciblé à $`p=0.8`$

Pour une paire $`(I_L,J_L)`$ tirée sous la Palm des événements de fusion
réalisés, supposons que son corridor collapsed aux rangs réels contienne
$`N_L`$ blocs complets tels que :

1. $`N_L\to\infty`$ en probabilité ;
2. chaque bloc possède un coefficient répliqué au plus
   $`\kappa<1`$ conditionnellement à son état de bord ;
3. la somme des erreurs de compression de bord et de découplage est $`o(1)`$.

Alors

```math
\mathbb E_{\text{Palm d'événement}}
\left[H_{\mathcal C}^{\mathrm{réel}}(I_L,J_L)^2\right]
\le
\mathbb E[\kappa^{N_L}]+o(1)
\longrightarrow0.
\qquad\text{(7.1)}
```

Combiné à l'annulation des racines, à la disparition sous-critique et à une
certification directe des blocs multiports aux rangs réalisés, (7.1)
interdirait la weak recovery à $`p=0.8`$, puis à tout $`p\le0.8`$ par
dégradation BSC des observations.

Le contenu utile est l'identification d'un bloc heat bath optimal. La
tensorisation scalaire ne réduit pas automatiquement le transfert multiport ;
celui-ci doit être calculé avec son état de bord réel.

## 8. Les trois verrous restants sur la grille triangulaire

### G1 — transfert direct aux rangs réalisés

Construire les transferts de masse et tordus d'une cellule multiport en
conservant chaque rang réalisé. Le théorème 4.1 peut servir pour une
coordonnée véritablement mono-bit, mais ne fournit aucune comparaison
automatique du bloc collapsed. Une domination cible-spécifique reste
admissible si elle est démontrée sous la vraie loi de bord.

### G2 — compression de l'état de bord

Construire un état fini $`Z_r`$ tel que le transfert d'un bloc dépende de
$`Z_r`$ et non de toute l'histoire ancestrale. Sur cactus, $`Z_r`$ est fini
exactement à bord libre : le fichier 21 montre qu'il se réduit à la parité
répliquée et calcule le noyau fermé. Sur une bande de largeur fixée, sa taille est
exponentielle dans la largeur mais reste calculable.

### G3 — abondance de blocs contractants

Montrer sous la Palm d'événement réelle qu'un nombre divergent de blocs a
soit une petite interface, soit un message de bord screené. La distance entre
$`i`$ et $`j`$ ne suffit pas : un corridor macroscopique peut contenir des
goulots pivotaux mais aussi de grandes coupes presque déterministes.

Le cactus résout désormais $`G2`$ grâce à ses articulations scalaires.
L'ordre recommandé est : test de dernière interaction, quotient de jauge
Markov-fermé ou bloc multi-update, puis seulement identification du déficit
à une signature locale du corridor. Une cellule T2 qui conserve le twist
dans son état complet a $`|U|=K`$ et ne peut pas fournir ce déficit.

## 9. Audit et contre-audit

| affirmation | statut | limite exacte |
|---|---|---|
| Le corridor collapsed est un heat bath valide | Établi | bloc pair-spécifique |
| Il est plus contractant qu'un sweep des mêmes nœuds | Établi en $`L^2`$ | ne compare pas les temps de calcul |
| Blackwell se tensorise dans le surrogate produit mono-bit | Établi | un seul caractère latent par bucket et indépendance conditionnelle des observations |
| Le même ordre vaut pour un corridor collapsed multiport | Faux en général | contre-exemple $`m=2`$ au même squelette et à la même taille |
| Un bucket critique domine un bucket tardif de taille différente | Faux en général | contre-lemme 4.2 : expériences incomparables à $`p=t=4/5`$ |
| La géométrie non marquée sous le seuil dépend de $`p`$ | Faux en coordonnée $`q`$ | établi sous la loi jointe annealed ; les ancêtres postcritiques restent à traiter |
| Le facteur $`(1-q_\triangle)^{m-1}`$ donne une queue Palm | Non démontré | il faut contrôler l'entropie des coupes et le biais de paire |
| Les parités du surrogate doivent être indépendantes | Faux | un prior corrélé arbitraire est permis, mais chaque bucket doit rester mono-bit |
| La formule produit (5.1) vaut toujours | Faux | elle exige le prior factorisé et l'absence d'état de bord partagé |
| Le cas LCA critique est le plus favorable sur un cactus | Établi dans le fichier 21 | aucun transfert automatique à la grille |
| La persistance cactus tend vers zéro à $`p=0.8`$ | Établi exactement | chaîne d'articulations, pas cycles chevauchants |
| Une infinité de blocs $`m=2`$ prouve $`p=0.8`$ | Seulement sous (7.1) | abondance et transfert de bord non démontrés sur la grille |
| La perte collapsed prouve la perte top-down | Pas nécessaire | le bloc collapsed fournit directement un couplage pairwise valide |
| La constante de Nishimori suit de Blackwell | Faux | Blackwell ordonne les niveaux mais ne contrôle pas l'abondance des blocs |

La voie privilégiée est donc désormais à deux étages : le corridor collapsed
pour la preuve d'impossibilité, et le sweep top-down comme dynamique
séquentielle de comparaison et diagnostic. Les deux utilisent les mêmes
facteurs hiérarchiques $`\Lambda_v`$ ; le premier marginalise leur feedback
dans un bloc exact au lieu de le supposer indépendant.

Le [certificat cactus](../results/hierarchical/21_CACTUS_COLLAPSED_CERTIFICATE.md) réalise le premier
cas non trivial de ce programme et déplace le prochain verrou vers une bande
triangulaire de largeur deux.
