# Premiers résultats sur le corridor à $`p_0=161/200`$

Cette note consigne le premier cycle d'exécution de la
[sous-feuille de route](../../archive/roadmaps/27_SUBROADMAP_CORRIDOR_P0805.md). Elle sépare les
identités exactes, les certificats finis et les diagnostics géométriques.
Elle ne revendique pas encore la borne
$`p_{\mathrm{WR}}\ge161/200`$.

> [!IMPORTANT]
> **Décision scientifique, corrigée par le fichier 29.** La piste
> hiérarchique n'est pas réfutée, mais le motif « bucket $`m=2`$ avec message
> neutre » et le déficit sur état fidèle sont fermés. Le
> [fichier 30](../../active/30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md) introduit
> les projections collapsed imbriquées ; le
> [programme 35](../../active/35_DISTANCE_ENTROPIE_ERGODICITE.md) fixe la voie
> active actuelle. La jauge T2 reste un contre-test d'une compression
> spéciale, pas le préalable principal.

> [!CAUTION]
> La criticalisation du corridor final annoncée dans la première version de
> cette note est réfutée sur une cellule multiport exacte. Les quantités
> $`q_v^{\mathrm{fav}}`$ ci-dessous sont des proxies algébriques, pas des
> majorations informationnelles. Voir
> [l'audit à froid](../29_AUDIT_FROID_PIVOT_RANGS_REELS.md).

## 1. Ce qui est désormais fixé exactement

### 1.1 Deux représentations de Palm à ne plus confondre

Pour une coupe candidate entre deux blocs $`A,B`$ présents juste avant le
rang $`q=q_p(\beta)`$, sa contribution à l'intensité LCA-Palm d'une paire
lointaine est

```math
u_p s_p(\beta)
m(A,B)N_\rho(A,B)\,d\beta
=
\frac{m(A,B)N_\rho(A,B)}{1-q}\,dq.
\qquad\text{(1.1)}
```

Dans un arbre de Kruskal déjà réalisé, la course des arêtes a déjà produit
le facteur $`m`$ et le hazard. Un nœud de fusion réalisé doit donc être
pondéré seulement par $`N_\rho`$. Les deux estimateurs sont équivalents
uniquement lorsqu'ils portent sur la même fine fenêtre de rang.

Le contre-audit fini vérifie exactement

```math
mN_\rho
=
(\text{taux événementiel }m)
\times
(\text{poids du nœud réalisé }N_\rho).
\qquad\text{(1.2)}
```

Pondérer de nouveau un nœud réalisé par $`mN_\rho`$ créerait le faux poids
$`m^2N_\rho`$.

Un second audit vérifie, environnement par environnement,

```math
\sum_{u\text{ réalisé}}
N_\rho(C_{u,1},C_{u,2})
=
\#\left\{
(i,j):d_L(i,j)\ge\rho L,
i\leftrightarrow j\text{ à }q_1
\right\},
\qquad
q_1=2p-1.
\qquad\text{(1.3)}
```

Chaque paire ordonnée connectée est ainsi comptée exactement une fois, au
nœud qui est son LCA.

### 1.2 Deux géométries distinctes

Le calcul conserve deux expériences séparées.

1. **Benchmark snapshot critique.** À $`q_c`$, les coupes candidates de la
   partition critique sont pondérées par leur intensité $`mN_\rho`$ et un
   LCA synthétique est ajouté au rang $`q_c`$.
2. **Corridor final réel.** Kruskal est exécuté jusqu'à
   $`q_1=2p_0-1=0.61`$ ; les nœuds réalisés sont pondérés par $`N_\rho`$ et
   les deux bras réels de la paire sont conservés.

Le premier objet change le squelette. Il constitue un benchmark d'intensité
critique, pas une domination de Blackwell du second. Sur le corridor final,
on calculait aussi le rang tronqué diagnostique

```math
q_v^{\mathrm{fav}}
=
\min(q_v,q_c),
\qquad
\text{squelette, tailles, incidences et attaches inchangés}.
\qquad\text{(1.4)}
```

Cette opération ne domine pas le transfert multiport réel. Elle sert
uniquement à définir l'ancien proxy de charge ; tout nouveau transfert doit
utiliser $`q_v`$.

## 2. Diagnostic géométrique du corridor réel

Pour une paire et son corridor $`\mathcal C_{ij}`$, définissons le proxy

```math
G_{M,J}(i,j)
=
\#\left\{
v\in\mathcal C_{ij}:
2\le m_v\le M,
m_v h_{p_0}(q_v^{\mathrm{fav}})^2\le J
\right\}.
\qquad\text{(2.1)}
```

Les buckets physiques comptés le long des deux bras sont deux à deux
edge-disjoint. Cette propriété ne signifie pas qu'ils sont séparés de
l'extérieur : ni les ports latéraux, ni le message $`B_v`$, ni le transfert
répliqué ne sont encore calculés.

Le diagnostic suivant utilise

```math
p_0=0.805,
\qquad
\rho=\frac14,
\qquad
M=8,
\qquad
J=1.
\qquad\text{(2.2)}
```

Comme $`q_v^{\mathrm{fav}}\le q_c`$ et
$`h_c^2=0.161994444381\ldots`$, le filtre de charge exclut en fait
$`m\ge7`$. Le choix $`M=8`$ laisse donc apparaître que le seuil $`J=1`$,
et non la troncature $`M`$, est actif sur les plus grandes tailles.

Les erreurs ci-dessous sont des jackknives delete-one-environment. Les
nœuds d'un même environnement ne sont jamais traités comme des observations
i.i.d.

| côté $`L`$ | environnements | coupes du corridor | buckets $`m=2`$ | proxy criticalisé algébrique $`G_{8,1}`$ |
|---:|---:|---:|---:|---:|
| 8 | 100 | $`12.526\pm0.140`$ | $`2.236\pm0.091`$ | $`5.786\pm0.135`$ |
| 12 | 50 | $`19.002\pm0.328`$ | $`2.929\pm0.143`$ | $`8.687\pm0.281`$ |
| 16 | 20 | $`26.071\pm0.606`$ | $`3.348\pm0.171`$ | $`12.067\pm0.495`$ |
| 24 | 6 | $`42.214\pm1.661`$ | $`4.487\pm0.556`$ | $`20.145\pm1.189`$ |

La dernière ligne est exploratoire, car elle ne contient que six
environnements. Sur cette plage, le proxy familial croît beaucoup plus vite
que le seul compte $`m=2`$. Cela donne un **go** au transfert sur une famille
de petites coupes et un **no-go** à une preuve reposant exclusivement sur
$`m=2`$.

Le benchmark snapshot critique donne aussi une croissance du proxy, de
$`6.427\pm0.269`$ à $`L=8`$ à $`14.857\pm1.097`$ à $`L=24`$. Ces nombres
ne doivent pas être comparés terme à terme à ceux du corridor final : les
deux colonnes proviennent de géométries et de lois différentes.

> [!CAUTION]
> La croissance de $`G_{8,1}`$ n'est pas encore un lemme d'abondance. Une
> coupe de faible charge peut être polarisée par ses ancêtres ou contournée
> par une route latérale. Le diagnostic sert à choisir l'état du transfert,
> pas à conclure que la corrélation spin--spin décroît.

## 3. Premier transfert répliqué certifié : cellule E1+

La première cellule explicite possède deux ports gauches $`L_0,L_1`$, deux
ports droits $`R_0,R_1`$ et les arêtes

```math
L_0R_0,
\qquad
L_1R_0,
\qquad
L_1R_1,
\qquad
R_0R_1.
\qquad\text{(3.1)}
```

Le triplet $`L_1,R_0,R_1`$ forme un triangle. Les quatre arêtes sont
conditionnées fermées au rang $`q_c`$. Dans la jauge plantée, leurs marques
résiduelles sont indépendantes, de paramètre

```math
s_c
=
\frac{p_0-q_c}{1-q_c}
=
0.701242667184\ldots.
\qquad\text{(3.2)}
```

Pour un environnement résiduel $`Z`$, soit $`K_Z`$ le heat bath conjoint
exact des deux ports droits conditionnellement aux ports gauches. Le noyau
répliqué correct est

```math
\mathbb E_Z[K_Z\otimes K_Z],
\qquad\text{(3.3)}
```

avec le **même** $`Z`$ dans les deux copies. Il se réduit, sous les deux
symétries de flip global, en un bloc de masse et un bloc tordu
$`\chi\otimes\chi`$ sur quatre états de parités.

Le bloc de masse est stochastique. Dans le secteur tordu, le poids uniforme
donne le certificat rationnel

```math
\|\mathscr U_{\mathrm{E1+}}\|_{\infty}
<
0.293993788341
<
0.3
\qquad
\text{à }p_0=\frac{161}{200}.
\qquad\text{(3.4)}
```

L'encadrement de $`q_c`$ est certifié par le changement de signe de
$`q^3-3q+1`$ entre deux rationnels décimaux consécutifs ; toutes les
opérations suivantes utilisent des `Fraction` et des intervalles sortants.
Une seconde implémentation énumère directement les environnements sur des
chaînes de profondeur zéro, un et deux.

À titre de diagnostic, la répétition indépendante de cette cellule neutre
donne

| profondeur | second moment |
|---:|---:|
| 1 | $`0.293993788340`$ |
| 2 | $`0.0649753038062`$ |
| 3 | $`0.0135754848472`$ |
| 10 | $`1.89285427006\times10^{-7}`$ |

Deux environnements indépendants entre les répliques donneraient à tort le
coefficient $`0.086432347583`$. Ce contre-factuel sous-estime fortement la
persistance et verrouille la convention d'environnement partagé.

### Portée exacte de (3.4)

Le résultat (3.4) est un certificat **E1+** pour une cellule neutre dont les
arêtes sont toutes fermées. Ce n'est pas encore le bloc E2/T2 de la feuille
de route. Il manque :

- l'arête gagnante, ou de façon équivalente sa marginalisation par le facteur
  Palm de fusion ;
- la partition ouverte du squelette et le statut pivotal ;
- les quatre $`\Lambda_v^{ab}`$ ancestraux ;
- la compatibilité de Kruskal entre cellules successives ;
- la loi LCA-Palm et les attaches en peigne.

Le coefficient (3.4) n'est donc ni une borne de weak recovery, ni un
certificat bidimensionnel.

## 4. No-go exact : potentiel extérieur non borné

Ajoutons au heat bath de la cellule le facteur extérieur

```math
\exp\left(
\frac{B}{2}(R_0+R_1)
\right).
\qquad\text{(4.1)}
```

Toutes les vraisemblances résiduelles sont strictement positives. Lorsque
$`B\to+\infty`$, le bord droit converge donc uniformément vers
$`(R_0,R_1)=(+1,+1)`$. Le second moment brut de transport de la parité
converge vers un :

```math
A_{\mathrm{E1+}}(B)
\longrightarrow
1.
\qquad\text{(4.2)}
```

Les valeurs finies illustrent cette obstruction.

| champ $`B`$ | second moment |
|---:|---:|
| 0 | $`0.293993788340`$ |
| 2 | $`0.726322292861`$ |
| 4 | $`0.938706909042`$ |
| 8 | $`0.998663483928`$ |
| 20 | $`0.999999991755`$ |

Il est donc impossible d'obtenir une contraction **absolue brute** uniforme
sur tous les potentiels extérieurs non bornés. Ce no-go n'exclut ni une norme
centrée, ni un état enrichi qui suit la polarisation, ni une contraction
annealed sous la loi réelle des messages.

La conséquence opérationnelle est nette : il ne faut pas tenter de prouver
$`|B|\le B_0`$ comme une hypothèse gratuite. Le potentiel extérieur doit
être une coordonnée du transfert, et les états très polarisés doivent être
payés par leur loi ou par un drift.

## 5. Composition réaliste : déficit tordu dépendant de l'état

Le no-go précédent suggère de remplacer la marge uniforme par une
représentation de Feynman--Kac. Après une normalisation de Doob commune,
supposons que le transfert positif levé d'un bloc soit

```math
\mathsf T_r(z,dz',d\epsilon),
\qquad
\epsilon\in\{-1,+1\}.
\qquad\text{(5.1)}
```

Définissons son noyau de masse et son secteur tordu par

```math
K_r(z,dz')
=
\sum_\epsilon
\mathsf T_r(z,dz',d\epsilon),
\qquad
U_r(z,dz')
=
\sum_\epsilon
\epsilon\,
\mathsf T_r(z,dz',d\epsilon).
\qquad\text{(5.2)}
```

Si $`K_r`$ est markovien, alors la mesure de variation totale vérifie
$`|U_r|\le K_r`$. On peut donc poser

```math
r_r(z,z')
=
\frac{d|U_r|}{dK_r}(z,z')
\in[0,1],
\qquad
d_r(z,z')
=
-\log r_r(z,z').
\qquad\text{(5.3)}
```

La convention est $`d_r=+\infty`$ lorsque $`r_r=0`$.

Pour toute suite inhomogène de blocs et toute fonction terminale bornée,
l'inégalité de variation totale donne, avec la convention
$`(U_rg)(z)=\int U_r(z,dz')g(z')`$,

```math
\left|
U_1\cdots U_Ng
\right|(z_0)
\le
\mathbb E_{z_0}^{K_1,\ldots,K_N}
\left[
\exp\left(
-\sum_{r=1}^N d_r(Z_{r-1},Z_r)
\right)
|g(Z_N)|
\right].
\qquad\text{(5.4)}
```

Cette forme possède trois avantages.

1. Les états polarisés peuvent avoir $`d_r\simeq0`$ sans invalider le
   théorème.
2. Les bons blocs n'ont pas besoin d'avoir une marge identique.
3. La cible géométrique devient directement une espérance exponentielle
   sous un noyau de masse positif.

Le lemme fini (5.4) est maintenant implémenté pour des suites inhomogènes et
des espaces d'états de dimensions variables. Deux audits en arithmétique
rationnelle comparent exactement la récurrence dynamique à l'énumération de
tous les chemins. Les transitions de masse nulle ont bien
$`U=0`$ et reçoivent la convention $`r=0`$.

Sur les puissances de la cellule E1+, la borne de Feynman--Kac améliore la
simple puissance du pire coefficient dès la profondeur deux.

| profondeur | valeur signée | enveloppe Feynman--Kac | borne uniforme |
|---:|---:|---:|---:|
| 1 | $`0.293993788340`$ | $`0.293993788340`$ | $`0.293993788340`$ |
| 2 | $`0.0649753038062`$ | $`0.0738919329503`$ | $`0.0864323475826`$ |
| 5 | $`0.000564617656372`$ | $`0.00113949947561`$ | $`0.00219629550382`$ |
| 10 | $`1.89285427006\times10^{-7}`$ | $`1.08695758136\times10^{-6}`$ | $`4.82371394009\times10^{-6}`$ |

Le transformé de Doob rétrograde du fichier 29 ferme maintenant la
normalisation abstraite en horizon fini. Le verrou est la preuve

```math
\mathbb E
\left[
\exp\left(
-\sum_{r=1}^{N_{ij}}
d_r(Z_{r-1},Z_r)
\right)
\right]
\longrightarrow0
\qquad
\text{sous la loi marquée du corridor réel}.
\qquad\text{(5.5)}
```

La formule (5.4) est donc établie et contre-auditée en dimension finie.
L'identification exacte de $`K_r,U_r`$ au heat bath collapsed, la
localisation des fonctions de Doob rétrogrades, le cocycle de bord et (5.5)
restent les étapes globales.

## 6. Cellule quotient T2-Kruskal : première couche construite

Le fichier 29 construit maintenant une cellule quotient conditionnelle avec
facteurs de fusion exacts, gagnantes marginalisées, attache ancestrale et
polarisation. Elle n'est pas encore embarquée dans une géométrie Palm réelle.
Elle suffit à réfuter la criticalisation uniforme sur l'état multiport
complet. L'enrichissement suivant doit ajouter exactement les objets
nécessaires au déficit composable.

### État minimal à conserver

```math
z
=
(\mathcal G,\Pi,\Psi,\xi_{\mathrm{front}},\mathfrak p),
\qquad\text{(6.1)}
```

où $`\mathcal G`$ contient le squelette et les groupes d'incidence,
$`\Pi`$ la partition signée des ports, $`\Psi`$ le potentiel extérieur
projectif, $`\xi_{\mathrm{front}}`$ les orientations encore nécessaires aux
facteurs futurs et $`\mathfrak p`$ le statut de fusion ou d'attache en
peigne. L'état ne doit pas conserver une orientation interne après sa
dernière interaction : si le twist est mesurable depuis la transition, le
no-go du fichier 29 donne $`|U|=K`$ et $`d=0`$.

### Cellule exigée

Le bloc T2-Kruskal enrichi doit contenir :

1. une fusion réelle, avec l'identité de l'arête gagnante oubliée et le
   facteur $`\Lambda e^{(1-\beta)\Lambda}`$ exact ;
2. au moins trois ports et une route latérale possible ;
3. une attache postcritique à rang variable, avec le rang critique comme
   colonne comparative seulement ;
4. les quatre valeurs $`\Lambda_v^{ab}`$ d'au moins un ancêtre ;
5. le même environnement pour les deux répliques ;
6. soit une jauge de ports dont la transition de $`\Psi`$ est Markov-fermée,
   soit plusieurs updates allant jusqu'à la dernière interaction d'une
   orientation ;
7. une énumération directe et une élimination dynamique indépendantes.

La première sortie utile n'est pas un unique rayon spectral. C'est la carte

```math
(\mathcal G,\Pi,\Psi,\mathfrak p)
\longmapsto
r(z,z')
\quad\text{ou}\quad
d(z,z'),
\qquad\text{(6.2)}
```

qui permet de relier exactement le transfert à un événement géométrique
observable dans le diagnostic Palm. Cette carte n'est utile que si plusieurs
signes sont effectivement agrégés dans une même transition quotient ; sur
l'état fidèle complet, elle doit retourner $`d=0`$.

## 7. Ordre de travail révisé

### Priorité 1 — identifier la dernière interaction

Pour chaque petite attache du corridor réel, calculer le premier ancêtre après
lequel ses deux enfants ne sont plus distingués par aucun facteur futur du
séparateur retenu. Si ce temps est typiquement macroscopique, un bloc local
exact ne pourra pas porter la contraction.

### Priorité 2 — fermer le quotient de frontière

Construire une jauge de ports avec transition exacte de $`\Psi`$, ou un bloc
multi-update éliminant une orientation à sa dernière interaction. Calculer
alors seulement $`K,U,d`$. La projection en un pas de déficit
$`0.008259\ldots`$ est non composable et ne constitue pas cette fermeture.

### Priorité 3 — enrichir le diagnostic géométrique

Pour chaque bucket au rang réel, enregistrer le nombre de ports, le type
d'attache, une signature locale des contournements, le message et le rang de
dernière interaction. Après fermeture de la priorité 2 seulement, remplacer
le simple compte de coupes par

```math
D_{ij}^{\mathrm{diag}}
=
\sum_{r}
d_{\mathrm{T2}}(z_{r-1},z_r).
\qquad\text{(7.1)}
```

### Priorité 4 — seulement alors prouver l'abondance

Le lemme planaire doit viser une transformée de Laplace du nombre d'outlets T2
protégés visités avec $`|\Psi|\le M`$, sans borner le maximum de $`\Psi`$ et
sans chercher la loi limite complète du dendrogramme. RSW, finite energy et
quasi-multiplicativité ne deviennent pertinents qu'une fois le quotient T2
précisément défini.

### Priorité 5 — certificat à $`p_0`$

Après fermeture de la composition et de l'abondance : rationaliser les
poids, contrôler la troncature de $`\Psi`$ et certifier les intervalles à
$`p_0=161/200`$. Tester $`p=0.81`$ seulement après cette clôture.

## 8. Portes go/no-go mises à jour

| test | go | no-go local |
|---|---|---|
| géométrie | beaucoup d'orientations quittent un séparateur local | dernière interaction typiquement macroscopique |
| bord | quotient Markov-fermé et retours suffisants de $`\Psi`$ dans un compact | projection non fermée ou potentiel sans retours modérés |
| transfert | déficit positif après élimination exacte du twist | $`\lvert U\rvert=K`$ parce que le twist reste mesurable |
| composition | Doob rétrograde exact et Feynman--Kac sans facteur exponentiel caché | rayons individuels inférieurs à un mais produit incontrôlé |
| seuil | marge certifiée dominant troncature et mauvais blocs à $`0.805`$ | optimisation numérique avant fermeture du théorème global |

## 9. Bilan de ce cycle

### Établi ou certifié

- conventions Palm sans double biais ;
- partition exacte des paires connectées par les LCA réalisés ;
- lemme de criticalisation mono-bit et contre-exemple multiport exact ;
- transfert répliqué à environnement partagé pour la cellule E1+ ;
- coefficient tordu E1+ strictement inférieur à $`0.3`$ par intervalles
  rationnels ;
- impossibilité d'une contraction absolue brute uniforme sur les potentiels
  extérieurs non bornés ;
- domination de Feynman--Kac pour toute suite finie de transferts positifs,
  avec normalisation rétrograde automatique ;
- cellule quotient T2-Kruskal à rangs paramétrés, gagnantes marginalisées et
  bord polarisé ; elle n'est pas échantillonnée sous la Palm réelle ;
- no-go exact $`|U|=K`$ pour un update dont l'état complet révèle le twist ;
- projection en un pas de coefficient $`0.991774\ldots`$, explicitement non
  Markov-fermée et non composable.

### Indiqué numériquement

- croissance du nombre de petites coupes edge-disjoint de charge bornée le
  long du corridor réel ;
- avantage net de la famille de petites tailles sur le seul motif $`m=2`$.

### Toujours ouvert avant toute borne de weak recovery nouvelle

- rang de dernière interaction sous la Palm réelle ;
- quotient de ports Markov-fermé ou bloc multi-update d'élimination ;
- transfert T2 composable avec fusion, pivotalité, ports et $`\Lambda`$
  ancestraux ;
- loi et drift du potentiel extérieur ;
- composition exacte jusqu'au heat bath collapsed ;
- abondance du déficit sous la loi LCA-Palm réelle ;
- passage à $`\mathbb E[A_{I_LJ_L}]\to0`$.

Le résultat le plus important de ce cycle n'est donc pas la valeur
$`0.29399\ldots`$ prise isolément. C'est le double filtre logique : le
Feynman--Kac fini est disponible, mais un état fidèle trop riche rend son
déficit trivial. La route restante passe par une élimination de dernière
interaction Markov-fermée, puis seulement par l'abondance annealed sur le
corridor réel.
