# Pivot : dissipation $`L^2`$ du secteur impair

> [!NOTE]
> **Socle actif, pas feuille de route canonique.** Les identités de
> dissipation de cette note restent utilisées. Leur assemblage actuel est le
> [programme distance–entropie](35_DISTANCE_ENTROPIE_ERGODICITE.md).

Cette note prend acte des deux no-go du
[fichier 29](../diagnostics/29_AUDIT_FROID_PIVOT_RANGS_REELS.md) et remplace la priorité
« déficit Feynman--Kac local sur un état T2 » par une cible qui reste fidèle à
la dynamique hiérarchique : la perte d'énergie du caractère de paire sous des
**projections collapsed imbriquées**.

> [!IMPORTANT]
> La route locale à état fidèle est fermée : si le twist est mesurable depuis
> une transition, $`|U|=K`$ et son déficit vaut zéro. Le pivot ci-dessous
> capte les cancellations après sommation dans $`L^2`$ ; il ne prétend pas
> construire une chaîne de bord de dimension fixe.

> [!CAUTION]
> Aucun nouveau seuil de weak recovery n'est obtenu ici. Le petit tore exact
> a déjà livré un avertissement : la dissipation est réelle, mais très
> concentrée. Une preuve multiscalaire ne commence que si la cellule à deux
> mises à jour explique cette concentration et fournit une marge énergétique
> robuste pour la fonction effectivement propagée.

## 1. Pourquoi le transfert borné est la mauvaise fermeture

Avec $`b`$ ports binaires, l'espace d'orientations modulo flip global est

```math
\Omega_b
=
\{\pm1\}^b/\{\text{flip global}\},
\qquad
|\Omega_b|=2^{b-1}.
\qquad\text{(1.1)}
```

Un potentiel extérieur positif général $`\Psi`$, modulo multiplication par
une constante, possède donc $`2^{b-1}-1`$ degrés de liberté. Une jauge exacte
peut transporter son orbite complète et le secteur impair associé, mais une
simple re-jauge ne supprime pas cette dimension. L'élimination générique crée
en outre des interactions de tous ordres entre les ports.

Cela ne prouve pas qu'aucune compression spéciale du modèle n'existe. Cela
montre en revanche que la fermeture tabulaire locale naturelle devient une
inférence à frontière croissante, dont le coût est gouverné par la largeur
d'élimination ; voir la synthèse de
[Peyrard et al.](https://arxiv.org/abs/1506.08544). Le diagnostic de dernière
incidence du fichier 29 renforce ce constat : éliminer les petites attaches
réduit l'état, mais ne libère presque jamais par ce seul certificat
l'orientation globale porteuse du twist aux tailles testées.

## 2. Cadre fini exact

Fixons $`(O,D,i,j)`$ avec $`i,j`$ dans une même racine finale et écrivons

```math
\pi_D(d\sigma)=\nu_O(d\sigma\mid D),
\qquad
f_{ij}(\sigma)=\sigma_i\sigma_j.
\qquad\text{(2.1)}
```

Soit $`\Omega_D`$ l'espace des orientations de la racine commune. Pour des
sous-groupes croissants de masques de flips hiérarchiques

```math
H_0=\{0\}
\subset H_1
\subset\cdots\subset H_K,
\qquad\text{(2.2)}
```

notons $`[\sigma]_{H_k}`$ l'orbite de $`\sigma`$ sous $`H_k`$ et posons

```math
\mathcal F_k
=
\sigma(O,D,[\sigma]_{H_k}),
\qquad
M_k
=
\mathbb E_{\pi_D}[f_{ij}\mid\mathcal F_k].
\qquad\text{(2.3)}
```

Les tribus $`\mathcal F_k`$ décroissent et
$`M_0=f_{ij}`$, donc $`\|M_0\|_2^2=1`$. De plus,

```math
M_k
=
\mathbb E_{\pi_D}[M_{k-1}\mid\mathcal F_k].
\qquad\text{(2.4)}
```

### Théorème 2.1 — dissipation collapsed, statut : établi

Pour tout $`1\le k\le K`$,

```math
\|M_{k-1}\|_2^2
-
\|M_k\|_2^2
=
\|M_{k-1}-M_k\|_2^2.
\qquad\text{(2.5)}
```

Si $`\|M_{k-1}\|_2>0`$, définissons la perte relative

```math
\alpha_k
=
\frac{\|M_{k-1}-M_k\|_2^2}{\|M_{k-1}\|_2^2}
\in[0,1].
\qquad\text{(2.6)}
```

Avec la convention qu'un facteur nul annule le produit suivant,

```math
\|M_K\|_2^2
=
\|M_0\|_2^2\prod_{k=1}^K(1-\alpha_k)
=
\prod_{k=1}^K(1-\alpha_k)
\le
\exp\!\left(-\sum_{k=1}^K\alpha_k\right).
\qquad\text{(2.7)}
```

#### Preuve

L'application $`g\mapsto\mathbb E[g\mid\mathcal F_k]`$ est la projection
orthogonale de $`L^2(\pi_D)`$ sur les fonctions $`\mathcal F_k`$-mesurables.
Comme $`M_k`$ est la projection de $`M_{k-1}`$, le théorème de Pythagore donne
(2.5). Diviser par la norme précédente et télescoper donne (2.7). La mesure
$`\pi_D`$ n'a pas besoin d'être invariante sous $`H_k`$ : l'espérance
conditionnelle est la moyenne **pondérée** sur chaque orbite, exactement
comme dans le diagnostic fini.

Cette identité n'utilise ni état Markov quotient, ni transformée de Doob, ni
valeur absolue entrée par entrée. Elle peut donc être stricte alors que le
transfert fidèle vérifie $`|U|=K`$.

Les racines finales distinctes sont traitées avant cette filtration : leur
recoloration indépendante annule exactement la parité. Dans (2.3),
$`\mathcal F_0=\sigma(O,D,\sigma)`$ contient toutes les orientations de la
racine commune, ce qui justifie bien $`M_0=f_{ij}`$.

### Corollaire 2.2 — retour à la corrélation postérieure

Les projections préservent les constantes et la moyenne. Par conséquent,

```math
\pi_D(f_{ij})^2
\le
\|M_K\|_2^2.
\qquad\text{(2.8)}
```

Après désintégration selon $`D`$, le théorème H donne donc toujours

```math
Q_L
\le
\mathbb E\!\left[\|M_K\|_2^2\right].
\qquad\text{(2.9)}
```

Une contraction de norme d'opérateur sur tout $`L^2`$ est impossible, car
les constantes sont fixées. Le seul objet pertinent est la perte relative le
long de $`f_{ij}`$ et des fonctions $`M_k`$ qu'il engendre réellement.

## 3. Critère multiscalaire suffisant

Découpons le corridor en $`K_L`$ blocs collapsed imbriqués, associés à des
annuli ou outlets épaissis aux rangs réalisés. Posons

```math
a_{k,L}
=
\mathbb E[M_k^2],
\qquad
\mathcal E_{k,L}
=
\mathbb E[(M_{k-1}-M_k)^2].
\qquad\text{(3.1)}
```

### Proposition 3.1 — critère annealed, statut : établi conditionnellement

Supposons qu'il existe des nombres $`\gamma_{k,L}\in[0,1]`$ et
$`\varepsilon_{k,L}\ge0`$ tels que

```math
\mathcal E_{k,L}
\ge
\gamma_{k,L}a_{k-1,L}
-
\varepsilon_{k,L}.
\qquad\text{(3.2)}
```

Alors

```math
a_{K_L,L}
\le
a_{0,L}
\prod_{k=1}^{K_L}(1-\gamma_{k,L})
+
\sum_{r=1}^{K_L}
\varepsilon_{r,L}
\prod_{k=r+1}^{K_L}(1-\gamma_{k,L}).
\qquad\text{(3.3)}
```

Il suffit donc que $`\sum_k\gamma_{k,L}\to+\infty`$ et que la somme d'erreur
amortie dans (3.3) tende vers zéro. Avec l'annulation exacte des racines
distinctes et le contrôle des classes sous-critiques, (2.9) implique alors
$`Q_L\to0`$.

La formulation géométrique utile doit être **pondérée par l'énergie**. Pour
des événements de bloc $`G_k`$, il faudrait obtenir

```math
\mathcal E_{k,L}
\ge
\delta\,
\mathbb E[\mathbf1_{G_k}M_{k-1}^2]
-
\varepsilon_{k,L},
\qquad\text{(3.4)}
```

puis

```math
\mathbb E[\mathbf1_{G_k}M_{k-1}^2]
\ge
\rho\,a_{k-1,L}
-
\eta_{k,L}.
\qquad\text{(3.5)}
```

Compter des gadgets sous la Palm ordinaire ne suffit pas : une petite classe
d'environnements très polarisés peut porter presque toute l'énergie. Une
alternative équivalente consiste à contrôler directement une transformée de
Laplace des pertes relatives sous la loi annealed réelle.

## 4. Blocs à viser, sans objectif analytique irréaliste

La première cible n'est ni un gap uniforme du Gibbs sampler, ni une
classification de la Palm critique. Elle comporte deux niveaux.

### D1 — cellule exacte à deux mises à jour

Construire une cellule réelle avec trois bits d'orientation : branche
principale, petite attache et sibling ancestral. Conserver les facteurs
$`\Lambda e^{(1-\beta)\Lambda}`$ aux deux rangs réalisés et un potentiel
extérieur projectif. Pour le caractère de paire, calculer exactement :

- la moyenne carrée et les normes avant, entre et après les deux projections ;
- les deux pertes relatives de (2.6), la seconde étant calculée sur la
  fonction réellement produite par la première ;
- le déficit Feynman--Kac fidèle, attendu nul, à côté de la dissipation
  $`L^2`$ ;
- la pire dissipation sur un compact explicite de potentiels atteignables ;
- la norme complète, qui doit rester égale à un à cause des constantes.

Une seule mise à jour serait trop faible : elle ne teste ni la propagation de
la fonction, ni l'incrément de dissipation. Les projections de cette
filtration sont emboîtées et commutent ; la non-commutation concerne les heat
baths séquentiels non emboîtés, pas D1.

Le no-go Feynman--Kac fidèle et l'invariance des constantes étant déjà
établis en général, le module D1 ci-dessous ne les recalcule pas. Il se
concentre sur les deux pertes propagées et sur la famille exacte de potentiels
extérieurs atteints.

**Résultat D1.** Une cellule exacte est implémentée dans
[`two_step_projective_l2_cell.py`](../computations/two_step_projective_l2_cell.py).
Elle provient d'un véritable environnement $`L=4`$ à $`p=0.805`$ et comporte
trois bits de tailles $`1,4,1`$ : branche endpoint, attache, puis sibling
ancestral. Les deux nœuds ont leurs rangs réalisés
$`0.19324`$ et $`0.20258`$. On obtient

```math
\|M_0\|_2^2=1,
\qquad
\|M_1\|_2^2=0.9691266,
\qquad
\|M_2\|_2^2=0.8368704,
```

donc

```math
\alpha_1=0.0308734,
\qquad
\alpha_2=0.1364694.
\qquad\text{(4.1)}
```

Le flip du nouveau sibling, appliqué directement à $`f_{ij}`$, a une perte
exactement nulle puisqu'il ne contient aucun endpoint. La seconde perte vient
donc bien de l'élargissement appliqué à $`M_1`$, et non d'une réinitialisation
du caractère.

Les 128 cosets extérieurs de masse positive donnent 128 potentiels projectifs
effectivement atteints. Sur les 64 potentiels strictement positifs, de masse
postérieure totale $`0.94805`$, le ratio énergétique agrégé de la seconde
perte vaut $`0.14420`$ ; la perte relative minimale vaut $`0.00303`$ et sa
médiane $`0.18048`$. Mais les 64 potentiels de bord font tomber la marge
minimale globale à **zéro**. L'erreur de factorisation projective est au plus
$`1.43\times10^{-14}`$ et les erreurs pythagoriciennes au plus
$`1.12\times10^{-16}`$.

Ce résultat ferme la question algébrique minimale : une seconde dissipation
peut survivre après propagation. Il ne ferme pas la question probabiliste.
Le witness a été choisi après un scan, sa paire est à distance $`1`$ et ses
rangs sont très précritiques ; il n'est ni un échantillon Palm non biaisé, ni
une preuve d'abondance. La prochaine porte hiérarchique n'est donc plus de
chercher une meilleure cellule, mais de vérifier si ces potentiels intérieurs
portent une fraction d'énergie non négligeable sous une loi de paire
macroscopique non sélectionnée.

### D2 — petit tore exact, diagnostic exploratoire effectué

L'énumération de $`\pi_D`$ est possible sur un tore $`L=4`$, avec une jauge
globale et les rangs réels. Une paire est tirée uniformément, indépendamment
du dendrogramme, parmi celles à distance au moins $`L/2`$ ; conditionner
ensuite sur sa connexion dans une racine finale ne rajoute aucun poids de
coupe. Les projections sont calculées par moyennes d'orbites pondérées et le
profil de $`\alpha_k`$ mesure :

- les blocs collapsed imbriqués ;
- la projection collapsed finale ;
- la moyenne postérieure, qui doit rester sous toutes ces persistances.

Deux ordres de sweep séquentiels peuvent servir de contre-audit ultérieur,
mais ne font pas partie du premier diagnostic collapsed.

Les incertitudes sont calculées par environnement entier. Ce calcul ne
fournit aucune extrapolation en $`L`$ et son ordre « des feuilles vers le
LCA » n'est pas canonique : changer la packetisation change le nombre
effectif de niveaux, même si la projection finale reste la même.

**Premier résultat.** Le diagnostic est implémenté dans
[`nested_projection_l2_diagnostic.py`](../computations/nested_projection_l2_diagnostic.py).
À $`p=0.805`$, avec la graine $`20260726`$, 94 paires sur 96 à distance
exactement $`2=L/2`$ appartiennent à une même racine finale. Conditionnellement
à cette connexion, on observe :

| quantité | moyenne |
|---|---:|
| persistance collapsed finale | $`0.82045\pm0.02732`$ |
| carré de la moyenne postérieure | $`0.79657`$ |
| nombre de nœuds du corridor | $`6.404`$ |
| perte absolue totale | $`0.17955`$ ; médiane $`0.04372`$ |
| nombre effectif de niveaux, sachant une perte non nulle | $`1.547`$ ; médiane $`1.396`$ |
| part du plus grand niveau, sachant une perte non nulle | $`0.798`$ ; médiane $`0.834`$ |
| seconde perte pré-LCA, ratio énergétique agrégé | $`0.02837`$ |
| seconde perte pré-LCA, médiane relative | $`1.22\times10^{-5}`$ |

L'identité (2.5) est recalculée à partir de la norme de $`M_{k-1}-M_k`$ ;
son erreur maximale vaut $`1.67\times10^{-15}`$. La seconde perte est
positive au seuil numérique $`10^{-14}`$ dans $`61.9\%`$ des 84 corridors
qui possèdent deux niveaux strictement pré-LCA, mais ce taux masque une
forte concentration : seulement $`26.2\%`$ dépassent $`1\%`$ de perte
relative et les $`10\%`$ de cas les plus dissipatifs portent $`78.0\%`$ de
la seconde perte absolue.

Le fait que $`82.0\%`$ de la perte poolée soit pré-LCA n'est pas, seul, une
preuve d'étalement : dans 54 cas sur 94, ajouter le LCA n'augmente même pas le
rang du groupe de flips, donc sa perte est structurellement nulle. D2 est
ainsi un **avertissement**, pas un feu vert. D1 doit déterminer si la queue
forte correspond à une classe de bords géométriquement certifiable ; sinon
l'accumulation brute des projections est une impasse probable.

### D1-pop — audit non sélectionné, signal critique

Le test de population demandé après D1--D2 est maintenant implémenté dans
[`two_step_l2_population_diagnostic.py`](../computations/two_step_l2_population_diagnostic.py).
Pour chaque dendrogramme $`L=4`$, il tire d'abord une paire à distance
maximale, indépendamment de la hiérarchie, puis énumère toutes ses cellules
strict-arm consécutives. Chaque postérieur et chaque potentiel extérieur
atteint sont calculés exactement.

Deux graines, avec 32 puis 64 répétitions à $`p=0.805`$, donnent 302 cellules
sur 96 paires connectées. Globalement, le ratio énergétique de seconde perte
vaut $`0.03722`$, sa médiane relative $`7.87\times10^{-5}`$ et la queue
supérieure porte l'essentiel de la dissipation. Une borne uniforme reste
fausse : $`65.2\%`$ des cellules ont un potentiel atteint de bord à marge
nulle.

La nouveauté est la localisation en rang. Dans la fenêtre

```math
|q_{\mathrm{sup}}-q_c|\le0.02,
\qquad\text{(4.2)}
```

14 cellules sur 302 portent seulement $`4.13\%`$ de l'énergie entrante,
mais $`34.12\%`$ de la seconde perte. Leur ratio énergétique agrégé est
$`0.30779`$. Avec une fenêtre de largeur $`0.05`$, 42 cellules portent
$`58.35\%`$ de la perte, pour un ratio $`0.16897`$.

Les deux graines restent variables : la fenêtre étroite contient 10 cellules
sur 102 dans la première et 4 sur 200 dans la seconde. Ce calcul ne démontre
donc ni une fréquence asymptotique, ni une contraction. Il indique néanmoins
que la queue dissipative n'est pas répartie arbitrairement : elle est
fortement enrichie près de $`q_c`$. C'est la seule raison scientifique de
continuer la branche hiérarchique.

### D3 — annuli dyadiques, gelé

L'audit non biaisé ne valide pas des annuli génériques : hors de la fenêtre
critique, la dissipation reste trop concentrée. D3 demeure donc gelé sous sa
forme initiale. La seule réactivation permise est la
[sous-feuille des cellules critiques](33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) :
des annuli near-critical définis par des probabilités de traversée, contenant
deux fusions consécutives et deux routes concurrentes. La cible reste
(3.4)--(3.5), sous pondération énergétique, sans demander l'indépendance des
annuli ni un nombre total de ports borné.

## 5. Angles et detectability lemma

Pour une projection $`P`$ et une fonction entrante $`g`$, la quantité
$`\alpha`$ est le carré du sinus de l'angle de $`g`$ à l'image de $`P`$.
Des angles ou valeurs singulières sont donc utiles sur une petite famille
cible-spécifique de fonctions entrantes.

Ils ne doivent pas être pris sur tout l'espace : deux images de projections
partagent toujours les constantes, et une direction presque invariante peut
tourner d'un bloc au suivant. De même, le detectability lemma contrôle
l'approximation de la projection sur l'intersection par un produit de
projecteurs ; il ne montre pas que $`f_{ij}`$ a une petite composante dans
cette intersection. C'est précisément la corrélation cherchée. Voir
[Anshu--Arad--Vidick](https://arxiv.org/abs/1602.01210). Cet outil n'est donc
pas prioritaire avant un lemme annulaire cible-spécifique.

## 6. Portes go/no-go

### Continuer

- l'audit non sélectionné montre déjà un enrichissement énergétique net dans
  $`|q-q_c|\le0.02`$ ;
- dériver une identité locale de variance pour ces cellules et certifier une
  marge sur un compact de potentiels intérieurs ;
- construire un événement critique blindé, mesurable avant le calcul du
  signal, qui soutient la pondération énergétique (3.5) ;
- prouver que le corridor rencontre un nombre divergent de telles cellules,
  sans conditionnement gratuit sur un bord modéré.

### Arrêter ou pivoter de nouveau

- la seconde perte s'annule dès que la première projection a propagé le
  caractère ;
- sa moyenne est portée par une queue de bords trop rare pour satisfaire
  (3.5) ;
- toute marge disparaît sur un compact de bords modérés atteignables ;
- la queue near-critical observée à $`L=4`$ disparaît sous un audit de volume
  ou une majoration analytique ;
- le nombre de cellules critiques énergétiquement actives reste tendu ;
- les seuls bons annuli exigent un événement supplémentaire trop rare ;
- l'énergie se concentre sur des bords extrêmes que les motifs locaux ne
  visitent pas.

## 7. Verdict

Le déficit Feynman--Kac additif à séparateur borné est une impasse
structurelle pour l'état fidèle et une impasse probable pour toute fermeture
locale de dimension bornée.
La dissipation $`L^2`$ globale n'est pas une reformulation cosmétique : elle
conserve les cancellations que l'enveloppe $`|U|`$ détruit et possède une
identité télescopique exacte. Elle reste néanmoins risquée, car le contrôle du
bord réapparaît sous la forme de la fonction $`M_{k-1}`$ réellement propagée.

D2 a falsifié l'image optimiste d'une dissipation régulièrement répartie. D1
montre que la propagation à deux niveaux est algébriquement possible. D1-pop
confirme ensuite la queue rare, mais découvre qu'elle est fortement enrichie
près de $`q_c`$. La voie large d'accumulation collapsed est donc une impasse
probable ; la voie étroite des cellules critiques n'en est pas encore une.
Elle est désormais formulée dans le [fichier
33](33_SOUS_FEUILLE_ROUTE_CELLULES_CRITIQUES_L2.md) autour d'un seul verrou :
une minoration annealed, pondérée par l'énergie, répétée sur un nombre
divergent d'outlets near-critical. Une borne à $`p=0.805`$ par cette route
reste ouverte et ne doit pas être annoncée avant ce lemme.

Pour l'objectif plus court « dépasser rigoureusement $`0.8`$ », le jalon
rationnel A0 du [fichier 11](../results/non_hierarchical/11_TRIANGLE_BLOCK_SDPI.md) est analytiquement
plus court. Son certificat less-noisy exact est donné dans le [fichier
31](../archive/certificates/31_CERTIFICAT_RATIONNEL_A0.md). Cette piste parallèle ne valide pas la
dynamique hiérarchique ; elle évite de faire dépendre le premier gain
quantitatif du verrou multiscalaire.
