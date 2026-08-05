# Sous-feuille de route : cellules critiques et dissipation $`L^2`$

> [!NOTE]
> **Sous-programme actif.** Cette note décrit la cellule locale à deux
> projections. La stratégie globale et ses portes go/no-go sont désormais
> centralisées dans le
> [programme distance–entropie](35_DISTANCE_ENTROPIE_ERGODICITE.md).

## 1. Décision scientifique

La voie hiérarchique n'est pas entièrement dans une impasse, mais sa version
large l'est probablement.

- La criticalisation globale d'un corridor multiport est fausse au sens de
  Blackwell.
- Un état local assez riche pour être Markov-fermé rend le twist mesurable et
  annule le déficit Feynman--Kac.
- Les projections collapsed à rangs arbitraires dissipent surtout dans une
  queue rare ; les accumuler uniformément sur des annuli génériques n'est pas
  un programme crédible.
- En revanche, l'audit non sélectionné montre que cette queue est fortement
  enrichie lorsque deux fusions consécutives ont lieu près de $`q_c`$.

La décision est donc : **continuer uniquement sur les cellules critiques à
deux fusions, avec une pondération par l'énergie entrante**. Toute autre
extension hiérarchique reste gelée jusqu'à la fermeture du lemme central de
ce fichier.

Cette route est plus difficile que le certificat de triangle A0 et ne doit
pas retarder la première borne quantitative. Son intérêt est différent :
elle peut expliquer pourquoi la dynamique hiérarchique détruit la
corrélation à travers plusieurs échelles et, à terme, conduire à un seuil
intrinsèque plutôt qu'à un canal de comparaison choisi à la main.

## 2. Ce que dit réellement le diagnostic

Le module
[`two_step_l2_population_diagnostic.py`](../computations/two_step_l2_population_diagnostic.py)
tire, indépendamment du dendrogramme, une paire à distance maximale sur le
tore $`L=4`$. Il énumère toutes les cellules formées par deux fusions
strict-arm consécutives avant le LCA et calcule exactement le postérieur
complet de chacune.

Deux graines indépendantes donnent :

| quantité | graine 20260729 | graine 20260730 | agrégat |
|---|---:|---:|---:|
| paires connectées | 32 | 64 | 96 |
| cellules admissibles | 102 | 200 | 302 |
| cellules dans $`\lvert q_{\mathrm{sup}}-q_c\rvert\le0.02`$ | 10 | 4 | 14 |
| ratio énergétique global | $`0.04597`$ | $`0.03265`$ | $`0.03722`$ |
| part d'énergie entrante dans la fenêtre | $`8.70\%`$ | $`1.74\%`$ | $`4.13\%`$ |
| part de seconde perte dans la fenêtre | $`56.6\%`$ | $`17.6\%`$ | $`34.1\%`$ |
| ratio énergétique dans la fenêtre | $`0.299`$ | $`0.330`$ | $`0.308`$ |

Sur l'agrégat, la fenêtre reçoit donc environ huit fois plus de perte que ne
le prédirait sa seule part d'énergie. Pour la fenêtre de largeur $`0.05`$,
42 cellules portent $`58.35\%`$ de la perte et ont un ratio énergétique
$`0.16897`$.

L'agrégat est un pool descriptif de cellules recouvrantes, sans erreur-type
par environnement. Ce signal est compatible avec l'idée que la fusion
near-critical est le lieu
où une nouvelle orientation de cluster perturbe le plus la fonction déjà
propagée. Il ne démontre rien en volume croissant : les cellules se
recouvrent, les deux graines sont variables et $`65.2\%`$ des cellules
agrégées possèdent encore un potentiel atteint de bord à marge nulle.

## 3. Dynamique à analyser

Conditionnons par le dendrogramme non marqué $`D`$ et par ses rangs réels.
La loi de Nishimori induit une mesure $`\pi_D`$ sur les orientations de ses
clusters. Pour une paire $`i,j`$, gardons seulement les deux bras allant des
feuilles à leur LCA et les siblings qui s'y attachent.

La dynamique pertinente est un heat bath hiérarchique exact :

1. aux feuilles, mettre à jour l'orientation d'un singleton ou du plus petit
   cluster encore actif, comme dans Glauber ;
2. à une fusion intermédiaire, rééchantillonner conjointement l'orientation
   du sous-arbre qui vient de s'attacher, conditionnellement au complément ;
3. à une racine finale, rééchantillonner l'orientation de toute la composante,
   comme dans Swendsen--Wang.

Pour la preuve, il ne faut pas étudier le temps de mélange d'un sweep
non commutatif. On utilise les projections collapsed emboîtées

```math
P_k g
=
\mathbb E_{\pi_D}[g\mid\mathcal F_k],
\qquad
M_k=P_k f_{ij},
\qquad
f_{ij}(\sigma)=\sigma_i\sigma_j.
\qquad\text{(3.1)}
```

Ici $`\mathcal F_k`$ oublie un ensemble croissant d'orientations depuis les
feuilles vers les racines. Chaque bloc est donc parfaitement mélangé au
moment où il est éliminé. L'identité

```math
\|M_{k-1}\|_2^2-\|M_k\|_2^2
=
\|M_{k-1}-M_k\|_2^2
\qquad\text{(3.2)}
```

remplace un théorème spectral global irréaliste. Elle respecte exactement le
caractère Glauber aux feuilles et Swendsen--Wang aux racines.

## 4. Fenêtres critiques sans exposants non disponibles

Il ne faut pas importer les exposants exacts de la percolation de sites sur
le réseau triangulaire. Pour une échelle $`r`$, définir plutôt deux rangs
$`q_r^-<q_c<q_r^+`$ par des probabilités de traversée fixées, par exemple

```math
\mathbb P_{q_r^-}(\text{traversée de }[0,2r]\times[0,r])=\eta,
\qquad
\mathbb P_{q_r^+}(\text{traversée de }[0,2r]\times[0,r])=1-\eta,
\qquad\text{(4.1)}
```

avec $`0<\eta<1/2`$ fixé. Les fenêtres
$`W_r=[q_r^-,q_r^+]`$ sont near-critical par définition et leur contrôle
peut reposer sur RSW, la monotonie et la longueur de corrélation, sans valeur
explicite de l'exposant $`\nu`$.

Dans un annulus $`A(r,2r)`$ rencontré par l'un des deux bras, appelons
**cellule critique blindée** une paire de fusions consécutives $`u\prec v`$
telle que :

- le rang supérieur $`q_v`$ appartient à $`W_r`$ ;
- deux routes concurrentes relient les interfaces intérieure et extérieure ;
- une séparation duale locale empêche le reste du corridor de déterminer
  immédiatement l'orientation du nouveau sibling ;
- la cellule est déclarée à partir de la géométrie et des rangs avant de
  regarder la valeur de $`f_{ij}`$ ou sa dissipation.

Cette définition doit rester assez souple pour autoriser des coupes de taille
non bornée. Exiger un nombre total de ports borné à toutes les échelles
recréerait l'impasse déjà identifiée.

## 5. Le théorème hiérarchique à viser

Soient $`r_1<\cdots<r_{K_L}`$ des échelles séparées, avec
$`K_L\asymp\log L`$, et $`G_k`$ l'événement qu'une cellule critique blindée
soit disponible à l'échelle $`r_k`$. Posons

```math
a_{k,L}=\mathbb E[M_k^2],
\qquad
\mathcal E_{k,L}=\mathbb E[(M_{k-1}-M_k)^2].
\qquad\text{(5.1)}
```

Le résultat central réaliste est le couple de lemmes suivant.

**Lemme local pondéré.** Il existe $`\delta(p)>0`$ et une erreur sommable
$`\varepsilon_{k,L}`$ tels que

```math
\mathcal E_{k,L}
\ge
\delta(p)\,
\mathbb E[\mathbf1_{G_k}M_{k-1}^2]
-\varepsilon_{k,L}.
\qquad\text{(5.2)}
```

**Lemme d'occupation énergétique.** Il existe $`\rho(p)>0`$ et une erreur
sommable $`\eta_{k,L}`$ tels que

```math
\mathbb E[\mathbf1_{G_k}M_{k-1}^2]
\ge
\rho(p)a_{k-1,L}-\eta_{k,L}.
\qquad\text{(5.3)}
```

Ces deux énoncés donnent la récursion
$`a_{k,L}\le(1-\delta\rho)\,a_{k-1,L}+(\delta\eta_{k,L}+\varepsilon_{k,L})`$,
d'où

```math
a_{K_L,L}
\le
\bigl(1-\delta(p)\rho(p)\bigr)^{K_L}a_{0,L}
+
\sum_{r=1}^{K_L}
\bigl(1-\delta(p)\rho(p)\bigr)^{K_L-r}
\bigl(\delta(p)\eta_{r,L}+\varepsilon_{r,L}\bigr).
\qquad\text{(5.4)}
```

La conclusion $`a_{K_L,L}\to0`$ exige que la **somme amortie** des
erreurs — le second terme de (5.4) — tende vers zéro en $L$ ; la seule
sommabilité en $k$ à $L$ fixé ne suffit pas. C'est la formulation (3.3)
du [fichier 30](30_PIVOT_DISSIPATION_L2_SECTEUR_IMPAIR.md), qui prévaut.

Après l'annulation exacte des paires appartenant à deux racines finales
distinctes et la décroissance des connexions uniformément sous-critiques, le
théorème pairwise du fichier
[03](../foundations/03_HIERARCHICAL_WEAK_RECOVERY.md) implique l'absence de weak recovery.

L'énoncé (5.3), et non le calcul local, est le vrai verrou. Une densité Palm
positive de cellules ne suffit pas : il faut les voir sous la mesure inclinée
par $`M_{k-1}^2`$.

## 6. Décomposition du verrou en tâches faisables

### H1 — formule locale de variance, faisable

Pour deux fusions consécutives, écrire la seconde perte sous la forme

```math
\mathbb E\!\left[
r_v(1-r_v)
\bigl(M_{k-1}^{(+)}-M_{k-1}^{(-)}\bigr)^2
\right],
\qquad\text{(6.1)}
```

où $`r_v`$ est la probabilité heat bath du nouveau sibling conditionnellement
au complément. Il faut ensuite isoler un compact
$`r_v\in[\epsilon,1-\epsilon]`$ et une condition de sensibilité relative de
$`M_{k-1}`$. C'est une identité de variance conditionnelle, pas un calcul de
rayon spectral.

### H2 — certificat local sur un compact, faisable

Paramétrer la cellule par les rangs $`q_u,q_v`$, les incidences relatives et
le rapport projectif du potentiel extérieur. Sur une famille critique
blindée explicite, certifier (5.2) par arithmétique d'intervalles. Il est
acceptable que la marge s'annule sur le bord du simplexe : la masse de ce
bord doit être traitée par H3, pas par une fausse borne uniforme.

### H3 — dépolarisation par blindage, difficile mais ciblé

Utiliser les deux routes concurrentes et la séparation duale pour montrer
que, sous la mesure inclinée par $`M_{k-1}^2`$, le log-rapport du sibling a
une probabilité uniformément positive de rester dans un intervalle borné.
Une preuve par exploration de l'annulus et propriété de Markov spatiale est
plus plausible qu'une classification de toute la Palm du LCA.

### H4 — outlets near-critical répétés, difficile mais standardisable

Construire $`c\log L`$ annuli séparés et révéler leurs variables de
l'intérieur vers l'extérieur. RSW near-critical et une version conditionnelle
de la séparation des bras doivent donner une chance uniforme de produire
$`G_k`$. Il n'est pas nécessaire de rendre les annuli indépendants : une
minoration conditionnelle avant révélation suffit pour itérer (5.2)--(5.3).

### H5 — transfert sous la Palm à deux points, dernier pont

Le LCA d'une paire réalisée est pondéré par $`N_\rho`$, tandis que le flux
pré-saut porte $`mN_\rho`$. Le transfert doit être écrit avec ces deux
conventions et montrer que le dévoilement des annuli ne réintroduit pas un
facteur $`m`$. Ce contrôle est plus étroit que la loi complète du corridor :
seuls les événements $`G_k`$ et leur énergie sont nécessaires.

### Appuis disponibles et limite du transfert

Cette architecture n'est pas inventée sans précédent. Damron--Sapozhnikov
obtiennent des contrôles de moments pour le nombre d'outlets par annulus, puis
des lois de grande échelle pour leur nombre dans une boîte en invasion
percolation bidimensionnelle : [outlets et IIC
multi-bras](https://arxiv.org/abs/0903.4496), [théorèmes limites pour les
outlets](https://arxiv.org/abs/1005.5696). Nolin organise la percolation
near-criticale par longueurs caractéristiques et comparabilité des événements
de bras, sans commencer par une limite de continuum : [near-critical
percolation](https://arxiv.org/abs/0711.4948). Enfin, la construction de la
MST plane par l'ensemble near-critical est développée par
Garban--Pete--Schramm : [near-critical et dynamique](https://arxiv.org/abs/1305.5526),
[MST et invasion](https://arxiv.org/abs/1309.0269).

Ces résultats justifient le vocabulaire d'outlets et l'espoir d'un nombre
logarithmique d'échelles actives. Ils ne prouvent pas H3--H5 : certains sont
formulés pour la percolation de sites triangulaire ou la percolation d'arêtes
carrée, aucun n'inclut l'inclinaison Nishimori par $`M_{k-1}^2`$, ni la Palm
d'une paire fixée, ni les marques de frontière du GSBM. La nouveauté requise
est exactement ce transfert énergétique, pas la théorie near-criticale de
base.

## 7. Ordre de travail et portes d'arrêt

1. Dériver H1 symboliquement et vérifier qu'il reproduit D1 cellule par
   cellule.
2. Choisir une seule géométrie blindée near-critical et fermer H2 sur un
   compact projectif explicite.
3. Tenter H3 sur un annulus unique, sans conditionnement LCA global.
4. Seulement si H3 réussit, établir H4 puis le transfert Palm H5.
5. Extraire enfin les constantes $`\delta(p),\rho(p)`$ et tester leur domaine
   en $`p`$.

Arrêter cette voie si l'un des faits suivants est démontré :

- même sous un blindage à deux routes, la mesure énergétique se concentre
  arbitrairement près de $`r_v=0`$ ou $`1`$ ;
- toute cellule qui satisfait H2 exige un événement multi-bras sommable sur
  les échelles ;
- le nombre de cellules énergétiquement actives reste tendu au lieu de
  diverger ;
- la Palm à deux points détruit toute minoration conditionnelle d'annulus.

## 8. Seuil que cette méthode pourrait identifier

Pour chaque $`p`$, définissons $`\delta(p)`$ comme la meilleure constante du
lemme local H2 et $`\rho(p)`$ comme la meilleure constante du lemme
d'occupation H3--H5. Le seuil naturel de la méthode serait

```math
p_{\mathrm{hier}}
=
\sup\{p:\delta(p)\rho(p)>0\}.
\qquad\text{(8.1)}
```

Il serait prématuré de lui attribuer une valeur numérique. Le petit tore dit
seulement que $`\delta(0.805)`$ n'est pas algébriquement forcée à zéro près de
$`q_c`$. Obtenir le seuil exact exigerait en plus de montrer que
$`\delta(p)\rho(p)=0`$ au-delà de la transition, ce qui est beaucoup plus
ambitieux que la borne d'impossibilité. La cible raisonnable est d'abord une
constante strictement positive sur un intervalle de $`p`$, puis son
optimisation.

## 9. Verdict opérationnel

La route n'est donc pas morte, mais elle ne dispose plus que d'une porte
étroite. Le signal critique observé justifie H1--H3. Il ne justifie ni un
théorème multiscalaire annoncé d'avance, ni une campagne de simulations plus
grandes, ni une nouvelle fermeture locale générique. Le prochain résultat
qui compterait est une version à un annulus de (5.2)--(5.3) ; sans elle, la
voie hiérarchique doit être classée comme impasse pour l'amélioration de la
borne.
