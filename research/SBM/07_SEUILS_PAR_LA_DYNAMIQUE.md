# La dynamique retrouve-t-elle les seuils du SBM classique ?

Cette note répond à la question de calibration préalable au GSBM
triangulaire :

> avant d'attaquer le Geometric SBM, vérifier si la dynamique
> hiérarchique **telle que définie** retrouve les seuils théoriques du SBM
> classique — weak recovery au niveau $\beta$ de percolation, exact
> recovery à $\beta=0$, où la dynamique dégénère en Glauber.

La réponse tient en un tableau, détaillé ensuite régime par régime.
« Retrouver » peut signifier trois choses de force croissante :
**(C)** une identité de calibration exacte, **(E)** un benchmark à
l'équilibre (stationnarité), **(D)** un théorème sur la dynamique
elle-même (nombre fini de sweeps, point de départ réaliste).

| régime | seuil théorique | niveau $\beta$ | (C) calibration | (E) équilibre | (D) dynamique |
|---|---|---|---|---|---|
| weak recovery | $`d\theta^2=1`$ (KS) | $`\beta_c^{\mathrm{geom}}`$ | **établie** : $`\beta_\chi=\beta_c^{\mathrm{geom}}\Leftrightarrow d\theta^2=1`$ | **établi** sur le broadcast (Jacobien $`\theta^2`$ + sandwich) ; port fini : no-go à KS | **ouvert** (porte SBM-DYN1) |
| almost exact | $`\lambda_n\to\infty`$ | $`\beta_{c,n}^{\mathrm{geom}}\to0`$ | **établie** : la coupe géométrique s'écrase sur Glauber ([05 §6](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md)) | **établi** : erreur par sommet $`\to0\Leftrightarrow\lambda_n\to\infty`$ (§5) | partiel : stabilité par sommet (§5) ; initialisation ouverte |
| exact | $`(\sqrt A-\sqrt B)^2=2`$ | $`\beta=0`$ (Glauber) | **établie** : $`\beta=0`$ séquentiel $=$ Glauber ([04 §4](04_DYNAMIQUE_HIERARCHIQUE.md)) | **établi ici** : stabilité de la vérité sous un sweep $\Leftrightarrow$ seuil (§4) | correction depuis un départ presque exact : littérature (graph splitting) ; Glauber seul : **ouvert** |

Aucun de ces énoncés ne constitue une preuve nouvelle des théorèmes de
seuil du SBM (Mossel–Neeman–Sly, Massoulié, Abbe–Bandeira–Hall) : la
colonne (E) vérifie que la dynamique est **calibrée juste**, c'est le
prérequis demandé avant tout transfert au GSBM.

## 1. Weak recovery au niveau de percolation

Quatre faits, tous détaillés dans [02](02_DEUX_DENDROGRAMMES_A_BETA_C.md),
[03](03_PREUVE_DU_SEUIL_WEAK_RECOVERY.md) et le
[pilote 37](../hierarchical-swendsen-wang/active/37_PILOTE_SBM_GIBBS_HIERARCHIQUE.md).

1. **Calibration exacte — établie.** Le temps informationnel
   $`\beta_\chi`$ ($`q_{\beta_\chi}=\theta^2`$) et la coupe géométrique
   $`\beta_c^{\mathrm{geom}}`$ ($`q=1/d`$) coïncident exactement au seuil
   de Kesten–Stigum :
   $`\beta_\chi=\beta_c^{\mathrm{geom}}\Longleftrightarrow d\theta^2=1`$.
2. **Jacobien répliqué — établi.** Deux répliques et deux dendrogrammes
   **indépendants**, chacun entièrement marginalisé, transmettent
   $`\theta^2`$ par arête dans le secteur overlap — pour **tout** niveau
   de coupe marginalisé. Le niveau $`\beta_c^{\mathrm{geom}}`$ est
   distingué par la géométrie (blocs critiques de degré moyen un), pas
   par l'information : c'est la calibration 1 qui le relie au seuil.
3. **Fermeture sur le broadcast — établie (littérature-dépendante).**
   La linéarisation $`d\theta^2`$ plus le sandwich classique
   (second moment EKPS / percolation d'information $\chi^2$) retrouvent
   le seuil $`d\theta^2=1`$, égalité comprise.
4. **SBM fini — obstruction identifiée.** Le port global recouple les
   racines ; à KS, la comparaison perturbative full-$D$ au broadcast est
   **réfutée** (deux géantes d'Erdős–Rényi anticorrélées,
   [39](../hierarchical-swendsen-wang/active/39_PORT_GLOBAL_SBM_RECOVERY.md)).
   Le candidat subcritique $`d\theta<1`$ reste ouvert.

Ce qui manque pour la colonne (D) : montrer qu'un **nombre fini de sweeps
hiérarchiques** à la coupe $`\beta_c^{\mathrm{geom}}`$ produit (ou détruit)
l'overlap au seuil KS. C'est la porte SBM-DYN1 (§6). Les deux no-gos du
pilote restent les garde-fous : dendrogramme figé $\Rightarrow d\theta$ ;
coupe partagée $\Rightarrow\theta^2+q(1-\theta)^2/(1-q)$.

## 2. Le régime divergent écrase la coupe sur Glauber

Le pont entre les deux extrémités de la question est déjà démontré dans
[05 §6](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md) : dès que
$`\lambda_n\to\infty`$,

```math
\beta_{c,n}^{\mathrm{geom}}
=
-\frac{\log(1-2/a_n)}{\log(a_n/b_n)}
\longrightarrow
0.
\qquad\text{(2.1)}
```

Le niveau de percolation **converge de lui-même vers l'extrémité
$`\beta=0`$** de la famille quand le signal diverge : demander « la
percolation pour la weak recovery, Glauber pour l'exact recovery » n'est
pas un changement de dynamique mais la limite continue de la même
calibration. À $`\beta=0`$ avec balayage séquentiel, la dynamique est
exactement le bain de Glauber de la postérieure
([04, (4.1)](04_DYNAMIQUE_HIERARCHIQUE.md)).

## 3. L'expérience locale à $`\beta=0`$

Notations de [05 §2](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md) :
$`p_n=a_n/n`$, $`q_n=b_n/n`$,
$`\rho_n=\sqrt{p_nq_n}+\sqrt{(1-p_n)(1-q_n)}`$. Sous l'a priori i.i.d.,
le bain de Glauber au sommet $v$, tous les autres labels étant à la
vérité $X_{-v}$, retire $`\sigma_v`$ selon la conditionnelle exacte ; sa
probabilité de flip est

```math
\varepsilon_v
=
\frac1{1+e^{\Delta_v}},
\qquad
\Delta_v
=
\sum_{u\ne v}X_uX_v
\left[
A_{uv}\log\frac{p_n}{q_n}
+(1-A_{uv})\log\frac{1-p_n}{1-q_n}
\right],
\qquad\text{(3.1)}
```

où $`A_{uv}`$ est l'indicatrice d'arête. $`\Delta_v`$ est le LLR exact de
l'expérience oracle de [05 §2](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md).

### Proposition 3.1 — identité d'affinité du sweep, statut : établi

Sous la vérité,

```math
\boxed{
\mathbb E\bigl[e^{-\Delta_v/2}\bigr]
=
\rho_n^{\,n-1},
}
\qquad\text{et}\qquad
\mathbb E[\varepsilon_v]
\le
\rho_n^{\,n-1}.
\qquad\text{(3.2)}
```

**Preuve.** Les $`n-1`$ arêtes potentielles incidentes sont
indépendantes conditionnellement aux labels. Pour un voisin de même
classe, le facteur vaut
$`p_n\sqrt{q_n/p_n}+(1-p_n)\sqrt{(1-q_n)/(1-p_n)}=\rho_n`$ ; pour une
classe opposée, $`q_n\sqrt{p_n/q_n}+(1-q_n)\sqrt{(1-p_n)/(1-q_n)}=\rho_n`$
également. Le produit donne $`\rho_n^{n-1}`$. La borne sur
$`\varepsilon_v`$ vient de $`1/(1+e^x)\le e^{-x/2}`$ pour tout
$`x\in\mathbb R`$. $\square$

L'identité est vérifiée à précision machine par
[sbm_glauber_stability_benchmark.py](../hierarchical-swendsen-wang/computations/sbm_glauber_stability_benchmark.py).

## 4. Exact recovery : stabilité de la vérité sous un sweep

Régime logarithmique $`a_n=A\log n`$, $`b_n=B\log n`$, $`A>B>0`$ ; alors
$`\rho_n^{n-1}=n^{-(\sqrt A-\sqrt B)^2/2+o(1)}`$
([05, (2.4)–(2.6)](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md)).

### Proposition 4.1 — stabilité au-dessus du seuil, statut : établi

Si $`(\sqrt A-\sqrt B)^2>2`$, un sweep séquentiel complet de la dynamique
à $`\beta=0`$ démarré de la vérité $X$ la laisse inchangée avec
probabilité $`1-n^{1-(\sqrt A-\sqrt B)^2/2+o(1)}\to1`$. Il en va de même
pour tout nombre $`T_n=n^{o(1)}`$ de sweeps.

**Preuve.** Soit $`\tau`$ le premier site flippé du balayage
$`v_1,\dots,v_n`$. Sur l'événement $`\{\tau=v_k\}`$, la configuration au
tour de $`v_k`$ est encore $X$, donc
$`\mathbb P(\tau=v_k)\le\mathbb E[\varepsilon_{v_k}]\le\rho_n^{n-1}`$
par la proposition 3.1. L'union sur les $n$ sites (et les $`T_n`$
sweeps) donne
$`\mathbb P(\text{au moins un flip})\le T_n\,n\,\rho_n^{n-1}
=T_n\,n^{1-(\sqrt A-\sqrt B)^2/2+o(1)}`$. $\square$

### Corollaire 4.2 — version Nishimori, statut : établi

Par l'identité de Nishimori, $(X,A)$ a la loi d'un couple
$(\sigma,A)$ avec $`\sigma\sim\mu_A`$ : une **réplique postérieure** est
donc aussi, avec la même probabilité, laissée invariante par le sweep.
La vérité est un point fixe typique de la dynamique de Glauber, au sens
fort, exactement quand l'exact recovery est possible.

### Proposition 4.3 — instabilité sous le seuil, statut : établi sous références standard

Si $`(\sqrt A-\sqrt B)^2<2`$ (et $`A,B>0`$), alors
$`\mathbb P(\Delta_v\le0)=n^{-(\sqrt A-\sqrt B)^2/2+o(1)}`$ (borne
inférieure de Chernoff–Cramér pour des sommes de Bernoulli
logarithmiques ; Abbe–Bandeira–Hall, lemmes de grandes déviations ;
Mossel–Neeman–Sly). Le nombre de sommets où le test local à la vérité
échoue a donc une espérance $`n^{1-(\sqrt A-\sqrt B)^2/2+o(1)}\to\infty`$,
et un argument de second moment standard (les $`\Delta_v`$ ne partagent
deux à deux qu'une arête) le rend $`\to\infty`$ en probabilité : la
vérité n'est **pas** stable, et l'exact recovery est impossible
([05, (4.5)](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md), théorème cité).

### Portée exacte de ces énoncés

La stabilité 4.1 est un benchmark de la dynamique **définie dans ce
dossier**, pas une preuve d'achievability : atteindre l'exact recovery
demande encore une initialisation presque exacte, puis une étape de
correction. La littérature ferme cette étape par graph splitting et
raffinement local (Abbe–Sandon ; Mossel–Neeman–Sly) ; la proposition 4.1
est exactement l'analyse de stabilité de ce raffinement, réécrite comme
propriété du sweep de Glauber. Que le **Glauber seul**, démarré d'une
configuration froide ou aléatoire, atteigne le seuil reste ouvert
(paysage vitreux ; porte SBM-DYN0).

## 5. Almost exact : le même mécanisme, par sommet

Dans tout régime divergent ($`a_n,b_n=o(n)`$), la proposition 3.1 donne
pour l'erreur moyenne par sommet

```math
\bar\varepsilon_n
\le
\rho_n^{\,n-1}
=
e^{-(1+o(1))(\sqrt{a_n}-\sqrt{b_n})^2/2},
\qquad\text{(5.1)}
```

donc

```math
\bar\varepsilon_n\longrightarrow0
\quad\Longleftrightarrow\quad
(\sqrt{a_n}-\sqrt{b_n})^2\to\infty
\quad\Longleftrightarrow\quad
\lambda_n\to\infty,
\qquad\text{(5.2)}
```

par l'encadrement $`\lambda_n\le(\sqrt{a_n}-\sqrt{b_n})^2\le2\lambda_n`$
([05, (3.3)](05_ALMOST_EXACT_ET_EXACT_RECOVERY.md)). L'équivalence de
droite avec la possibilité de l'almost exact recovery est le théorème
classique (Mossel–Neeman–Sly). Le sens direct de (5.2) — l'erreur par
sommet s'annule au-dessus du seuil — est démontré ici ; le sens réciproque
au niveau de la borne (5.1) demande la borne inférieure de Cramér comme
en 4.3. Un sweep démarré de la vérité flippe donc une **fraction**
$`o_{\mathbb P}(1)`$ des sommets exactement dans le régime almost exact
(Markov sur (5.1) pour le sens direct).

## 6. Portes dynamiques restantes

| porte | énoncé à établir | statut |
|---|---|---|
| SBM-DYN0 | Glauber ($`\beta=0`$) depuis un départ aléatoire atteint l'almost/exact recovery au seuil | ouvert (littérature : deux étapes avec graph splitting) |
| SBM-DYN1 | un nombre fini de sweeps hiérarchiques à $`\beta_c^{\mathrm{geom}}`$ voit le seuil KS $`d\theta^2=1`$ | ouvert ; équilibre calibré (§1), port fini : no-go à KS |
| SBM-DYN2 | interpolation : le sweep à la coupe $`\beta_{c,n}^{\mathrm{geom}}\to0`$ relie continûment SBM-DYN1 à la stabilité §4 | ouvert (formulé ici) |

Le transfert au GSBM triangulaire ne doit être entrepris qu'après SBM-DYN1
au minimum : c'est l'analogue exact de la porte TRI du
[statut canonique](../hierarchical-swendsen-wang/CURRENT_STATUS.md), et le
[problème central 42](../hierarchical-swendsen-wang/foundations/ancestral/42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md)
en est la version géométrique (chaîne des $`\Lambda_v`$ au-dessus du LCA).

## 7. Vérifications reproductibles

Le script
[sbm_glauber_stability_benchmark.py](../hierarchical-swendsen-wang/computations/sbm_glauber_stability_benchmark.py)
calcule **exactement** (énumération binomiale, sans Monte-Carlo) :

- l'identité d'affinité (3.2) à précision machine ;
- la probabilité de flip $`\mathbb E[\varepsilon_v]`$ et l'espérance de
  flips d'un sweep $`n\,\mathbb E[\varepsilon_v]`$ ;
- la localisation du croisement $`n\,\mathbb E[\varepsilon_v]\approx1`$
  autour de $`(\sqrt A-\sqrt B)^2=2`$ pour $n$ croissant ;
- la décroissance de $`\bar\varepsilon_n`$ avec $`\lambda_n`$.

Les tests unitaires sont dans
[test_sbm_glauber_stability_benchmark.py](../hierarchical-swendsen-wang/computations/test_sbm_glauber_stability_benchmark.py).
