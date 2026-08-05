# Calibration hiérarchique et preuve sur le broadcast

## 1. Énoncé du benchmark

Pour le SBM symétrique binaire à degré borné, posons

```math
d=\frac{a+b}{2},
\qquad
\theta=\frac{a-b}{a+b},
\qquad
\lambda=d\theta^2.
\qquad\text{(1.1)}
```

La transition de weak recovery est

```math
\boxed{\lambda_c=1}
\qquad\Longleftrightarrow\qquad
\boxed{(a-b)^2=2(a+b)}.
\qquad\text{(1.2)}
```

Dans la formulation binaire standard, aucune procédure ne garde un overlap
positif pour $\lambda\le1$ ; une procédure efficace existe pour
$\lambda>1$. Le but de cette note est de montrer exactement
où (1.2) apparaît dans le calcul hiérarchique, puis de distinguer cette
calibration de la preuve complète sur le graphe fini.

## 2. Le message Gibbs sur l'arbre

Considérons le broadcast binaire sur un arbre
$T\sim\mathrm{PGW}(d)$. Si $M_{v\to u}$ est la magnétisation transmise par
le sous-arbre enraciné en $v$, la récursion exacte est

```math
M_{u\to\mathrm{par}(u)}
=
\tanh\left[
\sum_{v\text{ enfant de }u}
\mathrm{atanh}\left(\theta M_{v\to u}\right)
\right].
\qquad\text{(2.1)}
```

Soit $M_t$ la magnétisation postérieure de la racine lorsque les spins à
profondeur $t$ sont révélés. La quantité de reconstruction est

```math
q_t=\mathbb E[M_t^2].
\qquad\text{(2.2)}
```

Elle est simultanément :

- le carré d'un message Gibbs ;
- une information $\chi^2$ de la racine vers le bord ;
- un overlap entre deux répliques ;
- le secteur quadratique naturel de la weak recovery.

## 3. Pourquoi deux dendrogrammes donnent $\theta^2$

La [note 02](02_DEUX_DENDROGRAMMES_A_BETA_C.md) a établi que, pour toute
coupe $\beta$ exactement marginalisée,

```math
\sum_{b_1,b_2}
\pi_{b_1}\pi_{b_2}c_{b_1}c_{b_2}
=
\theta^2.
\qquad\text{(3.1)}
```

Chaque branche transmet donc :

```math
\theta
\quad\text{dans une réplique},
\qquad
\theta\times\theta=\theta^2
\quad\text{dans le secteur overlap}.
\qquad\text{(3.2)}
```

Le nombre moyen de branches étant $d$, la linéarisation à l'origine vaut

```math
\boxed{
\mathcal L_{\mathrm{overlap}}=d\theta^2=\lambda.
}
\qquad\text{(3.3)}
```

La coupe à $\beta_c^{\mathrm{geom}}$ rend les blocs critiques et autorise
l'élimination
hiérarchique de (2.1). Le coefficient (3.3), lui, ne dépend pas de l'ordre
dans lequel les sommes exactes ont été effectuées.

## 4. La linéarisation seule ne constitue pas une preuve

Au voisinage de zéro, (2.1) suggère

```math
q_{t+1}
=
\lambda q_t+O(q_t^2).
\qquad\text{(4.1)}
```

Mais (4.1) ne contrôle ni les grands messages ni l'accumulation du reste.
La fermeture rigoureuse utilise deux bornes globales. Posons

```math
\ell_t(\lambda)
=
\left(
\sum_{s=0}^t\lambda^{-s}
\right)^{-1},
\qquad
r_0=1,
\qquad
r_{t+1}=1-e^{-\lambda r_t}.
\qquad\text{(4.2)}
```

Alors le broadcast binaire vérifie

```math
\boxed{
\ell_t(\lambda)\le q_t\le r_t(\lambda).
}
\qquad\text{(4.3)}
```

### Borne inférieure

Un estimateur linéaire de la racine est construit à partir du census du
bord. Le rapport

```math
\frac{
\mathbb E[\sigma_\rho Z_t]^2
}{
\mathbb E[Z_t^2]
}
\qquad\text{(4.4)}
```

est une borne inférieure à l'information optimale. Le calcul de second
moment sur le Galton--Watson donne $\ell_t(\lambda)$.

### Borne supérieure

La contraction $\chi^2$ d'une arête de corrélation $\theta$ vaut
$\theta^2$. L'information racine--bord est dominée par la probabilité qu'une
percolation auxiliaire, de nombre moyen d'enfants $d\theta^2=\lambda$,
atteigne la profondeur $t$. Pour un arbre $\mathrm{PGW}$, cette probabilité
satisfait précisément la récursion de $r_t$ dans (4.2).

Cette percolation d'information est un **outil de majoration**. Elle ne doit
pas être identifiée à la coupe physique de rétention
$q_{\beta_c^{\mathrm{geom}}}=1/d$.

## 5. Fermeture exacte sur le broadcast

Les deux côtés de (4.3) suffisent.

### Si $\lambda<1$

La récursion de $r_t$ est sous-critique et $r_t\to0$
exponentiellement. Donc $q_t\to0$.

### Si $\lambda=1$

On a

```math
\frac1{t+1}
\le
q_t
\le
r_t,
\qquad
r_t\sim\frac2t.
\qquad\text{(5.1)}
```

Ainsi $q_t\to0$ au point critique.

### Si $\lambda>1$

La borne inférieure donne

```math
\liminf_{t\to\infty}q_t
\ge
\frac{\lambda-1}{\lambda}
>
0.
\qquad\text{(5.2)}
```

Le broadcast est reconstructible.

On obtient donc

```math
\boxed{
\text{deux Gibbs exacts}
\ +\
\text{deux dendrogrammes indépendants marginalisés}
\ +\
\text{sandwich global}
\quad\Longrightarrow\quad
d\theta^2=1.
}
\qquad\text{(5.3)}
```

> [!NOTE]
> Le sandwich (4.3) n'est pas redémontré ici : la borne inférieure est le
> second moment exact du broadcast (Evans–Kenyon–Peres–Schulman) et la
> borne supérieure vient de la percolation d'information $\chi^2$
> (Abbe–Boix-Adserà) ; voir [REFERENCES](REFERENCES.md) pour les énoncés
> précis. Cette étape est donc **littérature-dépendante**, seule la
> calibration hiérarchique autour d'elle est propre au dossier.

## 6. Couper à $\beta_c$ et voir exactement le seuil

Plaçons sur la même horloge la géométrie et l'information. La coupe
géométrique d'une réplique et le temps de la percolation $\chi^2$ sont

```math
q_{\beta_c^{\mathrm{geom}}}=\frac1d,
\qquad
q_{\beta_\chi}=\theta^2.
\qquad\text{(6.1)}
```

Les deux dendrogrammes indépendants et leurs deux marginalisations donnent
d'abord $\theta\times\theta=\theta^2$. L'égalité de droite dans (6.1)
replace ensuite cette rétention sur l'horloge. Puisque $q_\beta$ est
strictement croissante,

```math
\boxed{
\beta_\chi=\beta_c^{\mathrm{geom}}
\quad\Longleftrightarrow\quad
\theta^2=\frac1d
\quad\Longleftrightarrow\quad
d\theta^2=1.
}
\qquad\text{(6.2)}
```

C'est le sens précis dans lequel le dendrogramme coupé à $\beta_c$ retrouve
**exactement** le seuil de weak recovery. Sous le seuil,
$\beta_\chi<\beta_c^{\mathrm{geom}}$ ; au-dessus,
$\beta_\chi>\beta_c^{\mathrm{geom}}$.

Cette lecture est un couplage de lois de percolation, pas un auxiliaire
révélé à l'estimateur. En particulier, l'intersection littérale des arêtes
ouvertes dans les deux dendrogrammes a une rétention $q_\beta^2$ : ce
n'est pas la coupe $q_{\beta_\chi}=\theta^2$. De même, couper ne signifie
jamais supprimer le corridor
$\beta_c^{\mathrm{geom}}<\beta\le1$ ; ses facteurs restent marginalisés
dans chacun des deux Gibbs complets.

## 7. Retour au SBM fini

Le calcul précédent est exact sur la limite locale de broadcast. Le
théorème du SBM fini ajoute deux ingrédients externes :

1. sous le seuil, une preuve de non-reconstruction/contiguïté montre
   qu'aucun algorithme n'a d'overlap non trivial ;
2. au-dessus du seuil, les méthodes spectrales non-backtracking ou le
   belief propagation construisent un estimateur corrélé.

Les résultats de Mossel--Neeman--Sly et de Massoulié identifient ainsi le
même seuil (1.2). Par le critère à deux répliques de la
[note 01](01_DU_CHAPITRE_11_AU_SBM.md#5-critère-informationnel-à-deux-répliques),
le théorème fini et le calcul de broadcast portent la même quantité
d'overlap.

### Ce que la hiérarchie explique exactement

- pourquoi la quantité pertinente est à deux répliques ;
- pourquoi une branche transmet $\theta^2$ et non $\theta$ ;
- pourquoi une coupe partagée déplace artificiellement le seuil ;
- comment couper à $\beta_c^{\mathrm{geom}}$ tout en gardant le Gibbs exact ;
- pourquoi la partie postcritique du dendrogramme ne peut pas être jetée.

### Ce qu'elle ne prouve pas encore seule

- une réduction du port global au broadcast : la route qui le déclarait
  négligeable est **réfutée** au seuil, où les deux grandes racines du
  full-$D$ ont des orientations asymptotiquement opposées ; le port doit
  être conservé exactement ou traité par une autre désintégration ;
- qu'un nombre explicite de sweeps du noyau hiérarchique contracte sous
  $\lambda=1$ ;
- que le calcul local demeure uniforme jusqu'aux profondeurs logarithmiques
  nécessaires sur le graphe avec cycles.

Le seuil est donc **retrouvé exactement comme calibration de la
représentation hiérarchique**, et non revendiqué comme une nouvelle preuve
autonome du théorème fini.

## 8. Deux faux raisonnements à écarter

### Figer tout un dendrogramme

Sur un arbre clairsemé, une fusion est souvent réalisée par une seule
arête. Conditionner par tout $D$ révèle alors sa parité. Le Gibbs à $D$ fixé
ne garde qu'un flip global par racine et voit le seuil $d\theta=1$, pas
$d\theta^2=1$.

### Appliquer seulement le Théorème 10 du chapitre 11

Appliqué littéralement au graphe potentiel complet du SBM, le graphe
Swendsen--Wang signé retient une arête avec probabilité
$(a-b)/n$. Sa percolation ne porte donc pas le seuil quadratique (1.2).
La preuve de weak recovery doit transporter la fiabilité du canal à travers
la géométrie ; la seule taille des composantes la perd.
