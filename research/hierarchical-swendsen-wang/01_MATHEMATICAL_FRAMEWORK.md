# Cadre mathématique de la dynamique hiérarchique

## 1. Postérieure et notations

On conserve les notations du chapitre 11. Pour une observation $O=(X,W)$, la postérieure est


$$
\mu_O(\sigma)
\propto
\mu_0(\sigma)
\exp\!\left[-\sum_{e\in E}|W_e|\,
\mathbf 1_{\{e\text{ non satisfaite par }\sigma\}}\right].
$$


Une arête $e=\{i,j\}$ est satisfaite si $`\sigma_i=\sigma_j`$ lorsque $`W_e>0`$, et si $`\sigma_i\ne\sigma_j`$ lorsque $`W_e<0`$.

Les résultats de weak recovery les plus nets de ce dossier sont d'abord formulés pour $K=2$, $`\sigma_i\in\{-1,+1\}`$, avec a priori i.i.d. uniforme. L'a priori général $`\mu_0`$ reste autorisé pour la dynamique, mais il doit apparaître dans les probabilités de recoloriage : on ne peut pas recolorer uniformément sans vérifier la balance vis-à-vis de $`\mu_0`$.

## 2. Horloges et filtration

Conditionnellement à $\sigma$, on pose indépendamment


$$
\xi_e\sim\operatorname{Exp}(|W_e|)
\quad\text{si }e\text{ est satisfaite},
\qquad
\xi_e=+\infty
\quad\text{sinon}.
$$


Une arête de poids nul est ignorée, ou de façon équivalente reçoit $`\xi_e=+\infty`$.

Pour $0\le t\le1$, soit $`\Pi_t`$ la partition en composantes connexes du graphe formé par les arêtes telles que $`\xi_e\le t`$. Alors


$$
\mathbb P(\xi_e\le1\mid e\text{ satisfaite})=1-e^{-|W_e|},
$$


donc $`\Pi_1`$ est exactement la partition gelée de Swendsen–Wang par arêtes.

Le dendrogramme $`D=(\Pi_t)_{0\le t\le1}`$ peut être calculé par Kruskal sur les poids $`\xi_e`$. La minimum spanning forest ne sert qu'à calculer la filtration : lorsqu'un nœud $u$ fusionne $`C_1`$ et $`C_2`$, les probabilités de flip utilisent **l'ensemble**


$$
E_u=\{\{i,j\}\in E:i\in C_1,\ j\in C_2\},
$$


et pas seulement l'arête choisie par Kruskal.

On note


$$
\Lambda_u(\sigma)
=
\sum_{e\in E_u}|W_e|\,
\mathbf 1_{\{e\text{ satisfaite par }\sigma\}},
\qquad
\beta_u
=
\min_{\substack{e\in E_u\\e\text{ satisfaite}}}\xi_e.
$$


Pour une coupe candidate fixée et $`\Lambda_u(\sigma)>0`$, le minimum brut est exponentiel de taux $`\Lambda_u(\sigma)`$. Sur l'intervalle $`[0,1]`$,


$$
\mathbb P(\beta_u\in d\beta\mid\sigma)=
\Lambda_u(\sigma)e^{-\beta\Lambda_u(\sigma)}\,d\beta,
\qquad 0\le\beta\le1.
$$


Il s'agit ici d'une sous-densité : la masse


$$
\mathbb P(\beta_u>1\mid\sigma)=e^{-\Lambda_u(\sigma)}
$$


correspond à l'absence de fusion avant la coupe. Conditionnellement à une fusion avant $1$, la densité est


$$
f(\beta_u\mid\beta_u\le1,\sigma)
=
\frac{\Lambda_u(\sigma)e^{-\beta_u\Lambda_u(\sigma)}}
{1-e^{-\Lambda_u(\sigma)}},
\qquad 0\le\beta_u\le1.
$$


Dans la loi globale du dendrogramme, on emploie les sous-densités, car l'existence ou non de chaque fusion fait partie de $D$.

## 3. Loi jointe de type Edwards–Sokal

La loi auxiliaire exacte peut s'écrire


$$
\nu_O(\sigma,D)
\propto
\mu_O(\sigma)\,\mathbb P(dD\mid\sigma).
$$


Avant simplification, la loi des horloges ordonnées prend la forme


$$
\mathbb P(dD\mid\sigma)
=
e^{-\Lambda_\infty(\sigma)}
\prod_{u\in D}
\Lambda_u(\sigma)e^{-\beta_u\Lambda_u(\sigma)}\,dD,
$$


où $`e^{-\Lambda_\infty(\sigma)}`$ regroupe les facteurs de survie des liens censurés à la coupe $1$. Ce terme est indispensable dans la dérivation, même s'il disparaît ensuite avec les facteurs correspondants de l'énergie.

Après annulation des facteurs correspondant à l'énergie, la conditionnelle utile est


$$
\boxed{
\nu_O(\sigma\mid D)
\propto
\mu_0(\sigma)
\prod_{u\in D}
\Lambda_u(\sigma)
e^{(1-\beta_u)\Lambda_u(\sigma)}.
}
$$


Le produit porte sur les fusions observées avant la coupe $1$, avec la convention précise adoptée dans la construction du dendrogramme. Cette formule, plutôt qu'un argument informel de recoloration, est le point de départ sûr pour prouver la stationnarité avec un a priori non uniforme.

## 4. Heat bath à un nœud

Soit $`u:C=C_1\mathbin{\dot\cup}C_2`$. Pour $a,b\in\{0,1\}$, $\sigma^{ab}$ désigne la configuration obtenue en flippant $`C_1`$ si $a=1$, et $`C_2`$ si $b=1$. Les quatre poids exacts sont


$$
q_u^{ab}
=
\mu_0(\sigma^{ab})
\prod_{v\succeq u}
\Lambda_v(\sigma^{ab})
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})},
$$


où $v\succeq u$ parcourt $u$ et ses ancêtres affectés par le flip. Le heat bath choisit


$$
\mathbb P((a,b)\mid\sigma,D,u)
=
\frac{q_u^{ab}}{\sum_{c,d\in\{0,1\}}q_u^{cd}}.
$$


La variable pertinente pour la transmission entre les deux fils est leur orientation relative. Ses log-odds exactes sont


$$
L_u
=
\log\frac{q_u^{00}+q_u^{11}}{q_u^{10}+q_u^{01}},
\qquad
\rho_u=\left|\tanh\frac{L_u}{2}\right|.
$$


$`\rho_u`$ est une première mesure locale de fiabilité. Elle n'est pas encore une capacité globale : les facteurs des ancêtres rendent les variables de nœuds dépendantes.

### Formule locale simplifiée

Posons


$$
T_u=\sum_{e\in E_u}|W_e|.
$$


Si l'a priori et les facteurs ancêtres sont neutres pour la comparaison des deux parités, alors $`T_u-\Lambda_u`$ devient le poids satisfait après inversion relative et


$$
L_u^{\mathrm{loc}}
=
\log\frac{\Lambda_u}{T_u-\Lambda_u}
+(1-\beta_u)(2\Lambda_u-T_u).
$$


Cette formule montre explicitement où intervient l'horloge $`\beta_u`$ de la première arête de fusion. Elle ne doit pas être employée lorsque les facteurs ancêtres ne s'annulent pas.

## 5. Stationnarité et programmes de flips

À $D$ fixé, chaque heat bath ci-dessus satisfait la balance détaillée pour $`\nu_O(\cdot\mid D)`$. Par composition, on peut donc choisir successivement un nombre arbitraire de nœuds — par exemple 1000 nœuds aléatoires, avec répétition — et appliquer les heat baths un à un.

Le choix des indices peut être déterministe ou aléatoire, dépendre de $D$ et de l'historique des indices déjà choisis, mais il doit être indépendant de la configuration courante $\sigma$. Une sélection adaptative dépendant de $\sigma$ ne conserve pas automatiquement la mesure ; elle exige une preuve de balance propre ou une correction de type Metropolis–Hastings.

Un pas complet est :

1. tirer $`D\sim\nu_O(\cdot\mid\sigma)`$ à l'aide des horloges ;
2. appliquer un programme $S$ de heat baths laissant $`\nu_O(\cdot\mid D)`$ invariante ;
3. oublier $D$.

Le noyau marginal ainsi obtenu laisse $`\mu_O`$ invariante.

Attention : faire mélanger uniquement la chaîne en $\sigma$ à $D$ fixé produit un tirage de $`\nu_O(\cdot\mid D)`$, pas une réplique indépendante de $`\mu_O`$. Pour obtenir une réplique postérieure indépendante, il faut faire mélanger la chaîne alternée qui rafraîchit aussi $D\mid\sigma$.

## 6. Deux extrémités de la hiérarchie

- **Orientations globales des arbres.** Sous a priori uniforme, le heat bath de l'orientation globale de chaque composante de $`\Pi_1`$ est uniforme ; on retrouve la recoloration Swendsen–Wang. Cette mise à jour à deux états est distincte du heat bath à quatre états d'un nœud interne supérieur.
- **Feuilles.** Le heat bath d'une feuille est une mise à jour mono-site de Glauber. Un noyau de Metropolis–Hastings mono-site satisfaisant la balance pour la même conditionnelle peut être substitué. Il vaut donc mieux écrire que la famille interpole entre Swendsen–Wang et les dynamiques locales heat bath/Metropolis, plutôt que d'identifier automatiquement une feuille à toute règle Metropolis–Hastings.

## 7. Points de vigilance mathématique

1. $D$ est une variable auxiliaire corrélée à une réplique postérieure ; elle n'est pas observée dans le problème statistique original.
2. Révéler $D$ à un estimateur est licite pour une preuve d'impossibilité par oracle. Une preuve de suffisance doit produire un estimateur à partir de $O=(X,W)$, ou justifier l'échantillonnage de $D$.
3. Les fusions sont sélectionnées par Kruskal. La loi du nombre d'arêtes satisfaites traversant une coupe sélectionnée n'est pas une loi binomiale naïve.
4. Les probabilités aux nœuds ne factorisent pas en général en canaux indépendants, à cause des ancêtres et des cycles du graphe initial.
5. Pour un a priori général, toute recoloration $`R_D`$ doit vérifier une balance explicite vis-à-vis de $`\mu_0`$. La mesure jointe ci-dessus est la voie recommandée.
