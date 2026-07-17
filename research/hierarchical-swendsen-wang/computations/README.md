# Calculs reproductibles

Ce sous-dossier doit contenir uniquement des calculs reproductibles liés aux énoncés du dossier de recherche. Aucun résultat numérique ne sera présenté comme preuve.

Les contrôles sans dépendance se lancent depuis la racine du dépôt avec :

```bash
python3 -m unittest discover \
  -s research/hierarchical-swendsen-wang/computations \
  -p 'test_*.py' -v
python3 research/hierarchical-swendsen-wang/computations/critical_band_thresholds.py
python3 research/hierarchical-swendsen-wang/computations/ancestral_lambda_chain.py
python3 research/hierarchical-swendsen-wang/computations/ancestral_lambda_estimation.py
python3 research/hierarchical-swendsen-wang/computations/critical_merger_oracle.py
python3 research/hierarchical-swendsen-wang/computations/triangle_block_sdpi.py
python3 research/hierarchical-swendsen-wang/computations/nishimori_hierarchical_entropy.py
python3 research/hierarchical-swendsen-wang/computations/ancestral_information_ledger.py
```

## Modules et extensions prévues

- `critical_band_thresholds.py` : vérification sans dépendance des trois seuils triangulaires et des temps $`\beta_c,t_\chi`$ ;
- `ancestral_lambda_chain.py` : calcul exact des quatre taux ancestraux, du message $`B_u`$ et de $`\mathbb E[\eta_u\mid\mathscr D]`$ sur un petit squelette homogène ;
- `test_ancestral_lambda_chain.py` : tests unitaires sans dépendance de la course conditionnelle et des identités quatre états ;
- `ancestral_lambda_estimation.py` : course pondérée exacte conditionnelle au squelette, moments des trois groupes et des quatre taux, certificat déterministe de queue ancestrale et transport vers la fiabilité LCA ;
- `test_ancestral_lambda_estimation.py` : contre-audit indépendant par énumération complète, réduction homogène fermée, bornes de contraste et constante de Lipschitz exacte $`2/(3\sqrt3)`$ de la fiabilité ;
- `critical_merger_oracle.py` : paramètres critiques fermés, somme finie $`\Gamma_m^c`$, borne exponentielle et limite $`m^{-1/2}`$ de l'oracle local ;
- `test_critical_merger_oracle.py` : contre-audits indépendants par LLR, expérience binaire symétrique, identité $`1/m`$ et limite gaussienne ;
- `triangle_block_sdpi.py` : profil $`c_q(t)`$ du canal de triangle, SDPI globale, audit scalaire et enveloppe multi-état conditionnelle ;
- `test_triangle_block_sdpi.py` : énumération directe du canal $`4\times8`$, information latérale, défaut de la droite naïve, matrice less-noisy $`3\times3`$ et conditions du candidat $`0.809909\ldots`$ ;
- `nishimori_hierarchical_entropy.py` : réduction exacte de l'équation (28) de Nishimori--Ohzeki à $`3h_2(p)-h_2((1+(2p-1)^3)/2)=1`$, racine unique, lois de bruit conditionnelles à quatre états et représentation par course exponentielle ;
- `test_nishimori_hierarchical_entropy.py` : contre-audit indépendant par énumération des huit mots de bruit, comparaison directe à l'équation publiée, monotonie de la balance et identité d'entropie du gagnant ;
- `ancestral_information_ledger.py` : bilan martingale générique d'une expérience binaire, bornes de logit, KL exacte d'une horloge exponentielle censurée, transport optimal $`2\times2`$ des quatre taux, information du dendrogramme d'une face et loi de Palm de la fusion critique ;
- `test_ancestral_information_ledger.py` : contre-audits du télescopage, des bornes de KL/transport, de l'énumération à huit états et des densités de Palm ;
- `exact_enumeration/` : énumération de toutes les configurations et horloges discrétisées sur petits graphes ;
- `symbolic/` : simplification des poids $`q_u^{ab}`$, contractions et seuils ;
- `triangular_grid/` : tores triangulaires, bandes et cactus ;
- `mixing/` : comparaison des programmes de nœuds ;
- `results/` : sorties brutes machine-readable, jamais uniquement des figures.

## Premiers tests unitaires

1. Vérifier que la marginale en $\sigma$ de $\nu(\sigma,D)$ est $\mu(\sigma)$.
2. Vérifier la balance détaillée des quatre états à chaque nœud.
3. Vérifier que 1000 heat baths successifs à $D$ fixé conservent $\nu(\cdot\mid D)$.
4. Vérifier que le cas racine donne
```math
    H(i,j)=\mathbf1_{\{i,j\text{ dans la même composante}\}}.
```
5. Comparer
```math
    Q_n,\qquad
    h_n(S),\qquad
    \lambda_{\max}(K_D^{\mathrm{info}})/n
```
    par énumération exacte.
6. Retrouver les trois constantes $0.673648\ldots$, $0.719224\ldots$ et $0.794659\ldots$ à partir des formules symboliques.
6 bis. Retrouver $`0.835805792367\ldots`$ comme unique racine supérieure de
$`H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3)=1`$ bit, puis vérifier algébriquement qu'il
s'agit exactement de l'équation (28), sans l'annoncer comme seuil démontré.
6 ter. Vérifier le bilan $L^2$ et entropique sur toute filtration binaire
finie, puis contre-auditer la formule

```math
D_{\mathrm{KL}}(P_{\lambda,\tau}\|P_{\mu,\tau})
=
(1-e^{-\lambda\tau})
\left(\log\frac\lambda\mu+\frac\mu\lambda-1\right)
```

par intégration indépendante. Sur le cactus de deux triangles, enregistrer
séparément corrélation réelle, fuite du sélecteur, incrément du LCA,
incrément du premier ancêtre et score terminal.
7. Pour chaque paire $i\ne j$, vérifier les identités LCA étendues aux racines distinctes
```math
    \mathbb E[m_{ij}^{\mathrm{LCA}}]
    =\mathbb E[(m_{ij}^{\mathrm{LCA}})^2]
    =\mathbb E[\eta_{ij}^{\mathrm{LCA}}]
    =A_{ij}^{(1)},
    \qquad
    c_{ij}^2\le A_{ij}^{(1)}.
```
8. Vérifier nœud par nœud
```math
    H_n^{\mathrm{LCA}}
    =\frac1{n^2}\mathbb E\left[
    n+2\sum_u|C_{u,1}||C_{u,2}|\eta_u
    \right]
    \le
    \frac1{n^2}\mathbb E\sum_R|R|^2.
```
9. Construire exactement la matrice du noyau $`K_{ij}^{\mathrm{LCA}}`$, contrôler que son spectre est dans $`[0,1]`$, puis comparer
```math
    A_{ij}^{(m)}
    =\langle f_{ij},K_{ij}^mf_{ij}\rangle
```
    à $`c_{ij}^2`$ pour $m=0,1,2,\ldots$.
10. Vérifier séparément la convention diagonale $`A_{ii}^{(m)}=1`$, puis les calibrations fermées : une arête, chemin, triangle isolé et cactus de triangles.
11. Au temps de fusion $t$, tester la loi conditionnelle
```math
    k=1+\mathrm{Bin}\left(m-1,
    \mathrm{logistic}(u_p(1-t))\right)
```

    conditionnellement au squelette non marqué, puis étudier séparément la loi géométrique du squelette choisi par Kruskal.
12. Pour chaque coupe $`\beta`$, vérifier exactement

```math
    Q_n
    \le
    S_n(\beta)+\mathcal M_n((\beta,1]),
```

    puis enregistrer séparément $`R_n^{>\beta}`$, $`\overline\eta_n^{>\beta}`$ et $`\kappa_n^{>\beta}`$.
13. Sur la grille triangulaire homogène, retrouver les trois seuils

```math
    p_{\mathrm{SW}}=0.673648\ldots,
    \qquad
    p_{\mathrm{info}}=0.794659\ldots,
    \qquad
    p_{\mathrm{pure}}=0.847296\ldots,
```

    et vérifier numériquement $`t_\chi(p_{\mathrm{info}})=\beta_c(p_{\mathrm{info}})`$.
14. Pour une coupe déterministe, simuler la décomposition

```math
    (R,S,U)
    \sim
    \mathrm{Mult}\left(
    m-1;
    h_p(t),
    \frac{1-h_p(t)}2,
    \frac{1-h_p(t)}2
    \right)
```

    et comparer les deux échelles $`m h_p(t)`$ et $`m h_p(t)^2`$.
15. Vérifier la calibration non oracle des bundles indépendants :

```math
    \gamma_m^{\mathrm{BSC}}
    =
    \sum_{k=0}^m
    \binom mk p^k(1-p)^{m-k}
    \tanh^2\left(\frac{u_p}{2}(2k-m)\right),
    \qquad
    \gamma_1^{\mathrm{BSC}}=(2p-1)^2.
```

16. Estimer la dérivée de $`\tau_{ij}(q_p(t))`$ par différences finies et la comparer au nombre moyen d'arêtes pivotales donné par la formule de Russo.
17. Sur des tores de diamètre $L$, tracer la mesure des temps de fusion dans la coordonnée proche-critique

```math
    \lambda
    =
    q_p'(\beta_c)L^{3/4}(t-\beta_c),
```

    en distinguant la masse géométrique, la masse pondérée par $`\eta_u`$ et le score signé.
18. Pour l'oracle local critique, vérifier indépendamment

```math
    s_c=\frac{p-q_c}{1-q_c},
    \qquad
    a_c=\log\frac{p-q_c}{1-p},
```

```math
    \Gamma_m^c(p_{\mathrm{SW}})=\frac1m,
```

    le rapport de vraisemblance $`P_+(K)/P_-(K)`$, la borne exponentielle et
    la limite

```math
    \Gamma_m^c\left(
    p_{\mathrm{SW}}+\frac{(1-q_c)\alpha}{2\sqrt m}
    \right)
    \longrightarrow
    \mathbb E[\tanh^2(\alpha Z+\alpha^2)].
```

19. Pour chaque bucket ancestral pondéré, vérifier indépendamment :

```math
\mathbb P(G_v=e\mid\mathscr S_u)
=\frac{w_es_{v,e}}{\sum_fw_fs_{v,f}},
```

les trois moyennes et covariances après marginalisation de $`G_v`$, leur
transport affine vers les quatre $`\Lambda_v^{ab}`$ et la borne

```math
|B_u-B_u^{(-I)}|\le\mathcal R_u(I).
```

Vérifier également, sur une grille déterministe et à proximité du point où la
dérivée est maximale,

```math
\left|
\tanh^2(L/2)-\tanh^2(\widetilde L/2)
\right|
\le
\frac{2}{3\sqrt3}|L-\widetilde L|
```

ainsi que $`\tanh^2(L/2)\le L^2/4`$.

20. Pour le canal de triangle, vérifier indépendamment

```math
c_q(1/4)=\eta_\triangle(q),
\qquad
c_q(1/2)=\frac{2q^2}{1+q^2},
```

le profil complet $`c_q(t)`$, l'échec des seules contraintes d'extrémité et
l'enveloppe affine. Contre-auditer $`P_\star`$ sur des a priori polarisés à
graine fixe en vérifiant tous les mineurs principaux de la matrice exacte,
sans confondre ce test fini avec une preuve. Toute sortie proche de
$`0.809909\ldots`$ doit conserver le statut « conditionnel au lemme
$`P_\star`$ ».

## Métadonnées minimales

Chaque expérience devra enregistrer :

- graphe et conditions de bord ;
- $`n,p,\mu_0`$ et graines ;
- règle exacte de génération de $D$ ;
- programme de nœuds $S$ ;
- burn-in, nombre de chaînes et diagnostics de mélange ;
- métrique utilisée : corrélation $`R_n`$, overlap $`\mathrm{ov}_n`$, $`Q_n`$ ou spectre de $`H_S`$.
- paire $(i,j)$, statut « LCA / racines distinctes / diagonale », puis, si le LCA existe, $`\beta_u`$, $`|E_u|`$, $`\Lambda_u`$, message $`B_u`$, $`m_u`$ et $`\eta_u`$.
- coupe critique utilisée, $`S_n(\beta)`$, masse de bande $`\mathcal M_n((\beta,1])`$ et facteurs $`R_n^{>\beta},\overline\eta_n^{>\beta},\kappa_n^{>\beta}`$ ;
- pour chaque coupe testée, multiplicité $m$, $`h_p(t)`$, nombre d'arêtes de bande et statut pivotal.

## Ordre de mise en œuvre

1. arbres de trois et quatre sommets ;
2. un triangle ;
3. deux triangles formant un cactus ;
4. chaînes de cactus de longueur croissante ;
5. bandes triangulaires ;
6. petits tores.

Cette progression doit détecter les erreurs de factorisation avant les simulations de grande taille.

## Protocole LCA sans burn-in annealed

Pour estimer la suite pair-spécifique avec $i\ne j$ :

1. générer $(\Sigma,O)$ et poser $\sigma^{(0)}=\Sigma$ ; par Nishimori, l'état initial est déjà stationnaire sous la moyenne annealed ;
2. à chaque itération, tirer de nouvelles horloges et construire le dendrogramme de Kruskal ;
3. si $`u_{ij}`$ existe, appliquer son heat bath exact avec les quatre $`q_u^{ab}`$ ; sinon, recolorer indépendamment les deux racines sous a priori uniforme ; puis oublier $D$ ;
4. estimer
```math
    A_{ij}^{(m)}
    =
    \mathbb E\left[
    \Sigma_i\Sigma_j\,
    \sigma_i^{(m)}\sigma_j^{(m)}
    \right].
```
La limite estime $`\mathbb E[c_{ij}(O)^2]`$. Pour $A^{(1)}$ agrégé sur toutes les paires, utiliser directement la somme sur les nœuds, sans énumérer les couples.

### Affecter tous les liens au bon nœud

Après construction de l'arbre de fusion, affecter chaque arête originale $e=\{x,y\}$ au LCA des feuilles $x,y$. Le bucket obtenu au nœud $u$ est exactement
```math
E_u=\{e:x\in C_{u,1},\ y\in C_{u,2}\},
```
y compris les arêtes de cycle. Cela garantit que $`\Lambda_u`$, $`T_u`$ et les flips utilisent tous les liens entre les deux fils, et pas seulement l'arête de la minimum spanning forest.
