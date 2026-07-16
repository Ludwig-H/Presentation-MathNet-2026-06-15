# Calculs et expériences à ajouter

Ce sous-dossier doit contenir uniquement des calculs reproductibles liés aux énoncés du dossier de recherche. Aucun résultat numérique ne sera présenté comme preuve.

## Structure prévue

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

   
$$
H(i,j)=\mathbf1_{\{i,j\text{ dans la même composante}\}}.
$$


5. Comparer

   
$$
Q_n,\qquad
   h_n(S),\qquad
   \lambda_{\max}(K_D^{\mathrm{info}})/n
$$


   par énumération exacte.
6. Retrouver les trois constantes $0.673648\ldots$, $0.719224\ldots$ et $0.794659\ldots$ à partir des formules symboliques.

## Métadonnées minimales

Chaque expérience devra enregistrer :

- graphe et conditions de bord ;
- $`n,p,\mu_0`$ et graines ;
- règle exacte de génération de $D$ ;
- programme de nœuds $S$ ;
- burn-in, nombre de chaînes et diagnostics de mélange ;
- métrique utilisée : corrélation $`R_n`$, overlap $`\operatorname{ov}_n`$, $`Q_n`$ ou spectre de $`H_S`$.

## Ordre de mise en œuvre

1. arbres de trois et quatre sommets ;
2. un triangle ;
3. deux triangles formant un cactus ;
4. bandes triangulaires ;
5. petits tores.

Cette progression doit détecter les erreurs de factorisation avant les simulations de grande taille.
