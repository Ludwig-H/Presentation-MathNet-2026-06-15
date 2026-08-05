# Oracles locaux et factorisés

Ces notes calculent exactement la corrélation au LCA ou dans un chemin
factorisé. Elles restent utiles comme benchmarks, mais leurs hypothèses ne
décrivent pas le corridor multiport réel.

> **Attention : « archivé » ne veut pas dire « périmé en bloc ».**
> Plusieurs énoncés de ces notes sont porteurs pour le problème central et
> sont désormais référencés depuis la
> [note 42](../../foundations/ancestral/42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md) :
>
> - [06](06_LCA_SPIN_CORRELATION.md) : définitions du heat bath au LCA,
>   identité $`L_u=B_u+\ell_u`$, borne $`Q_n\le H_n^{\mathrm{LCA}}`$
>   (conditionnelle à A1) — socle du dossier ancestral ;
> - [16](16_FLIP_PROBABILITIES_DESCENDANT_PATH.md) : lemme 4.1
>   (paramétrisation de Walsh complète), lemme 4.2
>   ($`J=J_{\mathrm{ext}}+\ell_u/2`$), lemme 6.1 (annulation exacte des
>   descendants — la justification structurelle du fait que seuls les
>   ancêtres entrent dans $`B_u`$) ;
> - [17](17_PATH_DECORRELATION_THRESHOLD.md) : proposition 7.1
>   (critère de contraction jointe par normes d'opérateurs tordus), seul
>   énoncé de la note valable hors PATH-FAC.
>
> Ce qui est réellement mort : PATH-FAC comme modèle de la vraie
> dynamique (contre-exemple exact, 16 §7) et le seuil
> $`p_{\mathrm{path}}(\alpha)`$, qui n'existe que sous un ansatz
> géométrique non démontré.

La [réduction pairwise](../../foundations/03_HIERARCHICAL_WEAK_RECOVERY.md)
reste la référence globale ; le programme actif remplace les oracles par
des blocs collapsed et une géométrie de paire
([38](../../active/38_DOUBLE_GEANTE_GIBBS_REPLIQUE.md),
[41](../../active/41_DESINTEGRATION_PALM_RESTE_SIGNE.md)).
