# Dossier ancestral : le problème central du programme

Ce dossier traite le problème qui gouverne toute la voie hiérarchique :

> pour une paire lointaine $`(i,j)`$ qui fusionne au nœud $`u`$ à un
> niveau $`\beta_u`$ proche du seuil de percolation $`\beta_c(p)`$,
> estimer les quatre taux $`\Lambda_v^{ab}`$ de **tous** les ancêtres
> $`v\succ u`$, donc le message ancestral $`B_u`$ et la persistance
> $`\eta_u=\tanh^2(L_u/2)`$ de la relation $`\sigma_i\sigma_j`$.

## Ordre de lecture

1. [42 — problème central](42_PROBLEME_CENTRAL_FUSION_CRITIQUE.md) :
   énoncé des deux volets (seuil de fusion, chaîne ancestrale), table de
   notations unique, identités structurelles, verrou G1–G3 et routes.
   **Commencer ici.**
2. [08 — chaîne exacte des taux](08_ANCESTRAL_LAMBDA_CHAIN.md) :
   décomposition des coupes ancestrales en trois groupes, quatre taux,
   réduction de Walsh $`(h_1,h_2,J)`$ et calcul de $`B_u`$.
3. [10 — estimation des taux](10_ANCESTRAL_LAMBDA_ESTIMATION.md) :
   loi exacte des marques sachant le squelette non marqué (course
   pondérée), moments, concentration, certificat de troncature
   $`R_u(I)`$, verrous G1–G3.
4. [14 — frontière critique](14_CRITICAL_COMPONENT_BOUNDARY.md) :
   annulation des arêtes internes, loi résiduelle des marques de
   frontière, biais de Palm de la paire lointaine, localisation gauche
   du LCA, certificat de majorité hiérarchique (théorème 7.1).

## Les trois résultats à retenir

- **Théorème de course conditionnelle** (10, th. 3.1 ; cas homogène
  08 §6) : sachant le squelette non marqué, les marques sont explicites —
  tout le biais restant est géométrique.
- **Certificat de troncature** $`|B_u-B_u^{(-I)}|\le R_u(I)`$ (10 §7) :
  seuls les ancêtres à déséquilibres non négligeables comptent.
- **Certificat de majorité hiérarchique** (14, th. 7.1) : des majorités
  conformes groupées à chaque ancêtre forcent la préférence de parité
  conforme — les ancêtres renforcent, ils n'annulent pas.

## Ce que ce dossier ne prouve pas

Ces notes ferment la partie « marques et fonctionnelles » du problème.
La loi du squelette groupé $`(m_{v,0},m_{v,1},m_{v,2},\beta_v)_{v\succ u}`$
sous la Palm de la paire critique (G1), la queue $`R_u`$ (G2), les coins
$`\ell_v\approx0`$ (G3) et les hypothèses CUT/ANC de la
[note 15](../../diagnostics/15_CRITICAL_GIANT_PAIR_FLIP.md) restent
ouvertes ; elles ne prouvent donc pas, à elles seules, une contraction du
corridor ni un seuil de weak recovery.
