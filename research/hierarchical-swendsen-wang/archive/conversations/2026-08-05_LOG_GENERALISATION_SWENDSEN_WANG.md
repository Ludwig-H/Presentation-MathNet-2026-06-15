**You:**

Prends connaissance du dossier dans mon repo GitHub Presentation-MathNet-2026-06-15/research/hierarchical-swendsen-wang/
Mon but est d'utiliser généralisation hiérarchique de Swendsen-Wang à l'aide d'horloges exponentielles pour montrer des résultats (de préférence exacts) sur la weak recovery.
Etant donné deux points i et j asymptotiquement infiniment éloignés, la meilleure façon d'utiliser cette dynamique est de regarder le seuil bêta où i et j fusionnent dans la même composante (l'on s'intéresse au seuil critique de weak recovery ; on peut donc se restreindre aux événements où cette fusion a lieu très proche du seuil critique).
Le tout est de bien savoir estimer les Lambda_v pour tout v au-dessus du noeud de fusion u dans la hiérarchie. Définis clairement le problème, cherche à le formaliser et des solutions. Mets le repo à jour pour qu'il soit propre et bien entretenu

* * *

**ChatGPT:**

Le dépôt est mis à jour sur la branche `agent/ancestral-lambda-chain` avec une PR brouillon #1.

Résultats principaux :

* Pour chaque ancêtre $v\succ u$, les quatre $\Lambda_v(\sigma^{ab})$ se déduisent exactement de trois groupes d’arêtes.
* Toute la chaîne ancestrale se réduit à quatre log-poids, puis aux trois scalaires $(h_1,h_2,J)$, avec $$B_u=2J+\log\cosh(h_1+h_2)-\log\cosh(h_1-h_2).$$
* Dans le GSBM homogène, conditionnellement au squelette de Kruskal non marqué, les buckets sont indépendants et $$K_{v,r}  
    =\mathbf1_{\{G_v=r\}}  
    +\mathrm{Bin}\!\left(m_{v,r}-\mathbf1_{\{G_v=r\}},s_p(\beta_v)\right).$$ Le biais restant porte donc sur la géométrie $(m_{v,r},\beta_v)$, pas sur les marques conditionnelles.
* Point important : le temps de fusion de deux sommets lointains n’est pas automatiquement concentré près de $\beta_c$. Une restriction à une bande infinitésimale demanderait un lemme de localisation supplémentaire. Il faut actuellement conserver toute la chaîne au-dessus du LCA.

Les livrables centraux sont , et .

Validation : 8 tests réussis, 15 fichiers Markdown vérifiés, seuils numériques reproduits et tree GitHub identique au tree local validé.

Je n’annonce pas encore de nouveau seuil planaire : la prochaine cible exacte réaliste est le cactus de triangles, puis les bandes de largeur fixée par matrices de transfert certifiées.

* * *

**You:**

Le cas le plus favorable pour la weak recovery est lorsque i et j lointains fusionnent au niveau du seuil beta associé à la percolation ; on peut donc se placer dans ce cas favorable pour mesurer le seuil de weak recovery. Pousse l'analyse effectivement plus loin dans cette direction. Et ne crée pas de branche supplémentaires : fusionne tout sur le main. Garde un repo propre et bien entretenu ! Sois extrêmement rigoureux mathématiquement et fais toujours des contre-audits pour les résultats difficiles.

* * *

**ChatGPT:**

Tout est fusionné sur `main` via la PR #1, commit `16682f0`.

Conclusion mathématique principale :

* Au temps critique, $$s_c=\frac{p-q_c}{1-q_c},\quad  
    h_c=\frac{2(p-p_{\mathrm{SW}})}{1-q_c},\quad  
    a_c=2\,\mathrm{artanh}(h_c).$$
* Pour l’oracle local $B_u=0$, $$\Gamma_m^c(p_{\mathrm{SW}})=\frac1m,  
    \qquad  
    \Gamma_m^c(p)\to1\quad(p>p_{\mathrm{SW}}).$$
* La transition finie vérifie $$p-p_{\mathrm{SW}}\asymp m^{-1/2},\qquad  
    \Gamma_m^c\to\mathbb E[\tanh^2(\alpha Z+\alpha^2)].$$
* Le contre-audit décisif est $$\mathcal C_{n,\delta}^c  
    \le S_n(\beta_c+\delta)-S_n(\beta_c).$$ Une fiabilité conditionnelle parfaite sur une classe de paires de masse $o(1)$ ne produit donc aucun overlap macroscopique.
* Ainsi, la frontière de cet oracle favorable est seulement $p_{\mathrm{SW}}=0.673648\ldots$, pas le seuil réel de weak recovery, déjà contraint par la borne d’information-percolation de référence (Abbe–Boix).
* Pour les ancêtres, les comptes sont encadrés entre les lois de paramètres $1/2$ et $s_c$, mais aucun ordre global sur $\eta_u$ n’est valide : les flips complémentent les comptes et $B_u$ peut annuler le LLR local.

Livrables principaux :

Validation : 19 tests, validateur mathématique, trois scripts reproductibles, audit aléatoire de 250 couples $(m,p)$, et deux remarques de revue GitHub corrigées. Le checkout local est propre, uniquement sur `main`, synchronisé avec `origin/main`. Aucune nouvelle branche n’a été créée ; seul le ref distant historique de la PR reste visible, le connecteur ne permettant pas sa suppression.

* * *

**You:**

Tu n'as pas compris du tout. Relis les slides consacrées à la dynamique hiérarchique beamer-presentation-reunion-2026-07-16/Presentation_2026-07-16_LouisHauseux_ReunionLouisNahuelKonstantin.pdf
La difficulté est d'estimer Les Lambda_v pour v au-dessus de u le LCA de i et j

* * *

**ChatGPT:**

Tu avais raison : le cœur du problème est bien la chaîne complète $v\succ u$, pas l’oracle local au nœud critique.

J’ai réorienté l’analyse autour du résultat exact suivant. Conditionnellement au squelette non marqué,

$$s_{v,e}  
=  
\frac{\pi_e e^{-w_e\beta_v}}  
{1-\pi_e+\pi_e e^{-w_e\beta_v}},  
\qquad  
\mathbb P(G_v=e\mid\mathscr S_u)  
=  
\frac{w_es_{v,e}}{\sum_fw_fs_{v,f}}.$$

Sachant la gagnante $G_v$, les autres marques sont des Bernoulli indépendantes. Cela donne exactement la loi jointe, la transformée génératrice, les moyennes et covariances de tous les $\Lambda_v^{ab}$.

Dans le cas homogène, la moyenne se ferme notamment en

$$\mathbb E[\Lambda_v^{ab}\mid\mathscr S_u]  
=  
\frac{u_p}{2}\left(m_v+\alpha_v A_v^{ab}\right),$$

avec

$$\alpha_v=h_v+\frac{1-h_v}{m_v},  
\qquad  
A_v^{ab}  
=  
m_{v,0}+(1-2a)m_{v,1}+(1-2b)m_{v,2}.$$

J’ai également établi un certificat déterministe de troncature :

$$|B_u-B_u^{(-I)}|\le \mathcal R_u(I),$$

où $\mathcal R_u(I)$ dépend explicitement des deux déséquilibres retournables, du plus petit des quatre taux et des temps $\beta_v$.

Le verrou restant est maintenant isolé proprement : sous la loi biaisée par une paire lointaine dont le LCA tombe dans la fenêtre critique, il faut établir :

* la convergence des premiers $(m_{v,0},m_{v,1},m_{v,2},\beta_v)$ ;
* la sommabilité de la queue $\mathcal R_u$ ;
* le contrôle des quatre taux proches de zéro.

La distinction entre loi du squelette et loi étendue « squelette × marques conditionnelles » est explicitée, tout comme le fait que l’expérience critique est vide avant $t=1$ lorsque $p<p_{\mathrm{SW}}$.

Fichiers principaux :

Contre-audits réussis : 32 tests, énumérations exhaustives pondérées, vérifications aléatoires des covariances, concentrations et bornes de queue, plus l’audit MathJax global.

Tout est directement sur `main`, sans nouvelle branche. `main` local et GitHub sont propres et synchronisés au commit `ac1d02c`.

* * *

**You:**

Continue l'exploration de cette voie. Peut-on obtenir une meilleure borne pour la weak recovery (par exemple sur le GSBM triangulaire)par cette voie ?

* * *

**ChatGPT:**

Oui, cette voie fournit un meilleur candidat, mais pas encore une borne inconditionnelle.

* L’approche scalaire par triangles est définitivement éliminée : sa SDPI globale vaut $$\frac{2q^2}{1+q^2},$$ et donne seulement $p<0{,}788675$, moins bien que la borne par arêtes $0{,}794659$. La valeur naïve $0{,}829491$ oubliait l’information latérale. C’est précisément le piège couvert par l’information-percolation multi-terminale de Polyanskiy–Wu et la baseline d’Abbe–Boix.
* En conservant les cinq états du triangle — vide, trois relations simples, plein — on obtient le candidat $$p_\star^{\mathrm{cond}}=0{,}809909289251919\ldots$$ où $q_\star=2p_\star-1$ est la racine dans $(0{,}61,0{,}63)$ de $$q^{10}+46q^8+45q^6-20q^4-12q^2+4=0.$$ Les paramètres vérifient $a_\star=e_\star$, donc se trouvent sur la surface critique de Chayes–Lei.
* Le verrou est maintenant très précis : prouver la positivité d’une matrice rationnelle $3\times3$, notée $P_\star$, pour les a priori ayant un atome $>1/2$. Ce lemme est exactement le critère less-noisy de Makur–Polyanskiy. Il est déjà démontré lorsque $\max_x\mu_x\le1/2$. Un contre-audit PSD reproductible n’a trouvé aucun contre-exemple dans le secteur polarisé, mais ce n’est pas une preuve.

Pour les $\Lambda_v$, la bonne quantité est donc la polarisation des quatre poids LCA :

$$M_u=\max_{a,b}\pi_u^{ab}.$$

Les estimations ancestrales peuvent certifier le secteur démontré via

$$\max\ell_u^{ab}-\min\ell_u^{ab}\le\log 3  
\quad\Longrightarrow\quad M_u\le\frac12.$$

Il reste toutefois à relier rigoureusement ces poids LCA à l’a priori du canal facteur, ou à fermer directement $P_\star$.

Tout est publié directement sur `main` dans le commit ea44785. Voir la et le . Les 43 tests et la validation des 18 fichiers mathématiques passent ; `main` est propre et synchronisé.

* * *

**You:**

non, mais je veux que tu utilises la dynamique hiérarchique !! je t'ai dit de montrer que les lemmes utiles our la weak recovery étaient de supposer le cas le plus favorable : i et j lointains sont dans un même arbre et se séparent au niveau du seuil de percolation

* * *

**ChatGPT:**

La voie est désormais correctement recentrée sur la dynamique hiérarchique.

Le résultat central est :

$$Q_L\le b_L+S_L(a_\varepsilon)  
+\Gamma_{L,\varepsilon}^{\mathrm{fav}}  
+\frac{2}{3\sqrt3}\zeta_{L,\varepsilon}  
+\delta_{L,\varepsilon},$$

où l’expérience favorable impose que $i,j$ soient lointains, dans le même arbre, et se séparent à la descente près de $\beta_c$. Le lemme HF formalise précisément que cet oracle critique domine les fusions postcritiques.

Surtout, sa fiabilité utilise bien

$$L_u=\ell_u^{\mathrm{crit}}+B_u,$$

avec $B_u$ calculé à partir des quatre $\Lambda_v^{ab}$ de chaque ancêtre $v\succ u$. La troncature après $K$ ancêtres vérifie exactement

$$|\eta_u-\eta_u^{(K)}|  
\le  
\min\!\left(1,\frac{2\mathcal R_u^{(>K)}}{3\sqrt3}\right).$$

Le contre-audit révèle aussi la difficulté réelle : à la borne information-percolation, l’oracle local $B_u=0$ a déjà une fiabilité $0{,}9978$ pour $m=64$. Une meilleure borne exige donc que le message ancestral compense effectivement le LLR local critique :

$$B_u\simeq-\ell_u^{\mathrm{crit}}.$$

Ce n’est pas encore démontré ; aucune nouvelle constante triangulaire n’est donc annoncée. Le canal de triangle du fichier 11 est maintenant explicitement relégué au rang de piste auxiliaire.

La nouvelle analyse est dans . Les README, roadmap et notes 04/07/09/10/11 ont été harmonisés.

Vérifications : 44 tests réussis, 19 fichiers Markdown validés, dépôt propre. Tout est publié directement sur `main`, sans nouvelle branche : commit 41b24d1.

* * *

**You:**

Regarde aussi s'il n'y aurait pas un moyen de réobtenir le seuil de la conjecture de Nishimori avec cette dynamique hiérarchique.

* * *

**ChatGPT:**

Oui — la dynamique hiérarchique permet de réobtenir exactement la constante conjecturée, mais pas encore de prouver qu’elle est le seuil de weak recovery.

Pour

$$r_p=\frac{1+(2p-1)^3}{2},$$

l’équation (28) de Nishimori–Ohzeki est exactement équivalente à

$$H(Z_1,Z_2,Z_3\mid Z_1Z_2Z_3)  
=3h_2(p)-h_2(r_p)=1\ \text{bit}.$$

Elle possède une unique racine dans $(1/2,1)$, encadrée par

$$0.835805792366<p_{\mathrm N}^{(0)}<0.835805792368.$$

Sur un triangle, les spins modulo flip global forment quatre états, en bijection avec les quatre mots de bruit compatibles avec le syndrome observé. Le heat bath hiérarchique collapsed des deux bits internes du dendrogramme, après marginalisation de $D$, est donc exactement une course de quatre horloges exponentielles. L’entropie ci-dessus est l’espérance de la surprise du gagnant.

Les contre-audits importants sont maintenant explicites :

* le heat bath d’un seul LCA n’est pas cette course de face : ses quatre flips contiennent deux copies liées par le flip global ;
* au LCA critique, les poids exacts restent $$q_u^{ab}=\mu_0(\sigma^{ab})  
    \prod_{v\succeq u}\Lambda_v^{ab}  
    e^{(1-\beta_v)\Lambda_v^{ab}},$$ donc tous les $\Lambda_v^{ab}$ ancestraux sont indispensables ;
* conditionner par le dendrogramme coûte $$H(\Sigma/\{\pm1\}\mid O)  
    =  
    H(\Sigma/\{\pm1\}\mid O,D)+I(\Sigma;D\mid O);$$
* une égalité d’entropie totale ne contrôle pas la weak recovery à grande distance.

La prochaine étape précise est de construire, sur cactus, les défauts autoduals $\widehat\Psi_0,\widehat\Psi_1$ sous la vraie loi biaisée du LCA critique, puis de les comparer à $\Gamma_1^{\mathrm{fav}}$ avec tous les ancêtres. C’est le pont NH1–NH3. La méthode améliorée par amas déplace déjà l’approximation vers $0.835985$, donc la hiérarchie complète pourrait corriger la conjecture initiale plutôt que la reproduire exactement. Nishimori–Ohzeki 2006, Ohzeki 2009.

Le dossier est documenté dans , avec le et son .

Tout est publié directement sur `main` jusqu’au commit `adb3af3`. Le dépôt est propre et synchronisé ; 50 tests et le validateur mathématique passent.

* * *

**You:**

Continue. Prends du recul et garde un repo propre et bien entretenu et pédagogique et bien présenté. N'hésite pas à faire du ménage.

* * *

**ChatGPT:**

Le dépôt est propre, aligné sur `main` et publié sans branche supplémentaire : commit 88e72c7.

Résultat mathématique principal :

$$\mathbb E[\eta_u]  
=  
\mathbb E[c_{ij}(O)^2]  
+\Delta_u^{\rm env}  
+\Delta_{u,0}^{(2)}  
+\sum_{k\ge1}\Delta_{u,k}^{(2)}.$$

C’est le bilan exact entre corrélation réelle, information du sélecteur LCA, fusion critique et ancêtres. Il montre notamment que les ancêtres ne peuvent pas diminuer en moyenne un score correctement marginalisé : l’annulation ponctuelle par $B_u$ n’est pas une stratégie de preuve.

Pour chaque ancêtre à coupe fixée, j’ai obtenu la brique quantitative

$$D_{\rm KL}(P_{\lambda,\tau}\Vert P_{\mu,\tau})  
=  
(1-e^{-\lambda\tau})  
\left(\log\frac{\lambda}{\mu}+\frac{\mu}{\lambda}-1\right),$$

puis une borne de transport $2\times2$ utilisant réellement les quatre taux  
$(\Lambda^{00},\Lambda^{11},\Lambda^{01},\Lambda^{10})$. La sélection géométrique du squelette reste explicitement séparée.

Le contre-audit triangulaire montre qu’au temps critique :

* à $p_{\rm info}=0.794659\ldots$, le dendrogramme grossier ajoute $0.0617163$ bit et un gain $L^2=0.0318349$ sur une paire ;
* à $p_{\rm N}^{(0)}=0.835805\ldots$, il ajoute encore $0.0427594$ bit et $L^2=0.0195231$ ;
* sous l’oracle de Palm critique, l’entropie de face résiduelle vaut $0.2515601$ bit, et non $1$ bit.

Conclusion prudente : aucune meilleure borne numérique n’est encore démontrée. Le raccourci identifiant l’équation de Nishimori à l’oracle critique est réfuté. Une amélioration au-delà de $0.794659\ldots$ reste possible seulement si le calcul sur cactus croissants montre que le score collapsed initial et les coûts géométriques/exponentiels sont sommables, puis si HF est établi.

Les éléments centraux sont dans et son . Les 58 tests passent et les 21 fichiers Markdown passent le validateur mathématique. Le prochain calcul décisif est maintenant clairement le bilan complet du cactus à deux triangles.

* * *

**You:**

Non mais tu n's rien compris à mon idée. Pour chaque p supérieur à 1/2 du GSBM, il y a un beta_c critique de percolation. L'idée est de dire que le cas le plus favorable pour l weak recovery est lorsque i et j loins fusionnent au niveau critique (pour récupérer le plus de liens Lambda et de meilleure qualité).
Le tout est de savoir la qualité de ces liens, c'est-à-dire si une majorité stricte est conforme à la ground truth, ce qui incitera la dynamique hiérarchique au niveau du noeud de fusion u à garder i et j avec un signe conforme. Or ces liens ramassés après beta_c sont l'ensemble des liens faux (non conformes à la ground truth) plus ceux dont l'horloge est entre beta _c et 1. formalise cette idée. Audite-la, puis contre-audite-la. Et nettoie le repo en insistant sur ce genre d'idées en rapport avec la dynamique hiérarchique.

* * *

**ChatGPT:**

C’est maintenant formalisé dans le sens exact de ton idée, directement sur `main`.

### Résultat mathématique

Pour un lien encore fermé au temps $t$,

$$(\text{vrai tardif},\text{vrai censuré},\text{faux})  
=  
\left(h_p(t),\frac{1-h_p(t)}2,\frac{1-h_p(t)}2\right),$$

avec

$$h_p(t)=\tanh\!\left(\frac{u_p(1-t)}2\right),  
\qquad  
s_p(t)=\frac{1+h_p(t)}2.$$

Ainsi, les vrais liens censurés après $1$ et les faux liens se compensent exactement. Tout l’excès de majorité conforme vient bien des vraies horloges dans $(t,1]$. Ce biais est maximal à $t=\beta_c$.

Au seuil triangulaire $q_c=2\sin(\pi/18)$, dont l’exactitude vient de la transformation étoile-triangle, les masses sont exactement

$$(q_c,\;2p-1-q_c,\;1-p,\;1-p)$$

pour, respectivement, les vrais précoces, vrais tardifs, vrais censurés et faux.

Deux seuils distincts apparaissent :

$$p_{\mathrm{SW}}=\frac{1+q_c}{2}=0.673648\ldots,  
\qquad  
p_{\mathrm{late}}=\frac{2+q_c}{3}=0.782432\ldots.$$

* Le paquet complet vu par $\Lambda$ a une majorité conforme stricte pour $p>p_{\mathrm{SW}}$.
* Les vrais tardifs battent à eux seuls tous les faux pour $p>p_{\mathrm{late}}$.

Le résultat nouveau le plus utile est le certificat hiérarchique exact :

> Sous a priori uniforme, si le bucket de fusion $u$ a une majorité conforme stricte et si, pour chaque ancêtre $v\succ u$, les deux groupes touchant respectivement $C_1$ et $C_2$ ont une majorité conforme non négative, alors
> 
> $$q_u^{00}+q_u^{11}>q_u^{10}+q_u^{01}.$$
> 
> La dynamique hiérarchique préfère donc strictement garder la parité conforme.

La preuve utilise la croissance et la convexité de $F_v(x)=xe^{(1-\beta_v)x}$, puis la stabilité par produit d’un cône de coefficients de Walsh positifs.

Conditionnellement au squelette et aux gagnantes, la probabilité de ce certificat est donnée exactement par un produit de queues binomiales. Sa concentration est gouvernée par

$$m_{v,r}h_p(\beta_v)^2,$$

avec un coût $\log H$ pour contrôler simultanément $H$ ancêtres.

### Contre-audit weak recovery et Nishimori

La majorité scalaire seule ne peut pas améliorer la borne actuelle :

$$p_{\mathrm{late}}=0.782432\ldots  
<  
p_{\mathrm{info}}=0.794659\ldots.$$

Elle devient donc favorable avant la zone encore ouverte. Une amélioration éventuelle doit venir des amplitudes des quatre $\Lambda_v^{ab}$, des tailles des groupes ancestraux, de leur dégradation lorsque $\beta_v>\beta_c$, et de la loi du squelette biaisée par la paire critique.

Au point conjectural de Nishimori,

$$h_c=0.4968797586\ldots,  
\qquad  
s_c=0.7484398793\ldots.$$

La majorité y est déjà fortement établie : elle ne peut pas, à elle seule, réobtenir $0.8358057923\ldots$. Il faut une identité sur les quatre poids collapsed, puis HF.

### Dépôt

* Note centrale :
* Calcul exact :
* Contre-audit :

L’ancien « ancestral information ledger » et ses deux scripts, hors de l’axe demandé, ont été supprimés ; ils restent récupérables dans l’historique Git. Les documents 04, 05, 08–10, 12, 13 et le README ont été réorientés.

Validation : 63 tests réussis, audit MathJax réussi, arbre de travail propre. Commit publié directement sur `main`, sans branche supplémentaire : `64f6e5d`.

* * *

**You:**

En fait, il ne faut pas prendre en compte les liens faux internes à des clusters. De plus, l'hyothèse i et j appartiennent à la composante géante au moment de la percolation a certainement une influence. Bref, pars de mon idée. Formalise la correctement. Analyse tous les verrous mathématiques. Audits. Puis contre-audits. Puis mise à jour GitHub propre.

* * *

**ChatGPT:**

Mise à jour publiée directement sur `main`, sans branche supplémentaire : commit c46fe49. Le dépôt local et GitHub sont synchronisés et propres.

Résultats principaux :

* Une arête fausse interne à un enfant de la fusion ne contribue à aucun $\Lambda_v$ pertinent : elle est invariante sous le flip global et s’annule exactement.
* Une arête traversant $C_1\mid C_2$ reste cependant pertinente, même si elle devient interne à la composante parente après la fusion.
* Conditionnellement à la partition complète $\Pi_t$, seules les arêtes de frontière ont la loi indépendante $$(\text{vraie tardive},\text{vraie censurée},\text{fausse})  
    =  
    \left(h_p(t),\frac{1-h_p(t)}2,\frac{1-h_p(t)}2\right).$$ Cette loi ne doit jamais être appliquée aux arêtes internes.
* Conditionner une paire lointaine à être connectée au seuil sélectionne sa composante avec un biais de Palm proportionnel au nombre de paires lointaines qu’elle contient, et non uniformément ou proportionnellement à sa seule taille.
* Exactement au seuil plan, il faut parler de composante critique macroscopique, pas de géante de densité positive. Cette distinction repose sur les résultats planaires de box-crossing et de criticité (Grimmett–Manolescu).
* L’appartenance à la même composante critique implique nécessairement $\beta_{ij}\le\beta_c$. Pour une paire macroscopiquement distante, la décroissance exponentielle sous-critique et RSW localisent néanmoins le LCA vers $\beta_c$ par la gauche (Duminil-Copin–Tassion).
* Aucun meilleur seuil rigoureux de weak recovery n’est encore obtenu : le diagnostic de majorité de frontière vaut déjà $p_{\partial,\mathrm{late}}=0.782432\ldots<p_{\mathrm{info}}$. L’amélioration éventuelle doit venir du critère quatre états complet, de la géométrie Palm des coupes, de toute la chaîne des $\Lambda_v$, et de la domination HF.

Fichiers principaux :

Contre-audits validés : 66 tests passent, notamment l’énumération exhaustive conditionnelle sur un triangle, l’annulation des arêtes internes et le biais exact de Palm. Les 21 fichiers Markdown passent également le vérificateur mathématique.

* * *

**You:**

Va lire mon manuscrit de thèse (le chapitre 11). Je cherche un moyen d'adapter le théorème principal qui permet de lier le couplage Swendsen-Wang généralisé avec des résultats pour la wek recovery. Je cherche à adapter ce théorème au cas de la dynamique hiérarchique. Evidemment, pour obtenir des résultats plus forts. Quelle serait la "meilleure" dynamique hiérarchique pour mesurer la weak recovery ?

* * *

**ChatGPT:**

J’ai relu la version actuelle du chapitre 11, puis les slides de dynamique hiérarchique.

Ma conclusion est nette : la meilleure dynamique pour mesurer la weak recovery n’est pas un sweep local nœud par nœud. C’est un **Swendsen–Wang hiérarchique complètement bloqué et collapsé**, qui rééchantillonne conjointement les orientations conditionnellement à une hiérarchie aussi peu informative que possible. Sa version calculable est un **heat bath bloqué du LCA et de tout son cône ancestral**.

## 1. Ce qu’il faut réellement généraliser

Le théorème principal du chapitre utilise deux ingrédients distincts :

1. toute transition laissant la postérieure invariante permet de remplacer la vérité par sa sortie, grâce à Nishimori ;
2. conditionnellement aux clusters gelés, leurs recoloriages sont indépendants et uniformes.

Le second point est celui qui échoue dans la hiérarchie : les facteurs ancestraux

$$F_v(\Lambda_v)=\Lambda_v e^{(1-\beta_v)\Lambda_v}$$

couplent les orientations des descendants. Une « composante géante ou non » ne contient donc plus assez d’information.

La bonne généralisation consiste à remplacer les blocs indépendants par la **matrice de corrélation conditionnelle de toute la hiérarchie**.

## 2. Théorème hiérarchique proposé

Je l’énonce pour $K=2$, le cas du GSBM triangulaire. Soit $O$ l’observation et $D$ le dendrogramme auxiliaire, sous la loi jointe exacte $\nu_O$. Posons

$$C_D  
=  
\mathbb E_{\nu_O}\!\left[\sigma\sigma^\top\mid O,D\right].$$

Cette matrice est positive, de diagonale $1$. Pour $\delta>0$, définissons la masse spectrale persistante

$$\vartheta_{n,\delta}(D)  
=  
\frac1n  
\min\left\{  
|S|:  
\lambda_{\max}\!\left(C_D[S^c,S^c]\right)  
\le \delta n  
\right\},$$

puis

$$\vartheta_D^{\max}  
=  
\inf\left\{  
\vartheta:  
\lim_{\delta\downarrow0}\limsup_{n\to\infty}  
\mathbb P\bigl(\vartheta_{n,\delta}(D)>\vartheta\bigr)=0  
\right\}.$$

Alors, pour tout estimateur $\tau_n$ et tout $\eta>0$,

$$\boxed{  
\mathbb P\left[  
\mathrm{ov}_n(\Sigma_n,\tau_n)  
\ge  
\frac12+\frac12\bigl(\vartheta_D^{\max}+\eta\bigr)  
\right]  
\longrightarrow0.  
}$$

### Preuve

On utilise le noyau bloqué

$$K_D^\star(\sigma,d\sigma')  
=  
\nu_O(d\sigma'\mid O,D).$$

Il est invariant, donc Nishimori permet de remplacer $\Sigma$ par $\sigma'$. Conditionnellement à $O,D,\tau$, choisissons $S$ réalisant la définition de $\vartheta_{n,\delta}$. Avec $T=S^c$,

$$\left|\frac1n\sum_i\tau_i\sigma'_i\right|  
\le  
\frac{|S|}{n}  
+  
\left|\frac1n\sum_{i\in T}\tau_i\sigma'_i\right|.$$

Or

$$\mathbb E\left[  
\left(\frac1n\sum_{i\in T}\tau_i\sigma'_i\right)^2  
\middle|O,D,\tau  
\right]  
=  
\frac{\tau_T^\top C_D[T,T]\tau_T}{n^2}  
\le\delta.$$

Tchebychev, puis $n\to\infty$ et $\delta\downarrow0$, donnent le résultat.

### Vérification Swendsen–Wang

Pour une recoloration indépendante des clusters classiques,

$$C_D=\bigoplus_C s_Cs_C^\top,$$

où $s_C\in\{-1,+1\}^{|C|}$ encode les relations figées dans $C$. Chaque bloc a pour unique valeur propre non nulle $|C|$. Ainsi, à l’arrondi près,

$$\vartheta_{n,\delta}(D)  
=  
\frac1n\sum_C (|C|-\delta n)_+.$$

Lorsque $\delta\downarrow0$, cela redonne exactement la masse macroscopique $\theta^{\max}$ du chapitre 11. Le nouveau théorème est donc une véritable extension, pas seulement une analogie.

Il peut être strictement plus fort : une composante géométrique de taille $n$ peut subsister, mais si ses orientations internes se décomposent en modes de taille $o(n)$, alors $\vartheta_D^{\max}=0$, malgré $\theta^{\max}=1$.

## 3. La quantité exacte pour la weak recovery

Pour la weak recovery à probabilité positive, l’objet naturel est

$$Q_n  
=  
\frac1{n^2}  
\mathbb E\,\mathrm{tr}(C_O^2),  
\qquad  
C_O=\mathbb E[\sigma\sigma^\top\mid O].$$

Dans le cas binaire uniforme,

$$\liminf Q_n>0  
\quad\Longleftrightarrow\quad  
\text{weak recovery}.$$

Comme $C_O=\mathbb E[C_D\mid O]$, Jensen donne

$$\boxed{  
Q_n  
\le  
B_n(D)  
:=  
\frac1{n^2}\mathbb E\,\mathrm{tr}(C_D^2).  
}$$

Cette quantité fournit la hiérarchie de bornes suivante :

$$Q_n  
\le B_n(D)  
\le H_n^{\mathrm{LCA}}  
\le  
\frac1{n^2}\mathbb E\sum_{R\text{ racine}}|R|^2.$$

La seconde inégalité vient du fait que le heat bath local au LCA conditionne sur davantage d’information que $D$ seul. C’est exactement l’inégalité de projection conditionnelle déjà sous-jacente au score LCA.

Ainsi :

* la percolation ne conserve que les tailles des racines ;
* le score LCA conserve la fiabilité locale $\eta_u$ ;
* $B_n(D)$ marginalise conjointement toutes les orientations et est plus fort.

## 4. Quelle dynamique est « meilleure » ?

À $D$ fixé, la réponse optimale dans la classe naturelle des noyaux de heat bath positifs et auto-adjoints est

$$\boxed{  
K_D^\star f  
=  
\mathbb E[f\mid O,D].  
}$$

En effet, pour tout noyau conditionnel positif $K_D$ laissant $\nu_O(\cdot\mid D)$ invariante,

$$\langle f,K_Df\rangle  
=  
\left\|\mathbb E[f\mid O,D]\right\|_2^2  
+  
\langle f-\mathbb E[f\mid O,D],  
K_D(f-\mathbb E[f\mid O,D])\rangle,$$

donc

$$\langle f,K_Df\rangle  
\ge  
\left\|\mathbb E[f\mid O,D]\right\|_2^2.$$

Le rééchantillonnage conditionnel complet atteint l’égalité simultanément pour toutes les observables. Il est donc le benchmark optimal.

Mais $K_D^\star$ peut être aussi difficile à simuler que la postérieure elle-même. La meilleure version exploitable est alors :

1. choisir une paire lointaine $(i,j)$, ou directement un nœud $u$ avec poids proportionnel au nombre de paires lointaines qu’il sépare ;
2. prendre $u=\mathrm{LCA}(i,j)$ ;
3. rééchantillonner conjointement la parité de $u$ et toutes les orientations de son cône ancestral ;
4. utiliser un random scan symétrique ou un sweep palindromique ;
5. rafraîchir $D$ entre les macro-étapes.

Le poids naturel du nœud est

$$N_u(R)  
=  
\#\{(i,j)\in C_{u,1}\times C_{u,2}:d(i,j)\ge R\},$$

et non un poids uniforme sur les nœuds. C’est exactement la mesure de Palm pertinente pour la weak recovery longue portée.

## 5. Le heat bath ancestral collapsé

Écrivons $r_v\in\{-1,+1\}$ pour les variables d’orientation hiérarchiques et

$$w_D(r)  
=  
\mu_0(\sigma(r))  
\prod_v F_v(\Lambda_v(\sigma(r))).$$

Pour le cône ancestral $A(u)$, posons

$$Z_u^\pm  
=  
\sum_{\substack{r_{A(u)}\\r_u=\pm1}}  
w_D(r).$$

Le log-rapport pertinent devient

$$L_u^{\mathrm{coll}}  
=  
\log\frac{Z_u^+}{Z_u^-},  
\qquad  
\eta_u^{\mathrm{coll}}  
=  
\tanh^2\!\left(\frac{L_u^{\mathrm{coll}}}{2}\right).$$

C’est plus fort que le heat bath local des slides, car ce dernier conditionne sur les orientations ancestrales courantes. Par inclusion des tribus,

$$\mathbb E\,\eta_u^{\mathrm{coll}}  
\le  
\mathbb E\,\eta_u^{\mathrm{LCA}}.$$

Surtout, cette formulation traite correctement les $\Lambda_v$ au-dessus de $u$ : ils ne sont ni ignorés ni remplacés par $B_u=0$, mais sommés dans un rapport de fonctions de partition.

Pour chaque ancêtre $v\succ u$, les statistiques minimales restent

$$\Lambda_v^{ab}  
=  
\frac12\left[  
T_{v,0}+T_{v,1}+T_{v,2}  
+X_{v,0}+(-1)^aX_{v,1}+(-1)^bX_{v,2}  
\right].$$

Aucune arête fausse interne à un enfant n’intervient. Seules les arêtes des frontières ancestrales courantes comptent.

## 6. Spécialisation au cas critique favorable

Ta réduction favorable doit être formulée en volume fini par

$$d(i,j)\ge \rho L,  
\qquad  
\beta_c-\varepsilon_L  
\le \beta_{ij}\le\beta_c,$$

puis $L\to\infty$ et $\varepsilon_L\downarrow0$. Une égalité exacte $\beta_{ij}=\beta_c$ a probabilité nulle.

Dans cette expérience :

* la composante est sélectionnée avec un biais proportionnel à son nombre de paires lointaines ;
* tous les ancêtres stricts $v\succ u$, jusqu’à la coupe $1$, restent présents ;
* conditionnellement au squelette, la qualité d’une marque de frontière à l’ancêtre $v$ vaut

$$s_p(\beta_v)  
=  
\mathrm{logistic}\!\left(u_p(1-\beta_v)\right).$$

Elle décroît vers $1/2$ pour les ancêtres tardifs. Le verrou n’est donc pas seulement une majorité au nœud critique : c’est le produit/transfert des messages signés sur toute la chaîne ancestrale.

La dynamique recommandée dans ce régime est donc le **heat bath du quotient critique, bloqué sur le cône ancestral du LCA**, avec les composantes sous-critiques contractées.

## 7. Contre-audit

Les restrictions suivantes sont indispensables :

* « Optimal » vaut dans la classe des noyaux de heat bath positifs. Des noyaux invariants avec autocorrélations négatives peuvent artificiellement battre un score à un pas.
* Un MAP ou une règle de majorité n’est pas un noyau invariant ; elle ne peut pas remplacer le heat bath dans le couplage de Nishimori.
* Retenir davantage d’information dans $D$ augmente généralement $B_n(D)$. Si $\bar D=g(D)$, $$B_n(\bar D)\le B_n(D).$$ À l’extrême, $D$ constant donne $B_n=Q_n$, mais cette « solution » est circulaire.
* Tronquer simplement $D$ à $\beta_c$ sans recalculer la loi conditionnelle est faux : le facteur résiduel des slides doit être intégré.
* Le LCA pair-spécifique fournit une borne sur $Q_n$, mais pas directement une unique dynamique globale pour l’overlap. Le théorème spectral ci-dessus résout cette différence.
* Conditionner sur la composante critique est un oracle. Pour en tirer une impossibilité globale, il faut encore le lemme de domination favorable HF ; pour une preuve de récupération, cet oracle est inutilisable tel quel.

## Verdict

La meilleure dynamique théorique est :

$$\boxed{  
\text{hiérarchie la plus grossière encore calculable}  
\;+\;  
\text{rééchantillonnage conditionnel conjoint de toutes ses orientations}.  
}$$

La meilleure approximation concrète est :

$$\boxed{  
\text{random scan pondéré par les paires lointaines}  
\;+\;  
\text{heat bath collapsé du LCA et de tout son cône ancestral}.  
}$$

Et la bonne observable n’est plus seulement $\beta_c$, $\theta^{\max}$, ou la majorité locale, mais

$$B_n(D)=\frac1{n^2}\mathbb E\,\mathrm{tr}(C_D^2)$$

pour la weak recovery qualitative, et $\vartheta_D^{\max}$ pour la fraction récupérable avec haute probabilité.

Cette voie peut réellement battre le théorème de percolation du chapitre 11, y compris en présence d’une composante géante. En revanche, un nouveau seuil numérique sur le GSBM triangulaire exige encore le contrôle Palm de la chaîne ancestrale et de son rapport de fonctions de partition ; aucun seuil amélioré n’est encore démontré par les calculs locaux seuls.

* * *

**You:**

Travaillons avec méthode. Ce qu'il nous faudrait pour commencer ce sont des lemmes pour p fixé et beta entre 0 et 1, et un cluster au niveau beta, savoir quelle est la proportion d'arêtes internes à un cluster fausses, et celle d'horloges non encore activées. Est-ce que cette proportion dépend de si on fait l'hypothèse d'être dans la composante géante ?

* * *

**ChatGPT:**

Oui — mais la réponse dépend du dénominateur. Le conditionnement « composante géante » ne change pas la répartition des marques parmi les arêtes internes encore fermées. Il change en revanche leur proportion parmi toutes les arêtes potentielles internes, via la densité d’arêtes déjà ouvertes.

Je travaille sous la loi annealed, dans la jauge de la ground truth.

## 1. Les quatre catégories d’une arête

Fixons $p\in(1/2,1)$, $q=1-p$,

$$u_p=\log\frac pq,  
\qquad  
a_\beta=p(1-e^{-u_p\beta}),  
\qquad  
b_\beta=pe^{-u_p\beta}.$$

Ainsi $a_\beta$ est le paramètre de la percolation ouverte au temps $\beta$, et

$$1-a_\beta=q+b_\beta.$$

Pour $\beta\in[0,1]$, une arête potentielle appartient exactement à l’une des catégories suivantes.

| Catégorie | Probabilité |
| --- | --- |
| conforme et déjà ouverte, $\xi_e\le\beta$ | $a_\beta$ |
| conforme, $\beta<\xi_e\le1$ | $b_\beta-q$ |
| conforme, $\xi_e>1$ | $q$ |
| fausse | $q$ |

En particulier, la masse totale des arêtes conformes dont l’horloge n’est pas encore activée à $\beta$ vaut $b_\beta$.

## 2. Lemme fondamental : loi résiduelle des arêtes fermées

Posons

$$A_e(\beta)  
=  
\mathbf 1_{\{Z_e=+1,\ \xi_e\le\beta\}}.$$

### Lemme 1 — marquage résiduel exact

Conditionnellement au graphe ouvert complet $A(\beta)$, les marques des arêtes telles que $A_e(\beta)=0$ sont indépendantes. Leur loi est

$$\boxed{  
\mathbb P(\text{fausse}\mid A_e=0)  
=  
f_\beta  
:=  
\frac{q}{q+b_\beta},  
}$$ $$\boxed{  
\mathbb P(\text{conforme mais non activée}\mid A_e=0)  
=  
s_\beta  
:=  
\frac{b_\beta}{q+b_\beta}.  
}$$

En distinguant la censure à $1$,

$$\left(  
\text{activation dans }(\beta,1],  
\text{censurée après }1,  
\text{fausse}  
\right)  
\sim  
\left(  
\frac{b_\beta-q}{q+b_\beta},  
\frac{q}{q+b_\beta},  
\frac{q}{q+b_\beta}  
\right).$$

### Conséquence importante

Le lemme reste inchangé après tout conditionnement mesurable par rapport à $A(\beta)$, notamment :

* $C$ est une composante de $\Pi_\beta$ ;
* $C$ est la plus grande composante ;
* $C$ est la composante géante ;
* deux sommets donnés appartiennent à $C$.

Donc l’hypothèse « être dans la géante » ne change pas $f_\beta$ ni $s_\beta$.

Elle change seulement le nombre et la géométrie des arêtes internes fermées.

## 3. Application à un cluster

Soit $C\in\Pi_\beta$. Il faut considérer toutes les arêtes potentielles

$$E(C)=\{\{x,y\}\in E:x,y\in C\}.$$

Posons

$$M_C=|E(C)|,  
\qquad  
O_C=\sum_{e\in E(C)}A_e(\beta),  
\qquad  
N_C=M_C-O_C.$$

Ici $O_C$ est le nombre d’arêtes internes déjà ouvertes et $N_C$ celui des arêtes internes encore fermées.

Notons :

* $F_C$ le nombre d’arêtes internes fausses ;
* $U_C$ le nombre d’arêtes conformes internes dont l’horloge est $>\beta$.

### Lemme 2 — décomposition interne conditionnelle

Conditionnellement à $A(\beta)$,

$$\boxed{  
F_C\sim\mathrm{Bin}(N_C,f_\beta),  
\qquad  
U_C=N_C-F_C.  
}$$

En particulier,

$$\mathbb E[F_C\mid A(\beta)]  
=  
N_C\frac{q}{1-a_\beta},$$ $$\mathbb E[U_C\mid A(\beta)]  
=  
N_C\frac{b_\beta}{1-a_\beta}.$$

Et, si $N_C>0$,

$$\mathbb P\left(  
\left|\frac{F_C}{N_C}-f_\beta\right|>\varepsilon  
\middle|A(\beta)  
\right)  
\le  
2e^{-2\varepsilon^2N_C}.$$

Ainsi, dès que $N_C\to\infty$,

$$\frac{F_C}{N_C}\longrightarrow f_\beta,  
\qquad  
\frac{U_C}{N_C}\longrightarrow s_\beta,$$

même si $C$ est conditionné à être la géante.

## 4. Proportions parmi toutes les arêtes internes

Posons la densité interne ouverte

$$\alpha_C(\beta)=\frac{O_C}{M_C}.$$

Alors, si $N_C\to\infty$,

$$\boxed{  
\frac{F_C}{M_C}  
=  
(1-\alpha_C(\beta))  
\frac{q}{1-a_\beta}  
+o_{\mathbb P}(1),  
}$$ $$\boxed{  
\frac{U_C}{M_C}  
=  
(1-\alpha_C(\beta))  
\frac{b_\beta}{1-a_\beta}  
+o_{\mathbb P}(1).  
}$$

C’est ici, et seulement ici, que l’hypothèse de composante géante intervient : elle modifie $\alpha_C(\beta)$.

Si l’on veut la proportion d’horloges non activées parmi les seules arêtes conformes internes, elle converge vers

$$\frac{(1-\alpha_C)s_\beta}  
{\alpha_C+(1-\alpha_C)s_\beta}.$$

Il n’existe donc pas de formule ne dépendant que de $p$ et $\beta$, sauf si l’on connaît aussi la limite de $\alpha_C$.

## 5. Conditionnement par un ensemble de sommets fixé

Fixons un ensemble $C\subset V$, avec

$$M=|E(C)|.$$

Conditionnellement à l’événement « $C$ est exactement une composante de $\Pi_\beta$ » :

* toutes les arêtes de la frontière de $C$ sont fermées ;
* les arêtes internes sont des Bernoulli$(a_\beta)$ conditionnées à rendre $G[C]$ connexe.

Notons

$$R_C(a)  
=  
\mathbb P_a\bigl(G[C]\text{ est connexe}\bigr)$$

le polynôme de fiabilité de $G[C]$.

### Lemme 3 — densité interne exacte sous connexité

On a

$$\boxed{  
\mathbb E[O_C\mid C\in\Pi_\beta]  
=  
a_\beta M  
+  
a_\beta(1-a_\beta)  
\frac{d}{da}\log R_C(a)\Big|_{a=a_\beta}.  
}$$

Par conséquent,

$$\boxed{  
\mathbb E\left[\frac{F_C}{M}\middle|C\in\Pi_\beta\right]  
=  
q\left[  
1-\frac{a_\beta}{M}  
\frac{d}{da}\log R_C(a)\Big|_{a=a_\beta}  
\right],  
}$$

et

$$\boxed{  
\mathbb E\left[\frac{U_C}{M}\middle|C\in\Pi_\beta\right]  
=  
b_\beta\left[  
1-\frac{a_\beta}{M}  
\frac{d}{da}\log R_C(a)\Big|_{a=a_\beta}  
\right].  
}$$

Comme la connexité est un événement croissant,

$$\frac{d}{da}R_C(a)\ge0.$$

Donc

$$\mathbb E\left[\frac{F_C}{M}\middle|C\in\Pi_\beta\right]\le q,  
\qquad  
\mathbb E\left[\frac{U_C}{M}\middle|C\in\Pi_\beta\right]\le b_\beta.$$

La connexité biaise les arêtes internes vers davantage d’ouvertures.

Un contre-exemple extrême est instructif : si $G[C]$ est un arbre, sa connexité impose que ses $M=|C|-1$ arêtes soient toutes ouvertes. Il n’y a alors aucune arête potentielle interne fausse ou tardive.

## 6. Que change précisément la composante géante ?

Sur un graphe transitif, supposons $a_\beta>p_c$. Pour une arête typique $e=\{x,y\}$, la densité ouverte parmi les arêtes potentielles internes à la composante infinie est naturellement

$$\alpha_\infty(a_\beta)  
=  
\mathbb P_{a_\beta}  
\left(  
A_e=1  
\mid  
x,y\in\mathcal C_\infty  
\right).$$

Sous les hypothèses ergodiques usuelles, les proportions dans la géante sont alors

$$\frac{F_{\mathrm{giant}}}{M_{\mathrm{giant}}}  
\longrightarrow  
(1-\alpha_\infty)  
\frac{q}{1-a_\beta},$$ $$\frac{U_{\mathrm{giant}}}{M_{\mathrm{giant}}}  
\longrightarrow  
(1-\alpha_\infty)  
\frac{b_\beta}{1-a_\beta}.$$

Par Harris–FKG,

$$\alpha_\infty(a_\beta)\ge a_\beta.$$

Donc la géante contient, parmi toutes ses arêtes potentielles internes, au plus la proportion non conditionnée $q$ d’arêtes fausses et $b_\beta$ d’horloges conformes non activées.

En revanche, $\alpha_\infty(a_\beta)$ n’a pas de formule fermée connue, même sur la grille triangulaire. C’est désormais un problème purement géométrique de percolation.

Sur la grille triangulaire :

$$a_\beta>q_c=2\sin(\pi/18)$$

est nécessaire pour parler d’une composante infinie. À $a_\beta=q_c$, il n’y a pas de géante de densité positive : il faut parler d’une composante critique macroscopique ou de la composante sélectionnée par une paire lointaine.

## 7. Réponse courte

* Parmi les arêtes ouvertes qui constituent le cluster : la proportion de fausses est exactement $0$.
* Parmi les arêtes internes encore fermées : les proportions sont universelles, $$\text{fausses}=\frac{1-p}{1-p+pe^{-u_p\beta}},  
    \qquad  
    \text{vraies non activées}=  
    \frac{pe^{-u_p\beta}}{1-p+pe^{-u_p\beta}},$$ indépendamment de l’hypothèse de géante.
* Parmi toutes les arêtes potentielles internes : elles dépendent de la géante par le seul facteur $1-\alpha_C(\beta)$.

Dernier contre-audit : ces arêtes internes n’entrent pas dans les $\Lambda_v$ situés au-dessus du cluster. Elles restent dans le même enfant à tous les ancêtres et s’annulent dans le heat bath. Pour les $\Lambda_v$, les quantités directement pertinentes restent les marques des frontières successives, dont la loi exacte est déjà établie dans la note sur les frontières critiques.

* * *

**You:**

Parfait merci. Maintenant faisons l'hypothèseque i et j (très espacés) sont dans la composante géante. Que le noeud u fusionnant C_1 et C_2 leurs clusters soit de niveau beta_u = beta_c (associé au seuil de percolation). Je voudrais un équivalent (quand ils s'éloignent infiniment) de la probabilité d'avoir le retournement de clusters pour C_i (0,0) ou (1,1) en fonction de p (supposé donc strictement supérieur à p_c le seuil de percolation). Mets à jour le GitHub avec ces lemmes et cette piste. Audit et contre-audit les résultats. Sois très rigoureux mathématiquement. Et conserve un repo GitHub bien entretenu !

* * *

**ChatGPT:**

Le résultat est publié directement sur `main`, sans branche supplémentaire : commit e713e33.

Le résultat principal est le suivant. Pour une coupure critique de taille $m$, sans message ancestral $B_u$,

$$\overline P_m^c(p)  
=\mathbb P\bigl((a,b)\in\{(0,0),(1,1)\}\bigr)  
=\frac{1+\Gamma_m^c(p)}2.$$

En posant

$$s_c(p)=\frac{p-q_c}{1-q_c},\qquad  
a_c(p)=\log\frac{p-q_c}{1-p},$$

on obtient :

$$\overline P_m^c(p_{\mathrm{SW}})  
=\frac12+\frac1{2m},  
\qquad  
p_{\mathrm{SW}}=\frac{1+q_c}{2}\simeq0.673648,$$

et, pour tout $p>p_{\mathrm{SW}}$,

$$1-\overline P_m^c(p)  
\sim  
\frac{C_{m\bmod 2}(p)}{\sqrt m}  
e^{-mI_c(p)},$$

avec

$$I_c(p)  
=  
D\!\left(\frac12\middle\|s_c(p)\right)  
=  
-\frac12\log(1-h_c(p)^2),$$

et deux préfacteurs explicites $C_0(p),C_1(p)$ selon la parité de $m$.

Pour la dynamique hiérarchique complète,

$$L_u=B_u+\ell^c_{M_L,K_L}(p).$$

Sous les hypothèses

$$M_L\longrightarrow\infty,  
\qquad  
\frac{B_L}{M_L}\longrightarrow0,$$

la probabilité paire tend vers $1$. Mais le contre-audit montre que la distance $d(i,j)\to\infty$ et l’appartenance à une composante macroscopique n’impliquent pas automatiquement $M_L\to\infty$ : deux grands amas peuvent fusionner par une seule arête pivotale.

Autres conclusions importantes :

* les arêtes internes, fausses ou vraies tardives, n’entrent pas dans les $\Lambda_v$ ancestraux ;
* conditionnellement au graphe ouvert, l’hypothèse de composante géante ne modifie pas la qualité des marques encore fermées ; elle modifie seulement leur géométrie et leur nombre ;
* à $\beta_c$, la grille infinie n’a pas de géante de densité positive : l’hypothèse doit être interprétée par exhaustion finie, fenêtre critique ou loi de Palm.

La démonstration complète, les audits et contre-audits sont dans 15_CRITICAL_GIANT_PAIR_FLIP.md. Le README, la roadmap et les scripts de calcul ont aussi été mis à jour. Validation finale : 75 tests réussis, compilation Python et vérification intégrale des formules Markdown.

* * *

**You:**

Tu n'as pas compris où je voulais en venir. Je voudrais aussi que tu ajoutes sur le repo un certain nombre de lemmes permettant de comprendre explicitement comment calculer les probabilités de retournement. Par exemple au niveau des racines, c'est du 1/2 1/2. Au niveau des feuilles, on retrouve Metropolis-Hastings. Et au niveau d'un noeud u ?
Je me demande s'il n'y aurait pas moyen de le voir aussi avec les descendants plutôt que les ascendants. Par exmple, pour i et j très éloignés, il suffirait d'observer leur chemin dans l'arbre couvrant. Avec des niveaux qui sont des exponentielles plus petites que beta. Comme tous les poids sont les mêmes dans le GSBM triangulaire, on peut certainement en déduire quelque chose asymptotiquement sur les probabilités que i et j soient conformes ou inversés par rapport à la ground truth. mêmes consignes : audit + contre audit et garde un repo pédagogique avec ces lemmes en plus sur la dynamique hiérarchique.

* * *

**ChatGPT:**

Tout est publié directement sur `main`, commit `aa19e7c`. Le dépôt local et GitHub sont propres et synchronisés.

Résultats principaux :

* Racine finale : retournement $1/2,1/2$ sous a priori uniforme.
* Feuille : la dynamique définie donne exactement Glauber/Barker, $$P_{\mathrm{HB}}=\frac{R}{1+R},$$ tandis que Metropolis–Hastings donne $\min(1,R)$. Même cible, probabilités différentes.
* Nœud interne : les quatre probabilités $p_u^{00},p_u^{01},p_u^{10},p_u^{11}$ sont explicites via les trois coordonnées de Walsh $(h_1,h_2,J)$. En particulier, $$P_u^{\rm pair}  
    =  
    \frac{e^J\cosh(h_1+h_2)}  
    {e^J\cosh(h_1+h_2)+e^{-J}\cosh(h_1-h_2)}.$$
* Sans message ancestral, pour une coupe homogène de taille $m$ au niveau $t$, $$K=1+\mathrm{Bin}(m-1,s_p(t)),  
    \qquad  
    \ell_{m,k}  
    =  
    \log\frac{k}{m-k}+u_p(1-t)(2k-m),$$ et $$\mathbb E p^{00}=\mathbb E p^{11}=\frac{1+\Gamma_m(t;p)}4,  
    \qquad  
    \mathbb E p^{01}=\mathbb E p^{10}=\frac{1-\Gamma_m(t;p)}4.$$ À $t=1$, les deux premières valent $1/4+1/(4m)$. Pour $t<1$ fixé et $m\to\infty$, elles tendent chacune vers $1/2$.

Pour la voie descendante, l’identité de parité sur les deux bras vers le LCA est exacte. Dans l’oracle factorisé PATH-FAC,

$$P_{ij}^{\rm même\ relation}  
=  
\frac12\left(1+\prod_{w\in\mathcal P(i,j)}  
\Gamma_{m_w}(t_w;p)\right),$$

et, après recoloration globale équitable,

$$P(\text{tous deux conformes})  
=  
P(\text{tous deux inversés})  
=  
\frac14\left(1+\prod_w\Gamma_w\right).$$

Pour le GSBM triangulaire critique, PATH-FAC conserve asymptotiquement la relation vraie si, pour un chemin de longueur $H$,

$$m_{\min}\ge  
\bigl(c_c(p)^{-1}+\varepsilon\bigr)\log H.$$

Le contre-audit montre cependant que cette factorisation n’est pas valable automatiquement dans la vraie dynamique : les décisions de heat bath sont dépendantes. La formulation exacte est désormais une récursion de transfert tordue sur les messages de frontière. C’est le verrou précis restant pour transformer cette piste en meilleure borne de weak recovery sur la grille triangulaire.

Fichiers principaux :

* Note mathématique complète
* Calculs reproductibles
* Contre-tests

Validation : 89 tests réussis et 23 fichiers Markdown validés mathématiquement.

* * *

**You:**

Continue d'explorer cette voie. Il faut repérer quand est-ce que (en fonction de p) asymptotiquement cette probabilité se rapproche de 1/2 et donc que i et j perdent leur corrélation.

* * *

**ChatGPT:**

Le critère de décorrélation est maintenant explicite et publié sur `main`.

Dans l’oracle PATH-FAC,

$$P_L(p)=\frac12(1+e^{-A_L(p)}),\qquad  
A_L(p)=-\sum_{w\in\mathcal P_L}\log\Gamma_{m_w}(t_w;p).$$

Donc $P_L\to 1/2$ exactement lorsque $A_L\to\infty$.

Résultats principaux :

* À $p<1$ fixé, un nombre divergent de coupes $2\le m_w\le M$ suffit à perdre la corrélation.
* Pour de grandes coupes hétérogènes, le critère devient $$\Phi_L^{\rm desc}(p)  
    =\sum_w m_w^{-1/2}e^{-m_w I(t_w;p)}.$$ Sa divergence donne $P_L\to1/2$, sa convergence vers zéro donne $P_L\to1$.
* Si $m\sim\alpha\log H$ au niveau critique, $$p_{\rm path}(\alpha)=  
    \frac{1+q_\triangle+(1-q_\triangle)\sqrt{1-e^{-2/\alpha}}}{2}.$$
* Seuls les descendants vérifiant $$\beta_c-t_w=O(m_w^{-1})$$ ont le même poids exponentiel qu’une fusion critique. Le fait que le seul LCA soit critique ne suffit donc pas.
* Pour $\alpha=7.053596\ldots$, calibré sur Nishimori, le seuil régulier vaut :
    
    | Niveau $t=\theta\beta_c$ | Seuil |
    | --- | --- |
    | $0$ | $0.7484399$ |
    | $\beta_c/2$ | $0.7983165$ |
    | $\beta_c$ | $0.8358058$ |
    

Ainsi, retrouver Nishimori par cette voie nécessiterait qu’une fraction géométriquement suffisante du chemin soit concentrée dans la bande critique $O(1/m)$. Sinon le seuil redescend rapidement vers la baseline $0.794659\ldots$.

Pour la dynamique jointe réelle, la factorisation des $\Gamma_w$ reste non démontrée. Le critère rigoureux obtenu est

$$|c_{ij}^{\rm joint}|\le\prod_r\kappa_r,$$

où les $\kappa_r$ sont les normes des opérateurs de transfert tordus. Il faut désormais établir une contraction de blocs et mesurer la loi jointe des $(m_w,t_w)$.

Le contre-audit de la littérature conclut — par inférence à partir de la portée des théorèmes existants — que ni la mesure pivotale, ni la limite du MST, ni les outlets d’invasion ne donnent cette loi de buckets sous le conditionnement de paire critique : mesure pivotale, MST planaire, outlets d’invasion.

Dépôt mis à jour :

* Commit GitHub : `b5fb51c`

Validation : 103 tests réussis, 24 fichiers Markdown audités, compilation Python réussie, `main` propre et synchronisé, sans branche supplémentaire.

* * *

**You:**

Prends du recul. Il n'y a donc pas moyen d'utiliser cette dynamique hiérarchique pour retrouver facilement le seuil de Nishimori pour la weak recovery. À tout le moins, battre les seuils précédents ?

* * *

**ChatGPT:**

Oui : retrouver « facilement » Nishimori avec le chemin hiérarchique paraît irréaliste. Mais battre la meilleure borne rigoureuse actuelle $0.794659\ldots$ reste crédible — avec une cible beaucoup plus raisonnable que $0.8358058\ldots$.

| Seuil | Statut |
| --- | --- |
| $0.673648\ldots$ | Swendsen–Wang |
| $0.719224\ldots$ | dynamique triangulaire |
| $0.794659\ldots$ | meilleure borne rigoureuse, information-percolation |
| $0.800000$ | cible rationnelle réaliste |
| $0.809909\ldots$ | candidat conditionnel multi-état |
| $0.8358058\ldots$ | Nishimori conjectural |

### Pourquoi la voie descendante n’est probablement pas la bonne première attaque

PATH-FAC demande simultanément :

1. la loi fine des $(m_w,t_w)$ sous le biais d’une paire critique ;
2. une quantité suffisante de niveaux dans la bande $\beta_c-t_w=O(1/m_w)$ ;
3. une comparaison entre le produit des $\Gamma_w$ et la dynamique jointe.

Même si les deux premières étapes réussissent, la troisième reste un vrai théorème de dépendance. Cette voie peut expliquer un seuil, mais elle ne semble pas la manière la plus courte de le prouver.

La voie ancestrale des $\Lambda_v$ est encore plus délicate pour une borne d’impossibilité : une grande coupe critique possède un LLR local très informatif pour tout $p>p_{\rm SW}$. Pour obtenir une contraction, le message ancestral devrait satisfaire approximativement

$$B_u\simeq-\ell_u^{\rm crit}.$$

Il faudrait donc démontrer une annulation structurée, alors que les ancêtres peuvent tout aussi bien renforcer le signal.

### La cible réaliste : $p=4/5$

Le canal multi-état des triangles donne une attaque beaucoup plus courte. À

$$p=\frac45,\qquad q=2p-1=\frac35,$$

l’enveloppe d’effacement devient entièrement rationnelle :

$$(a,s,e)=\frac1{1445}(441,162,518).$$

Elle est strictement sous-critique :

$$e-a=\frac{77}{1445}>0.$$

Si l’on prouve l’inégalité less-noisy correspondante pour tout a priori quatre états, alors on obtient rigoureusement

$$p<\frac45  
\quad\Longrightarrow\quad  
\text{pas de weak recovery}.$$

Ce serait une amélioration réelle sur $0.794659\ldots$.

Le secteur non polarisé $\max_x\mu_x\le1/2$ est déjà traité. Le verrou restant est la positivité d’une matrice rationnelle $3\times3$ lorsque un état possède une masse $>1/2$. Les audits numériques ne trouvent aucun contre-exemple, mais une certification par fractions de Bernstein, arithmétique d’intervalles ou sommes de carrés reste nécessaire. C’est nettement plus circonscrit que la géométrie complète du chemin de Kruskal.

### Où utiliser réellement la dynamique hiérarchique

Son meilleur rôle n’est probablement pas de multiplier naïvement des fiabilités locales, mais de contrôler la polarisation créée par l’information extérieure :

* sur les nœuds où $\max\pi_u^{ab}\le1/2$, la comparaison multi-état est déjà disponible ;
* sur les nœuds polarisés, les $\Lambda_v^{ab}$ permettent de quantifier exactement le défaut ;
* on peut alors payer explicitement la masse des mauvais nœuds dans une contraction de blocs.

Autrement dit, la hiérarchie pourrait fermer le secteur polarisé que la SDPI triangulaire laisse ouvert. C’est une synergie plus convaincante que PATH-FAC.

### Mon ordre de travail recommandé

1. Prouver d’abord la borne rationnelle $p_\star=4/5$.
2. Optimiser ensuite la certification jusqu’au candidat $0.8099092892\ldots$.
3. Construire une suite de blocs hiérarchiques de taille croissante — triangle, cactus, bande — et étudier si leurs seuils augmentent.
4. Ne revenir à Nishimori que si cette suite semble converger naturellement vers $0.8358058\ldots$.

Donc mon verdict est : Nishimori n’est pas à portée facile de la dynamique actuelle, mais une nouvelle borne autour de $0.80$, puis potentiellement $0.809909\ldots$, est une cible sérieuse et mathématiquement beaucoup plus propre. Le dossier pertinent est ; la difficulté proprement hiérarchique reste résumée dans .

* * *

**You:**

Est-ce qu'il y aurait déjà un bon moyen d'estimer asymptotiquement la probabilité que i et j soient assignés conformément à la ground truth lorsqu'il sont très éloignés et sous l'hypothèse qu'ils apparaissent dans la même composante au beta associé à la percolation.
On prendra pour dynamique : un parcours complet (top-down ou bottom-up) de la hiérarchie avec possible retournement sur chacun des noeuds.

* * *

**ChatGPT:**

Oui. Le bon objet asymptotique n’est pas une probabilité locale de retournement, mais un coefficient de corrélation obtenu par un produit d’opérateurs de transfert « tordus » le long des deux branches reliant $i$ et $j$ à leur LCA.

Posons

$$\mathcal F_R=\bigl\{d(i,j)=R,\; i\leftrightarrow j  
\text{ dans }G_{\beta_c}\bigr\},  
\qquad R\longrightarrow\infty,$$

et supposons la configuration initiale égale à la ground truth $\Sigma$.

### 1. La condition de connexion force effectivement le LCA vers $\beta_c$

Si $\beta_{ij}$ est le niveau où les composantes de $i$ et $j$ fusionnent, alors, pour tout $\varepsilon>0$,

$$\mathbb P\!\left(  
\beta_{ij}\leq \beta_c-\varepsilon  
\,\middle|\,\mathcal F_R  
\right)  
=  
\frac{  
\mathbb P_{q_p(\beta_c-\varepsilon)}(i\leftrightarrow j)  
}{  
\mathbb P_{q_c}(i\leftrightarrow j)  
}.$$

Le numérateur décroît exponentiellement en $R$, puisqu’il correspond à une percolation strictement sous-critique. Le dénominateur ne décroît que polynomialement, par RSW et les estimations multi-échelles critiques. Par conséquent,

$$\boxed{\quad  
\beta_{ij}\xrightarrow[R\to\infty]{\mathbb P(\,\cdot\,\mid\mathcal F_R)}  
\beta_c .  
\quad}$$

Ainsi, il n’est pas nécessaire d’imposer artificiellement $\beta_{ij}=\beta_c$ : la condition $i\leftrightarrow j$ à $\beta_c$ produit asymptotiquement cette localisation critique. Sous les exposants near-critical standards, on s’attend plus précisément à

$$\beta_c-\beta_{ij}=O_{\mathbb P}(R^{-3/4}),$$

mais le premier énoncé ne requiert pas cet exposant précis. Voir les résultats near-critical de Nolin et les résultats d’universalité isoradiale de Grimmett–Manolescu.

### 2. Identité exacte pour un parcours complet

Numérotons dans l’ordre du parcours les nœuds actualisés $v_1,\dots,v_H$. Soient

* $A_r\in\{0,1\}$, l’indicateur que le cluster actualisé en $v_r$ est retourné ;
* $\chi_r(i,j)=1$ si cette actualisation sépare $i$ de $j$, c’est-à-dire retourne exactement l’un des deux.

Alors, parcours top-down ou bottom-up,

$$\frac{\sigma_i'\sigma_j'}{\sigma_i\sigma_j}  
=  
\prod_{r=1}^{H}(-1)^{A_r\chi_r(i,j)}.$$

Seuls les nœuds appartenant aux deux branches

$$i\rightsquigarrow u  
\qquad\text{et}\qquad  
j\rightsquigarrow u,$$

où $u=\mathrm{LCA}(i,j)$, interviennent directement dans cette parité. Les autres sous-arbres interviennent indirectement dans les probabilités conditionnelles de retournement.

Définissons

$$C_R^{\mathsf S}(p)  
=  
\mathbb E\!\left[  
\prod_{r=1}^{H_R}(-1)^{A_r\chi_r(i,j)}  
\,\middle|\,\mathcal F_R  
\right],$$

où $\mathsf S\in\{\mathrm{BU},\mathrm{TD}\}$ désigne l’ordre du sweep. On obtient exactement

$$\boxed{  
\mathbb P\!\left(  
\sigma_i'\sigma_j'=\Sigma_i\Sigma_j  
\,\middle|\,\mathcal F_R  
\right)  
=\frac{1+C_R^{\mathsf S}(p)}2 .  
}$$

C’est la probabilité pertinente pour la weak recovery, définie modulo le retournement global.

Si un retournement équitable de la racine impose la symétrie globale, alors

$$\mathbb P(\sigma_i'=\Sigma_i,\sigma_j'=\Sigma_j\mid\mathcal F_R)  
=\frac{1+C_R^{\mathsf S}(p)}4.$$

Il faut donc distinguer soigneusement :

* « relation $i,j$ conforme » : $(1+C_R)/2$ ;
* « les deux labels littéralement conformes » : $(1+C_R)/4$.

### 3. Le meilleur estimateur asymptotique : transfert tordu sous loi Palm critique

Soit $X_r$ l’état frontière minimal avant l’actualisation $r$ : partitions, messages des branches latérales et décisions ancestrales encore pertinentes. Si

$$Q_r(x,a,dx')$$

est le noyau joint donnant le retournement $a$ et l’état suivant, définissons

$$(\mathcal T_r^{\mathsf S}f)(x)  
=  
\sum_{a=0}^{1}  
(-1)^{a\chi_r}  
\int f(x')Q_r^{\mathsf S}(x,a,dx').$$

Alors, conditionnellement à toute la hiérarchie,

$$c_{ij}^{\mathsf S}  
=  
\lambda^{\mathsf S}  
\mathcal T_1^{\mathsf S}\cdots  
\mathcal T_{H_R}^{\mathsf S}\mathbf 1,  
\qquad  
C_R^{\mathsf S}  
=  
\mathbb E[c_{ij}^{\mathsf S}\mid\mathcal F_R].$$

C’est exact. Cela fournit aussi la borne rigoureuse

$$\boxed{  
|C_R^{\mathsf S}|  
\leq  
\mathbb E\!\left[  
\prod_{k}\kappa_k^{\mathsf S}  
\,\middle|\,\mathcal F_R  
\right],  
}$$

où $\kappa_k^{\mathsf S}$ est la norme d’un opérateur tordu portant sur un bloc de plusieurs niveaux consécutifs.

Si, sous la loi conditionnelle $\mathcal F_R$, une proportion positive des blocs vérifie uniformément

$$\kappa_k^{\mathsf S}\leq1-\delta$$

pour un $\delta>0$, alors

$$C_R^{\mathsf S}\longrightarrow0,  
\qquad  
\mathbb P(\text{relation conforme}\mid\mathcal F_R)  
\longrightarrow\frac12.$$

C’est précisément le critère de perte de corrélation recherché.

Sous une limite stationnaire ergodique de l’environnement observé depuis le chemin critique, cela devient un exposant de Lyapunov :

$$\gamma_{\mathsf S}(p)  
=  
\lim_{H\to\infty}  
\frac1H  
\log  
\left\|  
\mathcal T_1^{\mathsf S}\cdots\mathcal T_H^{\mathsf S}  
\right\|\leq0.$$

Lorsque $\gamma_{\mathsf S}(p)<0$,

$$|c_{ij}^{\mathsf S}|  
=  
\exp\!\left\{  
H_R\gamma_{\mathsf S}(p)+o(H_R)  
\right\}.$$

Pour la moyenne annealed $C_R$, il faut cependant contrôler les chemins rares ; l’exposant quenched seul ne suffit pas.

### 4. Quelle direction de parcours choisir ?

Je commencerais par le bottom-up.

Il permet d’intégrer exactement les branches latérales en messages remontants, puis de combiner les deux branches au LCA. Sur un cactus de triangles, un graphe de largeur bornée ou une approximation par séparateurs finis, les $\mathcal T_r$ deviennent des matrices finies et le calcul est exact.

Le top-down reste parfaitement définissable, mais l’état doit conserver les décisions ancestrales susceptibles d’influencer les retournements descendants. Il n’existe pas, a priori, d’ordre stochastique universel disant que top-down conserve davantage ou moins la corrélation que bottom-up : il faut comparer leurs rayons spectraux ou leurs normes tordues.

### 5. Influence réelle de l’hypothèse de composante critique

La condition $\mathcal F_R$ a deux effets essentiels :

1. elle localise le LCA près de $\beta_c$ ;
2. elle remplace l’environnement typique par un environnement critique biaisé par la connexion de la paire.

Une IIC ordinaire peut décrire localement le voisinage d’un endpoint conditionné à se connecter loin — voir Járai — mais elle ne donne pas, à elle seule, la loi de toute la hiérarchie de Kruskal entre deux endpoints. L’objet manquant est plutôt une loi Palm à deux points du backbone critique, décorée par les horloges et les marques du GSBM.

En résumé, la voie exploitable est

$$\boxed{  
\text{conditionnement critique à deux points}  
\;+\;  
\text{transfert tordu par blocs}  
\;+\;  
\text{exposant spectral/Lyapunov}.  
}$$

Elle donne déjà une caractérisation exacte en volume fini et un critère asymptotique rigoureux. Le verrou principal n’est plus la formule de retournement : c’est l’estimation de la loi Palm des états-frontières le long des deux branches critiques. Les bases correspondantes sont exposées dans et .

* * *

**You:**

Voilà exactement où je voulais en venir, bravo ! explore maintenant à fond cette voie ! le but final étant de trouver un p supérieur aux valeurs déjà connues pour lequel la weak recovery ne serait pas possibl

* * *

**ChatGPT:**

J’ai poussé directement sur `main`, sans branche supplémentaire : commit 57bc640. Le dépôt est propre et synchronisé.

Le résultat central est le suivant. Pour un sweep hiérarchique $S$, en posant

$$H_S(i,j)=\mathbb E_S[\sigma_i\sigma_i'\sigma_j\sigma_j'\mid O,D,\sigma],$$

on a rigoureusement

$$\mathbb E\!\left[\frac{\lambda_{\max}(H_S)}n\right]  
\le  
\sqrt{\mathbb E[H_S(I_n,J_n)^2]}.$$

Il faut donc faire disparaître le second moment conditionnel, et non seulement montrer que la probabilité annealed de conformité tend vers $1/2$. Son calcul exact requiert deux sweeps indépendants dans le même environnement $(O,D,\sigma)$ : c’est le transfert tordu répliqué.

Sous le conditionnement $i\leftrightarrow j$ à $\beta_c$,

$$\mathbb P(\beta_{ij}\le\beta_c-\varepsilon\mid i\leftrightarrow j\text{ à }\beta_c)  
=  
\frac{\tau_{ij}(q_p(\beta_c-\varepsilon))}  
{\tau_{ij}(q_\triangle)}  
\longrightarrow0.$$

Cela formalise proprement ton cas favorable : le LCA se localise au seuil par la gauche. La preuve combine décroissance exponentielle sous-critique et estimations RSW critiques. Aizenman–Barsky, Köhler-Schindler–Tassion.

Verdict sur le seuil :

* $p=4/5=0{,}8$ est désormais le premier certificat précis à viser, strictement au-dessus de $p_{\mathrm{info}}=0{,}794659\ldots$.
* Une preuve à $p=0{,}8$ donnerait l’impossibilité pour tout $p\le0{,}8$ par dégradation BSC.
* Cette preuve n’est pas encore obtenue. Le dossier ne présente donc pas $0{,}8$ comme un nouveau seuil démontré.
* PATH-FAC se décorrèle très rapidement : à $L=64$, la corrélation diagnostique vaut environ $0{,}0053$ à $p=0{,}8$.
* Le contre-audit avec tous les $\Lambda_v^{ab}$ montre cependant que la vraie dynamique jointe reste très corrélée sur les petits tores : second moment autour de $0{,}82$–$0{,}87$ à $L=8,p=0{,}8$, même après plusieurs sweeps. Les ancêtres et l’état-frontière transportent donc un mode global absent de PATH-FAC.

Le verrou exact est maintenant HF-S2 : démontrer que le transfert répliqué de la paire critique favorable domine celui des fusions postcritiques, puis certifier un rayon spectral $\rho_{2,S}(0.8)<1$, avec contrôle des environnements rares.

Documents principaux :

* note mathématique complète
* diagnostic géométrique Palm
* sweep hiérarchique joint exact

Validation : 114 tests réussis, contrôle mathématique des 25 fichiers Markdown réussi, compilation Python réussie. La prochaine cible rationnelle est le transfert répliqué exact sur un cactus de triangles, puis un certificat par arithmétique d’intervalles à $p=4/5$.

* * *

**You:**

Oui, c'est exactement l'idée. continue l'exploration de l'arbre de recherche de façon méthodique. Par exemple, il faudra montrer que i et j liés à beta_c est en quelque sorte le cas le plus favorable (si par exemple, il n'étaient même pas connectés à beta=1, il seraient sweepés indépendamment 1/2 1/2).
Continue l'exploration pour déjà obtenir p=0,8

* * *

**ChatGPT:**

J’ai publié directement sur `main`, sans branche ni PR supplémentaire : commit `0d3eaa9`.

L’avancée centrale est une formalisation rigoureuse du « cas favorable à $\beta_c$ » :

* Si $i,j$ appartiennent à deux racines distinctes à $\beta=1$, alors $$H_{\rm TD}(i,j)=H_{\rm BU}(i,j)=0$$ exactement, pour les deux ordres de sweep.
* Un sweep est un produit de projections de heat bath et $$\mathbb E[H_S(i,j)^2\mid O,D]  
    =\|K_Sf_{ij}\|_2^2  
    =\langle f_{ij},K_S^*K_Sf_{ij}\rangle.$$ Cela fournit aussi une décomposition pythagoricienne exacte de la perte d’information. Ce cadre rejoint la théorie des Gibbs samplers comme projections alternées : Diaconis–Khare–Saloff-Coste.
* Pour un bucket de taille $m$, $$K\mid X=+\sim1+\mathrm{Bin}(m-1,s_p(t)),\qquad  
    K\mid X=-\sim\mathrm{Bin}(m-1,1-s_p(t)).$$ J’ai démontré que, pour tout $t\ge\beta_c$, $$\mathcal E_{m,s_p(\beta_c)}  
    \succeq_{\rm Blackwell}  
    \mathcal E_{m,s_p(t)}.$$ Ainsi le niveau critique est réellement le plus favorable sous la loi complète du bucket, à taille fixée et avec information latérale conditionnellement indépendante. La preuve passe par la domination des courbes ROC et le théorème de Blackwell.
* Le contre-audit reste essentiel : la monotonie pointwise est fausse à $p=0{,}8$. Un exemple admissible donne $$\eta_{\beta_c}\simeq7.4\,10^{-8},  
    \qquad  
    \eta_{0.8}\simeq0.0693.$$ Il ne contredit pas Blackwell, car il fige le même compte anti-aligné au lieu de moyenner sous les lois correctes.

À $p=0{,}8$, le bloc $m=2$ est entièrement explicite :

$$s_c=0.693582222752\ldots,  
\qquad  
\kappa_2(b)=s_c+(1-s_c)\tanh^2(b/2).$$

Entre $\beta_c$ et $t=1$, le canal tardif s’obtient par un effacement supplémentaire de probabilité

$$1-\frac{1/2}{s_c}=0.279104937240\ldots.$$

Les diagnostics rejouables donnent des contrastes critique–tardif positifs pour les deux ordres ; à $L=8$, ils ne sont toutefois qu’à environ une erreur standard de zéro. Cette réserve est explicitement documentée.

Il manque encore deux étapes avant de pouvoir annoncer rigoureusement l’impossibilité de weak recovery à $p=0{,}8$ :

1. transporter la domination de Blackwell à travers le changement de géométrie de Kruskal et la réutilisation dynamique des buckets ;
2. prouver l’accumulation d’un nombre divergent de blocs critiques contractants.

La route top-down est maintenant prioritaire, car elle préserve le mieux l’indépendance conditionnelle requise avant que les descendants ne réutilisent le bucket.

La note complète est ici : 19_FAVORABLE_SWEEP_PROJECTIONS.md. Validation finale : 126 tests réussis, 26 fichiers mathématiques contrôlés, dépôt local propre et aligné sur `main`.

* * *

**You:**

Parfait, merci, on est sur la bonne route. Continue l'exploration ! n'oublie pas de bien tenir à jour le repo'. Cette voie est maintenant cele privilégiée, n'hésite pas à nettoyer le repo' d'anciennes idées sans intérêt. Fais un repo pédagogique : la description de la dynamique, les idées intuitives, la formalisation, les résultats prouvés, ceux qu'il reste à montrer, etc.

* * *

**ChatGPT:**

Le dépôt est à jour, propre et directement fusionné sur `main`, sans branche ni PR : commit `3eedfed`.

Principaux progrès mathématiques :

* Le corridor collapsed est formalisé comme le heat bath conjoint pair-spécifique. Il est $L^2$-optimal parmi les sweeps utilisant les mêmes nœuds.
* La globalisation vers l’impossibilité de weak recovery est prouvée paire par paire, sans supposer une matrice de Gram commune.
* La domination critique/postcritique de Blackwell se tensorise exactement sur un corridor et reste valable pour un prior arbitrairement corrélé des parités.
* Contre-audit important : cette domination exige des tailles fixées. À $p=t=4/5$, le bucket critique $m=4$ et le bucket tardif $m=2$ sont rigoureusement incomparables. Les deux témoins sont certifiés par arithmétique rationnelle, sans décision fondée sur des flottants.
* Nouvelle séparation exacte : sous la loi jointe annealed de Nishimori, la géométrie non marquée sous le seuil critique, exprimée en coordonnée $q$, ne dépend pas de $p$. Le bruit entre seulement par $$s_c(p)=\frac{p-q_\triangle}{1-q_\triangle}.$$ Cette séparation ne couvre toutefois pas encore les ancêtres postcritiques.
* À $p=0{,}8$, les blocs $m=2$ neutres ont la contraction exacte $s_c=0.693582\ldots$, mais aucune impossibilité nouvelle n’est encore annoncée.

Le README pédagogique, la feuille de route et la nouvelle note centrale distinguent désormais clairement résultats établis, résultats conditionnels, diagnostics et conjectures. Les anciennes voies utiles sont conservées comme contre-audits, mais reléguées.

Validation finale : 27 fichiers mathématiques vérifiés, 135 tests réussis, worktree propre et `main` synchronisé.

Le prochain verrou prioritaire est maintenant très précis : construire le transfert collapsed exact sur un cactus de deux puis trois triangles, avec état de bord fini et double certification indépendante.

* * *

**You:**

Continue

* * *

**ChatGPT:**

J’ai publié directement sur `main` le certificat exact du cactus triangulaire : commit 5da9389.

Le résultat principal est, sous la densité de Palm fixant le LCA au rang $q$,

$$A_h^{\mathrm{LCA}}(p,q)  
=  
\kappa_{\mathrm{flux}}(p,q)  
\kappa_{\mathrm{conn}}(p,q)^{h-1},$$

avec

$$\kappa_{\mathrm{conn}}  
=  
\frac{1+(2p-1)q-q^2}{1+q-q^2},  
\qquad  
\kappa_{\mathrm{flux}}  
=  
\frac{1+(4p-2)q-3q^2}{1+2q-3q^2}.$$

Ces deux coefficients décroissent strictement avec $q$. Sur ce modèle, fusionner exactement au seuil critique est donc bien le cas postcritique le plus favorable.

À $p=0{,}8$,

$$q_\triangle=0.347296355334,\qquad  
\beta_c=0.410716539196,$$

et l’équivalent est même une identité exacte :

$$P_h^{\mathrm{conf}}-\frac12  
=  
0.446308680939\,  
e^{-0.120189290653h}.$$

Ainsi, pour $h=40$, la conformité Nishimori vaut $0.503645301931$.

Les contre-audits importants sont explicités :

* le coefficient local n’est pas $s_c$ : les histoires « arête directe d’abord » transmettent parfaitement la relation ;
* « connecté avant $q$ » et « LCA exactement à $q$ » donnent deux coefficients différents ;
* $q_\triangle$ n’est pas le seuil du cactus, qui est quasi unidimensionnel : le conditionnement a la masse rare $c(q_\triangle)^h$ ;
* il s’agit du heat bath collapsed pair-spécifique, pas encore d’un équivalent top-down/bottom-up ;
* ce résultat ne prouve donc pas encore l’impossibilité à $p=0{,}8$ sur la grille.

La note pédagogique complète, le calcul reproductible et ses tests sont intégrés. La feuille de route pointe maintenant vers une bande triangulaire de largeur deux, puis vers la vraie géométrie Palm bidimensionnelle.

Validation finale : 149 tests réussis, contrôle MathJax réussi, dépôt propre et synchronisé sur `main`.

* * *

**You:**

Continue. Il faut vraiment utiliser la hiérarchie du mieux possible : fusion au beta_c entre i et j. Vaut-il mieux seulement retourner les deux clusters juste sous cette fusion contenant respectivement i et j ? ou vaut-il mieux continuer les retrournements jusqu'aux feuilles i et j ?

* * *

**ChatGPT:**

Il vaut mieux poursuivre jusqu’aux feuilles pour obtenir l’obstruction de weak recovery la plus forte.

Attention : retourner simultanément les deux clusters $(1,1)$ ne change jamais $\sigma_i\sigma_j$. Le LCA doit utiliser le heat bath complet à quatre états ; seuls $(1,0)$ et $(0,1)$ inversent la relation.

La comparaison générale est exacte :

$$\|P_u f_{ij}\|_2^2  
=  
\|P_{\downarrow}f_{ij}\|_2^2  
+  
\|(P_u-P_{\downarrow})f_{ij}\|_2^2.$$

Ainsi :

* le LCA seul est le cas le plus favorable à la conservation de la relation ;
* le corridor descendant collapsed est toujours au moins aussi contractant, donc préférable pour prouver l’impossibilité.

Sur le cactus, sous Palm-LCA à $\beta_c$,

$$A_h^{\mathrm{LCA\ only}}=\kappa_{\mathrm{flux}},  
\qquad  
A_h^{\mathrm{full}}  
=  
\kappa_{\mathrm{flux}}\kappa_{\mathrm{conn}}^{h-1}.$$

À $p=0{,}8$ :

* LCA seul : conformité $0.895765368433$, indépendamment de la distance ;
* corridor complet, $h=40$ : conformité $0.503645301931$.

Donc le LCA critique fournit le bon bloc pivotal, mais toute l’atténuation asymptotique provient des deux bras descendants.

Pour un seul sweep séquentiel :

* bottom-up, avec le LCA mis à jour physiquement en dernier, ne peut pas augmenter la persistance du LCA seul ;
* top-down n’a pas cette garantie à cause de la non-commutation des projections ;
* le heat bath collapsed complet reste optimal. Des sweeps répétés convergent vers lui en volume fini.

Descendre ne suffit toutefois pas automatiquement : il faudra encore prouver une atténuation cumulée divergente ou un rayon spectral strictement inférieur à un sur la bande puis la grille.

La démonstration, les audits et le contre-exemple top-down exact sont dans 22_LCA_VS_FULL_HIERARCHY.md. Publication directe sur `main` : commit 07a6b1d. Les 151 tests passent et le dépôt est propre.

* * *

**You:**

Et quelle serait la meilleure stratégie pour montrer l'obstruction à la weak recovery ?

* * *

**ChatGPT:**

La meilleure stratégie est désormais :

$$\boxed{\text{corridor collapsed complet}  
+\text{criticalisation à squelette fixé}  
+\text{contraction annulaire multiscale}.}$$

Pour chaque paire $i,j$, on considère

$$A_{ij}=\|P_{ij}^{\downarrow}(\sigma_i\sigma_j)\|_2^2,  
\qquad  
Q_L\le \mathbb E[A_{I_LJ_L}].$$

Il suffit donc de prouver $\mathbb E[A_{I_LJ_L}]\to0$.

Les deux lemmes réellement décisifs seraient :

$$A_{ij}\le \kappa(p)^{N_L}+o(1),\qquad \kappa(p)<1,$$

où $N_L$ compte les blobs hiérarchiques screenés et ambigus traversés par le corridor, puis

$$\mathbb P(G_k\mid\mathcal F_{k-1},\mathrm{Palm}_{i\leftrightarrow j})  
\ge a(p)>0$$

sur $K_L\asymp\log L$ échelles. Ils donneraient

$$\mathbb E[A_{ij}\mid\mathrm{Palm}]  
\le  
\bigl(1-a(p)(1-\kappa(p))\bigr)^{K_L}+o(1)  
\longrightarrow0.$$

Le point subtil est la notion de « cas critique favorable » :

* À squelette et tailles de buckets fixés, remplacer les niveaux postcritiques par $\beta_c$ est rigoureusement favorable, par tensorisation de Blackwell.
* Remplacer tout le corridor tardif par une véritable géométrie critique n’est pas encore justifié : les tailles des buckets changent et l’ordre de Blackwell peut alors échouer.
* La voie la plus robuste consiste donc à criticaliser les canaux sur leur squelette réel, puis à prouver la contraction uniformément pour les corridors postcritiques.
* La domination géométrique par le cas $i,j$ reliés à $\beta_c$ reste un raccourci possible, mais nécessite un nouveau théorème.

Contre-audits essentiels :

* une arête pivotale isolée donne $m=1$, donc un canal parfait : il faut compter des blobs avec cycles ou coupes multi-arêtes, pas les pivots seuls ;
* une contraction locale ne suffit pas si les messages latéraux contournent le bloc ;
* RSW non conditionné ne suffit pas sous la loi Palm à deux points : il faut quasi-multiplicativité et séparation conditionnelle des interfaces, dans l’esprit de Garban–Pete–Schramm et des estimées proches-critiques de Nolin ;
* les exposants exacts de la percolation de sites triangulaire de Smirnov–Werner ne doivent pas être transférés sans preuve à notre percolation par arêtes.

À $p=0.8$, le cactus donne déjà $\kappa_{\rm conn}=0.8867525\ldots<1$. Le prochain objectif concret est donc le transfert répliqué exact d’une bande triangulaire de largeur deux, afin d’identifier le plus petit motif annulaire contractant. Ce ne sera pas encore la preuve bidimensionnelle, mais cela fixera précisément l’événement géométrique à chercher.

J’ai ajouté la synthèse complète dans , réorganisé le README, la feuille de route et la bibliographie, puis publié directement sur `main` : commit 50e8933.

Validation : 30 fichiers mathématiques contrôlés, 151 tests réussis, dépôt propre et synchronisé. Aucune impossibilité nouvelle à $p=0.8$ n’est encore annoncée.

* * *

**You:**

Explique-moi les idées de façon pédagogique et avec beaucoup plus de détails mathématiques

* * *

**ChatGPT:**

La bonne façon de comprendre l’idée est de séparer quatre objets :

1. la vérité plantée ;
2. une réplique tirée de la postérieure ;
3. le dendrogramme aléatoire construit par les horloges ;
4. le couplage hiérarchique utilisé pour mesurer ce qu’il reste d’information entre deux répliques.

La dynamique n’est pas seulement un algorithme : c’est surtout un outil pour construire un couplage postérieur dont la persistance donne une borne supérieure sur la weak recovery.

## 1. La quantité statistique qu’il faut annuler

Notons $X_i\in\{-1,+1\}$ la vérité et, pour une arête $e=\{i,j\}$,

$$Y_{ij}=X_iX_j Z_{ij},  
\qquad  
\mathbb P(Z_{ij}=+1)=p.$$

Une réplique postérieure est une configuration

$$\sigma\sim\mu_O,$$

où $O=(Y_e)_e$ représente les observations. On considère deux répliques indépendantes conditionnellement à $O$,

$$\sigma^{(1)},\sigma^{(2)}\overset{\mathrm{i.i.d.}}{\sim}\mu_O.$$

Leur overlap est

$$R_L(\sigma^{(1)},\sigma^{(2)})  
=  
\frac1{|V_L|}  
\sum_{i\in V_L}\sigma_i^{(1)}\sigma_i^{(2)}.$$

La quantité fondamentale est

$$Q_L  
=  
\mathbb E\left\langle  
R_L(\sigma^{(1)},\sigma^{(2)})^2  
\right\rangle.$$

En développant le carré,

$$Q_L  
=  
\frac1{|V_L|^2}  
\sum_{i,j}  
\mathbb E\left[  
\left\langle\sigma_i\sigma_j\right\rangle_{\mu_O}^2  
\right].  
\tag{1}$$

Autrement dit, si $I_L,J_L$ sont deux sommets uniformes indépendants,

$$Q_L  
=  
\mathbb E\left[  
c_{I_LJ_L}(O)^2  
\right],  
\qquad  
c_{ij}(O)  
=  
\mathbb E_{\mu_O}[\sigma_i\sigma_j].  
\tag{2}$$

Pourquoi est-ce la bonne quantité ?

* Si une weak recovery est possible, alors la postérieure conserve un overlap non nul avec la vérité.
* Par l’identité de Nishimori, cela implique qu’elle conserve aussi un overlap non nul entre deux répliques.
* Réciproquement, si $Q_L$ reste positif, tirer une réplique postérieure fournit déjà une procédure informationnelle de weak recovery.

Ainsi,

$$\boxed{  
Q_L\longrightarrow0  
\quad\Longrightarrow\quad  
\text{pas de weak recovery}.  
}  
\tag{3}$$

Notre but n’est donc pas seulement de montrer que la corrélation signée moyenne tend vers zéro. Il faut annuler son second moment.

## 2. Construction du dendrogramme hiérarchique

Une arête est satisfaite par $\sigma$ lorsque

$$Y_{ij}\sigma_i\sigma_j=+1.$$

Conditionnellement à $\sigma$, on lui associe l’horloge

$$\xi_e\sim\mathrm{Exp}(u_p)$$

si elle est satisfaite, et

$$\xi_e=+\infty$$

sinon, avec

$$u_p=\log\frac{p}{1-p}.$$

Au temps $t$, on ouvre les arêtes telles que $\xi_e\le t$. Les composantes fusionnent progressivement et définissent un dendrogramme $D$.

Un nœud interne

$$u:C_1\mathbin{\dot\cup}C_2=C$$

correspond à la fusion de deux composantes $C_1,C_2$, au temps $\beta_u$.

Il faut insister sur un point : le dendrogramme est non marqué. Il retient

* les deux composantes fusionnées ;
* le temps $\beta_u$ ;

mais pas l’identité de l’arête gagnante de Kruskal.

Par conséquent, au nœud $u$, toute la coupe

$$E_u  
=  
\bigl\{\{x,y\}\in E:x\in C_1,\ y\in C_2\bigr\}  
\tag{4}$$

doit être utilisée. On définit

$$\Lambda_u(\sigma)  
=  
\sum_{e\in E_u}  
|W_e|  
\mathbf 1_{\{e\text{ satisfaite par }\sigma\}}.  
\tag{5}$$

Conditionnellement à $\sigma$, le minimum des horloges satisfaites de cette coupe a pour densité

$$\Lambda_u(\sigma)  
e^{-\beta\Lambda_u(\sigma)}\,d\beta.$$

Après combinaison avec la densité postérieure, la conditionnelle exacte de $\sigma$ sachant le dendrogramme est

$$\boxed{  
\nu_O(\sigma\mid D)  
\propto  
\mu_0(\sigma)  
\prod_{v\in D}  
\Lambda_v(\sigma)  
e^{(1-\beta_v)\Lambda_v(\sigma)}.  
}  
\tag{6}$$

C’est cette formule qui explique pourquoi tous les $\Lambda_v$ ancestraux interviennent dans un flip interne.

## 3. Les quatre possibilités de flip à un nœud

Au nœud $u:C_1\dot\cup C_2$, notons $\sigma^{ab}$ la configuration obtenue en retournant $C_1$ si $a=1$ et $C_2$ si $b=1$.

Il y a quatre états :

$$(0,0),\quad(1,0),\quad(0,1),\quad(1,1).$$

Leurs poids exacts sont

$$q_u^{ab}  
=  
\mu_0(\sigma^{ab})  
\prod_{v\succeq u}  
\Lambda_v(\sigma^{ab})  
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.  
\tag{7}$$

Seuls $u$ et ses ancêtres apparaissent parce qu’un flip de $C_1$ ou $C_2$ ne change pas les relations à l’intérieur de ces clusters.

Pour une paire $i\in C_1,j\in C_2$,

$$\frac{\sigma_i^{ab}\sigma_j^{ab}}  
{\sigma_i\sigma_j}  
=  
(-1)^{a+b}.$$

Donc :

* $(0,0)$ et $(1,1)$ conservent la relation entre $i$ et $j$ ;
* $(1,0)$ et $(0,1)$ l’inversent.

La probabilité d’inversion est

$$\mathbb P_u(\mathrm{impair}\mid\sigma,D)  
=  
\frac{q_u^{10}+q_u^{01}}  
{q_u^{00}+q_u^{01}+q_u^{10}+q_u^{11}}.  
\tag{8}$$

Les log-odds de conservation sont

$$L_u  
=  
\log  
\frac{q_u^{00}+q_u^{11}}  
{q_u^{10}+q_u^{01}},  
\tag{9}$$

et la moyenne signée de la parité du flip vaut

$$\tanh\frac{L_u}{2}.$$

La fiabilité locale en second moment est donc

$$\eta_u  
=  
\tanh^2\frac{L_u}{2}.  
\tag{10}$$

## 4. Comment interviennent exactement les $\Lambda_v$ ancestraux ?

Introduisons la variable de satisfaction signée

$$S_e(\sigma)  
=  
\mathrm{sgn}(W_e)\sigma_x\sigma_y  
\in\{-1,+1\},$$

pour $e=\{x,y\}$. Alors

$$\mathbf 1_{\{e\text{ satisfaite}\}}  
=  
\frac{1+S_e(\sigma)}2.$$

Si $\chi_{ab}(x)$ vaut $-1$ lorsque $x$ appartient à un cluster retourné et $+1$ sinon, alors

$$S_e(\sigma^{ab})  
=  
S_e(\sigma)\chi_{ab}(x)\chi_{ab}(y).$$

D’où la formule explicite

$$\boxed{  
\Lambda_v(\sigma^{ab})  
=  
\sum_{e=\{x,y\}\in E_v}  
|W_e|  
\frac{  
1+  
S_e(\sigma)\chi_{ab}(x)\chi_{ab}(y)  
}{2}.  
}  
\tag{11}$$

Au nœud $u$, toutes les arêtes de $E_u$ ont un endpoint dans chaque enfant. Par conséquent,

$$\Lambda_u^{00}=\Lambda_u^{11}=\Lambda_u,$$

et, en notant

$$T_u=\sum_{e\in E_u}|W_e|,$$

on a

$$\Lambda_u^{10}  
=  
\Lambda_u^{01}  
=  
T_u-\Lambda_u.  
\tag{12}$$

Pour un ancêtre $v\succ u$, seule une partie des arêtes de $E_v$ voit exactement un endpoint retourné. On n’a donc plus la simple complémentation $\Lambda_v\mapsto T_v-\Lambda_v$.

C’est précisément la difficulté relevée dans les slides.

On peut isoler tous les ancêtres dans un message

$$B_u  
=  
\log  
\frac{  
A_u^{00}+A_u^{11}  
}{  
A_u^{10}+A_u^{01}  
},$$

où

$$A_u^{ab}  
=  
\mu_0(\sigma^{ab})  
\prod_{v\succ u}  
\Lambda_v(\sigma^{ab})  
e^{(1-\beta_v)\Lambda_v(\sigma^{ab})}.$$

Alors

$$\boxed{  
L_u  
=  
B_u  
+  
\log\frac{\Lambda_u}{T_u-\Lambda_u}  
+  
(1-\beta_u)(2\Lambda_u-T_u).  
}  
\tag{13}$$

Cette identité résume toute la difficulté :

* les deux derniers termes constituent l’information locale de la coupe ;
* $B_u$ est le message de tous les ancêtres.

### Cas de poids homogènes

Dans le GSBM triangulaire, toutes les arêtes ont le même poids $u_p$. Si $m_u=|E_u|$ et si $K_u$ arêtes de la coupe sont satisfaites, alors

$$\Lambda_u=u_pK_u,  
\qquad  
T_u=u_pm_u,$$

et

$$L_u  
=  
B_u  
+  
\log\frac{K_u}{m_u-K_u}  
+  
(1-\beta_u)u_p(2K_u-m_u).  
\tag{14}$$

Une majorité stricte $K_u>m_u/2$ donne bien une contribution locale positive. Mais elle ne suffit pas si $B_u$ est fortement négatif.

C’est pourquoi « la majorité des liens est conforme » est une bonne intuition locale, mais pas encore un théorème de weak recovery.

### Pourquoi les fausses arêtes internes ne comptent pas

Une arête interne à $C_1$ a ses deux endpoints retournés simultanément lors d’un flip de $C_1$. Sa satisfaction ne change donc pas. Elle se simplifie dans tous les rapports de poids.

Au niveau $u$, seules les arêtes de la coupe $E_u$ comptent. Les arêtes internes fausses peuvent influencer la construction antérieure des clusters, mais elles ne votent pas directement dans le flip relatif de $C_1$ contre $C_2$.

## 5. Pourquoi utiliser tout le corridor descendant ?

Soit

$$u=\mathrm{LCA}_D(i,j).$$

La relation $\sigma_i\sigma_j$ peut se décomposer comme le produit des orientations relatives rencontrées :

* sur le bras descendant de $i$ vers $u$ ;
* au nœud $u$ ;
* sur le bras descendant de $u$ vers $j$.

Schématiquement :

```
flowchart TB
    A["Ancêtres : messages Λᵥ"] --> U["LCA u, niveau βc"]
    U --> CI["Cluster contenant i"]
    U --> CJ["Cluster contenant j"]
    CI --> BI["Blocs descendants"]
    CJ --> BJ["Blocs descendants"]
    BI --> I["Feuille i"]
    BJ --> J["Feuille j"]
```

Mettre à jour uniquement le LCA ne touche qu’une orientation relative. Toute l’information stockée dans les orientations descendantes reste figée.

Soit $P_u$ le heat bath du LCA. Comme tout heat bath,

$$P_u g  
=  
\mathbb E_{\pi_D}[g\mid\mathcal G_u]$$

est une projection orthogonale dans $L^2(\pi_D)$.

Pour tout le corridor descendant, posons

$$\mathcal A_{ij}^{\downarrow}  
=  
\bigcap_{v\in\mathcal C_{ij}^{\downarrow}}\mathcal G_v,$$

et

$$P_{ij}^{\downarrow}g  
=  
\mathbb E_{\pi_D}[g\mid\mathcal A_{ij}^{\downarrow}].  
\tag{15}$$

Comme

$$\mathcal A_{ij}^{\downarrow}\subseteq\mathcal G_u,$$

on a l’identité de Pythagore

$$\boxed{  
\|P_uf_{ij}\|_2^2  
=  
\|P_{ij}^{\downarrow}f_{ij}\|_2^2  
+  
\|(P_u-P_{ij}^{\downarrow})f_{ij}\|_2^2.  
}  
\tag{16}$$

Donc

$$\boxed{  
\|P_{ij}^{\downarrow}f_{ij}\|_2^2  
\le  
\|P_uf_{ij}\|_2^2.  
}  
\tag{17}$$

Le corridor complet est toujours au moins aussi contractant que le LCA seul.

Plus généralement, si

$$K=P_{v_1}\cdots P_{v_M}$$

est un sweep utilisant les mêmes nœuds, alors

$$\boxed{  
\|Kf\|_2^2  
=  
\|P_{\mathcal C}f\|_2^2  
+  
\|K(I-P_{\mathcal C})f\|_2^2.  
}  
\tag{18}$$

Ainsi, le heat bath collapsed de tout le corridor est la dynamique la plus favorable pour démontrer une obstruction parmi les programmes utilisant ces coordonnées.

Cela ne signifie pas qu’il s’agit du meilleur algorithme pratique de MCMC. C’est un ordre portant sur la persistance $L^2$ d’une paire.

## 6. Pourquoi cette persistance borne-t-elle la weak recovery ?

Posons

$$f_{ij}(\sigma)=\sigma_i\sigma_j$$

et

$$A_{ij}(O,D)  
=  
\|P_{ij}^{\downarrow}f_{ij}\|_{L^2(\pi_D)}^2.  
\tag{19}$$

Cette quantité admet une interprétation à deux répliques. Conditionnellement à ce qui est conservé hors du corridor, tirons deux orientations indépendantes du corridor, $\sigma^{(1)}$ et $\sigma^{(2)}$. Alors

$$A_{ij}  
=  
\mathbb E\left[  
f_{ij}(\sigma^{(1)})  
f_{ij}(\sigma^{(2)})  
\mid O,D  
\right],$$

après moyenne sur l’état extérieur.

Autrement dit, $A_{ij}$ mesure ce qui reste prévisible de la relation $i,j$ lorsqu’on rééchantillonne exactement tout son corridor.

Pour relier cela à $Q_L$, écrivons

$$c_{ij}(O)  
=  
\mathbb E_{\mu_O}[f_{ij}].$$

La désintégration par $D$ donne

$$c_{ij}(O)  
=  
\mathbb E_{D\mid O}  
\left[  
\mathbb E_{\pi_D}[f_{ij}]  
\right].$$

Par Jensen,

$$c_{ij}(O)^2  
\le  
\mathbb E_{D\mid O}  
\left[  
\mathbb E_{\pi_D}[f_{ij}]^2  
\right].$$

Or

$$\mathbb E_{\pi_D}[f_{ij}]  
=  
\mathbb E_{\pi_D}[P_{ij}^{\downarrow}f_{ij}],$$

donc, encore par Jensen,

$$\mathbb E_{\pi_D}[f_{ij}]^2  
\le  
\mathbb E_{\pi_D}  
\left[  
(P_{ij}^{\downarrow}f_{ij})^2  
\right]  
=  
A_{ij}.$$

Ainsi,

$$\boxed{  
c_{ij}(O)^2  
\le  
\mathbb E[A_{ij}\mid O].  
}$$

Après moyenne sur une paire uniforme,

$$\boxed{  
Q_L  
\le  
\mathbb E[A_{I_LJ_L}].  
}  
\tag{20}$$

C’est le théorème maître : pour interdire la weak recovery, il suffit de construire, séparément pour chaque paire, un corridor collapsed dont la persistance tend vers zéro.

## 7. La coordonnée de percolation

Sous la loi annealed de Nishimori, une arête est satisfaite avec probabilité $p$. Elle est ouverte avant le temps $t$ avec probabilité

$$q_p(t)  
=  
p(1-e^{-u_pt}).  
\tag{21}$$

Sur la grille triangulaire, le seuil de percolation par arêtes est

$$q_c  
=  
2\sin\frac{\pi}{18}  
=  
0.3472963553\ldots$$

Le temps critique est donc défini par

$$q_p(\beta_c)=q_c,$$

c’est-à-dire

$$\boxed{  
\beta_c(p)  
=  
-\frac1{u_p}  
\log\left(1-\frac{q_c}{p}\right).  
}  
\tag{22}$$

À $p=0.8$,

$$u_p=\log4=1.38629436\ldots,$$

et

$$\beta_c(0.8)  
=  
0.4107165392\ldots$$

La probabilité d’ouverture avant le temps final $1$ est

$$q_p(1)=2p-1.$$

À $p=0.8$,

$$q_p(1)=0.6>q_c.$$

## 8. La qualité exacte des liens non encore activés

C’est ici que votre intuition sur les liens faux et les horloges après $\beta_c$ devient exacte.

À un temps $t$, une arête non encore ouverte peut être :

1. fausse/non satisfaite, auquel cas son horloge vaut $+\infty$ ;
2. satisfaite, mais avec une horloge supérieure à $t$.

Sa probabilité conditionnelle d’être satisfaite est

$$s_p(t)  
=  
\mathbb P(  
e\text{ satisfaite}  
\mid  
\xi_e>t  
).$$

Or

$$\mathbb P(e\text{ satisfaite},\xi_e>t)  
=  
p e^{-u_pt}  
=  
p-q_p(t),$$

et

$$\mathbb P(\xi_e>t)=1-q_p(t).$$

Donc

$$\boxed{  
s_p(t)  
=  
\frac{p-q_p(t)}{1-q_p(t)}.  
}  
\tag{23}$$

Au seuil,

$$\boxed{  
s_c(p)  
=  
s_p(\beta_c)  
=  
\frac{p-q_c}{1-q_c}.  
}  
\tag{24}$$

Pour $p=0.8$,

$$s_c(0.8)  
=  
0.6935822228\ldots$$

Cette valeur est inférieure à $p=0.8$ parce que le conditionnement « l’arête ne s’est pas encore ouverte » enrichit le pool en arêtes fausses.

De plus,

$$s_p(0)=p,  
\qquad  
s_p(1)=\frac12,$$

et $s_p(t)$ décroît avec $t$. Plus une fusion est tardive, moins les arêtes résiduelles distinguent les deux orientations.

C’est la première justification mathématique de l’idée :

> à taille et géométrie fixées, une fusion au seuil critique fournit des liens résiduels de meilleure qualité qu’une fusion postcritique.

### Entre $\beta_c$ et $1$

Le pool restant à $\beta_c$ est exactement l’union de :

* toutes les arêtes fausses, d’horloge $+\infty$ ;
* les arêtes satisfaites dont l’horloge dépasse $\beta_c$.

Parmi ces dernières :

* celles avec $\beta_c<\xi_e\le1$ pourront participer à une fusion ultérieure ;
* celles avec $\xi_e>1$ resteront censurées, comme les arêtes fausses.

Il faut ensuite conditionner par le fait qu’une fusion a effectivement lieu à un certain temps. Cette sélection Palm garantit au moins une arête gagnante satisfaite : on ne peut donc pas utiliser une simple loi binomiale indépendante.

## 9. Le canal exact d’un bucket

Supposons qu’une coupe possède $m$ arêtes candidates. Soit $X\in\{-1,+1\}$ la parité latente que l’on cherche à déterminer et $K$ le nombre d’arêtes favorisant l’orientation de référence.

Avec

$$s=s_p(t),$$

la loi exacte sous le biais de fusion est

$$\boxed{  
K\mid X=+1  
\sim  
1+\mathrm{Bin}(m-1,s),  
}  
\tag{25}$$

tandis que

$$\boxed{  
K\mid X=-1  
\sim  
\mathrm{Bin}(m-1,1-s).  
}  
\tag{26}$$

Le décalage $+1$ vient de l’arête gagnante latente : le dendrogramme révèle qu’une arête compatible avec la fusion a sonné, mais pas laquelle.

### Le cas $m=1$

Si $m=1$,

$$K=  
\begin{cases}  
1,&X=+1,\\  
0,&X=-1.  
\end{cases}$$

Le canal est parfait. C’est un point crucial :

$$\boxed{  
\text{une arête pivotale isolée ne contracte aucune information.}  
}$$

Il ne faut donc pas compter les pivots ou les arêtes rouges comme des blocs de décorrélation. Il faut trouver des blobs où plusieurs routes ou plusieurs arêtes candidates créent une vraie ambiguïté.

## 10. En quel sens le seuil critique est-il favorable ?

Il y a deux affirmations différentes.

### 10.1 Domination des canaux : établie

À taille $m$ fixée, si $t_1\le t_2$, alors

$$s_p(t_1)\ge s_p(t_2).$$

Le canal $\mathcal E_{m,s_p(t_1)}$ Blackwell-domine le canal $\mathcal E_{m,s_p(t_2)}$. Il existe donc un noyau stochastique $G$ tel que

$$P_{t_2}^{x}  
=  
P_{t_1}^{x}G,  
\qquad x\in\{-1,+1\}.  
\tag{27}$$

Le compte tardif peut être simulé en ajoutant du bruit au compte précoce.

Sur un corridor fixé, les noyaux se tensorisent. Pour toute fonction $F$ des parités, même si les parités ont un prior corrélé,

$$\mathbb E\left[  
\mathbb E(F\mid K^{\mathrm{late}})^2  
\right]  
\le  
\mathbb E\left[  
\mathbb E(F\mid K^{\mathrm{crit}})^2  
\right].  
\tag{28}$$

On peut donc remplacer, sur le même corridor,

$$t_r  
\longmapsto  
t_r^{\mathrm{fav}}  
=  
\min(t_r,\beta_c).  
\tag{29}$$

Cela donne un oracle plus informatif et donc une borne supérieure sur la persistance réelle.

### 10.2 Domination des géométries : ouverte

Une vraie fusion plus tardive ne possède généralement pas :

* les mêmes clusters ;
* les mêmes coupes ;
* les mêmes tailles de buckets ;
* le même nombre de niveaux descendants.

Or la domination de Blackwell n’est pas totale lorsque $m$ change. Le dépôt contient même un contre-exemple exact entre un bucket critique $m=4$ et un bucket tardif $m=2$.

On ne peut donc pas encore écrire sans preuve

$$A_L^{\mathrm{late}}  
\le  
A_L^{\mathrm{critical\ geometry}}.$$

La formulation rigoureuse est plutôt :

1. garder le squelette tardif réel ;
2. améliorer artificiellement ses canaux jusqu’au niveau critique ;
3. prouver une contraction uniforme sur ces squelettes réels.

C’est plus robuste qu’un couplage global entre deux dendrogrammes différents.

## 11. Pourquoi conditionner $i,j$ à être reliés à $\beta_c$ ?

Pour une paire lointaine, décomposons selon son rang de fusion $T_{ij}=q_p(\beta_{ij})$ :

$$\begin{aligned}  
\mathbb E[A_{ij}]  
={}&  
\mathbb E[A;T<q_c-\delta]  
+  
\mathbb E[A;|T-q_c|\le\delta]  
\\  
&+  
\mathbb E[A;T>q_c+\delta,\ T\le q_p(1)]  
+  
\mathbb E[A;T>q_p(1)].  
\end{aligned}  
\tag{30}$$

Les deux extrêmes sont contrôlés.

### Fusion sous-critique précoce

Pour une paire à distance macroscopique,

$$\mathbb P(T<q_c-\delta)\longrightarrow0$$

par décroissance sous-critique. Comme $A\le1$, cette contribution est $o(1)$.

### Racines distinctes

Si $T>q_p(1)$, les deux sommets appartiennent à deux composantes finales différentes. Leurs orientations globales sont rééchantillonnées indépendamment avec probabilité $1/2$. Donc

$$A_{ij}=0.  
\tag{31}$$

Il reste les fusions critiques et postcritiques.

Le scénario où $i,j$ sont déjà connectés dans une composante macroscopique au seuil est intuitivement le plus favorable à la conservation :

* ils sont reliés aussi tôt que possible sans événement sous-critique rare ;
* les canaux résiduels ont leur qualité maximale $s_c(p)$ ;
* les composantes critiques peuvent offrir plusieurs routes et messages ancestraux cohérents.

Mais ce conditionnement est rare pour deux sommets uniformes : au seuil, il n’existe pas de géante de densité positive. Il doit donc servir d’oracle majorant, pas remplacer silencieusement la moyenne globale.

Enfin, l’événement exact

$$T_{ij}=q_c$$

a probabilité nulle en volume fini. Il faut l’interpréter comme :

* un conditionnement $T_{ij}\in[q_c-\delta_L,q_c+\delta_L]$ ;
* ou une espérance de Palm par rapport à la densité du flux de fusions.

## 12. Pourquoi la distance peut aider : les blocs multiscales

Le LCA seul ne gagne pas automatiquement avec la distance. Sur le cactus triangulaire du dépôt,

$$A_h^{\mathrm{LCA\ only}}  
=  
\kappa_{\mathrm{flux}}(p,q),$$

indépendamment de la longueur $h$.

En revanche, le corridor complet donne

$$A_h^{\mathrm{full}}  
=  
\kappa_{\mathrm{flux}}(p,q)  
\kappa_{\mathrm{conn}}(p,q)^{h-1}.  
\tag{32}$$

À $p=0.8,q=q_c$,

$$\kappa_{\mathrm{conn}}  
=  
0.8867525669\ldots<1.$$

La corrélation disparaît donc exponentiellement sur ce modèle.

Sur la grille bidimensionnelle, la factorisation (32) est fausse en général, car :

* les cycles se chevauchent ;
* plusieurs routes contournent un bloc ;
* les messages ancestraux sont corrélés ;
* le nombre de ports d’un bloc peut grandir.

Il faut remplacer le produit scalaire par un transfert d’états de bord.

## 13. Le transfert répliqué

Pourquoi deux répliques ?

Parce que

$$A_{ij}  
=  
\left\|  
\mathbb E[f_{ij}\mid\mathcal A]  
\right\|_2^2$$

est un second moment. On le représente par deux rééchantillonnages indépendants des orientations internes, partageant :

* les mêmes observations ;
* le même dendrogramme ;
* le même état extérieur.

L’état de bord minimal d’un bloc doit contenir au moins

$$Z_k=  
\left(  
\Pi_k,  
R_k^{(1)},  
R_k^{(2)},  
P_k  
\right),  
\tag{33}$$

où :

* $\Pi_k$ est la partition de connectivité des ports ;
* $R_k^{(1)},R_k^{(2)}$ sont les parités dans les deux répliques ;
* $P_k$ enregistre le statut de fusion ou de pivotalité.

Le bloc définit un opérateur de transfert

$$\mathscr U_{p,k}$$

sur ces états.

La contraction pertinente est, schématiquement,

$$\eta_k  
=  
\sup_{g\perp1}  
\frac{\|\mathscr U_{p,k}g\|_2^2}  
{\|g\|_2^2}.  
\tag{34}$$

On veut trouver des blocs pour lesquels

$$\eta_k\le\kappa(p)<1.  
\tag{35}$$

Le secteur $g\perp1$ est le secteur informatif, notamment celui qui porte la parité entre les deux extrémités.

## 14. Qu’est-ce qu’un bon bloc annulaire ?

Pour $d(i,j)\asymp L$, prenons des annuli de rayons géométriques

$$\mathcal A_k  
=  
B(i,\lambda^{k+1})  
\setminus  
B(i,\lambda^k),  
\qquad  
1\le k\le K_L,$$

avec

$$K_L\asymp\log L.$$

Un bon bloc doit satisfaire quatre propriétés.

### 1. Nombre borné de ports

Seul un nombre borné d’interfaces macroscopiques doit traverser le bloc. Sinon l’état de bord devient trop grand.

### 2. Écrantage latéral

Les routes extérieures ne doivent pas pouvoir contourner arbitrairement le canal étudié. Des interfaces duales doivent canaliser la connexion vers les ports retenus.

### 3. Ambiguïté réelle

Le bloc doit contenir au moins deux routes ou arêtes candidates contribuant à une coupe hiérarchique. Un simple pivot $m=1$ ne convient pas.

### 4. Message extérieur contrôlé

Le message ancestral

$$B_k$$

ne doit pas être arbitrairement grand. Pour un bucket critique $m=2$, une borne $|B_k|\le B_0$ donne le coefficient local

$$\boxed{  
\kappa_2(B_0;p)  
=  
s_c(p)  
+  
(1-s_c(p))  
\tanh^2(B_0/2)  
<1.  
}  
\tag{36}$$

Cette formule explique exactement le rôle du screening :

* pour $B_0=0$, la contraction vaut $s_c(p)$ ;
* lorsque $B_0\to\infty$, elle tend vers $1$.

Ainsi, même un bon canal local ne décorrèle rien si l’extérieur connaît déjà presque parfaitement la parité.

## 15. Les deux lemmes qui termineraient la preuve

### Lemme de canal

Si $N_L$ est le nombre de bons blocs rencontrés par le corridor, montrer

$$\boxed{  
A_{ij}  
\le  
\kappa(p)^{N_L}  
+  
\varepsilon_L^{\mathrm{tr}},  
}  
\tag{37}$$

avec

$$\mathbb E[\varepsilon_L^{\mathrm{tr}}]\longrightarrow0.$$

L’erreur couvre, par exemple, la troncature d’un nombre exceptionnellement grand de ports.

Cette inégalité doit être conditionnelle aux géométries et états de bord. On n’a alors pas besoin d’indépendance entre les blocs : les normes d’opérateurs sont sous-multiplicatives.

### Lemme géométrique sous Palm

Pour des annuli suffisamment séparés, montrer

$$\boxed{  
\mathbb P(  
G_k  
\mid  
\mathcal F_{k-1},  
\mathrm{Palm}_{i\leftrightarrow j}  
)  
\ge a(p)>0.  
}  
\tag{38}$$

C’est beaucoup plus fort qu’une simple probabilité RSW non conditionnée. La borne doit rester vraie après révélation des autres échelles et sous le biais imposant la connexion lointaine.

Sous (38), il y a typiquement $c\log L$ bons blocs. Plus directement, à chaque échelle le facteur conditionnel moyen est au plus

$$(1-a(p))\times1+a(p)\times\kappa(p)  
=  
1-a(p)(1-\kappa(p)).$$

Ainsi,

$$\boxed{  
\mathbb E[A_{ij}\mid\mathrm{Palm}]  
\le  
\left(  
1-a(p)(1-\kappa(p))  
\right)^{K_L}  
+o(1).  
}  
\tag{39}$$

Comme $K_L\asymp\log_\lambda L$,

$$\left(  
1-a(1-\kappa)  
\right)^{K_L}  
=  
L^{-\gamma+o(1)},$$

où

$$\gamma  
=  
-\frac{  
\log(1-a(1-\kappa))  
}{  
\log\lambda  
}  
>0.  
\tag{40}$$

On obtient alors une disparition polynomiale de la corrélation.

## 16. Le critère spectral plus général

La formulation par bons blocs est un critère suffisant simple. Le critère potentiellement optimal utilise le produit complet des transferts :

$$\mathscr U_{p,K}\cdots\mathscr U_{p,1}.$$

On cherche un exposant de Lyapunov négatif,

$$\lambda_{\mathrm{hier}}(p)  
=  
\limsup_{K\to\infty}  
\frac1K  
\mathbb E  
\log  
\left\|  
\mathscr U_{p,K}\cdots\mathscr U_{p,1}  
\right\|  
<0.  
\tag{41}$$

Dans un modèle fini et stationnaire, cela devient souvent une condition de rayon spectral :

$$r(\mathscr U_p)<1.  
\tag{42}$$

C’est pourquoi la bande triangulaire de largeur deux est le prochain calcul naturel. Son état de bord est fini, mais elle possède déjà :

* des cycles chevauchants ;
* plusieurs ports ;
* des messages latéraux ;
* un transfert véritablement répliqué.

Si l’on certifie

$$r(\mathscr U_{0.8,2})<1,$$

on obtient un vrai certificat de canal. Il restera ensuite à prouver que les configurations correspondantes apparaissent avec une densité d’échelles positive dans la géométrie bidimensionnelle critique.

## 17. Pourquoi $p=0.8$ est un bon premier objectif

La borne d’information-percolation arête par arête utilise le coefficient

$$\eta_{\mathrm{edge}}  
=  
(2p-1)^2.$$

Elle interdit la weak recovery lorsque

$$(2p-1)^2<q_c.$$

Donc

$$p<  
\frac{1+\sqrt{q_c}}2  
=  
0.794659\ldots  
\tag{43}$$

À $p=0.8$,

$$(2p-1)^2=0.36,$$

alors que

$$q_c=0.347296\ldots$$

La borne arête par arête échoue de peu. Le but de la hiérarchie est précisément d’exploiter :

* les coupes multi-arêtes ;
* le fait que l’arête gagnante est marginalisée ;
* les cycles et blobs ambigus ;
* plusieurs niveaux descendants ;

pour obtenir une contraction de bloc plus forte que la contraction arête par arête.

Le cactus montre que cette amélioration est possible dans une géométrie contrôlée. La bande de largeur deux doit montrer qu’elle survit aux premiers cycles chevauchants. Le lemme annulaire doit ensuite la transporter à la grille.

## 18. Résumé de la stratégie de preuve

La chaîne logique complète serait :

$$\text{weak recovery}  
\Longrightarrow  
Q_L\nrightarrow0,$$

mais

$$Q_L  
\le  
\mathbb E[A_{I_LJ_L}].$$

On décompose les paires :

* fusion sous-critique lointaine : probabilité $o(1)$ ;
* racines distinctes : persistance exactement nulle ;
* fusion critique/postcritique : partie difficile.

Sur chaque corridor difficile :

1. on conserve son squelette et ses tailles ;
2. on remplace les canaux tardifs par leurs versions critiques, plus informatives ;
3. on applique le heat bath collapsed de tout le corridor descendant ;
4. on décompose le corridor en blocs annulaires ;
5. on montre une contraction spectrale sur les bons blocs ;
6. on montre qu’il existe $c\log L$ bons blocs sous la loi Palm.

Alors

$$\mathbb E[A_{I_LJ_L}]  
\longrightarrow0,$$

donc

$$Q_L\longrightarrow0,$$

et la weak recovery est impossible.

## Ce qui est établi et ce qui reste ouvert

Établi :

* la loi jointe exacte avec tous les $\Lambda_v$ ;
* les quatre probabilités de flip ;
* la réduction pairwise vers $Q_L$ ;
* l’optimalité $L^2$ du corridor collapsed ;
* la domination de Blackwell à squelette et tailles fixés ;
* la contraction exacte sur le cactus ;
* l’annulation pour deux racines distinctes.

Ouvert :

* le transfert exact de largeur deux ;
* la compression correcte de l’état de bord ;
* le screening uniforme des messages ancestraux ;
* l’abondance de blobs ambigus sous la loi Palm à deux points ;
* la comparaison géométrique critique/postcritique, si la version uniforme ne suffit pas.

La synthèse mathématique vérifiée sur `main` est , au commit 50e8933.

* * *

**You:**

Prenons un peu de recul et essayons de voir si des choses plus simples peuvent marcher. Mon intuition est la suivante, au début il y a plus d'arêtes avec vraies information que fausse. D'où le fait qu'à faible beta de petits clusters informatifs apparaissent. Mais à partir d'un certain beta, entre deux clusters, il commence à y avoir trop de fausses arêtes par rapport aux arêtes non encore sonnées. À partir de ce moment, il n'y a plus d'information, donc plus de corrélation, donc plus de weak recovery. Le cas le plus favorable étant de supposer que deux points i et j lointains apparaissent dans un même cluster le plus tôt possible (au beta_c de percolation).
Peux-tu formaliser ces idées. Et montrer ces résultats. Mets à jour le GitHub avec la présentation pédagogique qe tu as faite et dans ce sens.

* * *

**ChatGPT:**

L’intuition contient un mécanisme réel, mais le « moment où les fausses arêtes deviennent majoritaires » doit être reformulé. Le bilan exact montre qu’il n’y a pas de basculement brutal avant $\beta=1$. La voie simple la plus prometteuse est plutôt l’accumulation de buckets ambigus de taille bornée.

## 1. Bilan exact des arêtes au temps $\beta$

Dans la jauge de Nishimori, une arête vraie reçoit une horloge $\mathrm{Exp}(u_p)$, tandis qu’une fausse arête reçoit $+\infty$, avec

$$u_p=\log\frac p{1-p}.$$

Au temps $\beta$, les masses sont exactement :

| catégorie | masse |
| --- | --- |
| vraie déjà sonnée | $p(1-e^{-u_p\beta})$ |
| vraie qui sonnera dans $(\beta,1]$ | $p(e^{-u_p\beta}-e^{-u_p})$ |
| vraie censurée après $1$ | $pe^{-u_p}=1-p$ |
| fausse | $1-p$ |

La masse de toutes les vraies arêtes non encore sonnées vaut

$$pe^{-u_p\beta}.$$

Comme $1-p=pe^{-u_p}$,

$$pe^{-u_p\beta}-(1-p)  
=  
p(e^{-u_p\beta}-e^{-u_p})  
>0$$

pour tout $\beta<1$. Donc :

$$\boxed{  
\text{les fausses arêtes ne dépassent jamais toutes les vraies non activées avant }\beta=1.  
}$$

L’égalité n’a lieu qu’au bord terminal.

## 2. Le croisement plus précoce existe, mais ce n’est pas un seuil d’information

Si l’on compare seulement les vraies arêtes qui sonneront encore avant $1$ aux fausses, le croisement satisfait

$$p(e^{-u_p\beta}-e^{-u_p})=1-p.$$

Il a lieu à

$$\boxed{  
\beta_{\rm act}(p)  
=  
1-\frac{\log2}{u_p}.  
}$$

À $p=0.8$,

$$\beta_{\rm act}=0.5.$$

Mais à ce temps, conditionnellement au fait qu’une arête soit encore fermée, les catégories

$$(\text{vraie future},\text{vraie censurée},\text{fausse})$$

ont chacune probabilité $1/3$. L’arête est donc encore vraie avec probabilité $2/3$. Ce croisement n’est pas une perte totale d’information.

Au seuil critique triangulaire, ce diagnostic donne

$$p_{\partial,\rm act}  
=  
\frac{2+q_c}{3}  
=  
0.7824321184\ldots,$$

strictement sous la borne information-percolation

$$p_{\rm info}  
=  
\frac{1+\sqrt{q_c}}2  
=  
0.7946592758\ldots.$$

La comparaison naïve « futures vraies contre fausses » ne peut donc pas, à elle seule, battre les résultats connus.

## 3. La vraie qualité résiduelle

Conditionnellement au fait qu’une arête ne soit pas ouverte avant $\beta$, sa probabilité d’être vraie est

$$\boxed{  
s_p(\beta)  
=  
\frac{pe^{-u_p\beta}}  
{1-p+pe^{-u_p\beta}}  
=  
\frac1{1+e^{-u_p(1-\beta)}}.  
}$$

Son biais signé est

$$h_p(\beta)  
=  
2s_p(\beta)-1  
=  
\tanh\left(\frac{u_p(1-\beta)}2\right).$$

On a

$$s_p(0)=p,\qquad s_p(1)=\frac12,  
\qquad s_p'(\beta)<0.$$

La qualité décroît donc continûment, sans transition locale abrupte.

À $p=0.8$,

$$\beta_c=0.4107165392\ldots<\beta_{\rm act}=0.5,$$

et

$$s_c(0.8)  
=  
s_p(\beta_c)  
=  
0.6935822228\ldots.$$

Sur une frontière critique, les proportions conditionnelles sont environ :

$$\begin{aligned}  
\text{vraie future}&=0.3871644,\\  
\text{vraie censurée}&=0.3064178,\\  
\text{fausse}&=0.3064178.  
\end{aligned}$$

Le LCA critique est donc fortement informatif. Cela confirme qu’il constitue bien le cas favorable le plus sévère au niveau du canal.

## 4. Deux mécanismes simples de perte d’information

### Grande coupe très tardive

Pour un bucket de taille $m$, la fiabilité locale oracle est

$$\Gamma_m(\beta;p)  
=  
\mathbb E\left[  
\tanh^2\left(\frac{\ell_{m,K}(\beta;p)}2\right)  
\right].$$

Pour tout $\beta<1$ fixé,

$$\Gamma_m(\beta;p)\longrightarrow1  
\qquad(m\to\infty).$$

Une grande coupe critique est donc extrêmement informative, et non bruitée.

En revanche, au niveau terminal,

$$\boxed{  
\Gamma_m(1;p)=\frac1m.  
}$$

Ainsi une grande coupe fusionnant très près de $1$, avec message extérieur neutre, devient asymptotiquement équitable. La fenêtre exacte est

$$1-\beta\asymp m^{-1/2}.$$

Le paramètre local pertinent est

$$m\,h_p(\beta)^2,$$

plutôt que le seul bilan moyen vrai/faux.

### Buckets ambigus de taille deux

Pour $m=2$, le canal est exactement

$$\begin{array}{c|ccc}  
&K=0&K=1&K=2\\  
\hline  
X=+1&0&1-s&s\\  
X=-1&s&1-s&0.  
\end{array}$$

C’est un canal d’effacement :

* avec probabilité $s=s_p(\beta)$, la parité est révélée ;
* avec probabilité $1-s$, le compte $K=1$ n’apporte aucune information locale.

Sa fiabilité exacte est

$$\boxed{  
\Gamma_2(\beta;p)=s_p(\beta)\le p<1.  
}$$

Avec un message ancestral $B$,

$$\boxed{  
\kappa_2(B;\beta,p)  
=  
s_p(\beta)  
+  
(1-s_p(\beta))\tanh^2(B/2).  
}$$

Si $|B|\le B_0<\infty$,

$$\kappa_2(B;\beta,p)  
\le  
p+(1-p)\tanh^2(B_0/2)  
=:\bar\kappa<1.$$

## 5. Le nouveau théorème simple à viser

Supposons que, sous la loi Palm favorable où $i,j$ fusionnent dans la fenêtre critique, leur corridor descendant contienne $N_L$ blocs disjoints tels que :

* chaque bucket satisfait $2\le m\le M$ ;
* les routes latérales sont écrantées ;
* les messages ancestraux sont bornés ;
* chaque transfert répliqué contracte d’un facteur au plus $\kappa<1$ ;
* $N_L\to\infty$ en probabilité.

Alors

$$A_{ij}  
\le  
\kappa^{N_L}+o(1),$$

et donc

$$\mathbb E[A_{ij}]\longrightarrow0.$$

Après réduction favorable des paires postcritiques,

$$Q_L\le\mathbb E[A_{I_LJ_L}]\longrightarrow0,$$

ce qui interdit la weak recovery.

Dans le corridor factorisé composé uniquement de buckets $m=2$, cette formule est déjà exacte :

$$A_N  
=  
\prod_{r=1}^N s_p(\beta_r)  
\le p^N.$$

À $p=0.8$, pour dix buckets critiques neutres,

$$s_c^{10}=0.0257619974\ldots,$$

et pour quarante,

$$s_c^{40}=4.4047182\times10^{-7}.$$

## 6. Le verrou désormais isolé

La meilleure voie simple est donc :

$$\boxed{  
\text{montrer sous Palm critique que }  
N_L\to\infty  
\text{ buckets bornés et écrantés apparaissent.}  
}$$

Deux difficultés restent :

1. contrôler le message ancestral $B_u$ ;
2. empêcher les cycles et routes latérales de contourner la contraction.

Si ce screening simple fonctionne, le transfert complet de bande de largeur deux devient inutile. S’il échoue, la bande reste le plan B pour conserver tout l’état de bord.

La domination géométrique du cas critique sur tous les corridors tardifs reste également ouverte. Ce qui est déjà prouvé est la domination des canaux à squelette et tailles fixés.

## Mise à jour du dépôt

J’ai ajouté la présentation complète et les preuves dans .

J’ai également :

* réordonné la feuille de route pour tester d’abord les buckets bornés ;
* conservé le transfert de bande comme plan B ;
* ajouté les fonctions exactes pour $\beta_{\rm act}$ et les deux bilans de masse ;
* ajouté deux contre-tests indépendants ;
* mis à jour le README et la documentation des calculs.

Tout est publié directement sur `main` au commit 882221d.

Validation : 31 fichiers mathématiques contrôlés, 153 tests réussis, dépôt propre et synchronisé. Aucune nouvelle borne globale à $p=0.8$ n’est encore revendiquée.

* * *

**You:**

Non, je pense que tu te fourvoies et ne prends pas en compte la géométrie. À beta fixé, les arêtes vraies non encore sonnées sont peut-être majoritaires, mais à l'intérieur des clusters déjà formés. Ce qui compte est au niveau des coupes entre deux clusters. Là j'e pense qu'à partir d'un certain moment, ces coupes cessent d'être informatives.

* * *

**ChatGPT:**

Tu avais raison de rejeter le bilan global : le bon objet est la coupe entre deux clusters. Le calcul conditionnel apporte toutefois une nuance décisive.

Conditionnellement à la partition complète $\Pi_\beta=\pi$, les arêtes d’une coupe $E(C_1,C_2)$ sont indépendantes et vérifient

$$s_p(\beta)  
=\mathbb P(\text{arête conforme}\mid \text{fermée à }\beta)  
=\frac{pe^{-u_p\beta}}{1-p+pe^{-u_p\beta}},$$

avec biais

$$h_p(\beta)=2s_p(\beta)-1  
=\tanh\!\left(\frac{u_p(1-\beta)}2\right).$$

Les vraies arêtes déjà activées sont bien internes aux clusters, mais cela ne crée pas d’appauvrissement supplémentaire des arêtes fermées de frontière. L’hypothèse « $i,j$ appartiennent à la composante géante » modifie donc la géométrie des coupes — tailles, formes, ports latéraux — mais pas leur paramètre unitaire $s_p(\beta)$, une fois la partition et la taille fixées.

Pour une coupe instantanée de taille $m$,

$$K\mid X=+1\sim\mathrm{Bin}(m,s_p(\beta)),\qquad  
K\mid X=-1\sim\mathrm{Bin}(m,1-s_p(\beta)),$$

et

$$L^{\mathrm{snap}}=u_p(1-\beta)(2K-m).$$

La charge géométrique pertinente est

$$\boxed{\mathcal J=m h_p(\beta)^2}.$$

Asymptotiquement :

$$\mathcal J\to0\Rightarrow\text{information}\to0,\qquad  
\mathcal J\to\infty\Rightarrow\text{information}\to1.$$

Près de la censure,

$$h_p(\beta)\sim \frac{u_p}{2}(1-\beta),$$

donc la fenêtre de perte est

$$1-\beta\asymp\frac{2}{u_p\sqrt m}.$$

Il n’existe donc pas de seuil universel en $\beta$ : il dépend de la croissance géométrique de $m(\beta)$.

Pour une fusion réelle, il faut ajouter l’arête gagnante :

$$K\mid X=+1\sim1+\mathrm{Bin}(m-1,s_p(\beta)),$$ $$L^{\mathrm{merge}}  
=\log\frac{K}{m-K}+u_p(1-\beta)(2K-m).$$

À $\beta=1$, sa fiabilité locale vaut exactement $1/m$, et une coupe $m=1$ reste parfaite.

Enfin, la coupe LCA d’une paire lointaine n’est pas typique. Son taux de fusion est

$$r_\beta(C_1,C_2)=m(C_1,C_2)u_ps_p(\beta),$$

et la loi Palm de la paire la repondère encore par le nombre $N_\rho(C_1,C_2)$ de paires lointaines séparées par ses enfants. Le poids géométrique exact est donc

$$\boxed{m(C_1,C_2)N_\rho(C_1,C_2)}.$$

À $p=0{,}8$,

$$h_p(\beta_c)=0.3871644455\ldots,\qquad  
\mathcal J_{\beta_c}=0.1498963\ldots\,m.$$

Une grande coupe critique est donc fortement informative, et le biais Palm favorise justement ces grandes coupes. L’obstruction à la weak recovery doit plutôt venir de l’accumulation, dans le corridor complet, de coupes screenées de faible charge. Il reste à contrôler conjointement

$$(m_v,\beta_v,Z_v,B_v)_{v\in\mathcal C_{ij}},$$

ainsi que les quatre $\Lambda_v^{ab}$ pour $v\succ u$. Une majorité scalaire locale ne suffit pas à cause de leur non-linéarité et des routes latérales.

Le dépôt a été nettoyé et publié directement sur `main`, sans branche supplémentaire :

* Nouvelle note géométrique complète
* Commit `74d03e9`

Validation : 32 fichiers mathématiques contrôlés, 159 tests réussis, dépôt propre et synchronisé.

* * *

**You:**

Merci. Remets en forme le repo' GitHub et notamment le dossier consacré à cette recherche. Qu'il soit à la fois pédagogique et que cette piste de recherche soit indiquée come prioritaire et mise en avant (ainsi que les lemmes qui vont avec comme le fait que le cas le plus favorable soit lorsque i et j lointains se retrouvent dan sla même composante juste à percolation). Fais un beau repo GitHub !

* * *

**ChatGPT:**

Le dépôt a été entièrement remis en forme et publié directement sur `main`, sans branche supplémentaire.

Points principaux :

* Nouvelle page d’accueil du dépôt, avec navigation, objectif scientifique et état des résultats.
* Nouveau programme de recherche prioritaire, qui distingue rigoureusement :
    * le résultat établi : à squelette hiérarchique fixé, faire fusionner $i,j$ au voisinage de $\beta_c$ constitue l’expérience favorable critique ;
    * le résultat exact sur le cactus triangulaire ;
    * le verrou encore ouvert : démontrer la domination favorable après conditionnement par la géométrie Palm critique sur la grille triangulaire.
* Présentation pédagogique de la dynamique, avec intuition, définitions, lemmes, statuts, contre-audits et ordre de lecture.
* Mise en avant de l’événement favorable $$d(i,j)\asymp L,\qquad  
    \beta_{ij}\in[\beta_c-\varepsilon,\beta_c],$$ puis du corridor hiérarchique complet entre $i,j$ et leur LCA — pas seulement du retournement au nœud de fusion.
* Tableau canonique des lemmes établis, conditionnels ou ouverts, notamment pour les lois de coupe, $m h^2$, la loi Palm pondérée par $mN_\rho$, les $\Lambda_v$ ancestraux et la contraction des corridors.
* Cible $p=0.8$ clairement affichée, avec calibration numérique et explication de ce qui manque encore pour obtenir l’obstruction à la weak recovery.
* README spécifiques pour chaque présentation.
* Suppression de sept artefacts LaTeX générés et inutiles ; sources et PDF conservés.
* Ajout d’un validateur de liens et d’un workflow GitHub Actions.

Validation finale :

* 35 documents mathématiques contrôlés ;
* 206 liens relatifs vérifiés ;
* 159 tests réussis ;
* dépôt local propre et synchronisé avec `main` ;
* GitHub Actions réussi.

Commit publié : `ceec3dc`.