Outre les fonctions données, voici briévement celles utilisées que nous avons crée.

plot_channel permet de tracer les figures pour les fonctions h représentatives des canaux
pour plot les différents BER, cela se trouve dans BER_plot. Attention si iteration est trop élevé le plot est très long, mais plus précis ce qui évite certaines 'cassure' apparente dans les plots

correcteur_1err et 2err corrigent respectivement 1 ou 2 erreurs pour un code BCH grâce à la méthode des syndromes
bch_simu_1err et 2err calculent la division de modulo pour les codes BCH corrigeant 1 ou 2 erreurs
On utilise binary_addition pour simplifier certaines expressions (juste addition modulo 2)
toeplitzMat calcule la matrice de Toeplitz d'un vecteur (comme dans l'énoncé)

BERThres permet de calculer (seulement pour BPSK) les BER théoriques dans les cas sans codage ou avec BCH
BERThresCanal applique un equalizer threshold pour le calcul théorique avec canal
BERDFE fait l'equalizer DFE
BERZE fait l'equalize ZF

Nmin, GainSNR et Eb_calc permettent de calculer les différents éléments pour plot un BER théorique