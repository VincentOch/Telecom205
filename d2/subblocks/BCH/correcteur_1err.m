function correct = correcteur_1err(mot_code,E)
% bit_code : le message à corriger, de taille 31
%E : la matrice de syndrome
% correct : le mot corrigé
%ATTENTION, il faut construire la matrice de syndrome avant (afin d'améliorer le temps de calcul)

%On regarde si oui ou non le mot est présent dans la matrice de syndrome
%puis on corrige
mot_code_mod = bch_simu_1err(mot_code);
mot_code_transform = str2double(sprintf('%d',mot_code_mod));
if ismember(mot_code_transform,E)
    index = find(E == mot_code_transform);
    %comparaison avec le syndrome pour voir si l'erreur est à la position
    %donnée. 
    mot_code(index) = binary_addition(mot_code(index)+1);
    correct = mot_code(1:26);
    return
end

%Si rien de détecté
correct = mot_code(1:26);
end
