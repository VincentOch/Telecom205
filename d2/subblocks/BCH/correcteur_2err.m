function correct = correcteur_2err(mot_code,E1,E2,indexList)
% bit_code : le message à corriger, de taille 31
%E1 matrice de syndorme si une erreur
%E2 matrice de syndrome si 2 erreurs
%indeList la liste de correspondance des positions des erreurs avec leur
%position dans E2
% correct : le mot corrigé
%ATTENTION, il faut construire la matrice de syndrome avant (afin d'améliorer le temps de calcul)




%On regarde si oui ou non le mot est présent dans la matrice de syndrome
%puis on corrige
mot_code_mod = bch_simu_2err(mot_code);
mot_code_transform = str2double(sprintf('%d',mot_code_mod));
if ismember(mot_code_transform,E1)
    index = find(E1 == mot_code_transform);
    %comparaison avec le syndrome pour voir si l'erreur est à la position
    %donnée. 
    mot_code(index) = binary_addition(mot_code(index)+1);
    correct = mot_code(1:21);
    return
end

%construction des syndromes pour 2 erreus

if ismember(mot_code_transform,E2)
    index = find(E2 == mot_code_transform);
    ii = indexList(index,1);
    jj = indexList(index,2);
    mot_code(ii) = binary_addition(mot_code(ii)+1);
    mot_code(jj) = binary_addition(mot_code(jj)+1);
    correct = mot_code(1:21);
    return
end

correct = mot_code(1:21);
end



