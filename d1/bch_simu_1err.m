function syst = bch_simu_1err(input_sequence)
% input_sequence est de taille  Nc = 31. Input_sequence contient n-k =
% 5 zéros pour les derniers bits parmis ces 31, ceci correspon au déclage de x^(n-k).
% Donc si le message est [m] en entré il y aura [m, 0...0] avec n-k zéros
% syst est la représentation systématique du vecteur input_sequence 
register = zeros(1,5); %initialisation du registre

%Calcul de m modulo g
for ii=1:length(input_sequence)
    temp = register(5); % conservation du dernier registre
    register(2:end) = register(1:end-1); % décalage de tous les registres

    % sommations pour ajouter des erreurs
    register(3) = binary_addition([register(3), temp]);
    register(1) = binary_addition([input_sequence(ii), temp]);
end
%on renvoit avec mc à la fin
syst = flip(register); %la puissance max est à gauche, terme en position 1
end
