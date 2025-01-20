function syst = bch_simu_2err(input_sequence)
% input_sequence est de taille  Nc = 31. Input_sequence contient n-k =
% 10 zéros pour les derniers bits parmis ces 31, ceci correspon au déclage de x^(n-k).
% Donc si le message est [m] en entré il y aura [m, 0...0] avec n-k zéros
% syst est la représentation systématique du vecteur input_sequence 

register = zeros(1,10); %initialisation du registre

for ii=1:length(input_sequence)
    temp = register(10); % conservation du dernier registre
    register(2:end) = register(1:end-1); % décalage de tous les registres

    % sommations pour ajouter des erreurs
    register(10) = binary_addition([register(10), temp]);
    register(9) = binary_addition([register(9), temp]);
    register(7) = binary_addition([register(7), temp]);
    register(6) = binary_addition([register(6), temp]);
    register(4) = binary_addition([register(4), temp]);
    register(1) = binary_addition([input_sequence(ii), temp]); % ... input
end
%on renvoit avec mc à la fin
syst = flip(register); %la puissance max est à gauche, terme en position 1

end
