function BERexp = BERThres(N0, mod, M, iteration, nb_err)
%Renvoit le BER théorique calculé pour différentes valeurs de N0, pour le cas où on utilise un equalizer threshold. 
%On calcule le BER en moyennant sur plusieurs valeurs. 
%mod la modulation souhaitée et M le nombre de points de la constellation
%nb_err le nombre d'erreurs. Si vaut 0 c'est un cas sans codage, si 1 c'est
%un BCH 1 erreur et si 2 c'est un BCH 2 erreurs
%Ne marche que pour une BPSK

sprintf('Il y a %d SNR à calculer', length(N0))


BERexp = zeros(1,length(N0)); 
iter = iteration; % TODO: changer à 10^5
Nc = 31;

%cas sans condage
if (nb_err == 0)
    
    long_vect = 21; %La longueur des mots générés (mc)
    long = length(bits2symbols(zeros(1,long_vect),mod,M)); %longeur des mots s, soit m modulé par mod avec une constellation de M


    for jj=1:length(N0) %Calcul pour toutes les valeurs de SNR

        
        sprintf('calcul du BER correspondant au SNR numéro, cas sans codage : %d', jj)
    
        %On construit de nombreux mots de codes (iter mots de codes) pour
        %pouvoir moyenner
    
        Mc = randi([0,1], 1, iter*long_vect); %génération du vecteur d'entrée
        S = bits2symbols(Mc,mod,M); 
        W = sqrt(N0(jj)/2)*(randn(length(S),1)+1i*randn(length(S),1));
        Z = S + W;
        [S2, ~] = threshold_detector(Z, mod,M);
        M2 = symbols2bits(S2,mod,M);
    
        %Calcul du BER
        BERexp(jj) = sum(abs(Mc - M2))/length(Mc);
    end
    return
end


%Avec codage et correcteur 1 erreur
if (nb_err == 1)
    long_vect = 26; %La longueur des mots générés (mc)

    %Construction de la matrice d'erreur (à faire avant pour limiter les
    %calculs)
    E = zeros(1,Nc);
    for ll = 1:Nc
        err = zeros(1,Nc);
        err(ll) = 1;
        err_code = bch_simu_1err(err);
        E(ll) = str2double(sprintf('%d', err_code));
    end
    
    %calcul du ber et autre
    for jj=1:length(N0)
        I = 0;
        sprintf('calcul du BER correspondant au SNR numéro, cas 1 erreur : %d', jj)

        for ii=1:iter
            m = randi([0,1], 1, long_vect); %génération du vecteur d'entrée
            m_padded = [m zeros(1,5)];
            c = [m bch_simu_1err(m_padded)];
            s = bits2symbols(c,mod,M);
            w = sqrt(N0(jj)/2)*(randn(1,length(s))+1i*randn(1,length(s)));
            s = s+w';
            [s2, ~] = threshold_detector(s,mod,M);
            c2 = symbols2bits(s2,mod,M);
            m2 = correcteur_1err(c2,E);
            I = I + sum(abs(m - m2))/long_vect;
        end
        BERexp(jj) = I/iter;
    end
    return
end

%Cas 2 erreurs
if (nb_err == 2)
    long_vect = 21; 
    long = length(bits2symbols(zeros(1,long_vect),mod,M));
    
    %Génération de la atrice pour une erreur
    E1 = zeros(1,Nc);
    for iii = 1:Nc
        err = zeros(1,Nc);
        err(iii) = 1;
        E1(iii) = str2double(sprintf('%d',bch_simu_2err(err))); %on concatène les bits pour faire un seul entier et pouvoir comparer plus facilement
    end
    
    %génération de la matrice pour 2 erreurs
    E2 = zeros(1,Nc*(Nc-1)/2);
    indexList = zeros(2,Nc*(Nc-1)/2);
    u = 1;
    for kkk = 1:Nc-1
     for jjj = kkk+1:Nc
        err = zeros(1,Nc);
        err(kkk) = 1;
        err(jjj) = 1;
        E2(u) = str2double(sprintf('%d',bch_simu_2err(err)));
        indexList(u,1) = iii;
        indexList(u,2) = jjj;
        u=u+1;
     end
    end

    for jj=1:length(N0)
        I = 0;
        sprintf('calcul du BER correspondant au SNR numéro, cas 2 erreurs : %d', jj)

        for ii=1:iter
            m = randi([0,1], 1, long_vect); %génération du vecteur d'entrée
            m_padded = [m zeros(1,10)];
            c = [m bch_simu_2err(m_padded)];
            s = bits2symbols(c,mod,M);
            w = sqrt(N0(jj)/2)*(randn(1,length(s))+1i*randn(1,length(s)));
            s = s+w';
            [s2, ~] = threshold_detector(s,mod,M);
            c2 = symbols2bits(s2,mod,M);
            m2 = correcteur_2err(c2,E1,E2,indexList);   
            I = I + sum(abs(m - m2))/long_vect;
        end
        BERexp(jj) = I/iter;
    end
    return
end

end




