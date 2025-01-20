function throughput = ARQ(N0, mod, M, iteration, nb_err)
%Renvoit le BER théorique calculé pour différentes valeurs de N0, pour le cas où on utilise un equalizer threshold. 
%On calcule le BER en moyennant sur plusieurs valeurs. 
%mod la modulation souhaitée et M le nombre de points de la constellation
%nb_err le nombre d'erreurs. Si vaut 0 c'est un cas sans codage, si 1 c'est
%un BCH 1 erreur et si 2 c'est un BCH 2 erreurs
%Ne marche que pour une BPSK

seuil = 2; %ARQ avec 1 seule retransmission au plus

BERexp = zeros(1,length(N0)); 
throughput = zeros(1,length(N0));
iter = iteration; % TODO: changer à 10^5
Nc = 31;
%cas sans condage
if (nb_err == 0)
    
    long_vect = 31; %La longueur des mots générés (mc)
    Mc = randi([0,1], 1, iter*long_vect); %génération du vecteur d'entrée
    S = bits2symbols(Mc,mod,M); 

    for jj=1:length(N0) %Calcul pour toutes les valeurs de SNR

        for kk=1:seuil
            flag = kk;
            W = sqrt(N0(jj)/2)*(randn(length(S),1)+1i*randn(length(S),1));
            Z = S + W;
            [S2, ~] = threshold_detector(Z, mod,M);
            M2 = symbols2bits(S2,mod,M);
        
            %Calcul du BER
            BERexp(jj) = biterr(M2,Mc);
            if BERexp(jj) == 0
                if flag == 1
                    break
                end
            end
        end
        if BERexp(jj) == 0
            throughput(jj) = length(Mc)/(flag*length(S));
        end
        
    end
    return
end



%Cas 2 erreurs
if (nb_err == 2)
    long_vect = 21; 
    
    %Génération de la matrice pour une erreur
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

    C = [];
    Mc = [];
    flag = 2;
   
    for kk=1:iter
        m = randi([0,1], 1, long_vect); %génération du vecteur d'entrée
        m_padded = [m zeros(1,10)];
        c = [m bch_simu_2err(m_padded)];
        C = cat(2,C,c);
        Mc = cat(2,Mc,m);
    end
    S = bits2symbols(C,mod,M);

    for jj=1:length(N0)

        for ll=1:seuil

            I = 0;
            flag = ll;
            w = sqrt(N0(jj)/2)*(randn(1,length(S))+1i*randn(1,length(S)));
            s2 = S+w';
            [s2, ~] = threshold_detector(s2,mod,M);
            c2 = symbols2bits(s2,mod,M);

            for uu=1:iter
                m2 = correcteur_2err(c2((uu-1)*Nc+1:uu*Nc),E1,E2,indexList);  
                I = I + biterr(Mc((uu-1)*long_vect+1:uu*long_vect),m2);
            end

            if I == 0 
                if flag == 1
                    break
                end 
            end
        end
        if I == 0
            throughput(jj) = length(Mc)/(flag*length(S));
        end 
    end
    return
end





