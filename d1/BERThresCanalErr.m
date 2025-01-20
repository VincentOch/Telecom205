function BERexp = BERThresCanalErr(N0, mod, M , h, N,iteration,nb_err)
%Renvoit le BER théorique calculé pour différentes valeurs de N0, pour le
%cas où on utilise un equalizer DFE
%On calcule le BER en moyennant sur plusieurs valeurs.
%mod la modulation souhaitée et M le nombre de points de la constellation
%h correspond au filtre du canal utilisé et N la longueur des trames
%souhaités (également la taille de la matrice de Toeplitz H associé à h)

sprintf('Il y a %d SNR à calculer, via threshold', length(N0))
Nc = 31;

BERexp = zeros(1,length(N0));

if (nb_err == 0)
    %On va régler la taille pour qu'elle soit compatible avec toutes les
    %constellations
    a = 21; %La longueur voulue des mots générés (mc)
    m = log2(M);
    long_vect = m*floor(a/m);

    long = length(bits2symbols(zeros(1,long_vect),mod,M)); %longeur des mots s, soit m modulé par mod avec une constellation de M
    L = (length(h)-1)/2;

    %On modifie N pour avoir un N multiple de la taille des s, pour avoir un nombre entier de s quand on effectue z=Hs+w
    %On prend l'entier supérieur pour être au moins au dessus de 100
    mult = floor(N/long)+1;
    N = mult*long;

    H = toeplitzMat(h,L,N);

    %On divise par mult car on construit mult symboles s pour le calcul H*s
    iter = iteration/mult; % TODO: changer à 10^5

    parfor jj=1:length(N0)
        I = 0;
        sprintf('calcul du BER correspondant au SNR numéro, via Threshold, sans erreurs : %d', jj)
        for ii=1:iter

            W = sqrt(N0(jj)/2)*(randn(N,1)+1i*randn(N,1));
            Mc = [];
            S = [];

            %On fait des boucles au lieu de générer directement un vecteur de
            %la bonne taille afin de gérer les tailles des différentes
            %constellations
            for ll=1:mult
                mc = randi([0,1], 1,long_vect);
                S = cat(1,S,bits2symbols(mc,mod,M));
                Mc = cat(2,Mc,mc);
            end

            Z = H*S+W;

            %Calcul du BER
            [S2, ~] = threshold_detector(Z,mod,M);
            M2 = [];
            for kk=1:mult
                M2 = cat(2,M2,symbols2bits(S2((kk-1)*long+1:kk*long),mod,M));
            end
            I = I + sum(abs(Mc - M2))/length(Mc);
        end
        BERexp(jj) = I/(iter);
    end
end

if (nb_err == 1)
    %On va régler la taille pour qu'elle soit compatible avec toutes les
    %constellations
    a = 26; %La longueur voulue des mots générés (mc)
    
    m = log2(M);
    long_vect = m*floor(a/m);

    long = length(bits2symbols(zeros(1,Nc),mod,M)); %longeur des mots s, soit m modulé par mod avec une constellation de M
    L = (length(h)-1)/2;

    %On modifie N pour avoir un N multiple de la taille des s, pour avoir un nombre entier de s quand on effectue z=Hs+w
    %On prend l'entier supérieur pour être au moins au dessus de 100
    mult = floor(N/long)+1;
    N = mult*long;

    H = toeplitzMat(h,L,N);

    %On divise par mult car on construit mult symboles s pour le calcul H*s
    iter = iteration/mult; % TODO: changer à 10^5

    E = zeros(1,Nc);
    for ll = 1:Nc
        err = zeros(1,Nc);
        err(ll) = 1;
        err_code = bch_simu_1err(err);
        E(ll) = str2double(sprintf('%d', err_code));
    end

    parfor jj=1:length(N0)
        I = 0;
        sprintf('calcul du BER correspondant au SNR numéro, via Threshold, avec 1 erreur : %d', jj)
        for ii=1:iter

            W = sqrt(N0(jj)/2)*(randn(N,1)+1i*randn(N,1));
            Mc = [];
            S = [];

            %On fait des boucles au lieu de générer directement un vecteur de
            %la bonne taille afin de gérer les tailles des différentes
            %constellations
            for ll=1:mult
                mc = randi([0,1], 1,long_vect);
                m_padded = [mc zeros(1,5)];
                c = [mc bch_simu_1err(m_padded)];
                S = cat(1,S,bits2symbols(c,mod,M));
                Mc = cat(2,Mc,mc);
            end

            Z = H*S+W;

            %Calcul du BER
            [S2, ~] = threshold_detector(Z,mod,M);
            M2 = [];
            for kk=1:mult
                c2 = symbols2bits(S2((kk-1)*long+1:kk*long),mod,M);
                M2 = cat(2,M2,correcteur_1err(c2,E));
            end
            I = I + sum(abs(Mc - M2))/length(Mc);
        end
        BERexp(jj) = I/(iter);
    end
end

if (nb_err == 2)
    %On va régler la taille pour qu'elle soit compatible avec toutes les
    %constellations
    a = 21; %La longueur voulue des mots générés (mc)
    m = log2(M);
    long_vect = m*floor(a/m);

    long = length(bits2symbols(zeros(1,Nc),mod,M)); %longeur des mots s, soit m modulé par mod avec une constellation de M
    L = (length(h)-1)/2;

    %On modifie N pour avoir un N multiple de la taille des s, pour avoir un nombre entier de s quand on effectue z=Hs+w
    %On prend l'entier supérieur pour être au moins au dessus de 100
    mult = floor(N/long)+1;
    N = mult*long;

    H = toeplitzMat(h,L,N);

    %On divise par mult car on construit mult symboles s pour le calcul H*s
    iter = iteration/mult; % TODO: changer à 10^5

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

    parfor jj=1:length(N0)
        I = 0;
        sprintf('calcul du BER correspondant au SNR numéro, via Threshold, avec 2 erreurs : %d', jj)
        for ii=1:iter

            W = sqrt(N0(jj)/2)*(randn(N,1)+1i*randn(N,1));
            Mc = [];
            S = [];

            %On fait des boucles au lieu de générer directement un vecteur de
            %la bonne taille afin de gérer les tailles des différentes
            %constellations
            for ll=1:mult
                mc = randi([0,1], 1,long_vect);
                m_padded = [mc zeros(1,10)];
                c = [mc bch_simu_2err(m_padded)];
                S = cat(1,S,bits2symbols(c,mod,M));
                Mc = cat(2,Mc,mc);
            end

            Z = H*S+W;

            %Calcul du BER
            [S2, ~] = threshold_detector(Z,mod,M);
            M2 = [];
            for kk=1:mult
                c2 = symbols2bits(S2((kk-1)*long+1:kk*long),mod,M);
                M2 = cat(2,M2,correcteur_2err(c2,E1,E2,indexList));
            end
            I = I + sum(abs(Mc - M2))/length(Mc);
        end
        BERexp(jj) = I/(iter);
    end
end

end


