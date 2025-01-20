function BERexp = BERZF(N0, mod, M , h ,N,iteration)
%Renvoit le BER théorique calculé pour différentes valeurs de N0, pour le cas où on utilise un equalizer ZF.
%On calcule le BER en moyennant sur plusieurs valeurs.
%mod la modulation souhaitée et M le nombre de points de la constellation
%h correspond au filtre du canal utilisé et N la longueur des trames
%souhaités (également la taille de la matrice de Toeplitz H associé à h)

sprintf('Il y a %d SNR à calculer pour ZF', length(N0))


BERexp = zeros(1,length(N0));

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
P = inv(H); %inverse de H

%On divise par mult car on construit mult symboles s pour le calcul H*s
iter = iteration/mult; % TODO: changer à 10^5

for jj=1:length(N0)
    I = 0;
    sprintf('calcul du BER correspondant au SNR numéro, pour ZF : %d', jj)

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
        Sest = P*Z;
        
        %Calcul du BER
        [S2, ~] = threshold_detector(Sest,mod,M);
        M2 = [];
        for kk=1:mult
            M2 = cat(2,M2,symbols2bits(S2((kk-1)*long+1:kk*long),mod,M));
        end

        I = I + sum(abs(Mc - M2))/length(Mc);


    end
    BERexp(jj) = I/(iter);


end

end


