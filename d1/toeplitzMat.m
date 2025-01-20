function H = toeplitzMat(h,L,N)
%H la matrice de Toepliz composée de h(-L),...,h(0),...h(L) ou par
%hypothèse, h(-L) correspond au premier indice de h et h(L) au dernier. h
%est de longueur 2L+1. h est le filtre d'un canal

H = zeros(N,N);
for ii =1:N
    for jj =1:N
        index = L+1-ii+jj;
        if index >= 1 && index<= length(h)
            H(ii,jj) = h(index);
        end
    end
end
end

