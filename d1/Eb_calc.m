function Eb = Eb_calc(mod,M)
%Calcul le Eb de différentes constellations
%R le gain de codage si il y a codage (1 si pas de codage)

if (mod == 'PAM')
    Es = (M^2-1)/3;
    Eb = Es/(log2(M));
end

if (mod == 'PSK')
    Es = 1;
    Eb = Es/(log2(M));
end

if (mod == 'QAM')
    Es = 2*(M-1)/3;
    Eb = Es/(log2(M));
end

end
