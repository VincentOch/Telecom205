function gamma = GainSNR(mod,M)
%Calcul le gain de la constellation
if (mod == 'PAM')
    gamma = 6*log2(M)/(M^2-1);
end

if (mod == 'PSK')
    gamma = 2*log2(M)*sin(pi/M)^2;
end

if (mod == 'QAM')
    gamma = 3*log2(M)/(M-1);
end

end
