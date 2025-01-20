function  nmin = Nmin(mod,M)
%NMIN Summary of this function goes here
%   Detailed explanation goes here
if (mod == 'PAM')
    nmin = 2*(M-1)/M;
end

if (mod == 'PSK')
    nmin = min(2,M);
end

if (mod == 'QAM')
    nmin = 4/log2(M);
end

end
