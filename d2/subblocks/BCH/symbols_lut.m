%% [s] = symbols_lut(mod,M)
%%
%% Look Up Table for various constellations
%%
%% mod = constellation type
%% M = constellation size
%%
%% s = vector with all symbols of the considered constellation

%% Location: Telecom ParisTech
%% Author:  Ph. Ciblat <ciblat@telecom-paristech.fr>
%% Date : 16/10/2018

function [s] = symbols_lut(mod,M)

if (M-2^(floor(log2(M)))~= 0)
    disp ('Error : M is not a power of 2');
    s = 0;
    return;
end;


if(mod=='PAM')
    s = [0:(M-1)]*2-M+1;
    s = sqrt(log2(M))*s/sqrt(sum(([0:M-1]*2-M+1).^2)/M);
end

if(mod=='PSK')
    s = sqrt(log2(M))*exp(2*i*pi*([0:(M-1)]/M));
end

if(mod=='QAM')
    if (floor(sqrt(M))^2-M ==0)
        m=floor(sqrt(M));
        
        for jj=0:(m-1)
            s(1,[1:m]+jj*m) = [0:(m-1)]*2-m+1+1i*(jj*2-m+1);
        end
        
        s = sqrt(log2(M))*s/sqrt(2*sum(([0:m-1]*2-m+1).^2)/m);
    
    elseif (M==8)
        s = [(1+sqrt(3))*1i; 1+1i; 1+sqrt(3); 1-1i; -(1+sqrt(3)); -1+1i;-(1+sqrt(3))*1i; -1-1i].';
        s = s*sqrt(log2(M))./sqrt(mean(abs(s).^2)); % norm to average energy Es = 3Eb and Eb =1;
        
    else
        disp ('Error : M has to be square or equal to 8');
        s = 0;
        return;
    end
    
    
end

if(mod=='OOK')
    s = [0:(M-1)] * sqrt(M * log2(M) * 6 / ((M-1)*M*(2*M-1)));
end
