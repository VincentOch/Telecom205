%% [s]=bits2symbols(a,mod,M)
%%
%% a = bits sequence
%% mod = constellation type ('PAM','PSK','QAM')
%% M = constellation size
%%
%% s = symbols sequence (with Eb=1 when uncoded bits are sent) associated with the bits sequence
%%


%% Author : Philippe Ciblat (ciblat@telecom-paris.fr)
%% Location : Telecom Paris
%% Date : 22/07/2021


function [s]=bits2symbols(a,mod,M)

m=log2(M);

%% Size test
if ((length(a)-m*floor(length(a)/m))~=0)
    disp('Bits sequence length not adapted to constellation size')
    s=0;
    return;
end

%% Constellation PAM
if (mod=='PAM')
    V=reshape(a,m,floor(length(a)/m))';
    ipDec = bin2dec(num2str(V(:,[1:m])));
    ipGrayDec = bitxor(ipDec,floor(ipDec/2));
    
    s = ipGrayDec*2-M+1;
    s = sqrt(log2(M))*s/sqrt(sum(([0:M-1]*2-M+1).^2)/M);%%energy purpose
end

%% Constellation PSK
if(mod=='PSK')
    V=reshape(a,m,floor(length(a)/m))';
    ipDec = bin2dec(num2str(V(:,[1:m])));
    ipGrayDec = bitxor(ipDec,floor(ipDec/2));
    
    s = sqrt(log2(M))*exp(2*1i*pi*(ipGrayDec/M));
end


%% Constellation QAM
if(mod=='QAM')
    if (floor(sqrt(M))^2-M ==0)
        V=reshape(a,m,floor(length(a)/m))';
        
        % taking the first bits to code real axis
        ipDecRe = bin2dec(num2str(V(:,[1:m/2])));
        ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
        
        % taking the last bits to code imaginary axis
        ipDecIm = bin2dec(num2str(V(:,[m/2+1:m])));
        ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2));
        
        mm=floor(sqrt(M));
        
        s = (ipGrayDecRe*2-mm+1)+1i*(ipGrayDecIm*2-mm+1);
        s = sqrt(m)*s/sqrt(2*sum(([0:mm-1]*2-mm+1).^2)/mm);%% energy purpose
        
    elseif(M==8)
        % disp('Generating cross 8-QAM');
        constellation = [(1+sqrt(3))*1i; 1+1i; 1+sqrt(3); 1-1i; -(1+sqrt(3)); -1+1i;-(1+sqrt(3))*1i; -1-1i];
        a = reshape(a,m,[]).';
        idx = 1 + a*2.^(m-1:-1:0)';
        s = constellation(idx);
        s = s*sqrt(m)./sqrt(mean(abs(constellation).^2)); % norm to average energy Es = 3Eb and Eb=1;
    else
        
        disp ('Error: M has to be square or equal to 8');
        c = 0;
        return;
    end
end


end




