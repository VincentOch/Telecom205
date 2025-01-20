%% [a]=symbols2bits(s,mod,M)
%%
%%
%% s = symbols sequence
%% mod = constellation type ('PAM','PSK','QAM')
%% M = constellation size
%%
%% a = bits sequence (associated with the symbols sequence)
%%


%% Author : Philippe Ciblat (ciblat@telecom-paris.fr)
%% Location : Telecom Paris
%% Date : 22/07/2021

function [a]=symbols2bits(s,mod,M)

m=log2(M);

%% Demapping PAM
if (mod=='PAM')
    
    % Gray-Binary Look Up Table
    ipBin = [0:M-1];
    opGray = bitxor(ipBin,floor(ipBin/2));
    [tt ind] = sort(opGray);
    
    % PAM renormalization
    s=sqrt(sum(([0:M-1]*2-M+1).^2)/M)*s/sqrt(log2(M));
    c=(s+M-1)/2;
    
    % sorting Gray code elements to form the lookup table
    opDec = ind(c+1)-1;
    V=dec2bin([opDec,M-1]);
    U=V([1:length(s)] ,:);
    d=reshape(U',1,m*length(s));
    for ii=1:(m*length(s)) 
        a(ii)=str2num(d(ii));
    end
    
end

%% Demapping PSK
if (mod=='PSK')
    
    % Gray-Binary Look Up Table
    ipBin = [0:M-1];
    opGray = bitxor(ipBin,floor(ipBin/2));
    [tt ind] = sort(opGray);
    
    % PAM renormalization
    c=M*angle(s/sqrt(log2(M)))/(2*pi);
    c=round(c-M*floor(c/M));
    
    
    % sorting Gray code elements to form the lookup table
    opDec = ind(c+1)-1;
    V=dec2bin([opDec,M-1]);
    U=V([1:length(s)] ,:);
    d=reshape(U',1,m*length(s));
    for ii=1:(m*length(s)) a(ii)=str2num(d(ii));
    end
    
end

%% Demapping QAM
if(mod=='QAM')
    if (floor(sqrt(M))^2-M ==0)
        
    % Gray-Binary Look Up Table
    ipBin = [0:M-1];
    opGray = bitxor(ipBin,floor(ipBin/2));
    [tt ind] = sort(opGray);
    
    sRe=real(s);
    sIm=imag(s);
    
    % QAM renormalization
    mm=floor(sqrt(M));
    s=sqrt(2*sum(([0:mm-1]*2-mm+1).^2)/mm)*s/sqrt(log2(M));
    sRe=real(s);
    sIm=imag(s);
    
    cRe=(sRe+mm-1)/2;
    cIm=(sIm+mm-1)/2;
    
    % sorting Gray code elements to form the lookup table
    opDecRe = ind(cRe+1)-1;
    opDecIm = ind(cIm+1)-1;
    
    VRe=dec2bin([opDecRe,mm-1]);
    VIm=dec2bin([opDecIm,mm-1]);
    
    V=[VRe,VIm];
    
    U=V([1:length(s)] ,:);
    d=reshape(U',1,m*length(s));
    for ii=1:(m*length(s)) 
        a(ii)=str2num(d(ii));
    end
    
    elseif(M==8) % assumes cross 8-QAM
        %constellation = [(1+sqrt(3))*1i; 1+1i; 1+sqrt(3); 1-1i; -(1+sqrt(3)); -1+1i;-(1+sqrt(3))*1i; -1-1i];      
        [~,loc] = threshold_detector(s,mod,M) ;
        loc = loc -1;
        U = dec2bin(loc,m);
        b = reshape(U',1,length(loc)*m);
        for ii=1:(length(b)) 
            a(ii)=str2num(b(ii));
        end


    else
        disp ('Error: M has to be square or equal to 8');
        return;
    end
end





