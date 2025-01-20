%% [rej] = d4_perfs(mode,rx)


%% Rejection rate for deliverable D4
%%
%% mode = selected data rate for all the users ('40k','400k','4M','40M')
%% rx = selected receiver for all the users ('zf','dfe')

%% rej = rejection rate versus K (number of users) 

%% Location : Telecom Paris
%% Author : Philippe Ciblat <ciblat@telecom-paris.fr>
%% Date   : 09/06/2023



function [K,rej] = d4_perfs_students(mode,rx)

Ts = 1/(20e6);
MC = 1000; %%number of Monte-Carlo simulations
Pmax = 20; %(dB)
c = 3.0e8;
f0 = 2.4e9;
lambda = c/f0; 
B = 30e6;
N0 = -174; %dBm/Hz
Noise = 10*log10(10^(N0/10)*B);

if(strcmp(mode,'40k') == 1) 
    R = 40*10^(3); 
    Kmax = 2000; 
end

if(strcmp(mode,'400k') == 1) 
    R = 400*10^(3); 
    Kmax = 200; 
end

if(strcmp(mode,'4M') == 1) 
    R = 4*10^(6); 
    Kmax = 20;
end

if(strcmp(mode,'40M') == 1) 
    R = 40*10^(6);
    Kmax = 2; 
end


if(strcmp(rx,'zf') == 1)
    SNR_min_bpsk = [8.9 13.6 21];%minimum Es/N0 for each channel (3-length vector);
    SNR_min_8qam = [11.35 16.2 24];%minimum Es/N0 for each channel (3-length vector);
    SNR_min_16qam = [13.5 18.3 25.8];%minimum Es/N0 for each channel (3-length vector);
end

if(strcmp(rx,'dfe') == 1)
    SNR_min_bpsk = [8.62 11.14 13.6];%minimum Es/N0 for each channel (3-length vector);
    SNR_min_8qam = [11.56 13.85 16.5];%minimum Es/N0 for each channel (3-length vector);
    SNR_min_16qam = [13.1 15.7 18.7];%minimum Es/N0 for each channel (3-length vector);
end


K =  1:(floor(Kmax/10)+1):Kmax ;
rej = zeros(1,length(K));

for kk = 1:length(K)
    aux = 0;

    for mm = 1:MC
        %On génère entre 0 et 1, on mutliplie par 2 et on enlève 1 pour
        %avoir des valeurs entre -1 et 1. On converti ensuite en m
        xx = (-1 + 2*rand(1, K(kk)))*1e3;
        % uniformly-distributed x-axis in a square of semi length 1 (K(kk)
        % lenght vector
        yy = (-1 + 2*rand(1, K(kk)))*1e3; 
        % uniformly-distributed y-axis in a square of semi-length 1 (K(kk) 
        % length vector);

        d2 = xx.^2+yy.^2;%% square-distance from the origin 
        a2 = min(1, lambda^2./(4^2*pi^2.*d2));% square magnitude attenuation in Friis equation  
        a2_dB = 10*log10(a2);
        SNR_rx_max = a2_dB + Pmax - Noise; %calculate ??? : SNRmax at TX in dB.

        SNR_min = SNR_min_16qam; %%ZF or DFE/16QAM: channel 1, channel 2, channel 3

        if(R*Ts*K(kk)<3)  
            SNR_min = SNR_min_8qam; 
        end%%ZF or DFE/8QAM: channel 1, channel 2, channel 3

        if(R*Ts*K(kk)<1) 
            SNR_min = SNR_min_bpsk; 
        end %%ZF or DFE/BPSK: channel 1, channel 2, channel 3

        SNR_rx_min = SNR_min(randi(3,K(kk),1));

        aux = aux+sum(SNR_rx_max < SNR_rx_min); %J'ai remplacé le calcul de aux pour gagner un peu en compréhension

    end

rej(kk) =  aux/(K(kk)*MC);


end


%% Channel display
