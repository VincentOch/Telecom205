close all;
clear


addpath(genpath('../subblocks/'))
%%% Signal Generation %%%
adcSamplingRate    = 20e6;     % A sampling rate which is sufficiently high in order to be close to the continous time
% 







%%General%%%
Ts=1/adcSamplingRate; %sampling period 
Subsamp_fac=round(adcSamplingRate/adcSamplingRate);
N=2^13*Subsamp_fac;  %Number of signal points
BW=10e6; %Signal bandwidth
K=1.38e-23; %Boltzmann Constant
T=290;      %room temperature
t=0:Ts:(N-1)*Ts; %Time vetor  
f=0:adcSamplingRate/N:adcSamplingRate/2-adcSamplingRate/N; %Frequency vector
Rin=50;  %Matching impedance chosen equal to 50
AntennaNoise=randn(1,N)*sqrt(K*T*adcSamplingRate/2*Rin);

%%%%%%%%%Input signal %%%%%%%%%
fin=6.9015e6;    %Input sine frequency
Bin_in=round(fin/adcSamplingRate*N); %Determining the input bin
fin=Bin_in*adcSamplingRate/N;
Pin=-40;            %Pin in dBm
Ain=sqrt(10.^((Pin-30)/10)*2*Rin);

%%/Input Signal Generation%%/
basebandSig=Ain*sin(2*pi*fin*t+rand())+AntennaNoise;%Input signal
basebandSig=basebandSig';



%%% Analog to Digital Conversion %%%
nBitADC = 14;
delay   = 0; 
FullScaleADC    = 2;
% Perform conversion
basebandAnalog_adc = ADC(basebandSig,nBitADC,FullScaleADC/2,adcSamplingRate,delay,adcSamplingRate);



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Theo%%%%%%%%%
OSR_ADC=adcSamplingRate/(2*BW);
q=FullScaleADC./2^nBitADC;
NA_ADC=q^2/12/OSR_ADC/Rin; %We divide by R to convert to a power 
F_ADC=NA_ADC/(K*T*BW)+1;



NF_tot_theo=10*log10(F_ADC);

SNR_IN_theo=Pin-(10*log10(K*T*BW)+30);




%%%%%%%%%%%%%SNR Practical%%%%%%%%%
Bin_Limits=[1 round(BW/adcSamplingRate*N)];
Bin_Limits_ADC=[1 round(BW/adcSamplingRate*length(basebandAnalog_adc))];
[SNR_IN,PS]=perf_estim(basebandSig,Bin_in,5,Bin_Limits,0);

SNR_ADC=perf_estim(basebandAnalog_adc,Bin_in,5,Bin_Limits_ADC,0);
SNR_ADC_theo=SNR_IN_theo-NF_tot_theo;

NF_tot=SNR_IN-SNR_ADC;


disp(['The simulation NF is ',num2str(NF_tot),' dB'])  
disp(['The theoritical NF is ',num2str(NF_tot_theo),' dB'])  





 
 
 
figure(1)
subplot(2,1,1)
plot_spectrum(basebandSig*sqrt(1e3/Rin),1, adcSamplingRate,1);
xlabel('frequency(Hz)')
ylabel('PSD(dBm/bin)')
title('Signal at the filter input')

subplot(2,1,2)
plot_spectrum(basebandAnalog_adc*sqrt(1e3/Rin),1,adcSamplingRate,4);
xlabel('frequency(Hz)')
ylabel('PSD(dBm/bin)')
title('output of the ADC')
axis([0,adcSamplingRate/2,-200,1])

% 
% disp(['The maximum voltage at the filter output is ',num2str(max(Filter_out)),' V'])  
% disp(['The SNR at the ADC output is ',num2str(SNR_ADC),'dB'])

