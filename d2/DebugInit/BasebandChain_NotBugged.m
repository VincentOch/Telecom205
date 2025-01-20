close all;
clear


addpath(genpath('../subblocks/'))
%%% Signal Generation %%%
continuousTimeSamplingRate    = 1e9;     % A sampling rate which is sufficiently high in order to be close to the continous time
basebandSamplingRate_or       = 25e6;     % The sampling rate of the complex baseband signal ; Units : Samples/second


basebandSamplingRate          = continuousTimeSamplingRate/round(continuousTimeSamplingRate/basebandSamplingRate_or);
adcSamplingRate = basebandSamplingRate;






%%General%%%
Ts=1/continuousTimeSamplingRate; %sampling period 
Subsamp_fac=round(continuousTimeSamplingRate/adcSamplingRate);
N=2^17*Subsamp_fac;  %Number of signal points
BW=10e6; %Signal bandwidth
K=1.38e-23; %Boltzmann Constant
T=290;      %room temperature
t=0:Ts:(N-1)*Ts; %Time vetor  
f=0:continuousTimeSamplingRate/N:continuousTimeSamplingRate/2-continuousTimeSamplingRate/N; %Frequency vector
Rin=50;  %Matching impedance chosen equal to 50


%%%%%%%%% Error 1:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A multiplication by Rin was missing to convert from W to V^2
%AntennaNoise=randn(1,N)*sqrt(K*T*continuousTimeSamplingRate/2);
AntennaNoise=randn(1,N)*sqrt(K*T*continuousTimeSamplingRate/2*Rin);
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Input signal %%%%%%%%%
fin=6.9015e6;    %Input sine frequency
Bin_in=round(fin/continuousTimeSamplingRate*N); %Determining the input bin

%%%%%%%%% Error 2:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input signal was not in a FFT bin
fin=Bin_in*continuousTimeSamplingRate/N;
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pin=-60;            %Pin in dBm
Ain=sqrt(10.^((Pin-30)/10)*2*Rin);

%%/Input Signal Generation%%/
basebandSig=Ain*sin(2*pi*fin*t+rand())+AntennaNoise;%Input signal
basebandSig=basebandSig';





%%% Baseband Analog filter %%%
RXBB_Filt_NF    = 5;     %(in dB)
RXBB_Filt_Fcut  = 12.5e6;  % Filter TX BB Fcut 3dB Frequency
%%%%%%%%% Error 5:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Low pass noise folding was overlooked in the theoretical calculation
%RXBB_Filt_Order = 3;     % Filter TX BB Order
RXBB_Filt_Order = 10;     % Filter TX BB Order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Instanciate filter (due to numerical issue, the filter has to be instanciated as SOS)
[RXBB_Filt_z,RXBB_Filt_p,RXBB_Filt_k]=butter(RXBB_Filt_Order,RXBB_Filt_Fcut/(continuousTimeSamplingRate/2));
RXBB_Filt_sos = zp2sos(RXBB_Filt_z,RXBB_Filt_p,RXBB_Filt_k);
% Perform filtering
basebandAnalog_filtrx = basebandAnalogFilt(basebandSig,RXBB_Filt_sos,RXBB_Filt_NF,Rin,continuousTimeSamplingRate);


%%% Baseband Gain %%%
BBamp_Gain    = 30; % (dB)
BBamp_IIP3    = 100; % (dBm)
BBamp_NF      = 5; % (dB)
basebandAnalog_amp = BBamp(basebandAnalog_filtrx,BBamp_Gain,BBamp_NF,BBamp_IIP3,Rin,BW,continuousTimeSamplingRate);





%%% Analog to Digital Conversion %%%
nBitADC = 14;
delay   = 0; 
FullScaleADC    = 2;
% Perform conversion
%%%%%%%%% Error 4:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The ADC Full Scale was confused with its reference voltage
%basebandAnalog_adc = ADC(basebandAnalog_amp,nBitADC,FullScaleADC,adcSamplingRate,delay,continuousTimeSamplingRate);
basebandAnalog_adc = ADC(basebandAnalog_amp,nBitADC,FullScaleADC/2,adcSamplingRate,delay,continuousTimeSamplingRate);
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Theo%%%%%%%%%
F_Filter=10.^(RXBB_Filt_NF/10);
F_BBamp=10.^(BBamp_NF/10);

G_Filter=1;
G_VGA=10.^(BBamp_Gain/10); % The conversion is done with 10*log10 because we need a power gain 

OSR_ADC=adcSamplingRate/(2*BW);
q=FullScaleADC./2^nBitADC;
NA_ADC=q^2/12/OSR_ADC/Rin; %We divide by R to convert to a power 
F_ADC=NA_ADC/(K*T*BW)+1;

F_InToVGA=F_Filter+(F_BBamp-1)/G_Filter;
NF_InToVGA=10*log10(F_InToVGA);

F_tot_theo=F_Filter+(F_BBamp-1)/G_Filter+(F_ADC-1)/(G_Filter*G_VGA);
NF_tot_theo=10*log10(F_tot_theo);

SNR_IN_theo=Pin-(10*log10(K*T*BW)+30);

SNR_VGA_theo=SNR_IN_theo-NF_InToVGA;
SNR_ADC_theo=SNR_IN_theo-NF_tot_theo;
SNR_filter_theo=SNR_IN_theo-RXBB_Filt_NF;


%%%%%%%%%%%%%SNR Practical%%%%%%%%%
Bin_Limits=[1 round(BW/continuousTimeSamplingRate*N)];
Bin_Limits_ADC=[1 round(BW/adcSamplingRate*length(basebandAnalog_adc))];



%%%%%%%%% Error 2:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The bin position was wrongly transmitted to the SNR function 
% [SNR_IN,PS]=perf_estim(basebandSig,Bin_in,5,Bin_Limits,1);
% SNR_filter=perf_estim(basebandAnalog_filtrx,Bin_in,5,Bin_Limits,1);
% SNR_VGA=perf_estim(basebandAnalog_amp,Bin_in,5,Bin_Limits,1);
% SNR_ADC=perf_estim(basebandAnalog_adc,Bin_in,5,Bin_Limits_ADC,1);
[SNR_IN,PS]=perf_estim(basebandSig,Bin_in,5,Bin_Limits,0);
SNR_filter=perf_estim(basebandAnalog_filtrx,Bin_in,5,Bin_Limits,0);
SNR_VGA=perf_estim(basebandAnalog_amp,Bin_in,5,Bin_Limits,0);
SNR_ADC=perf_estim(basebandAnalog_adc,Bin_in,5,Bin_Limits_ADC,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




NF_tot=SNR_IN-SNR_ADC;


disp(['The simulation NF is ',num2str(NF_tot),' dB'])  
disp(['The theoritical NF is ',num2str(NF_tot_theo),' dB'])  



figure(2)
subplot(2,2,1)
plot(t*1e6,basebandSig);
xlabel('time (us)')
ylabel('Amplitude(v)')
xlim([0 2])
title('Signal at the filter input')

subplot(2,2,2)
plot(t*1e6,basebandAnalog_filtrx);
xlabel('time (us)')
ylabel('Amplitude(v)')
xlim([0 2])
title('Signal at the filter output')

subplot(2,2,3)
plot(t*1e6,basebandAnalog_amp);
xlabel('time (us)')
ylabel('Amplitude(v)')
xlim([0 2])
title('output of the VGA')

t_adc=0:1/adcSamplingRate:(length(basebandAnalog_adc)-1)/adcSamplingRate;
subplot(2,2,4)
plot(t_adc*1e6,basebandAnalog_adc)
xlabel('time (us)')
ylabel('Amplitude(v)')
xlim([0 2])
title('output of the ADC')


 
 
 
figure(1)
subplot(2,2,1)
plot_spectrum(basebandSig*sqrt(1e3/Rin),1, continuousTimeSamplingRate,1);
xlabel('frequency(Hz)')
ylabel('PSD(dBm/bin)')
title('Signal at the filter input')
axis([0,100e6,-200,1])
subplot(2,2,2)
plot_spectrum(basebandAnalog_filtrx*sqrt(1e3/Rin),1,continuousTimeSamplingRate,2);
xlabel('frequency(Hz)')
ylabel('PSD(dBm/bin)')
title('Signal at the filter output')
axis([0,100e6,-200,1])
subplot(2,2,3)
plot_spectrum(basebandAnalog_amp*sqrt(1e3/Rin),1,continuousTimeSamplingRate,3);
xlabel('frequency(Hz)')
ylabel('PSD(dBm/bin)')
title('output of the VGA')
axis([0,100e6,-200,1])
subplot(2,2,4)
plot_spectrum(basebandAnalog_adc*sqrt(1e3/Rin),1,adcSamplingRate,4);
xlabel('frequency(Hz)')
ylabel('PSD(dBm/bin)')
title('output of the ADC')
axis([0,adcSamplingRate/2,-200,1])

% 
% disp(['The maximum voltage at the filter output is ',num2str(max(Filter_out)),' V'])  
% disp(['The SNR at the ADC output is ',num2str(SNR_ADC),'dB'])

