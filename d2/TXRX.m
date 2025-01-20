% Rev: Dec. 2023, Germain
% 

close all;
clear;





K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
Rin = 50;       % Input impedance
N   = 2^19;     % Number of signal points ADC
BW  = 10e6;     % Signal bandwidth
Flo = 240e6;    % Central frequency      (divided by 10 to reduce the simulation time)

basebandSamplingRate          = 30e6;     % sampling frequency ADC
OSR_ct_emulation              = fix((19.98e9/100)/basebandSamplingRate); % sampling frequency Simulation (divided by 100 to reduce the simulation time)
continuousTimeSamplingRate    = OSR_ct_emulation*basebandSamplingRate;
Ts_Cont                       = 1/continuousTimeSamplingRate; % sampling period Simulation


N_Sim   = floor(N*continuousTimeSamplingRate/basebandSamplingRate); % Number of  Simulation points
t_Sim   = 0:Ts_Cont:(N_Sim-1)*Ts_Cont; % Time vector Simulation
f_Sim   = 0:continuousTimeSamplingRate/N_Sim:continuousTimeSamplingRate/2-continuousTimeSamplingRate/N_Sim; % Frequency vector Simulation



%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[InputAudio,fsAudio]    = audioread('lab_data/out.wav');
InputAudio              = InputAudio(1:N,1)/max(InputAudio(:,1));





%%%%%%%%%%%%%%%%% Transmistter Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DAC Specifications %%%
nBitDAC   = 12;
Vref_DAC  = 1;
dacType   = 'zoh';

%%%%%%% Filter %%%%%%%%%%%%
TXBB_Filt_Order   = 8;
TXBB_Filt_NF      = 2;
TXBB_Filt_Fcut    = 12e6;

%%%%%%%%%%%% PA %%%%%%%%%%
PA_NF   = 5;
PA_Gain = 14;
PA_IIP3 = 30;
 

%%%%%%%%%%%%%%%% Channel Attenuation %%%%%%%%%%%%

Channel_Attenuation       = 80;
Channel_Attenuation_lin   = 10^(-Channel_Attenuation/20);
AntennaNoise              =randn(1,N_Sim)*sqrt(K*T*continuousTimeSamplingRate/2*Rin);






%%%%%%%%%%%%%%%%% Receiver Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LNA %%%

LNA_Gain = 15;     % (dB)
LNA_IIP3 = -45;    % (dBm)
LNA_NF   = 2;      % (dB)



%%%%%%% Filter %%%%%%%%%%%%
RXBB_Filt_Order   = 10;
RXBB_Filt_NF      = 0;
RXBB_Filt_Fcut    = 12e6;

%%%%% ADC %%%%%%%
nBitADC   = 13;
Vref_ADC  = 1;
delay     = 0;
BB_gain   = 20;



TxOut       = TX(   InputAudio,...
                    InputAudio,...
                    Vref_DAC,...
                    nBitDAC,...
                    basebandSamplingRate,...
                    continuousTimeSamplingRate,...
                    dacType,...
                    TXBB_Filt_Order,...
                    TXBB_Filt_NF,...
                    TXBB_Filt_Fcut,...
                    Flo,...
                    PA_IIP3,...
                    PA_NF,...
                    PA_Gain);

rxSignal    = TxOut*Channel_Attenuation_lin+AntennaNoise';

[basebandAnalog_adc_I, basebandAnalog_adc_Q] = RX(  rxSignal,...
                                                    LNA_IIP3,...
                                                    LNA_NF,...
                                                    LNA_Gain,...
                                                    Flo,...
                                                    continuousTimeSamplingRate,...
                                                    basebandSamplingRate,...
                                                    nBitADC,...
                                                    Vref_ADC,...
                                                    delay,...
                                                    BB_gain);     



plot_spectrum(TxOut*sqrt(1e3/Rin),1,continuousTimeSamplingRate,1);
xlabel('frequency (Hz)')
ylabel('PSD (dBm/bin)')
title('TX Output')

plot_spectrum(basebandAnalog_adc_I*sqrt(1e3/Rin),2,basebandSamplingRate,1);
xlabel('frequency (Hz)')
ylabel('PSD (dBm/bin)')
title('ADC\_output')
% 
soundsc(basebandAnalog_adc_I,fsAudio);





