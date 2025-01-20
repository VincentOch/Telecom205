% Rev: Dec. 2023, Germain
% 
% This is the TXRX.m version  for TELECOM205

close all;
clear;





K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
Rin = 50;       % Input impedance
N   = 2^19;     % Number of signal points ADC
BW  = 10e6;     % Signal bandwidth
Flo = 240e6;    % Central frequency      (divided by 10 to reduce the simulation time)

basebandSamplingRate          = 30e6;     % sampling frequency ADC
OSR_ct_emulation              = fix((19.98e9/10)/basebandSamplingRate); % sampling frequency Simulation (divided by 10 to reduce the simulation time)
continuousTimeSamplingRate    = OSR_ct_emulation*basebandSamplingRate;
Ts_Cont                       = 1/continuousTimeSamplingRate; % sampling period Simulation


N_Sim   = floor(N*continuousTimeSamplingRate/basebandSamplingRate); % Number of  Simulation points
t_Sim   = 0:Ts_Cont:(N_Sim-1)*Ts_Cont; % Time vector Simulation
f_Sim   = 0:continuousTimeSamplingRate/N_Sim:continuousTimeSamplingRate/2-continuousTimeSamplingRate/N_Sim; % Frequency vector Simulation



%% %%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[InputAudio,fsAudio]    = audioread('lab_data/out.wav');
InputAudio              = InputAudio(1:N,1)/max(InputAudio(:,1));





%%%%%%%%%%%%%%%%% Transmistter Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DAC Specifications %%%
nBitDAC   = 12;
Vref_DAC  = 1;
dacType   = 'zoh';

%% %%%%% Filter %%%%%%%%%%%%
% Baseband FAKE Analog filter
TXBB_Filt_NF    = 3;    %(in dB)

% https://fr.mathworks.com/help/signal/ref/firpmord.html
TXBB_Filt_rp    = 0.1;         % Passband ripple in dB
TXBB_Filt_rs    = 40;          % Stopband ripple in dB
TXBB_Filt_fs    = continuousTimeSamplingRate; % Sampling frequency (Hz)
TXBB_Filt_f     = [15 20]*1e6;  % Cutoff frequencies (Hz)
TXBB_Filt_a     = [1 0];        % Desired amplitudes
% Convert the deviations to linear units. 
TXBB_Filt_dev   = [(10^(TXBB_Filt_rp/20)-1)/(10^(TXBB_Filt_rp/20)+1) 10^(-TXBB_Filt_rs/20)];
% Design the filter
[TXBB_Filt_n,TXBB_Filt_fo,TXBB_Filt_ao,TXBB_Filt_w] ...
                = firpmord(TXBB_Filt_f,TXBB_Filt_a,TXBB_Filt_dev,TXBB_Filt_fs);
disp('Designing the TX filter - This takes a while...')
TXBB_Filt       = firpm(TXBB_Filt_n,TXBB_Filt_fo,TXBB_Filt_ao,TXBB_Filt_w);

%% %%%%%%%%%% PA %%%%%%%%%%
% The power amplifier is chosen with amplifier_index
amplifier_index = 3;
[PA_Gain, PA_1dB_out, PA_IIP3, PA_NF, PA_PowerConsumption,Pmax_in] = choose_PA(amplifier_index);
 

%% %%%%%%%%%%%%%% Channel Attenuation %%%%%%%%%%%%

Distance = 1.4; % distance en mètres
Channel_Attenuation       = 20*log10(4*pi*Distance*Flo*10/3e8); % atténuation en espace libre
%Channel_Attenuation       = 80;
Channel_Attenuation_lin   = 10^(-Channel_Attenuation/20);
AntennaNoise              =randn(1,N_Sim)*sqrt(K*T*continuousTimeSamplingRate/2*Rin);






%% %%%%%%%%%%%%%%% Receiver Specs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % LNA %%%

LNA_Gain = 13.4;     % (dB)
LNA_IIP3 = 10;    % (dBm)
LNA_NF   = 4;      % (dB)



%% %%%%% Filter %%%%%%%%%%%%
%%% Baseband fake Analog filter %%%
RXBB_Filt_NF    = 3;     %(in dB)

RXBB_Filt_rp    = 0.1;         % Passband ripple in dB
RXBB_Filt_rs    = 40;          % Stopband ripple in dB
RXBB_Filt_fs    = continuousTimeSamplingRate; % Sampling frequency (Hz)
RXBB_Filt_f     = [15 20]*1e6;  % Cutoff frequencies (Hz)
RXBB_Filt_a     = [1 0];        % Desired amplitudes
% Convert the deviations to linear units. 
RXBB_Filt_dev   = [(10^(RXBB_Filt_rp/20)-1)/(10^(RXBB_Filt_rp/20)+1) 10^(-RXBB_Filt_rs/20)];
% Design the filter
[RXBB_Filt_n,RXBB_Filt_fo,RXBB_Filt_ao,RXBB_Filt_w] ...
                = firpmord(RXBB_Filt_f,RXBB_Filt_a,RXBB_Filt_dev,RXBB_Filt_fs);
disp('Designing the RX filter - This takes a while...')
RXBB_Filt       = firpm(RXBB_Filt_n,RXBB_Filt_fo,RXBB_Filt_ao,RXBB_Filt_w);

%% %%% ADC %%%%%%%
nBitADC   = 13;
Vref_ADC  = 1;
delay     = 0;
BB_gain   = 20;



TxOut       = TX_proj(  InputAudio,...
                        InputAudio,...
                        Vref_DAC,...
                        nBitDAC,...
                        basebandSamplingRate,...
                        continuousTimeSamplingRate,...
                        dacType,...
                        TXBB_Filt,...
                        TXBB_Filt_NF,...
                        Flo,...
                        PA_IIP3,...
                        PA_NF,...
                        PA_Gain,...
                        PA_PowerConsumption);

rxSignal    = TxOut*Channel_Attenuation_lin+AntennaNoise';

[basebandAnalog_adc_I, basebandAnalog_adc_Q] = RX_proj( rxSignal,...
                                                        LNA_IIP3,...
                                                        LNA_NF,...
                                                        LNA_Gain,...
                                                        RXBB_Filt,...
                                                        RXBB_Filt_NF,...
                                                        Flo,...
                                                        continuousTimeSamplingRate,...
                                                        basebandSamplingRate,...
                                                        nBitADC,...
                                                        Vref_ADC,...
                                                        delay,...
                                                        BB_gain);     

%% Plot des spectres en entrée et en sortie du système

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





