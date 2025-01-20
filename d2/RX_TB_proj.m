% This script demonstrates how to compose a complete RX chain using
% - LNA
% - Mixer
% - Fake Analog filter
% - ADC
% The RF signal is a thermally noised sine wave
% 
% This is the RX_TB.m version  for TELECOM205

% Rev: March 2023, Germain
% Rev: Dec. 2023, Germain

close all;
clear all;

%% Simulation parameters
continuousTimeSamplingRate = 19.98e9; % Operating Frequency to emulate the behavior of a continuous time system
Ts_Cont                    = 1/continuousTimeSamplingRate; % Continuous time sampling period
%N_Cont                     = 2^18;   % Number of signal points 
N                          = 2^13;       % Number of signal points (@ADC rate)

% System properties
BW  = 10e6;     % Signal bandwidth
K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
Rin = 50;       % Matching impedance chosen equal to 50

% RF parameters
Flo       = 2.4e9;
BW_rf     = 2*BW;

% ADC parameters
Fs_ADC        = 30e6;       % Sampling frequency ADC a.k.a basebandSamplingRate
Vref          = 1;          % Reference voltage of the ADC
Nbits_ADC     = 13;         % Number of bits for the ADC
BB_gain       = 20;         % dB
BB_gain_lin   = 10^(BB_gain/20); % in linear

% LNA
LNA_Gain = 13.4;     % (dB)
LNA_IIP3 = 10;    % (dBm)
LNA_NF   = ;    % (dB)

%% % Baseband fake Analog filter %%%
RXBB_Filt_NF    = 0;     %(in dB)

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

%% Simulation subproperties
N_Cont = N*continuousTimeSamplingRate/Fs_ADC; % Number of signal points (@continuous time rate)
t_Cont = (0:(N_Cont-1))*Ts_Cont; % Time vector (@continuous time rate)
Ts_ADC = 1/Fs_ADC;                                % Sampling period ADC
t_ADC  = 0:Ts_ADC:(N-1)*Ts_ADC;                   % Time vector Simulation

%% Input signal
fin_or       = 1.1e6;          % Input sine frequency
Bin_in       = round(fin_or/Fs_ADC*N); % Determining the input bin
fin          = Bin_in*Fs_ADC/N;        % 

Pin             = -90;  % Pin in dBm
Ain             = sqrt(10.^((Pin-30)/10)*2*Rin);
AntennaNoise    = randn(1,N_Cont)*sqrt(K*T*continuousTimeSamplingRate/2*Rin);
rxSignal        = Ain*sin(2*pi*(Flo+fin)*t_Cont+rand()*2*pi)+AntennaNoise; % Input signal



%%% RX section %%%%
%% Operate LNA
[rfLNASignal, LNAconsum] = rfLNA(rxSignal,LNA_Gain,LNA_NF,LNA_IIP3,Rin,continuousTimeSamplingRate/2);
fprintf("Consumption LNA = %fW\n", LNAconsum)

%% Mixing down to BB
[basebandAnalog_raw_I,basebandAnalog_raw_Q,downMixer_power] = downMixer(rfLNASignal,Flo,continuousTimeSamplingRate);
fprintf("Consumption down-mixer = %fW\n", downMixer_power)




%% Perform filtering
disp('Filtering with the RX filter - This takes a while...')
basebandAnalog_filtrx_I = basebandAnalogFiltFake(basebandAnalog_raw_I,RXBB_Filt,RXBB_Filt_NF,Rin,continuousTimeSamplingRate);
basebandAnalog_filtrx_Q = basebandAnalogFiltFake(basebandAnalog_raw_Q,RXBB_Filt,RXBB_Filt_NF,Rin,continuousTimeSamplingRate);

%% Baseband Gain
% BBamp_Gain    = 20; % (dB)
% BBamp_IIP3    = 40; % (dBm)
% BBamp_NF      = 10; % (dB)
% BBamp_band    = 10e6;% (MHz)
% [basebandAnalog_amp_I, BBampPower_I] = BBamp(basebandAnalog_filtrx_I,BBamp_Gain,BBamp_NF,BBamp_IIP3,Rin,BBamp_band,continuousTimeSamplingRate);
% [basebandAnalog_amp_Q, BBampPower_Q] = BBamp(basebandAnalog_filtrx_Q,BBamp_Gain,BBamp_NF,BBamp_IIP3,Rin,BBamp_band,continuousTimeSamplingRate);
% 
% fprintf("Consumption of Baseband Gain I = %fW and Q = %fW\n", BBampPower_I, BBampPower_Q)

%% Perform AD conversion
delay_ADC = 0; % A delay used to compensate for the delay caused by the filters
[basebandAnalog_adc_I, ADC_I_consum] = ADC(BB_gain_lin*basebandAnalog_filtrx_I,Nbits_ADC,Vref,Fs_ADC,delay_ADC,continuousTimeSamplingRate);
[basebandAnalog_adc_Q, ADC_Q_consum] = ADC(BB_gain_lin*basebandAnalog_filtrx_Q,Nbits_ADC,Vref,Fs_ADC,delay_ADC,continuousTimeSamplingRate);


fprintf("Consumption of ADC I = %fW and Q = %fW\n", ADC_I_consum, ADC_Q_consum)

%% IQ combination for complex baseband signals
basebandComplexDigital = complex(basebandAnalog_adc_I,basebandAnalog_adc_Q);

tot_consum_Rx = ADC_I_consum + ADC_Q_consum+LNAconsum+downMixer_power;
fprintf("Total consumption of the Rx = %fW",tot_consum_Rx)


%% % Plotting %%%
voltsq2mwatt    = 1e3/Rin; % Conversion factor from V^2 to milliWatt
window_number   = 1;
lineSpec_index  = 1;
fullband        = true;


subplot(2,2,1)
plot_spectrum(rfLNASignal*sqrt(voltsq2mwatt),window_number,continuousTimeSamplingRate,lineSpec_index);
title('spectrum output LNA')
xlabel('Frequency (Hz)')
ylabel('PSD (dBm/bin)')
subplot(2,2,2)

plot_spectrum(basebandAnalog_filtrx_I*sqrt(voltsq2mwatt),window_number,continuousTimeSamplingRate,lineSpec_index+1);
title('spectrum filter mixer')
xlabel('frequency(Hz)')
ylabel('PSD(dBm/bin)')
axis([0 100e6 -250 0])



subplot(2,2,3)
plot_spectrum(basebandAnalog_adc_I*sqrt(voltsq2mwatt),window_number,Fs_ADC,lineSpec_index+2);
title('spectrum ADC')
xlabel('Frequency (Hz)')
ylabel('PSD (dBm/bin)')


subplot(2,2,4)
plot_spectrum(basebandComplexDigital,window_number,Fs_ADC,lineSpec_index+3,fullband);
title('spectrum reconstructed complex output')
xlabel('Frequency (Hz)')
ylabel('PSD (dBm/bin)')


%% % SNR Computations %%%

% Compute bins with fftshifted psd (hence N/2 is the DC bin)
Bin_limits     = N/2+[round(-BW/Fs_ADC*N),round(BW/Fs_ADC*N)];
% Compute SNR
Bin_sig_shiftd = N/2+Bin_in;
bin_width      = 2;
SNR_out        = perf_estim(basebandComplexDigital,Bin_sig_shiftd,bin_width,Bin_limits,fullband);

disp(['The SNR at the output of the ADC is ',num2str(SNR_out), ' dB'])

