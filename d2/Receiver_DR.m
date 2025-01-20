% This script demonstrates how to use the : 
% - rfLNA() function
% - downMixer() function
% - filter() function (due to numerical issue, sosfilt() is used instead of filter())
% - ADC() function
% The input signal to the RX system is a thermally noised RF dualtone signal. 

% Rev: March 2022, Chadi+Germain
% Rev: March 2023, Germain
% Rev: Dec 2023, Germain

close all;
clear;

% Simulation parameters
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
Fs_ADC       = 30e6;       % Sampling frequency ADC
Vref         = 1;          % Reference voltage of the ADC
Nbits_ADC    = 10;         % Number of bits for the ADC

% Simulation subproperties
N_Cont = N*continuousTimeSamplingRate/Fs_ADC; % Number of signal points (@continuous time rate)
t_Cont = (0:(N_Cont-1))*Ts_Cont; % Time vector (@continuous time rate)
Ts_ADC = 1/Fs_ADC;                                % Sampling period ADC
t_ADC  = 0:Ts_ADC:(N-1)*Ts_ADC;                   % Time vector Simulation


% Input signal
fin_or       = [1.1e6,2.1e6];          % Input sine frequency
Bin_in       = round(fin_or/Fs_ADC*N); % Determining the input bin
fin          = Bin_in*Fs_ADC/N;        % 

Pin             = -40;  % Pin in dBm
Ain             = sqrt(10.^((Pin-30)/10)*2*Rin);
AntennaNoise    = randn(1,N_Cont)*sqrt(K*T*BW*Rin);

% Generate dualtone signal
in = 0;
for i = 1:length(fin)
  in = in+Ain*sin(2*pi*(Flo+fin(i))*t_Cont+rand()*2*pi); % Input signal
end
in        = in+AntennaNoise;

% LNA parameters
G_LNA        = 30;      % Gain of the first stage in dB
NF_LNA       = 4;       % Noise figure of the first stage in dB
IIP3_LNA     = -20;     % Set very high for Simplicity

% Baseband analog filter example 
% (cut-off frequency is "randomly" chosen just for the example)
Filter_order              = 7;
cutoff_freq               = 10*BW;
[filter_z,filter_p,filter_k] = butter(Filter_order,cutoff_freq/(continuousTimeSamplingRate/2));
% Instanciate filter (due to numerical issue, the filter has to be instanciated as SOS)
Filt_sos = zp2sos(filter_z,filter_p,filter_k);

% Please consider using FIRPM instead of BUTTER for results consistency

%%%% RX section %%%%

% pass the signal through the LNA
out_LNA     = rfLNA(in,G_LNA,NF_LNA,IIP3_LNA,Rin,continuousTimeSamplingRate/2); % output Signal of the amplifier

% downmix the signal
[outI,outQ] = downMixer(out_LNA,Flo,continuousTimeSamplingRate);

% Baseband low pass filter example (due to numerical issue, use sosfilt() instead of filter())
Out_Filter  = sosfilt(Filt_sos,outI);

% ADC
ADCdelay = 0;
Out_ADC  = ADC(Out_Filter,Nbits_ADC,Vref,Fs_ADC,ADCdelay,continuousTimeSamplingRate);


%%% Plotting %%%
voltsq2mwatt        = 1e3/Rin; % Conversion factor from V^2 to milliWatt
window_number       = 1;
lineSpec_index      = 1;

subplot(2,2,1)
plot_spectrum(out_LNA*sqrt(voltsq2mwatt),window_number,   ...
                  continuousTimeSamplingRate,lineSpec_index);
title('Spectrum output LNA')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')

subplot(2,2,2)
plot_spectrum(outI*sqrt(voltsq2mwatt),window_number,      ...
                  continuousTimeSamplingRate,lineSpec_index+1);
title('Spectrum output mixer')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')

subplot(2,2,3)
plot_spectrum(Out_Filter*sqrt(voltsq2mwatt),window_number,...
                  continuousTimeSamplingRate,lineSpec_index+2);
title('Spectrum output Filter')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')

subplot(2,2,4)
plot_spectrum(Out_ADC*sqrt(voltsq2mwatt),window_number,   ...
                  Fs_ADC,lineSpec_index+3);
title('Spectrum output ADC')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')


