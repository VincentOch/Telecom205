% This script demonstrates how to compose the begin section of the TX chain:
% - DAC
% - Fake Analog Filter
% - Mixer
% The signal is sampled at a baseband rate and then converted to a continuous time signal
% The signal is a sine wave without noise
% 
% This is the TX_Chain_woPA.m version  for TELECOM205

% Rev: March 2023, Germain
% Rev: Dec. 2023, Germain


close all;
clear;


% Simulation parameters
continuousTimeSamplingRate    = 19.98e9; % Operating Frequency to emulate the behavior of a continuous time system
Ts_Cont                       = 1/continuousTimeSamplingRate; % Continuous time sampling period
%N_Cont                       = 2^18;     % Number of signal points (@continuous time rate)
N                             = 2^13;     % Number of signal points (@baseband rate)

% System properties
BW  = 10e6;     % Signal bandwidth
K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
Rin = 50;       % Matching impedance chosen equal to 50

% RF parameters
Flo       = 2.4e9;
BW_rf     = 2*BW;

% DAC Specifications
Nbits_DAC   = 18;
Vref_DAC    = 1;
Fs_DAC      = 30e6; % DAC sampling frequency aka basebandSamplingRate
mode        = 'zoh';

%%% Baseband FAKE Analog filter %%%
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

% Simulation subproperties
Subsamp_fac   = round(continuousTimeSamplingRate/Fs_DAC);
N_Cont        = N*Subsamp_fac; % Number of signal points (@continuous time rate)
t_Cont        = (0:(N_Cont-1))*Ts_Cont; % Time vector (@continuous time rate)

Ts_BB         = 1/Fs_DAC;               % Baseband sampling period
t_BB          = 0:Ts_BB:(N-1)*Ts_BB;    % Baseband Time vector

% Input signal
fin_or          = 1.1e6; % Input sine frequency
Bin_in          = round(fin_or/Fs_DAC*N); % Determining the input bin
fin             = Bin_in*Fs_DAC/N;

Ain_I = 1; Ain_Q = 1;


%%% TX Section %%%
% Signal generation
phi_0             = rand()*pi;
Input_I           = Ain_I*sin(2*pi*fin*t_BB+phi_0)'; % Input signal I Channel
Input_Q           = Ain_Q*sin(2*pi*fin*t_BB+phi_0-pi/2)'; % Input signal Q Channel

% Operate DACs
DACOutputZOH_I    = DAC(Input_I,Nbits_DAC,Vref_DAC,mode,Fs_DAC,continuousTimeSamplingRate);
DACOutputZOH_Q    = DAC(Input_Q,Nbits_DAC,Vref_DAC,mode,Fs_DAC,continuousTimeSamplingRate);

% Filter DACs output
disp('Filtering the DAC output - This takes a while...')
Filtered_output_I = basebandAnalogFiltFake(DACOutputZOH_I,TXBB_Filt,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);
Filtered_output_Q = basebandAnalogFiltFake(DACOutputZOH_Q,TXBB_Filt,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);

% Mix IQ signals to RF
MixerOutput = upMixer(Filtered_output_I,Filtered_output_Q,Flo,continuousTimeSamplingRate);





%%% Plotting %%%
voltsq2mwatt    = 1e3/Rin; % Conversion factor from V^2 to milliWatt
window_number   = 1;
lineSpec_index  = 1;

plot_spectrum(MixerOutput*sqrt(voltsq2mwatt),window_number,continuousTimeSamplingRate,lineSpec_index);
xlabel('frequency (Hz)')
ylabel('PSD (dBm/bin)')
title('Mixer Output')

plot_spectrum(MixerOutput*sqrt(voltsq2mwatt),window_number+1,continuousTimeSamplingRate,lineSpec_index);
xlabel('frequency (Hz)')
ylabel('PSD (dBm/bin)')
axis([Flo-5*BW Flo+5*BW,-150 30])
title('Zoom Mixer Output')





