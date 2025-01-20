% This script demonstrates how to compose the begin section of the TX chain:
% - DAC
% - Analog Filter
% - Mixer
% The signal is sampled at a baseband rate and then converted to a continuous time signal
% The signal is a sine wave without noise

% Rev: March 2023, Germain


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

% Filter Specifications (Example)
TXBB_Filt_NF    = 3;          % Filter Noise Figure Filter
TXBB_Filt_Fcut  = 11.2479e6;  % Filter TXBB_Fil_Fcut 3dB Frequency
TXBB_Filt_Order = 6;          % Filter TXBB_Fil_Order

% Instanciate filter (due to numerical issue, the filter has to be instanciated as SOS)
[TXBB_Filt_z,TXBB_Filt_p,TXBB_Filt_k] = butter(TXBB_Filt_Order,TXBB_Filt_Fcut/(continuousTimeSamplingRate/2));
TXBB_Filt_sos = zp2sos(TXBB_Filt_z,TXBB_Filt_p,TXBB_Filt_k);

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
Filtered_output_I = basebandAnalogFilt(DACOutputZOH_I,TXBB_Filt_sos,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);
Filtered_output_Q = basebandAnalogFilt(DACOutputZOH_Q,TXBB_Filt_sos,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);

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





