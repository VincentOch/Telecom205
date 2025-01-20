% This script demonstrates the use of the DAC function
% - with a sampled sinusoidal input signal (@full scale)

% Rev: March 2023, Germain

close all;
clear;

% Simulation parameters
continuousTimeSamplingRate = 19.98e9; % Operating Frequency to emulate the behavior of a continuous time system
Ts_Cont                    = 1/continuousTimeSamplingRate; % Continuous time sampling period
N                          = 2^13;     % Number of signal points (@baseband rate)

% System properties
BW  = 10e6;     % Signal bandwidth (UNILATERAL)
K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
Rin = 50;       % Matching impedance chosen equal to 50

% DAC Specifications
Nbits_DAC               = 18;
Vref_DAC                = 1;
Fs_DAC                  = 30e6;
basebandSamplingRate    = Fs_DAC;
Ts_BB                   = 1/Fs_DAC; % Baseband sampling period

% Simulation subproperties
Subsamp_fdac = round(continuousTimeSamplingRate/basebandSamplingRate);
t_Cont       = (0:(N*Subsamp_fdac-1))*Ts_Cont; % Time vector (@continuous time rate)
t_BB         = (0:(N-1))*Ts_BB;                % Time vector (@baseband rate)


%%% Input Signal Specifications %%%%
fin_or    = 1e6;    % Input sine frequency
Bin_in    = round(fin_or/basebandSamplingRate*N); % Determining the input bin
fin       = Bin_in*Fs_DAC/N;
Ain       = 1;
% Input signal
Input     = Ain*sin(2*pi*fin*t_BB+rand()*2*pi);

mode              = 'impulse'; % "Zero Padding"
AnalogOutputZP    = DAC(Input',Nbits_DAC,Vref_DAC,mode,Fs_DAC,continuousTimeSamplingRate);

n_samples_plot=100;
figure(1)
stem(t_BB(1:n_samples_plot),Input(1:n_samples_plot));
hold all
plot(t_Cont(1:n_samples_plot*Subsamp_fdac),AnalogOutputZP(1:n_samples_plot*Subsamp_fdac),'linewidth',2);
hold off
xlabel('Time (s)')
ylabel('Amplitude (V)')
legend('Input (DT)','Output (emulated CT)')
title('DAC signal')








