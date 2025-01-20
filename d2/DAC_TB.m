% This script demonstrates the use of the DAC function
% - with a random colored noise input signal (set to full scale)
% - in both cases : ZOH and "impulse" mode

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

%%%% Input  Generation %%%%
% Let's generate a low pass noise signal
InputOrigin   = randn(1,N);             % Generate white noise
CutFil        = fir1(100,BW/(Fs_DAC/2));% Generate Low pass filter
Input         = filter(CutFil,1,InputOrigin); % Filter noise to make it band limited
Input         = Input/max(Input);       % Normalize

%%%%%%%%%%%%% Zero Order Hold %%%%%%%%%%%%%%%%%%%
mode              = 'zoh';
AnalogOutputZOH   = DAC(Input',Nbits_DAC,Vref_DAC,mode,Fs_DAC,continuousTimeSamplingRate);

%%%%%%%%%%%%% Zero Padding %%%%%%%%%%%%%%%%%%%
mode              = 'impulse';
AnalogOutputZP    = DAC(Input',Nbits_DAC,Vref_DAC,mode,Fs_DAC,continuousTimeSamplingRate);


%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
% Some custom parameters for the home made plot_spectrum function
window_plot_number = 1;
lineSpec_index     = 1;

figure(1)
subplot(2,2,1)
plot_spectrum(Input,window_plot_number, Fs_DAC,lineSpec_index);
xlabel('frequency (Hz)')
ylabel('PSD (dB/bin)')
title('Input Signal')
axis([0,100e6,-200,1])

subplot(2,2,2)
plot_spectrum(AnalogOutputZOH,window_plot_number,continuousTimeSamplingRate,lineSpec_index+1);  
xlabel('frequency (Hz)')
ylabel('PSD (dB/bin)')
title('DAC Output Zero Holder')
axis([0,100e6,-200,1])

subplot(2,2,3)
plot_spectrum(AnalogOutputZP,window_plot_number,continuousTimeSamplingRate,lineSpec_index+2);
xlabel('frequency (Hz)')
ylabel('PSD (dB/bin)')
title('DAC Output  Zero Padding')
axis([0,100e6,-200,1])

subplot(2,2,4)
plot_spectrum(AnalogOutputZP*(continuousTimeSamplingRate/Fs_DAC),window_plot_number,continuousTimeSamplingRate,lineSpec_index+2); 
hold on
plot_spectrum(AnalogOutputZOH,window_plot_number,continuousTimeSamplingRate,lineSpec_index+1);
plot_spectrum(Input,window_plot_number, Fs_DAC,lineSpec_index);
xlabel('frequency (Hz)')
ylabel('PSD (dB/bin)')
title('All signals (normalized)')
axis([0,100e6,-200,1])





