% This script demonstrates the usage a TX function that integrates the whole TX chain
% The baseband signal is dual tone IQ signal without noise
% 
% This file has no equivalent for TELECOM201/ICS905

% Rev: March 2023, Germain
% Rev: Dec. 2023, Germain

close all;
clear;


%% Simulation parameters
continuousTimeSamplingRate    = 19.98e9; % Operating Frequency to emulate the behavior of a continuous time system
Ts_Cont                       = 1/continuousTimeSamplingRate; % Continuous time sampling period
%N_Cont                       = 2^18;     % Number of signal points (@continuous time rate)
N                             = 2^13;     % Number of signal points (@baseband rate)

%% System properties
BW  = 10e6;     % Signal bandwidth
K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
Rin = 50;       % Matching impedance chosen equal to 50

%% RF parameters
Flo       = 2.4e9;
BW_rf     = 2*BW;

%% DAC Specifications
Nbits_DAC   = 18;
Vref_DAC    = 1;
Fs_DAC      = 30e6; % DAC sampling frequency aka basebandSamplingRate
mode        = 'zoh';

%% Baseband FAKE Analog filter %%%
Filt_NF    = 3;    %(in dB)

% https://fr.mathworks.com/help/signal/ref/firpmord.html
Filt_rp    = 0.1;         % Passband ripple in dB
Filt_rs    = 40;          % Stopband ripple in dB
Filt_fs    = continuousTimeSamplingRate; % Sampling frequency (Hz)
Filt_f     = [15 20]*1e6;  % Cutoff frequencies (Hz)
Filt_a     = [1 0];        % Desired amplitudes
% Convert the deviations to linear units. 
Filt_dev   = [(10^(Filt_rp/20)-1)/(10^(Filt_rp/20)+1) 10^(-Filt_rs/20)];
% Design the filter
[Filt_n,Filt_fo,Filt_ao,Filt_w] ...
                = firpmord(Filt_f,Filt_a,Filt_dev,Filt_fs);
disp('Designing the TX filter - This takes a while...')
Filt       = firpm(Filt_n,Filt_fo,Filt_ao,Filt_w);

%% Power Amplifier Specifications
amplifier_index = 5;
[PA_Gain, PA_1dB_out, PA_IIP3, PA_NF, PA_power,Pmax_in] = choose_PA(amplifier_index);

%% Simulation subproperties
Subsamp_fac   = round(continuousTimeSamplingRate/Fs_DAC);
N_Cont        = N*Subsamp_fac;          % Number of signal points (@continuous time rate)
t_Cont        = (0:(N_Cont-1))*Ts_Cont; % Time vector (@continuous time rate)

Ts_BB         = 1/Fs_DAC;               % Baseband sampling period
t_BB          = 0:Ts_BB:(N-1)*Ts_BB;    % Baseband Time vector

%% Input signal
fin_or    = [0.5e6, 0.7e6];         % Input sine frequency
Bin_in    = round(fin_or/Fs_DAC*N); % Determining the input bin
fin       = Bin_in*Fs_DAC/N;
Ain       = 0.5;


%% Generate baseband signal
phi_0     = rand()*pi;
Input_I   = Ain*sin(2*pi*fin(1)*t_BB+phi_0)      +Ain*sin(2*pi*fin(2)*t_BB+phi_0);     % Input signal I Channel
Input_Q   = Ain*sin(2*pi*fin(1)*t_BB+phi_0-pi/2) +Ain*sin(2*pi*fin(2)*t_BB+phi_0-pi/2);% Input signal Q Channel


%% %%%%%%% TX section %%%%%%%%%
TxOut = TX_proj(Input_I',Input_Q',...            % digital baseband signal
            Vref_DAC,Nbits_DAC,Fs_DAC,continuousTimeSamplingRate,mode,...  % DA conversion
            Filt,Filt_NF,...% Baseband filter
            Flo,...         % RF carrier frequency
            PA_IIP3,PA_NF,PA_Gain,PA_power,Pmax_in);% PA parameters


%% % ACPR %%%
voltsq2mwatt    = 1e3/Rin; % Conversion factor from V^2 to milliWatt
window_number   = 1;
lineSpec_index  = 1;

[acpr_left, acpr_right] = ACPR(BW_rf ,BW_rf + 10e6,Flo,TxOut*sqrt(voltsq2mwatt),continuousTimeSamplingRate);

%% plot
plot_spectrum(TxOut*sqrt(voltsq2mwatt),2,continuousTimeSamplingRate,1);
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')
title('TX Output')
axis([Flo-BW,Flo+BW,-150,30])
