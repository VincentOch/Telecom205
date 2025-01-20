clear;
close all;

Nsim    = 2^14;                         % Number of points of the simulation
Fs      = 30e6;                         % Sampling Frequency
Ts      = 1/Fs;                         % Sampling Period
t_sim   = 0:Ts:(Nsim-1)*Ts;             % Time Vector
t_sim   = reshape(t_sim,Nsim,1);        % Reshape Row vector to Column vector

fsig_or   = 1e6;                        % Input frequency (Hz)
sig_bin   = round(fsig_or/Fs*Nsim);     % Input bin of the signal frequency
fsig      = sig_bin*Fs/Nsim;            % re-adjusting in a signal bin

Amp   = 1;                              % Signal Amplitude
x     = Amp*sin(2*pi*fsig*t_sim);       % Input signal generation

sigma_noise   = 0.01;                   % Sigma of the noise (Ecart type)
x_noi         = x+sigma_noise*randn(size(t_sim));  % Adding noise to the signal
error         = x_noi - x;              % Calculating the error

PS_theo   = Amp.^2/2;                   % Theoretical Average Signal power
PN_theo   = sigma_noise.^2;             % Theoretical error or noise power
SNR_theo  = 10*log10(PS_theo/PN_theo);  % Theoretical SNR
disp(['The theoretical SNR is ', num2str(SNR_theo), ' dB'])

PS_time   = mean(x.^2);                 % Empirical time domain Signal power 
PN_time   = mean(error.^2);             % Empirical time domain error or noise power 
SNR_time  = 10*log10(PS_time/PN_time);  % Empirical time domain SNR 
disp(['The time domain SNR is ', num2str(SNR_time), ' dB'])

Nx_noi    = length(x_noi);
win       = blackman(Nx_noi,'periodic');   % Calculating the window
x_noiPSD  = abs(fft(x_noi(:).*win(:))).^2; % Calculating the PSD other windowed signal
x_noiPSD  = x_noiPSD/Nx_noi;               % Normalizing with respect to the length

figure()
plot(10*log10(x_noiPSD))
xlabel('bin index')
ylabel('PSD(dB/bin)') 
title('Plot of the raw FFT');

figure()
subplot(2,1,1)
plot(10*log10(x_noiPSD(1:Nsim/2)))
xlabel('bin index')
ylabel('PSD(dB/bin)') 
title('Plot of the FFT (positive frequencies only)');

subplot(2,1,2)
f_step  = Fs/Nsim;
f       = 0:f_step:(Nsim/2-1)*f_step;
plot(f/1e6,10*log10(x_noiPSD(1:Nsim/2)))
xlabel('Frequency (MHz)')
ylabel('PSD(dB/bin)') 


err_bin   = [1 , Nsim/2];   % beggining and end of the useful band (here, Nyquist band)
nintsig   = 2;              % Number of points around the useful signal to integrate
PS_freq   = sum(x_noiPSD(sig_bin+1-nintsig:sig_bin+1+nintsig)); % Empirical frequency domain Signal power 
PN_freq   = sum(x_noiPSD(1:Nsim/2))-PS_freq;                    % Empirical frequency domain Noise or error power
SNR_freq  = 10*log10(PS_freq/PN_freq);
disp(['The frequency domain SNR is ', num2str(SNR_freq), ' dB'])
