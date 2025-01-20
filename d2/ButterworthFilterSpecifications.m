close all
clear

continuousTimeSamplingRate = 900e6; %Frequency to emulate the behavior of a continuous time system
N       = 1e5;   % number of points
f1      = 10e6;  % Low frequency - Useful Bandwidth limit
f2      = 20e6;  % High Frequency - Stop Frequency
Amin    = 30;    % Min attenuation at f2 in dB
Amax    = 1;     % Max attenuation at f1 in dB

% Compute filter order
[Order,Omega_cut]   = buttord(2*pi*f1,2*pi*f2,Amax,Amin,'s');

% Compute filter coefficients
[b_butter,a_butter] = butter(Order,Omega_cut,'s');

% Compute frequency response
[butter_f_resp,butter_f_resp_omega] = freqs(b_butter,a_butter,N);

% Plot frequency response
semilogx(butter_f_resp_omega/(2*pi)/1e6,20*log10(abs(butter_f_resp)))
xlabel("Frequency (MHz)")
ylabel("Gain(dB)")
% axis([0.1,200,-200,1])
xlim([min(get(gca,'xlim')) 50])

disp(['The filter order is ',int2str(Order),' and the Fcut is ', num2str(Omega_cut/(2*pi)/1e6),' MHz'])
