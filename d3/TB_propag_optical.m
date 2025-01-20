clear
close all

%% Paramètres de simulation
Nb_points = 1000;
L = 10;  % Fiber length in km
alpha = 0.2; % Attenuation in dB/km
D = 17;  % Dispersion parameter in ps/nm/km
S = 0.09;  % Dispersion slope in ps/nm^2/km
c = 3.0e8;
T = sqrt(80*L*D*1e-3*(1550e-9)^2/c); %Pour un sinus cardinal
lambda = 1550;
Fsymb = 1/T;

%% Generation du signal
Imp_a = 111; %début de l'impulsion
input_signal = zeros(1,Nb_points);
input_signal(Imp_a:Imp_a+1) = 0.5 ;

% Paramétrage de l'émetteur
s = make_emlaser('P_opt_dBm',2);
[sig_mod, Ts, ~] = TX_optical_eml(input_signal,1/Fsymb,s);

%% Signal at the end of the optical fiber
sig_out = propag_optical(alpha,D,S,lambda,L, 2*Fsymb,sig_mod);

%% PLot
t = linspace(0, Nb_points/Fsymb, Nb_points);

fig = figure();
plot(t,10*log10(abs(sig_mod).^2), 'LineWidth',2, 'DisplayName',"Power of the input signal")
hold on
plot(t,10*log10(abs(sig_out).^2), 'LineWidth',2, 'DisplayName',"Power at the end of the optic fiber")
hold off
xlabel("time (s)")
ylabel("power (dB)")
legend()
title('Power of the signal before and after the optical fiber')
grid on

xlim([t(Imp_a-10) t(Imp_a+10)])


set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)