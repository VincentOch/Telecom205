% Simulation de propagation à 100km pour 10GHz avec diverses corrections
% pour assurer un BER en dessous de 10e-3
tic
clear
close all

%% Paramètres de simulation
% Specs du signal transmis
nbits_X = 1e6;
Fsymb = 10e9;

% Specs de l'émetteur
L = 100; %km
Tsymb = 1/Fsymb;

% Specs de la fibre
alpha = 0.2; % atténuation en dB/km
D = 17;  % Dispersion en ps/nm/km
S = 0.09;  % Pente de dispersion en ps/nm^2/km
lambda = 1550; % nm
c=3.e8;

% Nombre de points de puissance
%% émission
sig_bb_ook = round(rand(1, nbits_X)); % signal ook random dans {0;1}
power_dbm = 6;
consumption = 0;

% Paramètres RX/TX
% Modulation I (simple)
params_tx = make_emlaser('modulation', 'I', 'P_opt_dBm', power_dbm);
params_rx = make_photodetector('B_e', Fsymb); % adaptation BP photodiode

% émission laser avec modulation externe
[sig_opt_tx, Ts_opt, power_consumed_tx] = TX_optical_eml(sig_bb_ook, Tsymb, params_tx);
consumption = consumption+power_consumed_tx;

%% propagation

% pas de préam optique => bruit nul
noise_psd = 0; % Watts

% propagation dans la fibre
sig_fiber_out = propag_optical(alpha,D,S,lambda,L, 2*Fsymb,sig_opt_tx); %Fibre SSFM
% sig_fiber_out = propag_optical(0.5,-80,-0.5,lambda,17.53, 2*Fsymb,sig_fiber_out); %Fibre DCM

%% Réception
[sig_rx, Ts_out_rx, power_consumed_rx, snr_rx]= RX_photodetector(sig_fiber_out, Ts_opt, noise_psd, params_rx);
consumption = consumption+power_consumed_rx;

% Récupération des bits : démodulation à seuil (seuil = moyenne)
thres = mean(sig_rx);
sig_rx_demod = double(sig_rx>thres);

% Calcul du BER
[nb_errs, ratio_ber] = biterr(sig_bb_ook, sig_rx_demod);
fprintf("BER : %0.2f\n",ratio_ber)
fprintf("total consumption : %0.2f\n",consumption)

%% pour le plot (inutilisé à la fin)

% figure();
% subplot(2, 1, 1)
% semilogy(L_km, ratio_ber, 'DisplayName', 'BER (%)')
% title("Calcul du taux d'erreur en fonction de la distance de propagation, à puissance reçue constante")
% xlabel('distance L (km)')
% ylabel('BER')
% grid on
% 
% subplot(2, 1, 2)
% plot(L_km, nb_errs, 'DisplayName', 'Nombre d''erreurs')
% title("Nombre d'erreurs en fonction de la distance de propagation, à puissance reçue constante")
% xlabel('distance L (km)')
% grid on
% ylabel('Nombre d''erreurs')
% grid on

% saveas(fig,'output/ber_vs_power_100km_10G.png')
% saveas(fig,'output/ber_vs_power_100km_10G.fig')

toc