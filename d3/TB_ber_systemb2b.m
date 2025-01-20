% En utilisant TX (eml) et RX, étudier le BER d'un système back to back
% (sans la fibre). EML = modulation externe
% Modulation OOK, détection directe (Rb=2.5, 10Gbps)
% Trouver la puissance requise au photodétecteur pour avoir un BER
% inférieur à 10^-3 et 10^-6 (sans équalizer ni code correcteur)
tic
clear
close all

%% Paramètres de simulation
nbits_X = 1e8;
Fsymb = 2.5e9;
power_min_dbm = -23;
power_max_dbm = -20;
Tsymb = 1/Fsymb;
nb_pow = 20;

%% émission
sig_bb_ook = round(rand(1, nbits_X));

power_dbm = linspace(power_min_dbm, power_max_dbm, nb_pow);
nb_errs = zeros(1, length(power_dbm));
ratio_ber = zeros(1, length(power_dbm));

parfor ii=1:length(power_dbm)

    fprintf("Calcul pour puissance %0.2f\n", power_dbm(ii))

    % Paramètres RX/TX
    params_tx = make_emlaser('modulation', 'I', 'P_opt_dBm', power_dbm(ii));
    params_rx = make_photodetector('B_e', Fsymb);
    
    % émission laser avec modulation externe
    [sig_opt_tx, Ts_opt, power_consumed_tx] = TX_optical_eml(sig_bb_ook, Tsymb, params_tx);
    
    %% réception
    % Le système est back to back, donc pas de canal de propag optique
    
    % pas de préam optique => bruit nul
    noise_psd = 0; % Watts
    [sig_rx, Ts_out_rx, power_consumed_rx, snr_rx]= RX_photodetector(sig_opt_tx, Ts_opt, noise_psd, params_rx);
    
    % Récupération des bits : démodulation à seuil (seuil = moyenne)
    thres = mean(sig_rx);
    sig_rx_demod = double(sig_rx>thres);
    
    % Calcul du BER
    [nb_errs(ii), ratio_ber(ii)] = biterr(sig_bb_ook, sig_rx_demod);
end

%% plots
fig = figure();
subplot(2, 1, 1)
semilogy(power_dbm - 3, ratio_ber, 'DisplayName', 'BER (%)', 'LineWidth', 1.5)
title("Calcul du taux d'erreur en fonction de la puissance optique reçue")
xlabel('Puissance optique reçue (dBm)') % -3dB parce que c'est la moyenne (OOK equiprobable, la moitié des échatillons vaut 1)
ylabel('BER')
grid on

subplot(2, 1, 2)
plot(power_dbm - 3, nb_errs, 'DisplayName', 'Nombre d''erreurs', 'LineWidth', 1.5)
title("Nombre d'erreurs en fonction de la puissance optique reçue")
xlabel('Puissance optique reçue (dBm)')
grid on
ylabel('Nombre d''erreurs')
grid on

saveas(fig,'output/ber_vs_power_b2b_simple_2.5.png')
saveas(fig,'output/ber_vs_power_b2b_simple_2.5.fig')

toc
