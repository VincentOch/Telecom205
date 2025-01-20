% Montrer que le BER augmente après une certains longueur L_max, même en maintenant une
% puissance constante au niveau du récepteur.

tic
clear
close all

%% Paramètres de simulation
% Specs du signal transmis
nbits_X = 1e4;
Fsymb = 10e9;

% Specs de l'émetteur
L_min_km = 0;
L_max_km = 100;
nb_dist = 200;
Tsymb = 1/Fsymb;

% Specs de la fibre
alpha = 0.2; % atténuation en dB/km
D = 17;  % Dispersion en ps/nm/km
S = 0.09;  % Pente de dispersion en ps/nm^2/km
lambda = 1550; % nm
c=3.e8;

% Nombre de points de puissance
%% émission
sig_bb_ook = round(rand(1, nbits_X));
L_km = linspace(L_min_km, L_max_km, nb_dist);
power_dbm = -17 + alpha.*L_km;  % pré-calcul de la puissance à émettre pour avoir une puissance reçue constante
nb_errs = zeros(1, length(L_km));
ratio_ber = zeros(1, length(L_km));
moy = zeros(1, length(L_km));

for uu=1:10
    disp(uu)
    parfor ii=1:length(power_dbm)
        % fprintf("Calcul pour distance %0.2f\n", L_km(ii))
        % Paramètres RX/TX
        params_tx = make_emlaser('modulation', 'I', 'P_opt_dBm', power_dbm(ii));
        params_rx = make_photodetector('B_e', Fsymb);
    
        % émission laser avec modulation externe
        [sig_opt_tx, Ts_opt, power_consumed_tx] = TX_optical_eml(sig_bb_ook, Tsymb, params_tx);
    
        %% propagation
    
        % pas de préam optique => bruit nul
        noise_psd = 0; % Watts
    
        % propagation dans la fibre
        sig_fiber_out = propag_optical(alpha,D,S,lambda,L_km(ii), 2*Fsymb,sig_opt_tx);
    
        %% Réception
        [sig_rx, Ts_out_rx, power_consumed_rx, snr_rx]= RX_photodetector(sig_fiber_out, Ts_opt, noise_psd, params_rx);
        % disp(snr_rx);
        % Récupération des bits : démodulation à seuil (seuil = moyenne)
        thres = mean(sig_rx);
        sig_rx_demod = double(sig_rx>thres);
    
        % Calcul du BER
        [nb_errs(ii), ratio_ber(ii)] = biterr(sig_bb_ook, sig_rx_demod);
    end
    moy = moy+ratio_ber/100;
end
%% plots

% moy = zeros(length(ratio_ber));
% for jj=1:length(ratio_ber)
%     a = max(1,jj-2);
%     b = min(jj+2, length(ratio_ber));
%     moy(jj) = mean(ratio_ber(a:b));
% end


plot(L_km, moy, 'DisplayName', 'BER en fonction de la distance', LineWidth=1.5)
yline(10^-3, LineWidth=1.5, Color = [1 0 0], DisplayName='BER de 10^{-3}')
title("Calcul du taux d'erreur en fonction de la distance de propagation, à puissance reçue constante, 10Gbit/s")
xlabel('distance L (km)')
ylabel('BER')
grid on
legend(Location="best")
set(gca, 'FontSize', 18)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 19)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

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