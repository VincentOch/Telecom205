% Modéliser le système avec une fibre de 100km de long. Montrer que le BER
% augmente après une certains longueur L_max, même en maintenant une
% puissance constante au niveau du récepteur.
tic
clear
close all

%% Paramètres de simulation
% Specs du signal transmis
nbits_X = 1e6;
Fsymb = [2.5e9 10e9];

% Specs de l'émetteur
power_min_dbm = -7;
power_max_dbm = 3;
sig_bb_ook = round(rand(1, nbits_X));
% Specs de la fibre
L = 100; % Longeur de la fibre en km
alpha = 0.2; % atténuation en dB/km
D = 17;  % Dispersion en ps/nm/km
S = 0.09;  % Pente de dispersion en ps/nm^2/km
lambda = 1550; % nm
c=3.e8;

% Nombre de points de puissance
nb_pow = 20;
%% émission

power_dbm = linspace(power_min_dbm, power_max_dbm, nb_pow);
nb_errs = zeros(1, length(power_dbm));
ratio_ber = zeros(1, length(power_dbm));

for kk=1:2
    Tsymb = 1/(Fsymb(kk));
    
    
    
    parfor ii=1:length(power_dbm)
        fprintf("Calcul pour puissance %0.2f\n", power_dbm(ii))
        % Paramètres RX/TX
        params_tx = make_emlaser('modulation', 'I', 'P_opt_dBm', power_dbm(ii));
        params_rx = make_photodetector('B_e', Fsymb(kk));
        
        % émission laser avec modulation externe
        [sig_opt_tx, Ts_opt, power_consumed_tx] = TX_optical_eml(sig_bb_ook, Tsymb, params_tx);
        
        %% propagation
        
        % pas de préam optique => bruit nul
        noise_psd = 0; % Watts
        
        % propagation dans la fibre
        sig_fiber_out = propag_optical(alpha,D,S,lambda, L, 2*Fsymb(kk),sig_opt_tx);
    
        %% Réception
        [sig_rx, Ts_out_rx, power_consumed_rx, snr_rx]= RX_photodetector(sig_fiber_out, Ts_opt, noise_psd, params_rx);
        
        % Récupération des bits : démodulation à seuil (seuil = moyenne)
        thres = mean(sig_rx);
        sig_rx_demod = double(sig_rx>thres);
        
        % Calcul du BER
        [nb_errs(ii), ratio_ber(ii)] = biterr(sig_bb_ook, sig_rx_demod);
    end
    semilogy(power_dbm - 3 - 20, ratio_ber, 'DisplayName', sprintf('%0.1f Gbit/s',Fsymb(kk)*1e-9),'LineWidth',1.5)

    %% plots
    if kk==1
        hold on
    end

       

end
hold off
title("Taux d'erreur en fonction de la puissance optique reçue")
xlabel('Puissance optique reçue (dBm)')
ylabel('BER')
legend('Location','best')
grid on
% subplot(2, 1, 2)
% plot(power_dbm - 3, nb_errs, 'DisplayName', 'Nombre d''erreurs')
% title("Nombre d'erreurs en fonction de la puissance émise (10 GHz)")
% xlabel('Puissance du laser (dBm)')
% ylabel('Nombre d''erreurs')
% grid on
% 
% saveas(fig,'output/ber_vs_power_100km_10.png')
% saveas(fig,'output/ber_vs_power_100km_10.fig')


toc