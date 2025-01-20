% En utilisant TX (dml), étudier les performances d'une modulation directe
% DML = modulation directe
% Modulation OOK, détection directe (Rb=1Gbps)
% back to back et sur 20km

tic
clear
close all

%% Paramètres de simulation
nbits_X = 10;
Fsymb = 1e9;
Tsymb = 1/Fsymb;

% Specs de la fibre
L = 20; % Longeur de la fibre en km
alpha = 0.2; % atténuation en dB/km
D = 17;  % Dispersion en ps/nm/km
S = 0.09;  % Pente de dispersion en ps/nm^2/km
lambda = 1550; % nm
c=3.e8;

%% émission
sig_bb_ook_origin = round(rand(1, nbits_X)); %[0 1 0 1 0 0 0 0 1 0]; %used for the plot in the report

% Paramètres RX/TX
params_tx = make_laser_simple('v', c/lambda,'B_e', Fsymb); % propag dans une fibre
params_rx = make_photodetector('B_e', Fsymb); % adaptation BP photodiode
Is = params_tx.I_th;

bias_currents = [0e-3 30e-3 300e-3];

[sig_test, ~, ~] = TX_optical_dml(sig_bb_ook_origin, Tsymb, params_tx);
M = length(sig_test)/length(sig_bb_ook_origin);
N_samples = length(sig_bb_ook_origin)*M;

XX = zeros(1, N_samples);
for kk=1:M
    XX(1:M:end) = sig_bb_ook_origin;
end

power_fiber = zeros(3, N_samples);
power_sig_out = zeros(3, N_samples);

for ii=1:length(bias_currents)
    IDC = bias_currents(ii);
    sig_bb_ook = sig_bb_ook_origin*10*Is+Is+IDC; %+Is pour être sûr d'être au dessus du seuil pour les symboles valant zéro. On différencie bien le 0 et le 1
    % émission laser avec modulation externe
    [sig_opt_tx, Ts_opt, power_consumed_tx] = TX_optical_dml(sig_bb_ook, Tsymb, params_tx);
 
    %Dans le cas d'un signal B2B
    sig_opt_out = sig_opt_tx;
    
    % propagation dans la fibre Dans le cas propagation 20km
    sig_opt_out_fiber = propag_optical(alpha,D,S,lambda,L,2*Fsymb,sig_opt_tx);
    %% réception
    
    % pas de préam optique => bruit nul
    noise_psd = 0; % Watts
    power_sig_out(ii,:) = 10*log10(abs(sig_opt_out).^2);
    power_fiber(ii,:) = 10*log10(abs(sig_opt_out_fiber).^2);

    [sig_rx, Ts_out_rx, power_consumed_rx, snr_rx]= RX_photodetector(sig_opt_out, Ts_opt, noise_psd, params_rx);
    [sig_rx2, Ts_out_rx2, power_consumed_rx2, snr_rx2]= RX_photodetector(sig_opt_out_fiber, Ts_opt, noise_psd, params_rx);

end

%% Plot   
fig = figure();
subplot(3, 2, 1)
plot(power_sig_out(1,:), LineWidth= 1.5)
title(sprintf("courant de biais de %.2f mA, B2B", bias_currents(1)*1e3))
ylabel('puissance (dB)')
grid on
set(gca, 'FontSize', 13)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 15)
set(get(gca,'ylabel'), 'FontSize', 15)

subplot(3, 2, 2)
plot(power_fiber(1,:), LineWidth= 1.5)
title(sprintf("courant de biais de %.2f mA, propagation de 20km dans une fibre", bias_currents(1)*1e3))
ylabel('puissance (dB)')
grid on
set(gca, 'FontSize', 13)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 15)
set(get(gca,'ylabel'), 'FontSize', 15)

subplot(3, 2, 3)
plot(power_sig_out(2,:), LineWidth= 1.5)
title(sprintf("courant de biais de %.2f mA, B2B", bias_currents(2)*1e3))
ylabel('puissance (dB)')
grid on
set(gca, 'FontSize', 13)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 15)
set(get(gca,'ylabel'), 'FontSize', 15)

subplot(3, 2, 4)
plot(power_fiber(2,:), LineWidth= 1.5)
title(sprintf("courant de biais de %.2f mA, propagation de 20km dans une fibre", bias_currents(2)*1e3))
ylabel('puissance (dB)')
grid on
set(gca, 'FontSize', 13)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 15)
set(get(gca,'ylabel'), 'FontSize', 15)

subplot(3, 2, 5)
plot(power_sig_out(3,:), LineWidth= 1.5)
title(sprintf("courant de biais de %.2f mA, B2B", bias_currents(3)*1e3))
ylabel('puissance (dB)')
grid on
set(gca, 'FontSize', 13)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 15)
set(get(gca,'ylabel'), 'FontSize', 15)

subplot(3, 2, 6)
plot(power_fiber(3,:), LineWidth= 1.5)
title(sprintf("courant de biais de %.2f mA, propagation de 20km dans une fibre", bias_currents(3)*1e3))
ylabel('puissance (dB)')
grid on
set(gca, 'FontSize', 13)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 15)
set(get(gca,'ylabel'), 'FontSize', 15)

sgt = sgtitle("Puissance des différents symbole du signal en sortie du laser, avec différents courant de biais, en back to back ou sur une propagation de 20km");
sgt.FontSize = 17;


toc
