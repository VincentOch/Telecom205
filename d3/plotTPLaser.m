% Script pour afficher les courbes des mesures effectuées en TP

clear
close all

%% Q3
figure(3)
V = [0 0.69 0.7 0.72 0.72 0.73 0.74 0.75 0.74 0.74 0.75 0.76 0.77 0.78 0.78 0.79 0.79 0.8 0.8 0.8 0.81 0.81 0.81 0.82 0.82 0.82 0.83 0.83 0.83 0.83 0.84 0.84 0.84 0.85 0.85 0.86 0.866 0.86 0.86 0.87 0.88 0.9 0.9 0.91 0.92]; %V
I = [0 0.1 0.2 0.4 0.5 0.7 0.8 0.9 1 1.4 1.6 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 13 14 15 16 17 18 19 20 25 30 35 40 45 50]; %mA
plot(I,V,'LineWidth',2)
xlabel("Courant I (mA)")
ylabel("Tension (V)")
title("tension aux bornes du laser en fonction du courant de polarisation")

grid on
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)


%% Q6
figure(6)
Ipolar = [0 2 4 6 8 10 12 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 25 30 35 40]; %mA
P = [1.562e-3 0.229 0.540 1.007 1.748 2.986 5 8.365 9.416 10.85 12.26 14.2 16.22 18.87 21.6 28.06 71.5 0.145e3 0.207e3 0.28e3 0.97e3 1.625e3 2.253e3 2.866e3]*1e-3; %mW
plot(Ipolar, P,'LineWidth',2)
xlabel("Courant I (mA)")
ylabel("Puissance (mW)")
title("Puissance optique en sortie du laser en fonction du courant de polarisation")

grid on
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

%% Q8
figure(8)
I_plot = linspace(20,150, 100)*1e-3;
R = 2;
Is = 18.5e-3;
Vd = 0.83;  
eta_ext = 0.1283;
eta_wpe = eta_ext*(1-Is./I_plot)./(Vd+R*I_plot)*100; % en %
% ylim([0 10])
plot(I_plot, eta_wpe,'LineWidth',2)
xlabel("Courant I (mA)")
ylabel("rendement à la prise WPE (%)")
title("rendement à la prise WPE en fonction du courant de polarisation")


grid on
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

%% Q11
figure(11)
I_11 = [17 17.2 17.4 17.6 17.8 18 18.2 18.4 18.6 18.8 19 20 21 22 23 24 25 26 27 28 29 30 35 40 45 50 55 60 65 70 75 80];
lambda = (1554+[.3234 .315 .3108 .3024 .294 .2856 .273 .2688 .2646 .2604 .2604 .2604 .2646 .2688 .273 .2814 .294 .2982 .3024 .3066 .3108 .315 .336 .357 .3822 .399 .42 .4452 .4662 .4872 .5082 .5334])*1e-9;
plot(I_11, lambda*1e9, LineWidth=2)
xlabel("Courant I (mA)")
ylabel("longueur d'onde du pic central (nm)")
title('Dérive de la longueur d''onde d''émission en fonction du courant de polarisation')

grid on
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)
%% Q12
figure(12)
I_12 = [18.5 19 19.5 20 21 22 23 24 25 30 35 40 50 60 70 80];
Power_princi = [-37 -29.5 -17.54 -14.3 -11.08 -9.33 -8.1 -7.14 -6.34 -3.76 -2.2 -1.07 0.58 1.69 2.62 3.38]; %dBm
Power_second = [-41.96 -39.408 -38.712 -38.712 -39.176 -39.64 -40.104 -40.336 -40.8 -41.96 -42.696 -43.088 -43.872 -44.6 -45.128 -45.392]; %dBm
plot(I_12, -Power_second+Power_princi,LineWidth=2)
xlabel("Courant I (mA)")
ylabel('taux de compression (dB)')
title('Taux de suppression des modes secondaires en fonction du courant de polarisation')



grid on
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

%% Q13
figure(13)
c = 3e8;
n = 1;
f = c./(n*lambda);
fref = c/(n*1554.1624*1e-9);
delta = abs(fref-f)*1e-9;
plot(I_11(17:end), delta(17:end), LineWidth= 2);
xlabel('Courant I (mA)')
ylabel("delta fréquence (GHz)")
title("Dérive de la fréquence thermique en fonction du courant de polarisation et au dessus du seuil")

grid on
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

