clear
close all

%% cas de simulation
[K_40k_zf, rej_40k_zf] =  d4_perfs_students('40k','zf');
[K_40k_dfe, rej_40k_dfe] =  d4_perfs_students('40k','dfe');

[K_400k_zf, rej_400k_zf] =  d4_perfs_students('400k','zf');
[K_400k_dfe, rej_400k_dfe] =  d4_perfs_students('400k','dfe');

[K_4M_zf, rej_4M_zf] =  d4_perfs_students('4M','zf');
[K_4M_dfe, rej_4M_dfe] =  d4_perfs_students('4M','dfe');

[K_40M_zf, rej_40M_zf] =  d4_perfs_students('40M','zf');
[K_40M_dfe, rej_40M_dfe] =  d4_perfs_students('40M','dfe');



%% plots
% For 40kbits/s
fig = figure();
subplot(2,2,1)
plot(K_40k_zf, rej_40k_zf,'o--',"LineWidth",1.5,DisplayName='ZF')
hold on
plot(K_40k_dfe, rej_40k_dfe,'o--',"LineWidth",1.5,DisplayName='DFE')
hold off
grid on
title('rejection rate function of the number of user, Rl = 40 kbit/s')
xlabel('Number of users');
ylabel('Rejection rate');
legend('Location', 'best');

set(gca, 'FontSize', 18)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

%% For 400kbits/s
subplot(2,2,2)
plot(K_400k_zf, rej_400k_zf,'o--',"LineWidth",1.5,DisplayName='ZF')
hold on
plot(K_400k_dfe, rej_400k_dfe,'o--',"LineWidth",1.5,DisplayName='DFE')
hold off
grid on
title('rejection rate function of the number of user, Rl = 400 kbit/s')
xlabel('Number of users');
ylabel('Rejection rate');
legend('Location', 'best');

set(gca, 'FontSize', 18)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

%% For 4Mbits/s
subplot(2,2,3)
plot(K_4M_zf, rej_4M_zf,'o--',"LineWidth",1.5,DisplayName='ZF')
hold on
plot(K_4M_dfe, rej_4M_dfe,'o--',"LineWidth",1.5,DisplayName='DFE')
hold off
grid on
title('rejection rate function of the number of user, Rl = 4 Mbit/s')
xlabel('Number of users');
ylabel('Rejection rate');
legend('Location', 'best');

set(gca, 'FontSize', 18)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

%% For 40Mbits/s
subplot(2,2,4)
plot(K_40M_zf, rej_40M_zf,'o--',"LineWidth",1.5,DisplayName='ZF')
hold on
plot(K_40M_dfe, rej_40M_dfe, 'o--',"LineWidth",1.5,DisplayName='DFE')
hold off
grid on
title('rejection rate function of the number of user, Rl = 40 Mbit/s')
xlabel('Number of users');
ylabel('Rejection rate');
legend('Location', 'best');

set(gca, 'FontSize', 18)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 15)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)