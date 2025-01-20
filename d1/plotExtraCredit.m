load("workspaceChannel1ExtraCredit.mat")


semilogy(10*log10(Eb./N0),DFE0Err,'blue', DisplayName='BER DFE sans erreur', Marker='^');

hold on
semilogy(10*log10(Eb./N0),DFE1Err,'blue', DisplayName='BER DFE correcteur 1 erreur', Marker='o');
semilogy(gca,10*log10(Eb./N0),DFE2Err,'blue', DisplayName='BER DFE correcteur 2 erreurs', Marker='x');

semilogy(10*log10(Eb./N0),Thres0Err,'magenta', DisplayName='BER Threshold sans erreur', Marker='^');
semilogy(10*log10(Eb./N0),Thres1Err,'magenta', DisplayName='BER Threshold correcteur 1 erreur', Marker='o');
semilogy(gca,10*log10(Eb./N0),Thres2Err,'magenta', DisplayName='BER Threshold correcteur 2 erreurs', Marker='x');

semilogy(10*log10(Eb./N0),ZF0Err,'red', DisplayName='BER ZF sans erreur', Marker='^');
semilogy(10*log10(Eb./N0),ZF1Err,'red', DisplayName='BER ZF correcteur 1 erreur', Marker='o');
semilogy(gca,10*log10(Eb./N0),ZF2Err,'red', DisplayName='BER ZF correcteur 2 erreurs', Marker='x');

hold off
grid on
xlim([dBMin dBMax])
ylim([10^(-6) 1])

title('BER en fonction de Eb/N0 en échelle log')%Préciser le canal et la modulation
xlabel('Eb/N0 en dB');
ylabel('Valeur du BER');
legend

set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)
