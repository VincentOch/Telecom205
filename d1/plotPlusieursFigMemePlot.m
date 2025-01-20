%Supperposer différents fichiers fig

% Ouvrir les deux fichiers .fig
fig1 = openfig('fig\canal1 16QAM_3.fig');
fig2 = openfig('fig\canal1 PSK_3.fig');
fig3 = openfig('fig\canal1 8QAM_3.fig');
% Extraire les axes des figures
ax1 = gca(fig1);
ax2 = gca(fig2);
ax3 = gca(fig3);

% Extraire les objets courbes des deux axes
curves1 = findobj(ax1, 'Type', 'line');
curves2 = findobj(ax2, 'Type', 'line');
curves3 = findobj(ax3, 'Type', 'line');

% Créer une nouvelle figure
figure;
% Copier les courbes de la première figure
for i = 1:length(curves1)
    xData = get(curves1(i), 'XData');
    % xData = xData(xData>=0);
    yData = get(curves1(i), 'Ydata');
    % yData = yData(length(yData)-length(xData)+1:end);
    semilogy(xData, yData, 'DisplayName', sprintf('%s 16QAM',get(curves1(i), 'DisplayName')), ...
             'Color', get(curves1(i), 'Color'), ...
             'LineStyle', get(curves1(i), 'LineStyle'), ...
             'Marker', 'o', ...
             'LineWidth',1.25);  

    if i==1
        hold on
    end
end

% Copier les courbes de la deuxième figure
for i = 1:length(curves2)
    xData = get(curves2(i), 'XData');
    yData = get(curves2(i), 'YData');
    semilogy(xData, yData, 'DisplayName', sprintf('%s BPSK',get(curves1(i), 'DisplayName')), ...
             'Color', get(curves2(i), 'Color'), ...
             'LineStyle', get(curves2(i), 'LineStyle'), ...
             'LineWidth',1.25);      

end

for i = 1:length(curves3)
    xData = get(curves3(i), 'XData');
    yData = get(curves3(i), 'YData');
    semilogy(xData, yData, 'DisplayName', sprintf('%s 8QAM',get(curves3(i), 'DisplayName')), ...
             'Color', get(curves3(i), 'Color'), ...
             'LineStyle', get(curves3(i), 'LineStyle'), ...
             'LineWidth',1.25, 'Marker','^');      

end
hold off;
% Modifier les labels des courbes


% Afficher la légende
legend('Location','best')
xlim([0 20]);
grid on

% Ajuster les propriétés de l'axe si nécessaire
xlabel('Eb/N0 en dB');
ylabel('BER');
title('BER pour le canal 3, comparaison entre une BPSK, une 8QAM et une 16QAM');

% Garder les courbes visibles



set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)

% Fermer les figures initiales
close(fig1);
close(fig2);
close(fig3);
