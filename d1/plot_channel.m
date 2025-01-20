%Cette partie permet de tracer l'allure des h pour chaque channel

tau = [0,0.5,1,1.5,2];

%definition des channel
A_1 = [1,0.1,0.1,0.1,0.1];
A_2 = [1,0.8,0.6,0.4,0.2];
A_3 = [1,0.8,0.8,0.8,0.8];

Ts = 0.05e-6;
L = 6;
m = -L:1:L;

%calcul des filtres et normalisation
h1 = filtre_canal(m,A_1,tau,Ts,L);
h1_norm = h1/norm(h1);

h2 = filtre_canal(m,A_2,tau,Ts,L);
h2_norm = h2/norm(h2);

h3 = filtre_canal(m,A_3,tau,Ts,L);
h3_norm = h3/norm(h3);

%Tracé des figures
plot(m,h1_norm, 'red', DisplayName='h1 normalisé', Marker='o')
hold on
plot(m,h2_norm,'magenta', DisplayName='h2 normalisé', Marker='x');
plot(m,h3_norm,'blue', DisplayName='h3 normalisé', Marker='*');
hold off
grid on

title('h en fonction de m')
xlabel('valeurs m');
ylabel('Valeur du filtre h');
legend

set(gca, 'FontSize', 21)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)