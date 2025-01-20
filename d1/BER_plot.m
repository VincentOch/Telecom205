
figure();
%Variable à modifier pour tracer, pour la canal donné, l'allure du BER pour
%les 3 equalizers possibles

%ATTENTION l'approximation donne des fois des valeurs de -infini pour le
%BER en échelle semilog, ce qui n'est pas tracé par matlab. Pour affiner
%les résultats, il est possible d'augmenter iteration mais cela augmente
%grandement le temps de calcul

%Modulation souhaitée (parmis PSK, QAM et PAM)
mod = 'PSK';
M=2;

R=1; %26/31 si BCH 1 erreur et 21/31 si 2 erreurs. Rendement si codage présent (pas utilisé ici)

Nb_points = 50; %Nombre de points du tracé
canalSelect = 1; %Le canal choisit, entre 1,2 et 3
%le range en dB de l'échelle que l'on veut pour Eb/N0
%augmenter le nombre de point si l'échelle est plus large
dBMin = 0;
dBMax = 30;




%Définition des canaux
tau = [0,0.5,1,1.5,2];
A_1 = [1,0.1,0.1,0.1,0.1];
A_2 = [1,0.8,0.6,0.4,0.2];
A_3 = [1,0.8,0.8,0.8,0.8];
Ts = 0.05e-6;
L = 6;
m = -L:1:L;
N=100;

%% Choisit le canal

switch canalSelect
    case 1
        h = filtre_canal(m,A_1,tau,Ts,L);
    case 2
        h = filtre_canal(m,A_2,tau,Ts,L);
    case 3
        h = filtre_canal(m,A_3,tau,Ts,L);
end

h_norm = h/norm(h);


%% Création duvecteur de SNR
SNR = linspace(dBMin,dBMax,Nb_points);
%Calcul du N0 pour le bon SNR
Eb = 1;
N0 = Eb./(10.^(SNR/10));


%% Pour l'ARQ
N0new = 10.^(-linspace(dBMin,dBMax,Nb_points)/10);
EsBPSK = 1;
Es8QAM = 2*(8-1)/3;
Es16QAM = 2*(16-1)/3;
moyB = zeros(1,Nb_points);
moyBB = zeros(1,Nb_points);
moy8 = zeros(1,Nb_points);
moy8B = zeros(1,Nb_points);
moy16 = zeros(1,Nb_points);
moy16B = zeros(1,Nb_points);

moyenneTerme = 100;
for kk=1:moyenneTerme
    moyB = moyB + ARQ(N0new, 'PSK', 2, 1, 0)/moyenneTerme;
    moyBB = moyBB + ARQ(N0new, 'PSK', 2, 1, 2)/moyenneTerme;
    moy8 = moy8 + ARQ(N0new, 'QAM', 8, 3, 0)/moyenneTerme;
    moy8B = moy8B + ARQ(N0new, 'QAM', 8, 3, 2)/moyenneTerme;
    moy16 = moy16 + ARQ(N0new, 'QAM', 16, 4, 0)/moyenneTerme;
    moy16B = moy16B + ARQ(N0new, 'QAM', 16, 4, 2)/moyenneTerme;

end

%%

plot(10*log10(EsBPSK./N0new),moyB,'b--o', DisplayName='BPSK',LineWidth= 2);
hold on
plot(10*log10(EsBPSK./N0new),moyBB,'b--x', DisplayName='BPSK with BCH',LineWidth= 2);
plot(10*log10(Es8QAM./N0new),moy8,'r--o', DisplayName='8 QAM',LineWidth= 2);
plot(10*log10(Es8QAM./N0new),moy8B,'r--x', DisplayName='8 QAM with BCH',LineWidth= 2);
plot(10*log10(Es16QAM./N0new),moy16,'m--o', DisplayName='16 QAM',LineWidth= 2);
plot(10*log10(Es16QAM./N0new),moy16B,'m--x', DisplayName='16 QAM with BCH',LineWidth= 2);


hold off 

%% Pour tracer des courbes avec les canaux, hors de l'EC
%iteration = 10e5; %attention très long si élevé

%Calcul du BER théorique

% gamma = GainSNR(mod,M);
% nmin = 1;
% x = sqrt(gamma*Eb./N0/R);

%BERth = nmin*1/2*erfc(x/sqrt(2)); %we have NminQ(sqrt(gamma*Eb/N0))

% %BER dans les cas avec erreur, ne marche que pour la BPSK
% BERexp1err = BERThres(N0,mod,M,iteration,1);
% BERexp2err = BERThres(N0,mod,M,iteration,2);
% BERexp0err = BERThres(N0,mod,M,iteration,0);

%Calcul des BER
% BERexpThres = BERThresCanal(N0,mod,M,h_norm,N,iteration);
% BERexpZF = BERZF(N0,mod,M,h_norm,N,iteration);
% BERexpDFE = BERDFE(N0,mod,M,h_norm,N,iteration);
% 
% %semilogy(10*log10(Eb./N0),BERth, 'green', DisplayName='BER théorique', Marker='+')
% 
% hold on
% 
% %BER threshold, ZF et DFE
% semilogy(10*log10(Eb./N0),BERexpZF,'blue', DisplayName='BER expérimental avec ZF', Marker='.');
% semilogy(10*log10(Eb./N0),BERexpThres,'magenta', DisplayName='BER expérimental avec Threshold', Marker='x');
% semilogy(gca,10*log10(Eb./N0),BERexpDFE,'red', DisplayName='BER expérimental avec DFE', Marker='o');
% 
% %BER avec et sans codage, ne marche que pour la BPSK
% % 
% % semilogy(10*log10(Eb./N0),BERexp1err,'magenta', DisplayName='BER expérimental pour BCH 1 erreur', Marker='+');
% % semilogy(10*log10(Eb./N0),BERexp2err,'red', DisplayName='BER expérimental pour BCH 2 erreur', Marker='*');
% % semilogy(10*log10(Eb./N0),BERexp0err,'blue', DisplayName='BER expérimental sans codage', Marker='x');
% 
% hold off

%% Pour améliorer les courbes

grid on
xlim([dBMin dBMax])


title('throughput en fonction de Es/N0 en échelle log')%Préciser le canal et la modulation
xlabel('Es/N0 en dB');
ylabel('throughput');
legend

set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)