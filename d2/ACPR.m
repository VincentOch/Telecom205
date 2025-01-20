function [acpr_left, acpr_right, power_sig] = ACPR(largeurBandeUtile, largeurBandeRf, Fc, sig, Fs)
%largeurBandeRf : largeur de la bande si par exemple, 30 MHz ici, avec la
%garde
%largeurBandeUtile : largeur bande utile sans la garde (20 ici)
%Fc : fréquence centrale
%sig : le signal
%garde : bande de garde
%Fs : la fréquence d'échantillonage
garde = largeurBandeRf - largeurBandeUtile;

if garde < 0
    disp("Problème dans la bande de garde (négative). Vérifier arguments de la fonction ACPR")
    return
end

psd = plot_spectrum(sig, 5, Fs, 0); % le dernier zéro empêche le plot
len_sig = length(psd)*2;
freq_min = Fc - largeurBandeUtile/2;
freq_max = Fc + largeurBandeUtile/2;
freq_min_distor_right = garde + freq_max;
freq_max_distor_right = freq_min_distor_right + largeurBandeUtile;
freq_max_distor_left = freq_min - garde;
freq_min_distor_left = freq_max_distor_left - largeurBandeUtile;

bin_freq_min = fix(freq_min*len_sig/Fs);
bin_freq_max = fix(freq_max*len_sig/Fs);
bin_freq_min_distor_right = fix(freq_min_distor_right*len_sig/Fs);
bin_freq_max_distor_right = fix(freq_max_distor_right*len_sig/Fs);
bin_freq_min_distor_left = fix(freq_min_distor_left*len_sig/Fs);
bin_freq_max_distor_left = fix(freq_max_distor_left*len_sig/Fs);

power_sig = sum(psd(bin_freq_min:bin_freq_max));
power_right = sum(psd(bin_freq_min_distor_right:bin_freq_max_distor_right));
power_left = sum(psd(bin_freq_min_distor_left:bin_freq_max_distor_left));
acpr_right = power_right/power_sig;
acpr_left = power_left/power_sig;

power_sig = 10*log10(power_sig);
acpr_left = 10*log10(acpr_left) ;
acpr_right = 10*log10(acpr_right);

end

