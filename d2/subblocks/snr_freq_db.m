function [snr_result, signal_power] = snr_freq_db(signal, Fs, band, f0)
    % Détermine le SNR d'un signal temporel.
    % signal : Signal temporel complet
    % Fs : Fréquence d'échantillonnage du signal
    % band : bande de fréquences d'intérêt
    % f0 : fréquence d'intérêt du signal

    frequencies = linspace(0, Fs, length(signal));
   
    signal_PSD = plot_spectrum(signal,1,length(signal), 0);
    i_fmin = round(band(1)/Fs*length(signal)); % pour trouver le bon bin
    i_fmax = round(band(2)/Fs*length(signal));
    
    i_fmin = i_fmin + 1;
    i_fmax = i_fmax + 1;

    if i_fmax >= length(frequencies)
        i_fmax = i_fmax - 1;
    end

    [~,pos] = max(signal_PSD);%round(f0/Fs * length(signal));

    bins = pos-4:pos+4;
    bins = bins + 1 ; % Car matlab indexe à partir de 1
    signal_power = sum(signal_PSD(bins));
    
    noise = signal_PSD; 
    noise(bins) = 0;
    noise = noise(i_fmin:i_fmax);
    noise_power = sum(noise);
    snr_result = 10*log10(signal_power/noise_power);
end
