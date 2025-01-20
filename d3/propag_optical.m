function sig_out = propag_optical(alpha, D, S, lambda, L, Fe, sig_in)
    % Propagates the optical signal through a fiber with attenuation and chromatic dispersion
    % Inputs:
    %   alpha: Attenuation coefficient in dB/km
    %   D : Dispersion parameter in ps/nm/km
    %   S : Dispersion slope in ps/nm^2/km
    %   lambda : longueur d'onde en nm
    %   L: Length of the fiber in km
    %   Fe: Sampling frequency in Hz, equal to 2 times the maximum frequence.
    %   Pour une entrée en sinus cardinal de largeur T, c'est 2/T, pour un
    %   cosinus surélevé de facteur rho (rho entre 0 et 1), c'est (1+rho)*Rs
    %   sig_in: Input signal in the time domain
    % Output:
    %   sig_out: Output signal after propagation through the fiber
    
    D = D*1e-3;  % Dispersion parameter in s/m/km
    S = S*1e6;  % Dispersion slope in s/m^2/km
    lambda = lambda*1e-9; %m
    c=3.e8;
    beta2 = -D*lambda^2/(2*pi*c);
    beta3 = S*lambda^4/(4*pi^2*c^2);
    
    % Attenuation
    sig = sig_in * 10^(-alpha/20*L) ; % divisé par 10 pour passer en linéaire ; par 2 pour passer en amplitude
    
    % Chromatic dispersion
    sig_psd = fft(sig);
    sig_psd = fftshift(sig_psd);
    N = length(sig_psd);
    freq = linspace(-Fe/2, Fe/2, N);
    delta_w = 2 * pi * freq; % As we consider w0 = 0
    H = exp(-1i * (beta2 * delta_w.^2 / 2 + beta3 * delta_w.^3 / 6) * L);
    sig_psd_out = sig_psd .* H;
    
    % Inverse Fourier transform
    sig_psd_out = ifftshift(sig_psd_out);
    sig_out = ifft(sig_psd_out);
    
    end