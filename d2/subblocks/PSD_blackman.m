function y = PSD_blackman(x)
    % Calcule la densité de puissance spectrale d'un signal x fenêtré avec
    % Blackman.
    win = blackman(length(x), 'periodic'); % fenêtrage blackman pour la FFT
    X = fft(x(:).*win(:));
    y = abs(X).^2;
end