function simple_raw_spectrum(x)
% simple_plot_spectrum - A short function to plot the spectrum of a signal in natural matlab fft bin ordering
%
% Syntax:
%	simple_raw_spectrum(x)
%
% Inputs:
%	x  [1D array]   : the input signal as 1D vector
%
% Outputs:
%	none
%
% Example:
%	simple_raw_spectrum(sin(2*pi*1.5e6*(0:1/20e6:1)))
%
% See also:
%	
%
% Germain PHAM, C2S Telecom Paris, 2023

win           = blackman(length(x),'periodic');
x_windowed    = x(:).*win(:); % windowing
x_psd         = abs(fft(x_windowed)).^2;

% Define a vector that count the bins
Nx                    = length(x_psd);
bin_freq_val          = [1:Nx];

x_psd_db              = 10*log10(x_psd);

figure
% plot PSD with reordering and frequency axis
plot(bin_freq_val,x_psd_db)
xlabel('Bin');
ylabel('PSD (dB)');
title('PSD with Matlab natural fft bin ordering');


end