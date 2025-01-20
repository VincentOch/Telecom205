function simple_plot_spectrum(x,fs)
% simple_plot_spectrum - A short function to plot the spectrum of a signal
%
% Syntax:
%	simple_plot_spectrum(x,fs)
%
% Inputs:
%	x  [1D array]   : the input signal as 1D vector
%	fs [scalar]     : the sampling frequency in Hz, optional, default is 1 MHz
%
% Outputs:
%	none
%
% Example:
%	simple_plot_spectrum(sin(2*pi*1.5e6*(0:1/20e6:1)),20e6)
%
% See also:
%	
%
% Germain PHAM, C2S Telecom Paris, 2023

% Check the number of input arguments
% if the second argument is not provided, set the default value
if nargin < 2
    fs = 1e6;
end

win           = blackman(length(x),'periodic');
x_windowed    = x(:).*win(:); % windowing
x_psd         = abs(fft(x_windowed)).^2;

% Define a vector that count the bins
Nx                    = length(x_psd);
bin_freq_val          = [1:Nx];

% Define a vector that contains the frequency values
bin_freq_val_shift    = -(Nx-1)/2 : (Nx-1)/2;
freq_val_shift        = bin_freq_val_shift/Nx*fs;

x_psd_db              = 10*log10(x_psd);

figure
% plot PSD with reordering and frequency axis
plot(freq_val_shift/1e6,fftshift(x_psd_db))
xlabel('Frequency (MHz)');
ylabel('PSD (dB)');
title('PSD with reordering and frequency axis');


end