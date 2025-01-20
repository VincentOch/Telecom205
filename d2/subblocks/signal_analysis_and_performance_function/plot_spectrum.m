function  [OutputSpectrum] = plot_spectrum(data,win,fs,lineSpec_index,fullband,scaletype)
% plot_spectrum - plot a spectrum of a vector with a blackman windowing
%
% Syntax:   [spectre] = plot_spectrum(data, win ,fs,lineSpec_index,fullband,scaletype)
%
% Inputs:
%    data             - the input signal
%    win              - window number for display
%    fs               - Sampling frequency
%    lineSpec_index   - lineSpec index (shortcut for color+line style)
%    fullband         - if equal to 1, the spectrum is plot from -fs/2
%    to fs/2 otherwise it is plot between 0 and fs/2. The first
%    configuration is interesting if the signal is complex
%    scaletype        - 'lin' : plots the spectrum with linear horizontal scale
%                       'log' : plots the spectrum with logarithmic horizontal scale
%    
%
% Outputs:
%    OutputSpectrum - Output Spectrum
%
% Other m-files required: black.m
% Subfunctions: none
% MAT-files required: none
%
% 
% Author: Germain PHAM, Chadi JABBOUR
% C2S, COMELEC, Telecom Paris, Palaiseau, France
% email address: dpham@telecom-paris.fr
% Website: https://c2s.telecom-paristech.fr/TODO
% Feb. 2021, March 2022
%------------- BEGIN CODE --------------


if nargin<6
    if win<50
        scaletype = 'lin';
    else
        scaletype = 'log';
    end
end

if nargin<5
    fullband = 0;
end

if nargin<4
    color='b';
    linewidth=2;
   
else
     switch lineSpec_index
         case 0
           color=0;
      
         case 1
             color='b';
             linewidth=2;
         case 2
             color='r';
              linewidth=2;
         case 3
             color='g';
             linewidth=2;
         case 4
             color='k';
             linewidth=2;
                 
         case 5
             color='y';
              linewidth=2;
         
         case 6
             color='m';
              linewidth=2;
         
         case 7
             color='c';
                linewidth=2;
         
         case 8
             color='--b';
             linewidth=2;
         
         case 9
             color='--r';
              linewidth=2;
              
         case 10
             color='--g';
              linewidth=2;
         
         case 11
             color='b';
             linewidth=1;
         case 12
             color='r';
              linewidth=1;
         case 13
             color='g';
             linewidth=1;
         case 14
             color='black';
             linewidth=1;
                 
         case 15
             color='y';
              linewidth=1;
         
         case 16
             color='m';
              linewidth=1;
         
         case 17
             color='c';
                linewidth=1;
         
         case 18
             color='--b';
             linewidth=1;
         
         case 19
             color='--r';
              linewidth=1;
              
         case 20
             color='--g';
              linewidth=1;
          
         otherwise 
            color='b';
     end
end
             
             
    
[nb_lines,nb_col]=size(data);
if (nb_lines>nb_col) % This part turn a column vector to a line vector
    data=data.'; 
end 

N = length(data);

% Blackman window
w = blackman(length(data));

% Window Normalization 
U = (1/length(w))*sum(abs(w).^2);
w = sqrt(1/U)*w/length(w);

% Spectrum computation
OutputSpectrum = abs(fft(data.*w')).^2;

if fullband
    f              = -fs/2:fs/N:fs/2-fs/N;
    OutputSpectrum = fftshift(OutputSpectrum);
else
    f              = 0:fs/N:fs/2-fs/N;
    OutputSpectrum = OutputSpectrum(1:N/2)*2; 
    % On tronque le vecteur à Fs/2 et 
    % On double la PSD pour lire directement la puissance d'une sinusoide
    % d'un seul côté du spectre.  
end

if color
    figure(win)
    switch scaletype
    case 'lin'
        plot(f, 10*log10(OutputSpectrum),color,'linewidth',linewidth)
    case 'log'
        semilogx(f,10*log10(OutputSpectrum),color,'linewidth',linewidth)
    end
    
end

xlabel('f (Hz)');
ylabel('PSD (dB/bin)');
