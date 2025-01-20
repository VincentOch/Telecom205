function [SNDR,PS]= perf_estim(input,Bin_In,n,Bin_Limits,fullband )
% perf_estim - Calculates the SNDR and power for an input signal wave
%
% Syntax:  [SNDR,PS]= perf_estim(input,Bin_In,n,Bin_Limits,fullband)
%
% Inputs:
%   input     : (time domain) (vector) signal 
%   Bin_In    : bin of the input signal, 
%               if  k=0  the function detects the maximum and considers it as the input 
%   n         : number of points around the central bin on which the power is integrated
%   Bin_Limits: vector of 2 elements which specifies the frequency bin limits of 
%               the (noise) integration band 
%   fullband  : binary flag : 
%                - when equal to 1, the spectrum is plot from -fs/2 to fs/2
%                - otherwise, it is plot between 0 and fs/2
%               The first case is interesting if the signal is complex
%
% Outputs:
%    SNDR : the ratio of the power integrated around Bin_In and the rest of the
%           spectrum defined by Bin_Limits. It is the signal to noise (and distorsion)
%           ratio when Bin_In corresponds to the signal bin and Bin_Limits corresponds
%           to the useful band. 
%    PS   : the power at the bin Bin_In (+/- the window bins). It is the signal power
%           when Bin_In corresponds to the signal bin
%
% Other m-files required: black.m
% Subfunctions: none
% MAT-files required: none
%
% See also: ADC
% Author: Germain PHAM, Chadi JABBOUR
% C2S, COMELEC, Telecom Paris, Palaiseau, France
% email address: dpham@telecom-paris.fr
% Feb. 2021


    
if nargin<5
  fullband = 0;
end

% reshape the vector as lign vector if it is a column vector
size_input= size(input);
if size_input(1)>1 
       input=input.';
end 

% Blackman window
w = blackman(length(input));

% Window Normalization 
U = (1/length(w))*sum(abs(w).^2);
if fullband
  w = sqrt(1/U)*w/length(w);
else
  w = sqrt(2/U)*w/length(w);
end

% calculating the spectrum
OutputSpectrum = abs(fft(input.*w')).^2;    

if fullband
  OutputSpectrum = fftshift(OutputSpectrum);
end
plot(OutputSpectrum)

[~, Bin_In] = max(OutputSpectrum);
PS    = sum(OutputSpectrum(Bin_In-n:Bin_In+n));  %signal power calculation
PND   = sum(OutputSpectrum(Bin_Limits(1):Bin_Limits(2)))-PS+(2*n+1)*(OutputSpectrum(Bin_In-n-1)+OutputSpectrum(Bin_In+n+1))/2;          %noise and distortions calculations
  % The term (2*n+1) ..../2 was added to compensante for the noise inside the signal bin. We assume that the noise floor is equal to the noise before and after the signal 
SNDR  = 10*log10(PS/PND);


end