function [OutputI,OutputQ,powerConsumption] = downMixer(InputRF,Flo,continuousTimeFs)
    %downMixer - Implements an AMS/RF downmixer
    %   For sake of simplicity, in the project, mixers are ideal
    %
    % Syntax:  [OutputI,OutputQ] = downMixer(InputRF,Flo,continuousTimeFs)
    %
    % Inputs:
    %    InputRF          - RF signal
    %    Flo              - Frequency of the local oscillator (in Hz)
    %    continuousTimeFs - Sampling frequency emulating the continuous time (in Hz)
    %
    % Outputs:
    %    OutputI - Ouput analog (baseband) signal - I path
    %    OutputQ - Ouput analog (baseband) signal - Q path
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: upMixer
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020
    %------------- BEGIN CODE --------------

    sigLen       = length(InputRF);
    timeInstants = (0:(sigLen-1))/continuousTimeFs;

    loI          = cos(2*pi*Flo*timeInstants(:));
    loQ          = cos(2*pi*Flo*timeInstants(:)+pi/2); % == -sin(2*pi*Flo*timeInstants(:))
    
    % Quadrature Down Mixing
    OutputI = InputRF(:) .* loI;
    OutputQ = InputRF(:) .* loQ;

    % In practice, mixers have a limited bandwidth ; 
    % For sake of simplicity, this effect is not implemented here    
    % [b_butter,a_butter]=butter(10,0.9);
    % MixerOut=Butt_Filter(MixerOut_nofil,b_butter,a_butter,0,50,Fs_Cont);

    powerConsumption = 5e-12 * Flo;
    %------------- END OF CODE --------------

