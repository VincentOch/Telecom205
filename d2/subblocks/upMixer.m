function [rfMixerOut,powerConsumption] = upMixer(InputI,InputQ,Flo,continuousTimeFs)
    %upMixer - Implements an AMS/RF upmixer
    %   For sake of simplicity, in the project, mixers are ideal
    %
    % Syntax:  rfMixerOut = upMixer(InputI,InputQ,Flo,continuousTimeFs)
    %
    % Inputs:
    %    InputI           - Input analog (baseband) signal - I path
    %    InputQ           - Input analog (baseband) signal - Q path
    %    Flo              - Frequency of the local oscillator (in Hz)
    %    continuousTimeFs - Sampling frequency emulating the continuous time (in Hz)
    %
    % Outputs:
    %    rfMixerOut - Output RF signal
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: downMixer
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020
    %------------- BEGIN CODE --------------

    sigLen       = length([InputI(:) InputQ(:)]);
    timeInstants = (0:(sigLen-1))/continuousTimeFs;

    loI          = cos(2*pi*Flo*timeInstants(:));
    loQ          = cos(2*pi*Flo*timeInstants(:)+pi/2);  % == -sin(2*pi*Flo*timeInstants(:))
    
    % Quadrature Mixing
    rfMixerOut = (InputI.*loI + InputQ.*loQ);

    % In practice, mixers have a limited bandwidth ; 
    % For sake of simplicity, this effect is not implemented here
    % [b_butter,a_butter]=butter(10,0.9);
    % rfMixerOut=Butt_Filter(rfMixerOut_nofil,b_butter,a_butter,0,50,Fs_Cont);

    powerConsumption = 5e-12 * Flo;
    %------------- END OF CODE --------------

