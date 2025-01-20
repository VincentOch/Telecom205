function [analogFilteredOutput] = basebandAnalogFilt(input,FiltSOS,NF,Rin,continuousTimeFs)
    %basebandAnalogFilt - Implements a baseband analog filter
    %   due to numerical issues, the filter has to be provided as second-order sections (SOS)
    %   Please refer to help zp2sos for example for more information
    %
    % Syntax:  [analogFilteredOutput] = basebandAnalogFilt(input,FiltSOS,NF,Rin,continuousTimeFs)
    %
    % Inputs:
    %    input            - Input baseband signal
    %    FiltSOS          - The filter description as second-order sections
    %    NF               - Noise factor of the filter (in dB)
    %    Rin              - Input impedance of the filter (in ohm)
    %    continuousTimeFs - Sampling frequency emulating the continuous time (in Hz)
    %
    % Outputs:
    %    analogFilteredOutput - analog filtered output signal
    %
    % Other m-files required: none
    % Subfunctions: sosfilt
    % MAT-files required: none
    %
    % See also: zp2sos
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020, March 2022
    %------------- BEGIN CODE --------------
    
    if size(FiltSOS,2)~=6
        error('TELECOM205:basebandAnalogFilt:inputType','Argument "FiltSOS" must be an SOS filter');
    end

    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        k  = physical_constant("boltzmann constant");
    else
        k  = physconst('Boltzmann');
    end
    
    T  = 290; % room temperature
          
    
    AntennaNoise = randn(size(input,1),1)*sqrt(k*T*continuousTimeFs/2*Rin); 
    AddedNoise   = sqrt(10.^(NF/10)-1) * AntennaNoise;

    analogFilteredOutput = sosfilt(FiltSOS,input+AddedNoise); 

    %------------- END OF CODE --------------
