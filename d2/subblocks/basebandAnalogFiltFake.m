function [analogFilteredOutput] = basebandAnalogFiltFake(input,hfir,NF,Rin,continuousTimeFs)
    %basebandAnalogFiltFake - Implements a baseband fake analog filter
    %   due to numerical issues, the filter has to be provided as second-order sections (SOS)
    %   Please refer to help zp2sos for example for more information
    %
    % Syntax:  [analogFilteredOutput] = basebandAnalogFiltFake(input,hfir,NF,Rin,continuousTimeFs)
    %
    % Inputs:
    %    input            - Input baseband signal
    %    hfir             - The filter description as second-order sections
    %    NF               - Noise factor of the filter (in dB)
    %    Rin              - Input impedance of the filter (in ohm)
    %    continuousTimeFs - Sampling frequency emulating the continuous time (in Hz)
    %
    % Outputs:
    %    analogFilteredOutput - analog filtered output signal
    %
    % Other m-files required: none
    % Subfunctions: 
    % MAT-files required: none
    %
    % See also: 
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020, March 2022, Dec. 2023
    %------------- BEGIN CODE --------------
    

    k  = physconst('Boltzmann');
    T  = 290; % room temperature
          
    
    AntennaNoise = randn(size(input,1),1)*sqrt(k*T*continuousTimeFs/2*Rin); 
    AddedNoise   = sqrt(10.^(NF/10)-1) * AntennaNoise;

    analogFilteredOutput = filtfilt(hfir,1,input+AddedNoise); 

    %------------- END OF CODE --------------
