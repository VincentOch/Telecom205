function [AnalogOutput,powerConsumption] = DAC(input,Nbits,Vref,DACType,basebandFs,continuousTimeFs)
    %DAC - Implements an AMS Digital to Analog converter 
    %
    % Syntax:  [AnalogOutput] = DAC(input,Nbits,Vref,DACType,basebandFs,continuousTimeFs)
    %
    % Inputs:
    %    input            - the input digital signal
    %    Nbits            - number of bits of the DAC
    %    Vref             - reference voltage (in V)
    %    DACType          - can be 'zoh', 'impulse'
    %    basebandFs       - Sampling rate of the DAC (in Hz)
    %    continuousTimeFs - Sampling frequency emulating the continuous time (in Hz)
    %
    % Outputs:
    %    AnalogOutput - Analog Output of the DAC, updated at the rate of basebandFs
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: ADC
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020
    %------------- BEGIN CODE --------------
    
    % TODO : check imaginary part

    ctOSR = continuousTimeFs/basebandFs;
    
    powerConsumption = 1e-13 * basebandFs * 2^Nbits;
    
    if rem(continuousTimeFs,basebandFs)~=0
        error('TELECOM205:DAC:ctOSR','"basebandFs" should divide "continuousTimeFs"')
    end
    sizeIn = size(input);
    if length(input)~=sizeIn(1)
        error('TELECOM205:DAC:inputType','Argument "input" should be a column-wise matrix or vector')
    end

    % Quantization
    DigInput=floor((input+Vref)*2^(Nbits-1)/Vref);  % Digital Input of the DAC
    DigInput(DigInput<0)         = 0;               % Clipping Down
    DigInput(DigInput>2^Nbits-1) = 2^Nbits-1;       % Clipping Up
    AnalogOutput_dt = ( DigInput-2^(Nbits-1))/2^(Nbits-1)*Vref+Vref/2^Nbits; 
    
    % Continuous time conversion
    AnalogOutput = zeros(sizeIn(1)*ctOSR,sizeIn(2));
    switch lower(DACType)
    case 'zoh'
        %AnalogOutput(:,:) = repelem(AnalogOutput_dt, ctOSR,1);
        AnalogOutput(:) = repmat(AnalogOutput_dt', ctOSR,1);
    case 'impulse'
        AnalogOutput(1:ctOSR:end,:) = AnalogOutput_dt;
    end

    


    %------------- END OF CODE --------------
