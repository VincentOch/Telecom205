function [DigOutput,powerConsumption] = ADC(input,Nbits,Vref,adcSamplingRate,delay,continuousTimeFs)
    %ADC - Implements an AMS Analog to Digital converter 
    %
    % Syntax:  [DigOutput] = ADC(input,Nbits,Vref,adcSamplingRate,delay,continuousTimeFs)
    %
    % Inputs:
    %    input              - the input analog signal
    %    Nbits              - number of bits of the ADC
    %    Vref               - reference voltage (in V)
    %    adcSamplingRate    - Sampling rate of the ADC (in Hz)
    %    delay              - the delay should be used to removed unwanted
    %                         cotinuous time transient effects ; is given in number of 
    %                         samples of the input (@continuousTimeFs)
    %    continuousTimeFs   - Sampling frequency emulating the continuous time (in Hz)
    %
    % Outputs:
    %    DigOutput - Digital output 
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: DAC
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020, Apr. 2020
    %------------- BEGIN CODE --------------

    % TODO : check imaginary part
    
    ctOSR = continuousTimeFs/adcSamplingRate;
    if rem(continuousTimeFs,adcSamplingRate)~=0
        error('TELECOM205:ADC:ctOSR','"adcSamplingRate" should divide "continuousTimeFs"')
    end
    sizeIn = size(input);
    if length(input)~=sizeIn(1)
        error('TELECOM205:ADC:inputType','Argument "input" should be a column-wise matrix or vector')
    end

    % Sampling
    AnaInput_DT = input((1+delay):ctOSR:end,:);

    % Quantization
    quantizedInput=floor((AnaInput_DT+Vref)*2^(Nbits-1)/Vref);  % Quantizing the sampled data
    quantizedInput(quantizedInput<0)         = 0;               % Clipping Down
    quantizedInput(quantizedInput>2^Nbits-1) = 2^Nbits-1;       % Clipping Up
    DigOutput = (quantizedInput-2^(Nbits-1))/2^(Nbits-1)*Vref+Vref/2^Nbits; 

    powerConsumption = 1e-13 * adcSamplingRate * 2^Nbits;
 
    %------------- END OF CODE --------------
