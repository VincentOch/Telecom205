function [BBOutput,powerConsumption] = BBamp(input,Gain,NF,IIP3,R,BW_noise,continuousTimeSamplingRate)
    %BBamp - Simulates a BB amplifier with a noise figure and an IIP3
    %
    % Syntax:  [lnaOutput] = BBamp(input,Gain,NF,IIP3,R,BW_noise)
    %
    % Inputs:
    %    input    - the input signal
    %    Gain     - the gain (in dB)
    %    NF       - the noise figure (in dB)
    %    IIP3     - the input third-order intercept point (in dB)
    %    R        - the input impedance (in ohm) (should be matched)
    %    BW_noise - the noise bandwidth usually fixed at the Nyquist frequency of the simulation (in Hz)
    %    continuousTimeSamplingRate - Frequency to emulate the continuous
    %    time behavior
    
    %
    % Outputs:
    %    BBOutput - the (noisy and distorted) output amplified signal of
    %               the BBamplifier
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: rfPA
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020, Apr. 2020, Dec. 2020, March 2022, March 2023
    %------------- BEGIN CODE --------------

    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        k  = physical_constant("boltzmann constant");
    else
        k  = physconst('Boltzmann');
    end
    T  = 290; % room temperature

    % Noise
    AntennaNoise = randn(size(input))*sqrt(k*T*continuousTimeSamplingRate/2*R);
    NoiseFilter=fir1(10000,BW_noise/(continuousTimeSamplingRate/2));
    AntennaNoise_filtered=filter(NoiseFilter,1,AntennaNoise);
    
    
    AddedNoise   = sqrt(10.^(NF/10)-1)*AntennaNoise_filtered;

    % Nonlinearity effect
    alpha3 = -4*10.^(Gain/20)/(10.^((IIP3-30+10*log10(2*R))/10));

    % Output
    inputReferedNoisySig = input+AddedNoise;
    BBOutput = 10.^(Gain/20)*inputReferedNoisySig + alpha3*(inputReferedNoisySig.^3);
    
    powerConsumption = 1e-10 * BW_noise;

    %------------- END OF CODE --------------
    
