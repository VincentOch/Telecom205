function [lnaOutput,powerConsumption] = rfLNA(input,Gain,NF,IIP3,R,BW_noise)
    %rfLNA - Simulates a low noise amplifier with a noise figure and an IIP3
    %
    % Syntax:  [lnaOutput] = rfLNA(input,Gain,NF,IIP3,R,BW_noise)
    %
    % Inputs:
    %    input    - the input signal
    %    Gain     - the (power) gain of the LNA (in dB)
    %    NF       - the noise figure of the LNA (in dB)
    %    IIP3     - the input third-order intercept point (in dBm)
    %    R        - the input impedance of the LNA (in ohm) (should be matched)
    %    BW_noise - the noise bandwidth usually fixed at the Nyquist frequency of the simulation (in Hz)
    %
    % Outputs:
    %    lnaOutput - the (noisy and distorted) output amplified signal of the LNA
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
    % Feb. 2020, Apr. 2020, Dec. 2020
    %------------- BEGIN CODE --------------

    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        k  = physical_constant("boltzmann constant");
    else
        k  = physconst('Boltzmann');
    end
    
    T  = 290; % room temperature

    % Noise
    AntennaNoise = randn(size(input))*sqrt(k*T*BW_noise*R);
    AddedNoise   = sqrt(10.^(NF/10)-1)*AntennaNoise;

    % Nonlinearity effect
    alpha3 = -4*10.^(Gain/20)/(10.^((IIP3-30+10*log10(2*R))/10));

    % Output
    inputReferedNoisySig = input+AddedNoise;
    lnaOutput = 10.^(Gain/20)*inputReferedNoisySig + alpha3*(inputReferedNoisySig.^3);
    
    powerConsumption = 1e-2 / (10^(NF/10) - 1);

    %------------- END OF CODE --------------
    
