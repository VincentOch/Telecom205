% completeTxRx_proj - Script that performs the simulation of the complete communication chain
%
%   This file has no equivalent for TELECOM201/ICS905
%
%   The purpose of this script is to demonstrate a complete wireless
%   communication chain.
%   As provided, the chain is not optimized ; ONE OF THE TARGET OF THE
%   PROJECT IS TO OPTIMIZE THE PARAMETERS OF THE SYSTEM SUBBLOCKS IN ORDER
%   TO MEET THE REQUIREMENTS.
%   Please pay attention to the fact that FIR filters in this script cause 
%   transient phenomena that have not been compensated. 
%   IT IS YOUR DUTY TO FIND THE EXPRESSION OF THE TOTAL DELAY IN "THE
%   ANALOG TO DIGITAL CONVERSION" SECTION TO MATCH THE TRANSMIT AND RECEIVE
%   DATA. (no hardcoding please...)
%   
%   Hopefully, the script has been written to be self-explanatory. 
%
%   In this project, we use the Quadriga Channel Model ("QuaDRiGa") for
%   generating realistic radio channel impulse responses.
%   It is an open source project whose license can be found in the
%   directories of the project. In order to reduce the footprint in terms
%   of disk space, we do not distribute the documentation (PDF) of the
%   QuaDRiGa toolbox included in the original archive distributed on the
%   official site. Please consult the download page and download the
%   archive to retrieve the documentation.
%   https://quadriga-channel-model.de/
%
% Other m-files required:   ./subblocks/*.m
% Subfunctions:             ./subblocks/*.m
% MAT-files required: none
%
% Author: Germain PHAM, Chadi JABBOUR
% C2S, COMELEC, Telecom Paris, Palaiseau, France
% email address: dpham@telecom-paris.fr
% Website: https://c2s.telecom-paristech.fr/TODO
% Feb. 2020, Apr. 2020, Mar. 2022, Dec. 2023
%------------- BEGIN CODE --------------

addpath(genpath('./subblocks/'))
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Transmitter                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% % Signal Generation %%%
continuousTimeSamplingRate    = 20e9;     % A sampling rate which is sufficiently high in order to be close to the continous time
basebandSamplingRate_or       = 30e6;     % The sampling rate of the complex baseband signal ; Units : Samples/second
                                          % in this project it MUST BE a multiple of symbolRate

basebandSamplingRate          = continuousTimeSamplingRate/round(continuousTimeSamplingRate/basebandSamplingRate_or);


%% % Signal Characteristics %%%
symbolRate              = 15e6;  % The raw symbol rate : the raw complex QAM symbols are sent at this rate ; Units : Symbols/second
basebandOverSampling    = round(basebandSamplingRate/symbolRate);
NSamples_BB             = 1e3;   % Signal length (after RRC filter)

% Signal frequencies
freqVin_or1   = 7.12e6;
freqVin_or2   = 6.12e6;
% Place the sine waves in an FFT bin
freqVin1      = round(freqVin_or1/basebandSamplingRate*NSamples_BB)...
                  *basebandSamplingRate/NSamples_BB; 
freqVin2      = round(freqVin_or2/basebandSamplingRate*NSamples_BB)...
                  *basebandSamplingRate/NSamples_BB; 

% Time vector of the simulation
t = 0:1/basebandSamplingRate:(NSamples_BB-1)/basebandSamplingRate;

%% % Baseband (digital) shaping filter %%%
rollOff     = 0.25; % (for RRC filter, single sided output BW is (1+beta)*Rsymb/2 )
symbolSpan  = 25;   % This parameter is related to both the filter length and the attenuation of the stop band
% Instanciate filter
basebandRRC = rcosdesign(rollOff,symbolSpan,basebandOverSampling,'sqrt'); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Possible values of test_type are:
%%%%%                               'onetone' for a one-tone sine 
%%%%%                               'twotone' for a two-tone sine 
%%%%%                               'mod'     for a modulated QAM 16 signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              

test_type='mod';
switch test_type
   case 'onetone'
      %%% one tone signal%%%
      basebandSig = exp(1j*2*pi*freqVin1*t);
   case 'twotone'
      %%% two tone signal%%%
      basebandSig = exp(1j*2*pi*freqVin1*t)+exp(1j*2*pi*freqVin2*t);
   case 'mod'
      %%% Modulated signal
      modSize       = 16; % Modulation order for 16QAM
      nQAMSymbols   = round(NSamples_BB/basebandOverSampling); % Number of QAM symbols to be generated
      inSig         = randi([0 modSize-1],nQAMSymbols,1);      % generate symbols as integer
      % Perform modulation : convert integer symbols to complex symbols
      if isOctave
         qamSig        = qammod(inSig,modSize);
         qamSig        = qamSig/sqrt(mean(abs(qamSig).^2));
      else % Matlab
         qamSig        = qammod(inSig,modSize,'UnitAveragePower',true);
      end
      

      % Apply filter with upsampling to basebandSamplingRate 
      basebandSig   = resample(qamSig,basebandOverSampling,1,basebandRRC);
      % Resample (compared to upfirdn) generates a signal that is exactly the 
      % length we can predict without having to compensate for the delay introduced by the filter
      % https://groups.google.com/d/msg/comp.soft-sys.matlab/UGLNR9vFqhM/c56ZlfUlhhcJ

   otherwise
      %%% one tone signal%%%
      basebandSig = exp(1j*2*pi*freqVin1*t);
end
    
%% System properties
BW  = 10e6;     % Signal bandwidth
BW_rf     = 2*BW;
K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature

%% % IQ separation for real baseband signals %%%
[basebandDigital_I_unorm,basebandDigital_Q_unorm] = complx2cart(basebandSig(:));

%%% Digital to Analog Conversion %%%
nBitDAC = 12;           % Number of bits of the DAC
Vref    = 1;            % Voltage reference of the DAC
dacType = 'zoh';        % DAC type ; can be 'zoh' or 'impulse'

% Normalize signal for conversion
% Must use same scale factor for both wave (Take max of both)
normalize_factor    = max( max(abs(basebandDigital_I_unorm)),...
                           max(abs(basebandDigital_Q_unorm)));
basebandDigital_I   = basebandDigital_I_unorm/normalize_factor*Vref;
basebandDigital_Q   = basebandDigital_Q_unorm/normalize_factor*Vref;

%% Perform conversion
[basebandAnalog_dac_I, basebandAnalog_dac_I_power] = DAC(basebandDigital_I,nBitDAC,Vref,dacType,basebandSamplingRate,continuousTimeSamplingRate);
[basebandAnalog_dac_Q, basebandAnalog_dac_Q_power] = DAC(basebandDigital_Q,nBitDAC,Vref,dacType,basebandSamplingRate,continuousTimeSamplingRate);

fprintf("\tConsumption of DAC I = %fW and Q = %fW\n", basebandAnalog_dac_I_power, basebandAnalog_dac_Q_power)


%% % Baseband FAKE Analog filter %%%
Rin             = 50;    % Input impedance of the filter
TXBB_Filt_NF    = 5;    %(in dB)

% https://fr.mathworks.com/help/signal/ref/firpmord.html
TXBB_Filt_rp    = 0.1;         % Passband ripple in dB
TXBB_Filt_rs    = 40;          % Stopband ripple in dB
TXBB_Filt_fs    = continuousTimeSamplingRate; % Sampling frequency (Hz)
TXBB_Filt_f     = [15 20]*1e6;  % Cutoff frequencies (Hz)
TXBB_Filt_a     = [1 0];        % Desired amplitudes
% Convert the deviations to linear units. 
TXBB_Filt_dev   = [(10^(TXBB_Filt_rp/20)-1)/(10^(TXBB_Filt_rp/20)+1) 10^(-TXBB_Filt_rs/20)];
% Design the filter
[TXBB_Filt_n,TXBB_Filt_fo,TXBB_Filt_ao,TXBB_Filt_w] ...
                = firpmord(TXBB_Filt_f,TXBB_Filt_a,TXBB_Filt_dev,TXBB_Filt_fs);
disp('Designing the TX filter - This takes a while...')

TXBB_Filt       = firpm(TXBB_Filt_n,TXBB_Filt_fo,TXBB_Filt_ao,TXBB_Filt_w);

%% Perform filtering
disp('Filtering with the TX filter - This takes a while...')
basebandAnalog_filt_I = basebandAnalogFiltFake(basebandAnalog_dac_I,TXBB_Filt,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);
basebandAnalog_filt_Q = basebandAnalogFiltFake(basebandAnalog_dac_Q,TXBB_Filt,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);

%% % Mixing up to RF %%%
Flo      = 2.4e9; % Local Oscillator Frequency
[rfSignal, upMixer_power] = upMixer(basebandAnalog_filt_I,basebandAnalog_filt_Q,Flo,continuousTimeSamplingRate);
fprintf("\tConsumption up-mixer = %fW\n", upMixer_power)

%% % RF Amplification %%%
amplifier_index = 5;
[PA_Gain, PA_1dB_out, PA_IIP3, PA_NF, PA_power,Pmax_in] = choose_PA(amplifier_index);

%1. calcul de la puissance du signal dans la bande de 20MHz autour de Flo
voltsq2mwatt    = 1e3/Rin; % Conversion factor from V^2 to milliWatt


freq_min = Flo - 20e6/2;
freq_max = Flo + 20e6/2;
psd = plot_spectrum(rfSignal*sqrt(voltsq2mwatt), 6, continuousTimeSamplingRate, 0); % le dernier zéro empêche le plot
bin_freq_min = fix(freq_min*length(psd)*2/continuousTimeSamplingRate);
bin_freq_max = fix(freq_max*length(psd)*2/continuousTimeSamplingRate);
pow_sig = 10*log10(sum(psd(bin_freq_min:bin_freq_max)));

%2. calcul du gain correctif à appliquer
gain_correcteur = Pmax_in - pow_sig;

%3. on applique le gain négatif
rfSignal_corr = rfSignal * 10^(gain_correcteur/20);

rfPASignal = rfPA(rfSignal_corr,PA_Gain,PA_NF,PA_IIP3,Rin,continuousTimeSamplingRate/2);
fprintf("\tConsumption power amplifier = %fW\n", PA_power)

total_tx_power = basebandAnalog_dac_I_power + basebandAnalog_dac_Q_power + upMixer_power + PA_power;

fprintf("\tTotal consumption TX = %f W\n", total_tx_power)

%% % ACPR %%%
if strcmp(test_type, 'mod')

    [acpr_left, acpr_right,power_sig] = ACPR(BW_rf ,BW_rf + 10e6,Flo,rfPASignal*sqrt(voltsq2mwatt),continuousTimeSamplingRate);
    fprintf("ACPR = %f dB\n", acpr_left)
    fprintf("Power of the signal = %f dB\n", power_sig)

end

%% % Baseband fake Analog filter %%%
RXBB_Filt_NF    = 5;     %(in dB)

RXBB_Filt_rp    = 0.1;         % Passband ripple in dB
RXBB_Filt_rs    = 40;          % Stopband ripple in dB
RXBB_Filt_fs    = continuousTimeSamplingRate; % Sampling frequency (Hz)
RXBB_Filt_f     = [15 20]*1e6;  % Cutoff frequencies (Hz)
RXBB_Filt_a     = [1 0];        % Desired amplitudes
% Convert the deviations to linear units. 
RXBB_Filt_dev   = [(10^(RXBB_Filt_rp/20)-1)/(10^(RXBB_Filt_rp/20)+1) 10^(-RXBB_Filt_rs/20)];
% Design the filter
[RXBB_Filt_n,RXBB_Filt_fo,RXBB_Filt_ao,RXBB_Filt_w] ...
                = firpmord(RXBB_Filt_f,RXBB_Filt_a,RXBB_Filt_dev,RXBB_Filt_fs);
disp('Designing the RX filter - This takes a while...')
RXBB_Filt       = firpm(RXBB_Filt_n,RXBB_Filt_fo,RXBB_Filt_ao,RXBB_Filt_w);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Channel                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_points = 100;
dist = linspace(1.4,1400,nb_points);
BER = zeros(1,nb_points);
for ii = 1:nb_points
    fprintf("itération n° %f",ii)
    %%% Channel %%%
    carrierFreq        = Flo; % Center frequency of the transmission
    c                  = 3e8; % speed of light in vacuum
    distance           = dist(ii); % Distance between Basestation and UE : [1.4,1.4e3] metres
    % Amplitude Attenuation in free space
    ChannelAttenuation = (c/carrierFreq./(4*pi*distance));
    
    fprintf("Channel attenuation = %f dB\n", 20*log10(ChannelAttenuation))
    
    rxSignal           = rfPASignal*ChannelAttenuation;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                     Receiver                        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% % LNA %%%
    LNA_Gain = 15;   % (dB) maximum autorisé
    LNA_IIP3 = 0;  % (dBm) valeur typique
    LNA_NF   = 4;    % (dB) valeur typique 
    
    
    [rfLNASignal, LNAconsum] = rfLNA(rxSignal,LNA_Gain,LNA_NF,LNA_IIP3,Rin,continuousTimeSamplingRate/2);
    fprintf("\tConsumption LNA = %fW\n", LNAconsum)
    
    %% % Mixing down to BB %%%
    [basebandAnalog_raw_I,basebandAnalog_raw_Q,downMixer_power] = downMixer(rfLNASignal,Flo,continuousTimeSamplingRate);
    fprintf("\tConsumption down-mixer = %fW\n", downMixer_power)
    
        
    
    
    %% Perform filtering
    disp('Filtering with the RX filter - This takes a while...')
    basebandAnalog_filtrx_I = basebandAnalogFiltFake(basebandAnalog_raw_I,RXBB_Filt,RXBB_Filt_NF,Rin,continuousTimeSamplingRate);
    basebandAnalog_filtrx_Q = basebandAnalogFiltFake(basebandAnalog_raw_Q,RXBB_Filt,RXBB_Filt_NF,Rin,continuousTimeSamplingRate);
    
    
    %% % Baseband Gain %%%
    BBamp_Gain    = 20; % (dB)
    BBamp_IIP3    = 40; % (dBm)
    BBamp_NF      = 10; % (dB)
    BBamp_band    = 10e6;% (MHz)
    [basebandAnalog_amp_I, BBampPower_I] = BBamp(basebandAnalog_filtrx_I,BBamp_Gain,BBamp_NF,BBamp_IIP3,Rin,BBamp_band,continuousTimeSamplingRate);
    [basebandAnalog_amp_Q, BBampPower_Q] = BBamp(basebandAnalog_filtrx_Q,BBamp_Gain,BBamp_NF,BBamp_IIP3,Rin,BBamp_band,continuousTimeSamplingRate);
    
    fprintf("\tConsumption of Baseband Gain I = %fW and Q = %f W\n", BBampPower_I, BBampPower_Q)
    
    
    
    %% % Analog to Digital Conversion %%%
    nBitADC = 12;
    delay   = (RXBB_Filt_n+TXBB_Filt_n)/2; % WARNING : non trivial value !!! to be thoroughly analyzed
    adcSamplingRate = basebandSamplingRate;
    % Perform conversion
    [basebandAnalog_adc_I, ADC_I_consum] = ADC(basebandAnalog_amp_I,nBitADC,Vref,adcSamplingRate,delay,continuousTimeSamplingRate);
    [basebandAnalog_adc_Q, ADC_Q_consum] = ADC(basebandAnalog_amp_Q,nBitADC,Vref,adcSamplingRate,delay,continuousTimeSamplingRate);
    
    fprintf("\tConsumption of ADC I = %fW and Q = %fW\n", ADC_I_consum, ADC_Q_consum)
    
    tot_consum_Rx = ADC_I_consum + ADC_Q_consum+BBampPower_Q+BBampPower_I+LNAconsum+downMixer_power;
    fprintf("\tTotal consumption RX = %f W\n", tot_consum_Rx)
    
    %% % IQ combination for complex baseband signals %%%
    basebandComplexDigital                = complex(basebandAnalog_adc_I,basebandAnalog_adc_Q);
    
    % RX RRC and downsampling (reverse effect of resample(qamSig...) )
    % WARNING : this downsampling may create unexpected sampling effects due to butterworth filtering and phase distortion
    %           please check signals integrity before and after this step
    basebandComplexDigital_fir            = resample(basebandComplexDigital,1,basebandOverSampling,basebandRRC);
    
    % Normalize received symbols to UnitAveragePower (see qammod(inSig)) 
    % Trully effective when noise and distortions are not too large compared to the useful signal
    basebandComplexDigital_fir            = basebandComplexDigital_fir / sqrt(var(basebandComplexDigital_fir));
    
    % Coarse truncation of the transient parts. 
    % This should be optimized with respect to the filters delays. 
    basebandComplexDigital_fir_truncated  = basebandComplexDigital_fir;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             Plot section            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isOctave
        demod_sig = qamdemod(basebandComplexDigital_fir_truncated,modSize);
    else % Matlab
        demod_sig = qamdemod(basebandComplexDigital_fir_truncated,modSize,'UnitAveragePower',true);
    end
    
    ecart = length(inSig)-length(demod_sig);
    delay = abs(finddelay(demod_sig,inSig));
    inSig2 = inSig(delay+1:end-(ecart-delay));
    [~, ber] = biterr(demod_sig,inSig2);
    BER(ii)=ber;

end
%%

plot(dist,BER);

xlabel('Distance(m)')
ylabel('BER')
title("BER for a modulated signal function of distance")

set(gca, 'FontSize', 20)
set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)



    