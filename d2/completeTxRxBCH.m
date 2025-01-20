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
Nc = 31;
symbolRate              = 15e6;  % The raw symbol rate : the raw complex QAM symbols are sent at this rate ; Units : Symbols/second
basebandOverSampling    = round(basebandSamplingRate/symbolRate);
NSamples_BB             = 3e3;   % Signal length (after RRC filter)

% génération de la matrice de correction pour une erreur
E1 = zeros(1,Nc);
for iii = 1:Nc
    err = zeros(1,Nc);
    err(iii) = 1;
    E1(iii) = str2double(sprintf('%d',bch_simu_2err(err))); %on concatène les bits pour faire un seul entier et pouvoir comparer plus facilement
end

%génération de la matrice de correction pour 2 erreurs
E2 = zeros(1,Nc*(Nc-1)/2);
indexList = zeros(2,Nc*(Nc-1)/2);
u = 1;
for kkk = 1:Nc-1
 for jjj = kkk+1:Nc
    err = zeros(1,Nc);
    err(kkk) = 1;
    err(jjj) = 1;
    E2(u) = str2double(sprintf('%d',bch_simu_2err(err)));
    indexList(u,1) = iii;
    indexList(u,2) = jjj;
    u=u+1;
 end
end

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
mod = "PSK";
  %%% Modulated signal
modSize       = 2; % Modulation order for BPSK

%variables pour générer des mots de longueur 21 sans soucis de taille après
%codage
m = log2(modSize);
a=21;
long_vect = m*floor(a/m);
long = length(bits2symbols(zeros(1,Nc),mod,modSize));
mult = floor(NSamples_BB/basebandOverSampling/long)+1;
Nsamples_round = mult*long;

%Création des vecteurs d'entrée
Mc = [];
Cin = [];

%% génération des vecteurs : taille 21 puis concaténation avec la forme
%systématique

%Ajout d'un vecteur de 1 à l'entrée (pour repérer les décalages plus
%facilement)
m_start = 1 + zeros(1,long_vect);
Mc = cat(2,Mc,m_start);
m_padded_start = [m_start zeros(1,10)];
c = [m_start bch_simu_2err(m_padded_start)];
Cin = cat(2,Cin,c);

for jj=1:mult-2
  mc = randi([0,1], 1,long_vect);
  Mc = cat(2,Mc,mc);
  m_padded = [mc zeros(1,10)];
  c = [mc bch_simu_2err(m_padded)];
  Cin = cat(2,Cin,c);
end

%Ajout d'un vecteur de 1 à la fin pour repérer les décalages
m_end = 1 + zeros(1,long_vect);
Mc = cat(2,Mc,m_end);
m_padded_end = [m_end zeros(1,10)];
c = [m_end bch_simu_2err(m_padded_end)];
Cin = cat(2,Cin,c);

%Modulation du signal
 S = pskmod(Cin,modSize);
 S = S/sqrt(mean(abs(S).^2));

 S = reshape(S,[],1);
  

  % Apply filter with upsampling to basebandSamplingRate 
  basebandSig   = resample(S,basebandOverSampling,1,basebandRRC);
  % Resample (compared to upfirdn) generates a signal that is exactly the 
  % length we can predict without having to compensate for the delay introduced by the filter
  % https://groups.google.com/d/msg/comp.soft-sys.matlab/UGLNR9vFqhM/c56ZlfUlhhcJ

 
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Channel                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Channel %%%
carrierFreq        = Flo; % Center frequency of the transmission
c                  = 3e8; % speed of light in vacuum
distance           = 1.400; % Distance between Basestation and UE : [1.4,1.4e3] metres
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
%%
% Coarse truncation of the transient parts. 
% This should be optimized with respect to the filters delays. 
basebandComplexDigital_fir_truncated = basebandComplexDigital_fir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Plot section            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [S2, ~] = threshold_detector(basebandComplexDigital_fir_truncated,'PSK',modSize);
% S2 = round(real(S2));

M2 = [];

%démodulation du signal
c2 = pskdemod(basebandComplexDigital_fir_truncated,modSize);
c2 = reshape(c2,1,[]);

%calcul du delay et alignement des signaux en fonction du delay
inSig = Cin(32:end);
delayed = finddelay(c2,inSig);
if delayed<0
    outSig = c2(-delayed+1:end);
end

if delayed>=0
    ecart = abs(length(c2)-length(inSig));
    inSig = inSig(delayed+1:end-max(0,(ecart-delayed)));
    outSig = c2(1:end-max(0,delayed-ecart));
end

total_sym = floor(length(outSig)/long);

%Correction du signal
for uu=1:total_sym
        M2 = cat(2,M2,correcteur_2err(outSig((uu-1)*long+1:uu*long),E1,E2,indexList));
end

% 
% ecart = abs(length(M2)-length(Mc));
% delayed = finddelay(M2,Mc);
% if delayed>=0
%     inSig = Mc(delayed+1:end-max(0,(ecart-delayed)));
%     outSig = M2(1:end-max(0,delayed-ecart));
% end

[~, BER] = biterr(Mc(22:end),M2);

fprintf('Taux d''erreur binaire (BER) : %f\n', BER);


%%

window_number       = 1;
lineSpec_index      = 1;
fullband_spectrum   = true;

plot_spectrum(basebandSig,window_number,...
               adcSamplingRate,lineSpec_index,fullband_spectrum);
title('TX Digital Complex recombined signal')

window_number       = window_number+1;
plot_spectrum(basebandComplexDigital,window_number,...
               adcSamplingRate,lineSpec_index,fullband_spectrum);
title('Receiver complex recombined output')

window_number       = window_number+1;
fullband_spectrum   = false;
plot_spectrum(rfPASignal,window_number,...
               continuousTimeSamplingRate,lineSpec_index,fullband_spectrum);
title('PA spectrum')

%% 
if strcmp(test_type, 'mod')
   figure()
   subplot(1,2,1)
   plot(S,'d', LineWidth=1.5)
   set(gca, 'FontSize', 20)

set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)
   title('constellation at TX')
   xlim([-2 2])
   ylim([-0.2 0.2])
   subplot(1,2,2)
   plot(basebandComplexDigital_fir_truncated,'d',LineWidth=1.5)
   set(gca, 'FontSize', 20)

set(gca, 'FontName', 'default')
set(get(gca,'Title'), 'FontSize', 21)
set(get(gca,'xlabel'), 'FontSize', 21)
set(get(gca,'ylabel'), 'FontSize', 21)
   title('constellation at RX')
   xlim([-2 2])
   ylim([-0.2 0.2])
   % WARNING : the received constellation has almost no sense until filters
   %           delays have been thoroughly analyzed and compensated

elseif strcmp(test_type, 'onetone')
    fprintf("SNR en sortie du RX : %f dB\n", snr_freq_db(basebandComplexDigital,basebandSamplingRate, [0 10e6], Flo))
end



%------------- END OF CODE --------------
    