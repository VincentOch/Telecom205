function rfPASignal = TX_proj(basebandDigital_I,basebandDigital_Q,Vref, nBitDAC,basebandSamplingRate,continuousTimeSamplingRate,dacType,TXBB_Filt,TXBB_Filt_NF,Flo,PA_IIP3,PA_NF,PA_Gain,PA_power,Pmax_in)
    %TX_proj - A shortcut function to emulate the TX chain - TELECOM205 version
    %   Calls the functions DAC, basebandAnalogFiltFake, upMixer and rfPA
    %
    % Syntax:  rfPASignal = TX_proj(basebandDigital_I,basebandDigital_Q,Vref, nBitDAC,basebandSamplingRate,continuousTimeSamplingRate,dacType,TXBB_Filt,TXBB_Filt_NF,Flo,PA_IIP3,PA_NF,PA_Gain,PA_power,Pmax_in)
    %
    % Inputs:
    %    basebandDigital_I          - baseband analog equivalent I signal before DAC
    %    basebandDigital_Q          - baseband analog equivalent Q signal before DAC
    %    Vref                       - reference voltage of the DAC
    %    nBitDAC                    - number of bits of the DAC
    %    basebandSamplingRate       - baseband sampling rate (Hz)
    %    continuousTimeSamplingRate - simulation step rate (Hz)
    %    dacType                    - DAC type (string)
    %    TXBB_Filt                  - TX baseband filter impulse response ; FIR ONLY
    %    TXBB_Filt_NF               - TX baseband filter NF (dB)
    %    Flo                        - LO frequency (Hz)
    %    PA_IIP3                    - PA IIP3 (dBm)
    %    PA_NF                      - PA NF (dB)
    %    PA_Gain                    - PA Gain (dB)
    %
    % Outputs:
    %    rfPASignal                 - RF PA signal
    %
    % Other m-files required: DAC, basebandAnalogFiltFake, upMixer, rfPA
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: none
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Dec. 2023
    %------------- BEGIN CODE --------------


Rin = 50;

% Perform conversion
[basebandAnalog_dac_I, basebandAnalog_dac_I_power] = DAC(basebandDigital_I,nBitDAC,Vref,dacType,basebandSamplingRate,continuousTimeSamplingRate);
[basebandAnalog_dac_Q, basebandAnalog_dac_Q_power] = DAC(basebandDigital_Q,nBitDAC,Vref,dacType,basebandSamplingRate,continuousTimeSamplingRate);

fprintf("Consumption of DAC I = %fW and Q = %fW\n", basebandAnalog_dac_I_power, basebandAnalog_dac_Q_power)

% Perform filtering
disp('Filtering the DAC output - This takes a while...')
basebandAnalog_filt_I = basebandAnalogFiltFake(basebandAnalog_dac_I,TXBB_Filt,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);
basebandAnalog_filt_Q = basebandAnalogFiltFake(basebandAnalog_dac_Q,TXBB_Filt,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);

%% Mixing up to RF %%%
[rfSignal, upMixer_power] = upMixer(basebandAnalog_filt_I,basebandAnalog_filt_Q,Flo,continuousTimeSamplingRate);
fprintf("Consumption up-mixer = %fW\n", upMixer_power)


%% RF Amplification %%%
% Il faut atténuer le signal pour que la puissance soit égale à Pmax_in
%1. calcul de la puissance du signal dans la bande de 20MHz autour de Flo
freq_min = Flo - 30e6/2;
freq_max = Flo + 30e6/2;
%psd = abs(fft(rfSignal)).^2;
psd = plot_spectrum(rfSignal*sqrt(1e3/50), 4, continuousTimeSamplingRate);
bin_freq_min = fix(freq_min*length(psd)*2/continuousTimeSamplingRate);
bin_freq_max = fix(freq_max*length(psd)*2/continuousTimeSamplingRate);
pow_signal = sum(psd(bin_freq_min:bin_freq_max));


%2. calcul du gain correctif à appliquer
gain_correcteur = Pmax_in - pow_signal;

%3. on applique le gain négatif
rfSignal_corr = rfSignal * 10^(gain_correcteur/20);

rfPASignal = rfPA(rfSignal_corr,PA_Gain,PA_NF,PA_IIP3,Rin,continuousTimeSamplingRate/2);
fprintf("Consumption power amplifier = %fW\n", PA_power)

% On ignore une partie des composants dans le calcul de la puissance
% (comme les filtres)
total_tx_power = basebandAnalog_dac_I_power + basebandAnalog_dac_Q_power + upMixer_power + PA_power;

fprintf("Total consumption TX = %fW\n", total_tx_power)