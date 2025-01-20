% Permet de trouver le bon ENOB par le test
% Utilisation : 
% - il faut lancer une première fois completeTxRx_proj
% - choisir le nombre de bits de l'ADC (en-tête de la boucle parfor)

parfor nBitADC=10:16
    %% % LNA %%%
    LNA_Gain = 15;   % (dB) maximum autorisé
    LNA_IIP3 = 100;  % (dBm) à déterminer !
    LNA_NF   = 4;    % (dB) valeur typique
    %%
    [rfLNASignal, LNAconsum] = rfLNA(rxSignal,LNA_Gain,LNA_NF,LNA_IIP3,Rin,continuousTimeSamplingRate/2);
    fprintf("\tConsumption LNA = %fW\n", LNAconsum)

    %% % Mixing down to BB %%%
    [basebandAnalog_raw_I,basebandAnalog_raw_Q,downMixer_power] = downMixer(rfLNASignal,Flo,continuousTimeSamplingRate);
    %fprintf("\tConsumption down-mixer = %fW\n", downMixer_power)


    %% % Baseband fake Analog filter %%%
    RXBB_Filt_NF    = 0;     %(in dB)

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
    %disp('Designing the RX filter - This takes a while...')
    RXBB_Filt       = firpm(RXBB_Filt_n,RXBB_Filt_fo,RXBB_Filt_ao,RXBB_Filt_w);

    %% Perform filtering
    %disp('Filtering with the RX filter - This takes a while...')
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
delay   = (RXBB_Filt_n+TXBB_Filt_n)/2; % WARNING : non trivial value !!! to be thoroughly analyzed
adcSamplingRate = basebandSamplingRate;
% Perform conversion
[basebandAnalog_adc_I, ADC_I_consum] = ADC(basebandAnalog_amp_I,nBitADC,Vref,adcSamplingRate,delay,continuousTimeSamplingRate);
[basebandAnalog_adc_Q, ADC_Q_consum] = ADC(basebandAnalog_amp_Q,nBitADC,Vref,adcSamplingRate,delay,continuousTimeSamplingRate);

fprintf("\tENOB %d Consumption of ADC I = %fW and Q = %fW\n", nBitADC, ADC_I_consum, ADC_Q_consum)

tot_consum_Rx = ADC_I_consum + ADC_Q_consum+BBampPower_Q+BBampPower_I+LNAconsum+downMixer_power;

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
basebandComplexDigital_fir_truncated  = basebandComplexDigital_fir(10:end-10);

    fprintf("ENOB %d SNR en sortie du RX : %f dB\n",nBitADC, snr_freq_db(basebandComplexDigital,basebandSamplingRate, [0 10e6], Flo))
    
end