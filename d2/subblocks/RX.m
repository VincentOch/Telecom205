function [basebandAnalog_adc_I, basebandAnalog_adc_Q] = RX(rxSignal,LNA_IIP3,LNA_NF,LNA_Gain,Flo,continuousTimeSamplingRate,adcSamplingRate,nBitADC,Vref,delay,BB_gain)
%RX - A shortcut function to emulate the RX chain - TELECOM201/ICS905 version
%   Calls the functions rfLNA, downMixer, basebandAnalogFilt and ADC
%
% Syntax:  [basebandAnalog_adc_I, basebandAnalog_adc_Q]=RX(rxSignal,LNA_IIP3,LNA_NF,LNA_Gain,Flo,continuousTimeSamplingRate,adcSamplingRate,nBitADC,Vref,delay,BB_gain)
%
% Inputs:
%    rxSignal                   - RX signal at antenna port
%    LNA_IIP3                   - LNA IIP3 (dBm)
%    LNA_NF                     - LNA NF (dB)
%    LNA_Gain                   - LNA Gain (dB)
%    Flo                        - LO frequency (Hz)
%    continuousTimeSamplingRate - simulation step rate (Hz)
%    adcSamplingRate            - ADC sampling rate (Hz)
%    nBitADC                    - number of bits of the ADC
%    Vref                       - reference voltage of the ADC
%    delay                      - delay to compensate for the filter group delay
%    BB_gain                    - baseband gain (dB)
%
% Outputs:
%    basebandAnalog_adc_I       - baseband analog equivalent I signal after ADC
%    basebandAnalog_adc_Q       - baseband analog equivalent Q signal after ADC
%
% Other m-files required: rfLNA, downMixer, basebandAnalogFilt, ADC
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


Rin=50;

rfLNASignal = rfLNA(rxSignal,LNA_Gain,LNA_NF,LNA_IIP3,Rin,continuousTimeSamplingRate/2);

%%% Mixing down to BB %%%
[basebandAnalog_raw_I,basebandAnalog_raw_Q] = downMixer(rfLNASignal,Flo,continuousTimeSamplingRate);

%%% Baseband Analog filter %%%
RXBB_Filt_NF    = 0;    %(in dB)
RXBB_Filt_Fcut  = 15e6;  % Filter TX BB Fcut 3dB Frequency
RXBB_Filt_Order = 6;    % Filter TX BB Order
% Instanciate filter (due to numerical issue, the filter has to be instanciated as SOS)
[RXBB_Filt_z,RXBB_Filt_p,RXBB_Filt_k]=butter(RXBB_Filt_Order,RXBB_Filt_Fcut/(continuousTimeSamplingRate/2));
RXBB_Filt_sos = zp2sos(RXBB_Filt_z,RXBB_Filt_p,RXBB_Filt_k);
% Perform filtering
basebandAnalog_filtrx_I = basebandAnalogFilt(basebandAnalog_raw_I,RXBB_Filt_sos,RXBB_Filt_NF,Rin,continuousTimeSamplingRate);
basebandAnalog_filtrx_Q = basebandAnalogFilt(basebandAnalog_raw_Q,RXBB_Filt_sos,RXBB_Filt_NF,Rin,continuousTimeSamplingRate);

BB_gain_lin=10^(BB_gain/20);
% Perform conversion
basebandAnalog_adc_I = ADC(BB_gain_lin*basebandAnalog_filtrx_I,nBitADC,Vref,adcSamplingRate,delay,continuousTimeSamplingRate);
basebandAnalog_adc_Q = ADC(BB_gain_lin*basebandAnalog_filtrx_Q,nBitADC,Vref,adcSamplingRate,delay,continuousTimeSamplingRate);


