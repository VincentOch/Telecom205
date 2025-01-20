function rfPASignal=TX(basebandDigital_I,basebandDigital_Q,Vref, nBitDAC,basebandSamplingRate,continuousTimeSamplingRate,dacType,TXBB_Filt_Order,TXBB_Filt_NF,TXBB_Filt_Fcut,Flo,PA_IIP3,PA_NF,PA_Gain)
%TX - A shortcut function to emulate the TX chain - TELECOM201/ICS905 version
%   Calls the functions DAC, basebandAnalogFiltFake, upMixer and rfPA
%
% Syntax:  rfPASignal = TX(basebandDigital_I,basebandDigital_Q,Vref, nBitDAC,basebandSamplingRate,continuousTimeSamplingRate,dacType,TXBB_Filt_Order,TXBB_Filt_NF,TXBB_Filt_Fcut,Flo,PA_IIP3,PA_NF,PA_Gain)
%
% Inputs:
%    basebandDigital_I          - baseband analog equivalent I signal before DAC
%    basebandDigital_Q          - baseband analog equivalent Q signal before DAC
%    Vref                       - reference voltage of the DAC
%    nBitDAC                    - number of bits of the DAC
%    basebandSamplingRate       - baseband sampling rate (Hz)
%    continuousTimeSamplingRate - simulation step rate (Hz)
%    dacType                    - DAC type (string)
%    TXBB_Filt_Order            - TX baseband filter order
%    TXBB_Filt_NF               - TX baseband filter NF (dB)
%    TXBB_Filt_Fcut             - TX baseband filter cut-off frequency (Hz)
%    Flo                        - LO frequency (Hz)
%    PA_IIP3                    - PA IIP3 (dBm)
%    PA_NF                      - PA NF (dB)
%    PA_Gain                    - PA Gain (dB)
%
% Outputs:
%    rfPASignal                 - RF PA signal
%
% Other m-files required: DAC, basebandAnalogFilt, upMixer, rfPA
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


% Perform conversion
basebandAnalog_dac_I = DAC(basebandDigital_I,nBitDAC,Vref,dacType,basebandSamplingRate,continuousTimeSamplingRate);
basebandAnalog_dac_Q = DAC(basebandDigital_Q,nBitDAC,Vref,dacType,basebandSamplingRate,continuousTimeSamplingRate);


% Instanciate filter (due to numerical issue, the filter has to be instanciated as SOS)
[TXBB_Filt_z,TXBB_Filt_p,TXBB_Filt_k]=butter(TXBB_Filt_Order,TXBB_Filt_Fcut/(continuousTimeSamplingRate/2));
TXBB_Filt_sos = zp2sos(TXBB_Filt_z,TXBB_Filt_p,TXBB_Filt_k);
% Perform filtering
basebandAnalog_filt_I = basebandAnalogFilt(basebandAnalog_dac_I,TXBB_Filt_sos,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);
basebandAnalog_filt_Q = basebandAnalogFilt(basebandAnalog_dac_Q,TXBB_Filt_sos,TXBB_Filt_NF,Rin,continuousTimeSamplingRate);

%%% Mixing up to RF %%%
rfSignal = upMixer(basebandAnalog_filt_I,basebandAnalog_filt_Q,Flo,continuousTimeSamplingRate);

%%% RF Amplification %%%
rfPASignal = rfPA(rfSignal,PA_Gain,PA_NF,PA_IIP3,Rin,continuousTimeSamplingRate/2);