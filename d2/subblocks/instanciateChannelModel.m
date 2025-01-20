function [firRFchannel] = instanciateChannelModel(distance,Fc,Fs)
    %instanciateChannelModel - generate a FIR channel filter
    % 
    % The generation is made using the Quadriga Channel Model ("QuaDRiGa") for 
    % realistic radio channel impulse responses.
    % For further information, please browse the download page and get the archive 
    % to retrieve the documentation : 
    %   https://quadriga-channel-model.de/
    % The scripts outputs a simplified version of the actual impulse response 
    % provided by QuaDRiGa. It uses only the real part of the coefficients. 
    %
    % Syntax:  [firRFchannel] = instanciateChannelModel(distance,Fc,Fs)
    %
    % Inputs:
    %    distance   - Distance between the TX and RX chain (in m)
    %    Fc         - Center frequency (in Hz)
    %    Fs         - Sampling frequency of the simulation (in Hz)
    %
    % Outputs:
    %    firRFchannel - The sampled FIR impulse response of the channel
    %
    % Example: 
    %    [firRFchannel] = instanciateChannelModel(1e3,2.14e9,20e9)
    %
    % Other m-files required: ./channelModel/*
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Feb. 2020
    %------------- BEGIN CODE --------------
    
    s = qd_simulation_parameters;                           % Set up simulation parameters
    s.show_progress_bars = 0;                               % Disable progress bars
    s.center_frequency = Fc;                                % Set center frequency
    l = qd_layout(s);                                       % Create new QuaDRiGa layout

    % l.set_scenario('3GPP_38.901_Indoor_NLOS');
    % % For typical indoor deployments such as WiFi or femto-cells covering carrier 
    % % frequencies from 500 MHz to 100 GHz. 
    % % The BS antenna height is fixed to 3 meters and the MT antenna height to 1 meter. 
    l.set_scenario('WINNER_Indoor_A1_NLOS');
    % WINNER Indoor Hotspot
    % For typical indoor deployments such as WiFi or femto-cells.


    l.no_tx = 1;                                            % Number of Tx
    l.tx_position = [ 0 ; 0; 3 ];                           % Position the Tx to (0,0) at 3m height

    l.no_rx = 1;                                            % Number of Rx
    l.rx_position = [ 0 ; distance ; 3 ];                   % Position the Rx to (0,distance) at 3m height

    l.tx_array = qd_arrayant('omni') ;
    l.rx_array = qd_arrayant('omni') ;

    cb = l.init_builder;                                    % Initialize builder
    gen_parameters( cb );                                   % Generate small-scale-fading
    c  = get_channels( cb );                                % Get channel coefficients
    % disp( char( c.name ) )                                % Show the names of the channels

    % Convert to FIR filter, sample the impulse response, 
    % use only real coefficients for sake of simplicity
    selectNpaths  = size(c.coeff,3);
    channelCoeffs = squeeze(c.coeff(1,1,1: selectNpaths,1));
    channelDelays = squeeze(c.delay(1,1,1: selectNpaths,1));

    channelLen    = fix(max(channelDelays)*Fs)+1;
    channelbins   = round(channelDelays*Fs)+1;

    firRFchannel              = zeros(1,channelLen);
    firRFchannel(channelbins) = real(channelCoeffs);


    %------------- END OF CODE --------------
    