function [Gain,P1dB_out,IPP3,NF,PowerConsumption,Pmax_in] = choose_PA(pa_index)
    % The Power Amplifier index is a integer as follows:
    % 1. ZX60-V63+
    % 2. ZX60-V62+
    % 3. ZHL-42   
    % 4. RFLUPA05M06G
    % 5. ADL5606
    
    Pout = 20; % Output power must be 20 dBm

    switch pa_index
        case 1 % 
            Gain = 20;
            IPP3 = 30 - Gain;
            NF = 3.7;
            P1dB_out = 17.8;
            PowerConsumption = 5 * 69e-3;
        case 2
            Gain = 15.4;
            IPP3 = 33 - Gain;
            NF = 5.1;
            P1dB_out = 18.5;
            PowerConsumption = 5*82e-3;
        case 3
            Gain = 33;
            IPP3 = 38 - Gain;
            NF = 7.55;
            P1dB_out = 28.3;
            PowerConsumption = 15 * 880e-3;
        case 4
            Gain = 33;
            IPP3 = 40.5 - Gain;
            NF = 3;
            P1dB_out = 30;
            PowerConsumption = 12 * 280e-3;
        case 5
            Gain = 24.3;
            IPP3 = 45.5 - Gain;
            NF = 4.7;
            P1dB_out = 30.8;
            PowerConsumption = 5 * 362e-3;
    end

    if P1dB_out > Pout % c'est possible de sortir un signal assez puissant (20dBm)
        Pmax_in = Pout - Gain; % Pmax_in est la puissance max qu'il faut envoyer au PA pour qu'il sorte 20dBm
    else
        fprintf("Output power will be below requirements (%f), P1dB_(out) = %f dBm", Pout, P1dB_out)
        Pmax_in = P1dB_out - Gain + 1; %  on considère la compression à 1dB
    end
end

