function [S_out, Ts_out, powerOfBlock] = TX_optical_eml(X, Ts, params)
  % Outputs samples of optical field from an externally-modulated laser.
  %
  % Requires package (pkg load miscellaneous).
  %
  % Arguments:
  % * X (vector):		symbols to transmit.
  % * Ts (seconds):		sample period (temps d'échantillonage du vecteur
  % d'entrée)
  % * params (struct):		hardware parameters, from make_emlaser().
  %
  % Return values:
  % * S (vector):		samples of the optical field, in units of √W.
  % * Ts_out (seconds):	sampling period.
  % * powerOfBlock (W):	estimation of the power consumed by the hardware.
  %
  % If params.modulation is 'I', only the intensity is modulated;
  % if params.modulation is 'I+Q', intensity and phase are modulated.
  % 
  % X is a vector of symbols.  Each returned sample in S_out is proportional
  % to the corresponding symbol (in I+Q mode) or its magnitude (in I mode).
  % A symbol value of 1 yields an optical field with power params.P_opt_dBm.
  % If a symbol's magnitude is greater than 1, it will be clipped to 1.
  %
  % Ts is the symbol rate in seconds.  While the resulting samples do not
  % depend on it per se, it is used in estimating the power consumed:
  % a fixed base for laser and modulator from params.consumption_laser_base
  % and params.consumption_modulator_base, + a component proportional to
  % the laser emitted power and inverse of params.consumption_laser_efficiency,
  % + a component proportional to params.consumption_modulator_slope2 and
  % the square of the symbol rate (taken to be the operating frequency).
  % Modulator power is counted twice in I+Q mode.

  % Calculate the amplitude corresponding to symbol value 1.
  % Keep the same sample rate.
  laser_power = 10.^((params.P_opt_dBm - 30) / 10);
  amplitude_scale = sqrt(laser_power);
  Ts_out = Ts;

  % Separate symbols along magnitude and phase, clip magnitude,
  % reconstitute and scale amplitude.
  S_out = min(abs(X),1); %clip(abs(X), 1); % clip is an octave-specific function
  if (strcmp(params.modulation, 'I+Q'))
    % Use phase in I+Q mode only.
    S_out = S_out .* exp(1i * angle(X));
  end
  S_out = S_out * amplitude_scale;

  % Power of block: use laser base consumption and wall-plug efficiency,
  % and modulator base consumption and a component in the square of the
  % operating frequency, ×2 in I+Q mode.
  powerOfBlock = params.consumption_modulator_base ...
		 + params.consumption_modulator_slope2 / (Ts .^ 2);
  if (strcmp(params.modulation, 'I+Q'))
    powerOfBlock = powerOfBlock * 2;
  end
  powerOfBlock = powerOfBlock + params.consumption_laser_base ...
		 + laser_power / params.consumption_laser_efficiency;
end

