function [S, Ts, powerOfBlock,SNR_elec] = RX_photodetector(S_opt, Ts_opt, N_opt, params)
  % Output samples of (noisy) electrical signal converted from optical.
  %
  % Requires package "miscellaneous" (pkg load miscellaneous).
  %
  % Arguments:
  % * S_opt (vector):		optical field samples in units of √W.
  % * Ts_opt (seconds):	sample period.
  % * N_opt (W/Hz):		optical noise PSD.
  % * params (struct):		see make_photodetector().
  %
  % S is the vector of photocurrent samples converted from S_opt,
  % with a photoreceiver sensitivity of params.sensitivity, and
  % noise added according to the optical noise given
  % and the thermal noise given in params.
  %
  % Ts is the output sample rate in seconds; currently always equal to Ts_opt.
  %
  % powerOfBlock is an estimation of the hardware power consumption.
  % Currently a constant given by params.P_consumption.
  %
  % SNR_elec is the electrical SNR, in dB, computed after photodetection

  Ts = Ts_opt;                            % Same input/output sample rates.
  I = abs(S_opt).^2 * params.sensitivity; % Photocurrent: sensitivity * Popt.

  % Local helper values.
  N_samples = length(S_opt);
  fs = 1 / Ts;

  % Limit the bandwidth using a first-order filter (modeling capacitive effects),
  % except that we don't bother if the signal's own bandwidth is less than 10% B_e.
  if fs >= 0.1 * params.B_e
    %% Signal may be too fast compared to the photoreceiver's bandwidth B_e,
    % filter it.
    F = (0:N_samples-1) * fs / N_samples;   % Frequency vector.
    F = reshape(F, size(I));
    F(ceil(N_samples/2)+1:end) = F(ceil(N_samples/2)+1:end) - fs;

    %% Make a raised-cosine filter with a 0.5 roll-off around B_e
    % in the frequency domain.
    F_temp = F - max(min(F,0.5*params.B_e),-0.5*params.B_e);
    F_limit = max(min(F_temp,params.B_e),-params.B_e);
    RCos_filter = 0.5*(1+cos(F_limit/params.B_e*pi));

    %% Filter the signal in the frequency domain.
    S = ifft(RCos_filter .* fft(I));
    S = real(S); % Should already be real, but squash rounding errors.
  else
    S = I; % Signal already slow, don't bother filtering it.
  end

  %% Calculate noise variance in terms of photocurrent.
  % * thermal noise: N_th × electrical bandwidth / R_load.
  % * noise-noise beat: 2 × polarization modes × (sensitivity × noise)² ×(B_opt-Be/2)B_e.
  % * signal-noise beat: 4 × sensitivity² × optical power × noise × electrical bandwidth.
  % (Use already-calculated I = sensitivity * optical power.)
  sigma2 = params.N_th * params.B_e / params.R_load ...
	   + 2 * 2 * ((params.sensitivity * N_opt) ^ 2) ...
	     * (params.B_opt - params.B_e/2) * params.B_e ...
	   + (4 * params.sensitivity * N_opt * params.B_e) * I;
  SNR_elec = 10*log10(mean(abs(S).^2)./mean(sigma2));
  S = S + sqrt(sigma2) .* randn(size(S));
  
  %% Power of block: currently a constant.
  powerOfBlock = params.P_consumption;
end
