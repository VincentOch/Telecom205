function [S_out, Ts_out, powerOfBlock] = TX_optical_dml(X, Ts, laser)
  % Outputs samples of optical field from a directly-modulated laser.
  %
  % Requires package (pkg load miscellaneous).
  %
  % Arguments:
  % * X (vector):	symbols to transmit (current values, positive, A).
  % * Ts (seconds):	symbol period.
  % * laser (struct):	hardware parameters, from make_laser_simple().
  %
  % Return values:
  % * S_out (vector):		optical field values (√W).
  % * Ts_out (seconds):	sampling period.
  % * powerOfBlock (W):	estimation of the power consumed by the hardware.
  %
  % Simplified rate equations are of the form:
  % dN/dt = I/q - N/tau_e - A·(N - N_transp)·P
  % dP/dt = beta_sp·N/tau_e + (A·(N - N_transp) - 1/tau_p)·P
  % where N and P are the populations (not densities) of carriers and photons
  % in the active region; I the current; tau_e and tau_p are the carrier and
  % photon lifetimes; N_transp the number of carriers at transparency; A the
  % global differential gain; and beta_sp the efficiency of spontaneous emission
  % into the waveguide.
  %
  % Note that not only these equations are a gross oversimplification;
  % solving them Euler-style is notoriously unreliable.  This is a
  % simplified version for teaching; do not use for real-life calculations.
  %
  % Power consumption is estimated from the average output power.

  % Change the sample rate to accomodate laser bandwidth: oversample by
  % a factor M to cover at least 10 times B_e and 10 times 1/tau_p.
  M = ceil(max(10 * laser.B_e, 10 / laser.tau_p) * Ts);

  N_samples = length(X) * M;
  Ts_out = Ts / M;
  fs_out = 1 / Ts_out;
  if N_samples > 1e5
    warning('Oversampling factor %d results in %d samples, expect long computation', ...
	    M, N_samples);
  end
  XX = zeros(1, N_samples);
  for ii = 1:M
    XX(ii:M:end) = X;
  end

  % Prepare an order-1 RC filter to limit the input bandwidth.
  F = (0:N_samples-1) * fs_out / N_samples; % Frequency vector.
  F(ceil(N_samples/2)+1:end) = F(ceil(N_samples/2)+1:end) - fs_out;
  RC_filter = 1 ./ (1 + 1i * F / laser.B_e);

  % Do the filtering.  Take the magnitude of the filtered input signal to be
  % the current actually flowing through the laser.
  I = abs(ifft(RC_filter .* fft(XX)));

  % Prepare vectors for the results: N and P.
  N = zeros(1, N_samples);
  P = zeros(1, N_samples);
  phi = zeros(1, N_samples);

  % Local helper values, including well-known physical constants.
  c = 3e8;
  h = 6.626e-34;
  q = 1.6e-19;
  Ntr_Nth = laser.N_transp / laser.N_th;
  inv_tau_p = 1 / laser.tau_p;

  % Initialize the current value of N and P using quasi-static solution
  % of rate equations with first value of X.
  epsilon = 0.5/(1-laser.beta_sp) * ...
	    (1 - I(1)/laser.I_th - laser.beta_sp * (2 - Ntr_Nth) ...
	     + sqrt((I(1) / laser.I_th - 1 ...
		     + laser.beta_sp * (2 - Ntr_Nth)) .^ 2 ...
		    + 4 * laser.beta_sp * (1 - laser.beta_sp) * (1 - Ntr_Nth)));
  N_current = laser.N_th * (1 - epsilon);
  P_current = laser.beta_sp / laser.A / laser.tau_e * (1 ./ epsilon - 1);
  phi_current = 0; % Arbitrary phase reference.

  % Iterate over samples, solving rate equations Euler-style.
  for ii = 1:N_samples
    N_tau_e = N_current / laser.tau_e;
    A_N_m_N_transp = laser.A * (N_current - laser.N_transp);
    dN_dt = I(ii) / q - N_tau_e - A_N_m_N_transp * P_current;
    dP_dt = laser.beta_sp * N_tau_e + (A_N_m_N_transp - inv_tau_p) * P_current;
    dphi_dt = laser.alpha_H * dP_dt / (P_current + 1);
    N_current = N_current + dN_dt * Ts_out;
    P_current = P_current + dP_dt * Ts_out;
    phi_current = phi_current + dphi_dt * Ts_out;

    % Before storing values, articifially ensure N and P are positive
    % and phi within [0, 2pi] to somewhat avoid diverging.
    % N_current = clip(N_current, [0, inf]);
    N_current(N_current<0) =0;
    N_current(N_current>inf) = inf;
    % P_current = clip(P_current, [0, inf]);
    P_current(P_current<0) =0;
    P_current(P_current>inf) = inf;
    
    phi_current = mod(phi_current, 2 * pi);
    N(ii) = N_current;
    P(ii) = P_current;
    phi(ii) = phi_current;
  end

  % Output power emitted from front facet: half lineic photon density ×
  % speed of light in material × facet transmittivity × photon energy.
  P_out = 0.5 * P / laser.L * c / laser.n_g * (1-laser.R_front) ...
	  * h * laser.v;

  % Output signal: optical field.
  S_out = sqrt(P_out) .* exp(1i * phi);

  % Estimated consumption: base consumption + average output power / efficiency.
  powerOfBlock = laser.consumption_laser_base ...
		 + sum(P_out) / length(P_out) ...
		   / laser.consumption_laser_efficiency;
end
