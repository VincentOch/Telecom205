%
% Returns parameters of a simple semiconductor laser modeled by basic rate equations.
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
% Note that A and N_transp must be calculated as a function of the volume V
% of the laser's active region: N_transp = V·n_transp (with n_transp being
% the carrier density at transparency); and A = Gamma·g_N / V (with g_N
% being the material differential gain, and Gamma the confinement factor).
%
% Arguments (as a series of 'parameter_name', parameter_value, ...) specify:
% * v (mandatory):	Optical frequency (Hz).
% * n_g:		Group index (default 4).
% * B_e:		Electrical current bandwidth (Hz, default 2.5 GHz).
% * L:			Laser cavity length (m, default 300e-6 m).
% * s:			Cross-section of the active region (m², default 0.2e-12).
% * R_front:		Power reflectivity of emitting facet (default 0.3).
% * R_back:		Power reflectivity of non-emitting facet (default 1).
% * tau_e:		Carrier lifetime (s, default 2 ns).
% * loss:		Material loss (1/m, default 40/cm).
% * beta_sp:		Spontaneous emission efficiency (default 1e-4).
% * n_transp:		Carrier density at transparency (1/m³, default 1e24).
% * gain_diff:		Differential gain (Gamma·g_N, m³/s, default 1e-12).
% * alpha_H:		Linewidth enhancement factor (default 5).
% * consumption_laser_base (W):		(default 0)
% * consumption_laser_efficiency:		(default 0.01)
%			power consumption parameters.
%

function l = make_laser_simple(varargin)
  %% Make a parameters struct according to arguments and default values.
  l = make_param_struct(varargin,...
			{{'v', '(mandatory)', 'real', 0}, ...
			 {'n_g', 4, 'real', 1}, ...
			 {'B_e', 2.5e9, 'real', 0}, ...
			 {'L', 300e-6, 'real', 0}, ...
			 {'s', 0.2e-12, 'real', 0}, ...
			 {'R_front', 0.3, 'real', 0, 1}, ...
			 {'R_back', 1, 'real', 0, 1}, ...
			 {'tau_e', 2e-9, 'real', 0}, ...
			 {'loss', 40e2, 'real', 0}, ...
			 {'beta_sp', 1e-4, 'real', 0, 1}, ...
			 {'n_transp', 1e24, 'real', 0}, ...
			 {'gain_diff', 1e-12, 'real'}, ...
			 {'alpha_H', 5, 'real'}, ...
			 {'consumption_laser_base', 0, 'real', 0}, ...
			 {'consumption_laser_efficiency', 0.01, 'real', 0, 1} ...
			});

  %% Calculate helper values: wavelength, active region volume, global
  %  differential gain and carrier population at transparency.
  c = 3e8;
  q = 1.6e-19;
  l.lambda = c / l.v;
  l.V = l.L * l.s;
  l.A = l.gain_diff / l.V;
  l.N_transp = l.n_transp * l.V;

  %% Calculate photon lifetime, threshold current and carrier density
  %  at threshold.
  l.tau_p = l.n_g / c / (l.loss - log(l.R_front * l.R_back) / 2 / l.L);
  l.N_th = l.N_transp + 1 / l.A / l.tau_p;
  l.I_th = l.N_th * q / l.tau_e;
end

