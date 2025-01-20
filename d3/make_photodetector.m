%
% Returns parameters of a photodetector.
%
% Arguments (as a series of 'parameter_name', parameter_value, ...) specify:
% * sensitivity (A/W):	photocurrent per optical power (default 1 A/W).
% * B_e (Hz):		electrical bandwidth in baseband (default 10 GHz).
% * B_opt (Hz):	optical bandwidth around carrier (default 125 GHz).
% * R_load (Ω):	impedance of the output load (default 50 Ω).
% * N_th (W/Hz):	thermal noise PSD (default 2e-20 W/Hz: T_eq ≈ 400 K).
% * P_consumption (W):	overall power consumption (default 100 mW).
%

function s = make_photodetector(varargin)
  s = make_param_struct(varargin,...
			{{'sensitivity', 1, 'real'}, ...
			 {'B_e', 1e10, 'real', 0}, ...
			 {'B_opt', 125e9, 'real', 0}, ...
			 {'R_load', 50, 'real', 0}, ...
			 {'N_th', 2e-20, 'real', 0}, ...
			 {'P_consumption', 0.1, 'real', 0}});
end
