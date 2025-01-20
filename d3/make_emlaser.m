% 
% Returns parameters of an externally-modulated ideal laser.
% 
% Arguments (as a series of 'parameter_name', parameter_value, ...) specify:
% * modulation:		'I' or 'I+Q' (default 'I').
% * P_opt_dBm (dBm):		power of the laser being modulated (mandatory).
% * consumption_laser_base (W):		(default 0)
% * consumption_laser_efficiency:		(default 0.01)
% * consumption_modulator_base (W):		(default 0)
% * consumption_modulator_slope2 (W/Hz²):	(default 100 mW/GHz²)
% 				power consumption parameters for the laser
% 				and modulator, see function TX_optical_eml.


function s = make_emlaser(varargin)
  s = make_param_struct(varargin,...
			{{'modulation', 'I', 'char'}, ...
			 {'P_opt_dBm', '(mandatory)', 'real'}, ...
			 {'consumption_laser_base', 0, 'real', 0}, ...
			 {'consumption_laser_efficiency', 0.01, 'real', 0, 1}, ...
			 {'consumption_modulator_base', 0, 'real', 0}, ...
			 {'consumption_modulator_slope2', 0.1e-18, 'real', 0} ...
			});
end


