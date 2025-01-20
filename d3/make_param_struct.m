%
% Returns a struct containing the given parameters, checking their values.
%
% Arguments:
% * parameters:		cell array containing the parameters to set;
% * constraints:	cell array defining the allowed parameters and values.
%
% The argument "parameters" is of the form {name1, value1, name2, value2, ...}
% Each name (string) is a field name to set to the corresponding value
% in the resulting struct.
%
% The argument "constraints" lists what is allowed.  Each element of
% "constraints" must be a cell array that gives a parameter name and
% constraints on its value as explained below; or a string that is a
% parameter name with no constraints.  A string s is processed as if
% the single-element array {s} was given.
%
% For each array c in "constraints", the elements of c should be:
%
% * c{1}: a string, the name of the parameter concerned.  The special
%         name '...' indicates that parameters may be given with
%         other names than those given in the constraints list.
%
% * c{2}: the default value to assign to the parameter if not given
%         in the "parameters" array; alternatively, special values:
%         '(mandatory)' means that the parameter must be given;
%         '(empty)', or non-existent c{2}, means that the parameter
%         will not be assigned a default value.
%
% * c{3}: a string, representing the type of the value allowed:
%         'int' means integer;
%         'real' means real number;
%         'scalar' means any scalar;
%         'char' means character array or string;
%         if c{3} does not exist, then any type is allowed.
%
% * c{4}: minimum allowed value, if any, and if the type is 'int' or 'real'.
%
% * c{5}: maximum allowed value, likewise.
%
% Examples:
%
% make_param_struct({'a', 1, 'b', [2, 3], 'c', 'hello'}, {'...'})
% ans = scalar structure containing the fields: a = 1, b = [2, 3], c = 'hello'
%
% make_param_struct({'a', 1, 'b', [2, 3], 'c', 'hello'}, {})
% error: Unknown field name 'a' given
%
% make_param_struct({'a', 1, 'b', [2, 3]}, {'a', 'b', {'d', 4}})
% ans = scalar structure containing the fields: a = 1, b = [2, 3], d = 4
%
% make_param_struct({'a', 1, 'b', [2, 3]}, {'a', 'b', {'d', '(mandatory)'}})
% error: Must give mandatory parameter 'd'
%
function S = make_param_struct(parameters, constraints)
  % Check the types of our arguments.
  parameters_error_msg = ...
  'parameters must be {''name1'', value1, ''name2'', value2, ...}';
  constraints_error_msg = ...
  'constraints must be {{''name1'', ...}, {''name2'', ...}, ...}';
  if (~ iscell(parameters))
    error(parameters_error_msg);
  elseif mod(numel(parameters), 2)
    error(parameters_error_msg);
  elseif (~ iscell(constraints))
    error(constraints_error_msg);
  end

  % Store possible types as fields of a struct.  The associated value
  % is 1 if the given type may have a range.
  possible_types = struct('int', 1, 'real', 1, 'scalar', 0, 'char', 0);

  % Store allowed field names and constraints into a struct
  % + a special variable to denote that other field names are allowed.
  other_field_names = 0;
  fields = struct();
  for c = constraints
    c = c{1}; % Non-intuitive: iteration over cell array yields cell array.

    % Check this constraint's type, convert isolated strings
    % to single-element arrays.
    if ischar(c)
      c = {c};
    elseif (~ iscell(c)) || (numel(c) < 1) || (~ ischar(c{1}))
      error(constraints_error_msg);
    end

    % Store the constraint into fields, unless it is the special value
    % indicating that other field names are allowed.
    if strcmp(c{1}, '...')
      other_field_names = 1;
    else
      fields.(c{1}) = c;
    end
  end

  % Make the result struct, store the parameters, checking the values.
  S = struct();
  for ii = 1:2:numel(parameters)
    field_name = parameters{ii};
    value = parameters{ii + 1};

    if isfield(fields, field_name)
      c = fields.(field_name);
      l = numel(c);

      % Check the type given in c{3} if any.
      if (l >= 3)
	t = c{3};
	type_error_msg = ...
	sprintf('Value ''%s'' for field ''%s'' must be type ''%s''', ...
		value, field_name, t);
	if (~ isfield(possible_types, t))
	  error('Unknown field type ''%s''', disp(t(1:end-1)));
	elseif ((strcmp(t, 'char') && (~ ischar(value)))...
		|| ((~ strcmp(t, 'char')) && (~ isscalar(value)))...
		|| (strcmp(t, 'real') && (~ isreal(value)))...
		|| (strcmp(t, 'int') && ((~ isreal(value))...
					 || (mod(value, 1) ~= 0)))...
	       )
	  error(type_error_msg);
	end
      end

      % Check the range given in c{4} and c{5} if any, and if type allows it.
      % Testing for l >= 4 first ensures that t has already been checked.
      if (l >= 4 && possible_types.(t)...
	  && (value < c{4} || (l >= 5 && value > c{5})))
	error('Given value ''%f'' out of range', value);
      end
    elseif (~ other_field_names)
      error('Unknown field name ''%s'' given', disp(field_name(1:end-1)));
    end

    % We have satisfied the constraints, now store the value.
    S.(field_name) = value;
  end

  % Now that we have stored all the values given, check for default values
  % and mandatory parameters.
  all_fields = fieldnames(fields);
  for ii = 1:numel(all_fields)
    f = all_fields{ii};
    if (~ isfield(S, f)) && (numel(fields.(f)) >= 2)
      % The field named f has not been given but has a default action.
      default_value = fields.(f){2};
      if (ischar(default_value) && strcmp(default_value, '(empty)'))
	; % Do nothing, empty default value.
      elseif (ischar(default_value) && strcmp(default_value, '(mandatory)'))
	error('Must give mandatory parameter ''%s''', f);
      else
	% Store default value.
	S.(f) = default_value;
    end
  end
end

