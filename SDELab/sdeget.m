function opts = sdeget(options,name,default,noErrorCheck)
%SDEGET	Get SDE OPTIONS parameters.
%   VAL = SDEGET(OPTIONS,'NAME') extracts the value of the named property from
%   integrator options structure OPTIONS, returning an empty matrix if the
%   property value is not specified in OPTIONS. It is sufficient to type only
%   the leading characters that uniquely identify the property. Case is ignored
%   for property names. The empty array, [], is a valid OPTIONS argument.
%   
%   VAL = SDEGET(OPTIONS,'NAME',DEFAULT) extracts the named property as above,
%   but returns VAL = DEFAULT if the named property is not specified in OPTIONS.
%   For example
%   
%       val = sdeget(opts1,'SDEType','Ito');
%   
%   returns val = 'Ito' if the SDEType property is not specified in opts1.
%   
%   See also SDESET, SDE_EULER, SDE_MILSTEIN.

%   SDEGET is based on an updating of version 1.37.4.5 of Matlab's ODEGET

%   Andrew D. Horchler, adh9@case.edu, 10-28-10
%   Revision: 1.0, 3-29-12


% Undocumented usage for fast access with no error checking
if nargin == 4 && strcmp(noErrorCheck,'flag')
	opts = getknownfield(options,name,default);
	return;
end

if nargin < 2
	error('SDELab:sdeget:NotEnoughInputs','Not enough input arguments.');
end
if ~isempty(options) && ~isa(options,'struct')
	error(  'SDELab:sdeget:Arg1NotSDESETStruct',...
            'First argument must be an options structure created with SDESET.');
end

if nargin < 3
	default = [];
end
if isempty(options)
	opts = default;
	return;
end

Names = {   'SDEType'
            'DFFUN'
            'DGFUN'
            'RandFUN'
            'RandSeed'
            'Antithetic'
            'AdditiveNoise'
            'DiagonalNoise'
            'ConstFFUN'
            'ConstGFUN'
            'ConstDGFUN'
            'NonNegative'
        };

j = find(strncmpi(name,Names,length(name)));
if isempty(j)           % if no matches
	error(  'SDELab:sdeget:InvalidPropertyName',...
           ['Unrecognized property name ''%s''.  See SDESET for '...
            'possibilities.'],name);
elseif length(j) > 1	% if more than one match
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = [Names{j(1)} cell2mat(strcat({', '},Names(j(2:end)))')];
        error(  'SDELab:sdeget:AmbiguousPropertyName',...
               ['Ambiguous property name abbreviation ''%s'' (' msg ').'],name);
    end
end
if any(strcmp(fieldnames(options),Names{j}))
	opts = options.(Names{j});
	if isempty(opts)
        opts = default;
	end
else
	opts = default;
end


function v = getknownfield(s,f,d)
%GETKNOWNFIELD	Get field f from struct s, or else yield default d.

if isfield(s,f)	% s could be empty.
	v = s.(f);
    if isempty(v)
        v = d;
    end
else
	v = d;
end