function opts = sdeget(options,name,default,noErrorCheck)
%SDEGET	 Get SDE OPTIONS parameters.
%   VAL = SDEGET(OPTIONS,'NAME') extracts the value of the named property from
%   integrator options structure OPTIONS, returning an empty matrix if the
%   property value is not specified in OPTIONS. It is sufficient to type only
%   the leading characters that uniquely identify the property. Case is ignored
%   for property names. The empty array, [], is a valid OPTIONS argument.
%   
%   VAL = SDEGET(OPTIONS,'NAME',DEFAULT) extracts the named property as above,
%   but returns VAL = DEFAULT if the named property is not specified in OPTIONS.
%
%   Example:
%       % SDEGET returns 'Ito' if the SDEType property is not specified in opts
%       opts = sdeset('RandSeed',1);
%       val = sdeget(opts,'SDEType','Ito')
%       opts = sdeset(opts,'SDEType','Stratonovich');
%       val = sdeget(opts,'SDEType','Ito')
%   
%   See also:
%       SDESET, SDEPLOT, SDE_EULER, SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU

%   SDEGET is based on an updating of version 1.37.4.5 of Matlab's ODEGET

%   Andrew D. Horchler, horchler @ gmail . com, 10-28-10
%   Revision: 1.3, 6-21-15


% Undocumented usage for fast access with no error checking
if nargin == 4 && strcmp(noErrorCheck,'flag')
	opts = getknownfield(options,name,default);
	return;
end

if nargin < 2
	error('SDETools:sdeget:NotEnoughInputs','Not enough input arguments.');
end
if ~isempty(options) && ~isa(options,'struct')
	error('SDETools:sdeget:Arg1NotSDESETStruct',...
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
            'DGFUN'
            'RandSeed'
            'Antithetic'
            'RandFUN'
            'DiagonalNoise'
            'ConstFFUN'
            'ConstGFUN'
            'ConstDGFUN'
            'Refine'
            'MaxStep'
            'NonNegative'
            'NonNegativeFUN'
            'EventsFUN'
            'OutputFUN'
            'OutputYSelect'
            'OutputWSelect'
            'Diagnostics'
            'Stats'
        };

j = find(strncmpi(name,Names,length(name)));
if isempty(j)           % if no matches
	error('SDETools:sdeget:InvalidPropertyName',...
         ['Unrecognized property name ''%s''.  See SDESET for '...
          'possibilities.'],name);
elseif length(j) > 1	% if more than one match
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = [Names{j(1)} cell2mat(strcat({', '},Names(j(2:end)))')];
        error('SDETools:sdeget:AmbiguousPropertyName',...
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


function v=getknownfield(s,f,d)
%GETKNOWNFIELD	Get field f from struct s, or else yield default d.

if isfield(s,f)	% s could be empty.
	v = s.(f);
    if isempty(v)
        v = d;
    end
else
	v = d;
end