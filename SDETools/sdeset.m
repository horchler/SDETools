function options = sdeset(varargin)
%SDESET  Create/alter SDE OPTIONS structure.
%   OPTIONS = SDESET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the specified
%   values. Any unspecified properties have default values. It is sufficient to
%   type only the leading characters that uniquely identify the property. Case
%   is ignored for property names. 
%   
%   OPTIONS = SDESET(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%   
%   OPTIONS = SDESET(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties overwrite
%   corresponding old properties. 
%   
%   SDESET with no input arguments displays all property names and their
%   possible values.
%   
%SDESET PROPERTIES
%   
%SDEType - Type of SDE  [ {Stratonovich} | Ito ]
%   This property specifies the form of the SDE to be solved, and therefore if
%   particular SDE integration scheme is based upon the Stratonovich or Ito
%   stochastic integral. The choice of Ito or Stratonovich depends on the form
%   of the modelled noise process. Ito is appropriate for approximating discrete
%   noise processes, while Stratonovich is suited for approximating continuous
%   noise processes. The Ito and Stratonovich forms are equivalent for additive
%   noise, i.e., the stochastic function, g(t,y), is not a function of the state
%   variables, g(t,y) = g(t).
%
%DGFUN - Derivative of the stochastic function  [ function_handle | matrix ]
%   Set this property to a function handle to specify the derivative of the
%   stochastic function, dg(t,y)/dy, to use derivative methods if available
%   (SDE_MILSTEN only). If this property is not set, derivative-free methods are
%   used. If a floating-point matrix is specified instead of a function handle,
%   the ConstGFUN property is implicitly assumed. If all values are equal to
%   zero, then the ConstDGFUN property is also assumed.
%   
%RandSeed - Create random stream and set seed  [ 0 <= integer < 2^32 ]
%   Create a random number stream separate from Matlab's default stream with the
%   specified seed value. If a seed is not specified, Matlab's current default
%   stream is used.
%
%Antithetic - Use antithetic variates for Wiener increments  [ yes | {no} ]
%   Set this property to 'yes' to use the antithetic variates variance reduction
%   method to calculate the normal variates for the Wiener increments. Half as
%   many normal variates are generated and the remainder are the negative of
%   these. If Matlab's default random number stream is used (i.e., RandSeed not
%   specified) the stream's antithetic property is reset to it's previous state
%   after use.
%   
%RandFUN - Random variates  [ RandStream object | function_handle | matrix ]
%   Specify a random number stream by setting this property to a RandStream
%   object created via Matlab's RandStream function. To specify an alternative
%   random number generator function, instead of using Matlab's random number
%   generators, set this property to a function handle. The corresponding
%   function must take two arguments, M and N, and output an M-by-N matrix of
%   normal variates. If a LENGTH(TSPAN)-by-D floating-point matrix, W, of
%   integrated Wiener increments is specified, these are used directly. The
%   RandSeed and Antithetic properties are ignored if the RandFUN property is
%   specified.
%
%DiagonalNoise - Dimension of stochastic function  [ {yes} | no ]
%   In the default diagonal (uncorrelated) noise situation each dimension of the
%   stochastic function is only perturbed by a single component of the
%   D-dimensional Wiener process. GFUN (and DGFUN, if specified) must return a
%   scalar or a column vector of length N = LENGTH(Y0). If a scalar is returned,
%   this value is used across all N dimensions. Set this property to 'no' to use
%   the more general, but more computationaly intense, correlated noise case
%   where each dimension of the stochastic function can be disturbed by multiple
%   components of the Wiener process. In this case, the stochastic function
%   argument, GFUN (and DGFUN, if specified), must return a N-by-D matrix.
%
%ConstFFUN - Deterministic function is constant  [ yes | {no} ]
%   Set this property to 'yes' if the determinstic function, f(t,y), of the SDE
%   is constant and therefore not a function of time or state. The function is
%   then evaluated only once by the integration routine, improving performance.
%   An empty vector, [], is equivalent to 'no'.
%
%ConstGFUN - Stochastic function is constant  [ yes | {no} ]
%   Set this property to 'yes' in the case of time-independent additive noise,
%   i.e., if the stochastic function, g(t,y), of the SDE is constant and
%   therefore not a function of time or state. The function is then evaluated
%   only once by the integration routine, improving performance. An empty
%   vector, [], is equivalent to 'no'.
%
%ConstDGFUN - Derivative of stochastic function is constant  [ yes | {no} ]
%   Set this property to 'yes' if the derivative of th stochastic function,
%   dg(t,y)/dy, of the SDE is constant and therefore not a function of time or
%   state (SDE_MILSTEIN only). The function is then evaluated only once by the
%   integration routine, improving performance.
%
%MaxStep - Upper bound on step size  [ scalar >= 0 ]
%   If any step size in TSPAN is greater than MaxStep, the interval will be
%   equally subdivided until the step size is less than or equal to MaxStep.
%   These sub-steps and their associated Wiener increments are not output.
%
%NonNegative - Non-negative components  [ yes | {no} | vector of indices ]
%   Set to 'yes' to specify that all components of the solution, including the
%   initial conditions, are to be kept non-negative by taking the absolute value
%   after each step. A vector of indices specifies individual components of the
%   solution vector to be kept non-negative. An empty vector, [], is equivalent
%   to 'no'. This option does not apply to the stochastic process functions,
%   which are based on analytic solutions.
%
%EventsFUN - Locate multiple zero-crossing events  [ function_handle ]
%   Set this property to a function handle in order to specify an events
%   function. The corresponding function must take at least two inputs and
%   output three vectors: [Value,IsTerminal,Direction] = Events(T,Y). The
%   scalar input T is the current integration time and the vector Y is the
%   current state. For the i-th event, Value(i) is the value of the
%   zero-crossing function and IsTerminal(i) = 1 specifies that integration is
%   to terminate at at a zero or to continue if IsTerminal(i) = 0. If
%   Direction(i) = 1, only zeros where Value(i) is increasing are found, if
%   Direction(i) = -1, only zeros where Value(i) is decreasing are found,
%   otherwise if Direction(i) = 0, all zeros are found. If Direction is set to
%   the empty matrix, [], all zeros are found for all events. Direction and
%   IsTerminal may also be scalars.
%
%OutputFUN - Function called by solver after each time-step  [ function_handle ]
%   Set this property to a function handle in order to specify a function that
%   will be called upon completion of each time-step. By default, the output
%   function must take three arguments: the current integration time, T, a
%   vector of solution components, Y, and a status string. The OutputYSelect
%   option can be used to pass only subsets (or none) of the solution components
%   to the output function. The OutputWSelect option can be enabled to pass
%   integrated Wiener increments, W, as an optional fourth argument to the
%   output function. See SDEPLOT for further output function requirements.
%
%OutputYSelect - Y Output selection indices  [ {yes} | no | vector of integers ]
%   A vector of indices specifies a subset of solution vector components of the
%   solution vector to be passed to OutputFUN. If 'no' is specified, and the
%   default OutputWSelect is used, only the current integration time is passed
%   to OutputFUN. An empty vector, [], is equivalent to 'no'.
%
%OutputWSelect - W Output selection indices  [ yes | {no} | vector of integers ]
%   Set to 'yes' to specify that all integrated Wiener increment components, W,
%   are passed to OutputFUN. A vector of indices specifies a subset of
%   integrated Wiener increment components to be passed to OutputFUN. An empty
%   vector, [], is equivalent to 'no'.
%   
%   See also:
%       SDEGET, SDEPLOT, SDE_EULER, SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU,
%       FUNCTION_HANDLE, RANDSTREAM, ODEPLOT

%   SDESET is based on an updating of version 1.46.4.10 of Matlab's ODESET.

%   Andrew D. Horchler, horchler @ gmail . com, 10-27-10
%   Revision: 1.3, 6-8-16


options = struct(	'SDEType',          [],...
                    'DFFUN',            [],...
                    'DGFUN',            [],...
                    'RandSeed',         [],...
                    'Antithetic',       [],...
                    'RandFUN',          [],...
                    'DiagonalNoise',	[],...
                    'ConstFFUN',        [],...
                    'ConstGFUN',        [],...
                    'ConstDGFUN',       [],...
                    'MaxStep',          [],...
                    'NonNegative',      [],...
                    'EventsFUN',        [],...
                    'OutputFUN',        [],...
                    'OutputYSelect',    [],...
                    'OutputWSelect',    []...
                );

Values = {	'{Stratonovich} |  Ito '
            'function_handle | vector'
            'function_handle | matrix'
            '0 <= integer < 2^32'
            ' yes  | {no} '
            'RandStream object | function_handle | matrix'
            '{yes} |  no '
            ' yes  | {no}'
            ' yes  | {no}'
            ' yes  | {no}'
            'scalar >= 0'
            ' yes  | {no} | vector '
            'function_handle'
            'function_handle'
            '{yes} |  no  | vector '
            ' yes  | {no} | vector '
         };

Names = fieldnames(options);
m = length(Names);

% Print out possible values of properties in the form of a struct
if nargin == 0 && nargout == 0
    len = cellfun(@length,Names);
    blanks = max(len)-len+4;
    sp = ' ';
    for i = m:-1:1
        out{i} = [sp(ones(1,blanks(i))) Names{i} ': [ ' Values{i} ' ]\n'];
    end
    fprintf(1,[out{:} '\n']);
    clear options;
    return;
end

% Combine all leading options structures opt1, opt2,... in sdeset(opt1,opt2,...)
i = 1;
while i <= nargin
	arg = varargin{i};
	if ischar(arg)      % arg is an option name
        break;
	end
	if ~isempty(arg)	% [] is a valid options argument
        if ~isa(arg,'struct')
            error('SDETools:sdeset:NoPropertyNameOrStruct',...
                 ['Expected argument %d to be a string property name or an '...
                  'options structure created with SDESET.'],i);
        end
        for j = 1:m
            Name = Names{j};
            if any(strcmp(fieldnames(arg),Name))
                options.(Name) = arg.(Name);
            end
        end
	end
	i = i+1;
end

% A finite state machine to parse name-value pairs
if rem(nargin-i+1,2) ~= 0
	error('SDETools:sdeset:ArgNameValueMismatch',...
          'Arguments must occur in name-value pairs.');
end
while i <= nargin
	arg = varargin{i};
    if ~ischar(arg)
        error('SDETools:sdeset:NoPropertyName',...
              'Expected argument %d to be a string property name.',i);
    end
    j = find(strncmpi(arg,Names,length(arg)));
    if isempty(j)           % If no matches
        error('SDETools:sdeset:InvalidPropertyName',...
             ['Unrecognized property name ''%s''.  See SDESET for '...
              'possibilities.'],name);
    elseif length(j) > 1	% If more than one match
        k = find(strcmpi(arg,Names));
        if length(k) == 1
            j = k;
        else
            msg = [Names{j(1)} cell2mat(strcat({', '},Names(j(2:end)))')];
            error('SDETools:sdeset:AmbiguousPropertyName',...
                 ['Ambiguous property name abbreviation ''%s'' (' msg ').'],...
                 arg);
        end
    end
    options.(Names{j}) = varargin{i+1};
    i = i+2;
end