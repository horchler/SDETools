function sde_euler_benchmark(tests,N,Tol)
%SDE_EULER_BENCHMARK  Performance tests for SDE_EULER solver function.
%   SDE_EULER_BENCHMARK will perform all tests. See code (type 'edit
%   SDE_EULER_BENCHMARK' in the command window) for test details or to add
%   additional tests. The SDETools Toolbox must be on the Matlab path or in the
%   same directory as this function.
%
%   SDE_EULER_BENCHMARK(TESTS) where TESTS = [N1 N2 ... NM] is a vector of
%   specific test numbers (indices, 1 ... N) to run. Tests can be listed in any
%   order and will be run in the order specified.
%
%   SDE_EULER_BENCHMARK(TESTS,N)
%   SDE_EULER_BENCHMARK(TESTS,N,TOL)
%
%   Not part of the the SDETools Toolbox; used only for development.
%   
%   See also: SDE_EULER, SDEARGUMENTS, SDE_EULER_UNITTEST, SDE_EULER_VALIDATE,
%       SDE_MILSTEIN_BENCHMARK

%   Andrew D. Horchler, adh9 @ case . edu, Created 1-2-11
%   Revision: 1.2, 5-3-13


% Make sure toolbox on path, otherwise ensure we're in right location and add it
if ~isempty(strfind(path,'SDETools'))
    if exist('sde_euler','file') ~= 2
        error('SDETools:sde_euler_benchmark:FunctionNotFound',...
             ['The SDETools Toolbox is appears to be on the Matlab path, '...
              'but the SDE_EULER solver function cannot be found.']);
    end
else
    if exist('../SDETools','dir') ~= 7
        error('SDETools:sde_euler_benchmark:ToolboxNotFound',...
             ['The SDETools Toolbox is not be on the Matlab path and the '...
              'root directory of the of the toolbox, SDETools, is not in '...
              'the same directory as the ''test'' directory.']);
    end
    addpath('../SDETools');
    Cleanup = onCleanup(@()rmpath('../SDETools'));  % Make sure path is reset
    if exist('sde_euler','file') ~= 2
        error('SDETools:sde_euler_benchmark:FunctionNotFoundAfterAdd',...
             ['The SDETools Toolbox was added to the Matlab path, but the '...
              'SDE_EULER solver function cannot be found.']);
    end
end

% Validate input argument if it exists
if nargin >= 1
    if isempty(tests) || ~isnumeric(tests) || ~all(isfinite(tests))
        error('SDETools:sde_euler_benchmark:InvalidArg1','Invalid argument 1.');
    end
    if any(tests < 1)
        error('SDETools:sde_euler_benchmark:NotAnIndex',...
              'Tests are numbered as indices, from 1 to N.');
    end
    tests = floor(tests);
    RunTests = true;
else
    RunTests = false;
end

if nargin >= 2
    if isempty(N) || ~isnumeric(N) || ~all(isfinite(N))
        error('SDETools:sde_euler_benchmark:InvalidArg2','Invalid argument 2.');
    end
    N = floor(N);
else
    N = 1e3;
end

if nargin == 3
    if isempty(Tol) || length(Tol) ~= 1 || ~isnumeric(Tol) || ...
            ~all(isfinite(Tol)) || Tol <= 0
        error('SDETools:sde_euler_benchmark:InvalidArg2','Invalid argument 2.');
    end
    Tol = floor(Tol);
else
    Tol = 1e-3;
end

M = 16;
fi = cell(M,1);
gi = cell(M,1);
tspani = cell(M,1);
y0i = cell(M,1);
optsi = cell(M,1);


% no FOR loop cases

% scalar ICs

i = 1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = [];
oneout{i} = false;
description{i} = 'Scalar, constant drift and diffusion, variable step-size';

i = i+1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = [];
oneout{i} = false;
description{i} = 'Scalar, constant drift and diffusion, fixed step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Scalar, constant drift and diffusion functions, variable step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Scalar, constant drift and diffusion functions, fixed step-size';

i = i+1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no');
oneout{i} = false;
description{i} = 'Scalar, constant drift, constant non-diagonal diffusion, variable step-size';

i = i+1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no');
oneout{i} = false;
description{i} = 'Scalar, constant drift, constant non-diagonal diffusion, fixed step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Scalar, constant drift function, constant non-diagonal diffusion function, variable step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Scalar, constant drift function, constant non-diagonal diffusion function, fixed step-size';


% vector ICs

i = i+1;
fi{i} = ones(1000,1);
gi{i} = ones(1000,1);
tspani{i} = 0:0.001:1;
y0i{i}(1,1000) = 0;
optsi{i} = [];
oneout{i} = false;
description{i} = 'Vector, constant drift and diffusion, variable step-size';

i = i+1;
fi{i} = ones(1000,1);
gi{i} = ones(1000,1);
tspani{i} = 0:0.5:500;
y0i{i}(1,1000) = 0;
optsi{i} = [];
oneout{i} = false;
description{i} = 'Vector, constant drift and diffusion, fixed step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i}(1,1000) = 0;
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Vector, constant drift and diffusion functions, variable step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i}(1,1000) = 0;
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Vector, constant drift and diffusion functions, fixed step-size';


i = i+1;
fi{i} = ones(1000,1);
gi{i} = 1;
tspani{i} = 0:0.001:1;
y0i{i}(1,1000) = 0;
optsi{i} = sdeset('Diagonal','no');
oneout{i} = false;
description{i} = 'Vector, constant drift, constant non-diagonal diffusion, variable step-size';

i = i+1;
fi{i} = ones(1000,1);
gi{i} = 1;
tspani{i} = 0:0.5:500;
y0i{i}(1,1000) = 0;
optsi{i} = sdeset('Diagonal','no');
oneout{i} = false;
description{i} = 'Vector, constant drift, constant non-diagonal diffusion, fixed step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x(1)+1;
tspani{i} = 0:0.001:1;
y0i{i}(1,1000) = 0;
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Vector, constant drift function, constant non-diagonal diffusion function, variable step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x(1)+1;
tspani{i} = 0:0.5:500;
y0i{i}(1,1000) = 0;
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
oneout{i} = false;
description{i} = 'Vector, constant drift function, constant non-diagonal diffusion function, fixed step-size';



i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = [];
oneout{i} = false;
description{i} = 'Scalar, diagonal diffusion function, variable step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = [];
oneout{i} = false;
description{i} = 'Scalar, diagonal diffusion function, fixed step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('SDEType','Ito');
oneout{i} = false;
description{i} = 'Scalar, Ito type, diagonal diffusion function, variable step-size';

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('SDEType','Ito');
oneout{i} = false;
description{i} = 'Scalar, Ito type, diagonal diffusion function, fixed step-size';


% Memory allocation
i = i+1;
fi{i} = ones(500000,1);
gi{i} = ones(500000,1);
tspani{i} = 0:0.5:5;
y0i{i}(1,500000) = 0;
optsi{i} = [];
oneout{i} = false;
description{i} = 'Memory allocation: vector, constant drift and diffusion, fixed step-size';


m = length(fi);
if RunTests
    ts = tests(tests<=m);
else
    ts = 1:m;
end

Stream = RandStream('mt19937ar','Seed',0);

sp = ' ';
t(1,N) = 0;
lts = length(ts);
mt(1,lts) = 0;
st(1,lts) = 0;
ost = ceil(log10((lts+1)));
for k=1:lts
    i = ts(k);
    
    % get arguments
    f = fi{i};
    g = gi{i};
    tspan = tspani{i};
    y0 = y0i{i};
    opts = optsi{i};
    
    % Ensure that all tests start out with same random variates
    try
        RandStream.setGlobalStream(Stream);
    catch                                   	%#ok<CTCH>
        RandStream.setDefaultStream(Stream);	%#ok<SETRS>
    end
    
    % Warm up function before timing
    if oneout{i}
        if isempty(opts)
            y = sde_euler(f,g,tspan,y0);        %#ok<*NASGU>
        else
            y = sde_euler(f,g,tspan,y0,opts);
        end
    else
        if isempty(opts)
            [y,w] = sde_euler(f,g,tspan,y0);	%#ok<*ASGLU>
        else
            [y,w] = sde_euler(f,g,tspan,y0,opts);
        end
    end
    
    j = 0;
    while true  % Timing loop
        j = j+1;
        
        % Time an iteration
        if oneout{i}
            if isempty(opts)
                tic
                y = sde_euler(f,g,tspan,y0);
                t(j) = toc;
            else
                tic
                y = sde_euler(f,g,tspan,y0,opts);
                t(j) = toc;
            end
        else
            if isempty(opts)
                tic
                [y,w] = sde_euler(f,g,tspan,y0);
                t(j) = toc;
            else
                tic
                [y,w] = sde_euler(f,g,tspan,y0,opts);
                t(j) = toc;
            end
        end
        
        % Check if times are sufficient
        if N < 20 && N >= 10
            if (j >= 10 && sum(t(1:j)) > 0.1 && std(t(1:j)) < Tol) || j > N
                break
            end
        elseif N < 10
            if j > N
                break
            end
        else
            if (j >= 10 && sum(t(1:j)) > 0.1 && std(t(1:j)) < Tol) || ...
                    (j >= 20 && sum(t(1:j)) > 2) || j > N || sum(t(1:j)) > 5
                break
            end
        end
    end
    mt(k) = mean(t(1:j));
	st(k) = std(t(1:j));
    fprintf(1,[description{k} '\n']);
    fprintf(1,['Test %s%d took an average of %1.3e sec. (%1.3e '...
               'sec./time-step) over %d runs.\n\n'],...
              sp(ones(1,ost-ceil(log10((i+1))))),i,mt(k),mt(k)/length(tspan),j);
end
fprintf(1,'Total mean time for all %d tests: %4.4f +/- %4.4f sec.\n',...
	lts,sum(mt),mean(st));