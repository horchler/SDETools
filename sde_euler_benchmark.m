function sde_euler_benchmark(tests,N,Tol)
%SDE_EULER_BENCHMARK  Suite of performance tests for SDE_EULER solver function.
%   SDE_EULER_BENCHMARK will perform all tests. See code (type 'edit
%   SDE_EULER_BENCHMARK' in the command window) for test details or to add
%   additional tests. The SDELab Toolbox must be on the Matlab path or in the
%   same directory as this function.
%
%   SDE_EULER_BENCHMARK(TESTS) where TESTS = [N1 N2 ... NM] is a vector of
%   specific test numbers (indices, 1 ... N) to run. Tests can be listed in any
%   order and will be run in the order specified.
%
%   SDE_EULER_BENCHMARK(TESTS,N)
%%   SDE_EULER_BENCHMARK(TESTS,N,Tol)
%
%   Not part of the the SDELab Toolbox; used only for development.
%   
%   See also: SDE_EULER, SDEARGUMENTS, SDE_EULER_UNITTEST, SDE_EULER_VALIDATE.

%   Andrew D. Horchler, adh9@case.edu, Created 1-2-11
%   Revision: 1.0, 1-2-12


% make sure toolbox on path, otherwise ensure we're in right location and add it
if strfind(path,'SDELab')
    if exist('sde_euler','file') ~= 2
        error(  'SDELab:sde_euler_benchmark:FunctionNotFound',...
               ['The SDELab Toolbox is appears to be on the Matlab path, '...
                'but the SDE_EULER solver function cannot be found.']);
    end
    PathAdded = false;
else
    if exist('SDELab','dir') ~= 7
        error(  'SDELab:sde_euler_benchmark:ToolboxNotFound',...
               ['The SDELab Toolbox is not be on the Matlab path and the '...
                'root directory of the of the toolbox, SDELab, is in the '...
                'same directory as this function.']);
    end
    addpath SDELab
    if exist('sde_euler','file') ~= 2
        rmpath SDELab
        error(  'SDELab:sde_euler_benchmark:FunctionNotFoundAfterAdd',...
               ['The SDELab Toolbox was added to the Matlab path, but the '...
                'SDE_EULER solver function cannot be found.']);
    end
    PathAdded = true;   % we'll reset path at end
end

% validate input argument if it exists
if nargin >= 1
    if isempty(tests) || ~isnumeric(tests) || ~all(isfinite(tests))
        error('SDELab:sde_euler_benchmark:InvalidArg1','Invalid argument 1.');
    end
    if any(tests < 1)
        error(  'SDELab:sde_euler_benchmark:NotAnIndex',...
                'Tests are numbered as indices, from 1 to N.');
    end
    tests = floor(tests);
    RunTests = true;
else
    RunTests = false;
end

if nargin >= 2
    if isempty(N) || ~isnumeric(N) || ~all(isfinite(N))
        error('SDELab:sde_euler_benchmark:InvalidArg2','Invalid argument 2.');
    end
    N = floor(N);
else
    N = 1e3;
end

if nargin == 3
    if isempty(Tol) || length(Tol) ~= 1 || ~isnumeric(Tol) || ...
            ~all(isfinite(Tol)) || Tol <= 0
        error('SDELab:sde_euler_benchmark:InvalidArg2','Invalid argument 2.');
    end
    Tol = floor(Tol);
else
    Tol = 1e3;
end

M = 16;
fi = cell(M,1);
gi = cell(M,1);
tspani = cell(M,1);
y0i = cell(M,1);
optsi = cell(M,1);
paramsi = cell(M,1);


% no FOR loop cases

% scalar ICs

i = 1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = [];
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = [];
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = 1;
gi{i} = 1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;


% vector ICs

i = i+1;
fi{i} = ones(1000,1);
gi{i} = ones(1000,1);
tspani{i} = 0:0.001:1;
y0i{i} = zeros(1,1000);
optsi{i} = [];
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = ones(1000,1);
gi{i} = ones(1000,1);
tspani{i} = 0:0.5:500;
y0i{i} = zeros(1,1000);
optsi{i} = [];
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = zeros(1,1000);
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = zeros(1,1000);
optsi{i} = sdeset('ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;



i = i+1;
fi{i} = ones(1000,1);
gi{i} = 1;
tspani{i} = 0:0.001:1;
y0i{i} = zeros(1,1000);
optsi{i} = sdeset('Diagonal','no');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = ones(1000,1);
gi{i} = 1;
tspani{i} = 0:0.5:500;
y0i{i} = zeros(1,1000);
optsi{i} = sdeset('Diagonal','no');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x(1)+1;
tspani{i} = 0:0.001:1;
y0i{i} = zeros(1,1000);
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x(1)+1;
tspani{i} = 0:0.5:500;
y0i{i} = zeros(1,1000);
optsi{i} = sdeset('Diagonal','no','ConstFFUN','yes','ConstGFUN','yes');
paramsi{i} = [];
oneout{i} = false;







i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = [];
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = [];
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.001:1;
y0i{i} = 0;
optsi{i} = sdeset('SDEType','Ito');
paramsi{i} = [];
oneout{i} = false;

i = i+1;
fi{i} = @(t,x)x+1;
gi{i} = @(t,x)x+1;
tspani{i} = 0:0.5:500;
y0i{i} = 0;
optsi{i} = sdeset('SDEType','Ito');
paramsi{i} = [];
oneout{i} = false;


m = length(fi);
if RunTests
    ts = tests(tests<=m);
else
    ts = 1:m;
end

sp = ' ';
t = zeros(1,N);
lts = length(ts);
mt = zeros(1,lts);
st = zeros(1,lts);
ost = ceil(log10((lts+1)));
for k=1:lts
    i = ts(k);
    
    % get arguments
    f = fi{i};
    g = gi{i};
    tspan = tspani{i};
    y0 = y0i{i};
    opts = optsi{i};
    params = paramsi{i};
    
    % ensure that all tests start out with same random variates
	RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));
    
    % warm up function before timing
    if oneout{i}
        if isempty(opts) && isempty(params)
            y = sde_euler(f,g,tspan,y0);    %#ok<*NASGU>
        elseif ~isempty(opts) && isempty(params)
            y = sde_euler(f,g,tspan,y0,opts);
        elseif isempty(opts) && ~isempty(params)
            y = sde_euler(f,g,tspan,y0,[],params);
        else
            y = sde_euler(f,g,tspan,y0,opts,params);
        end
    else
        if isempty(opts) && isempty(params)
            [y w] = sde_euler(f,g,tspan,y0); %#ok<*ASGLU>
        elseif ~isempty(opts) && isempty(params)
            [y w] = sde_euler(f,g,tspan,y0,opts);
        elseif isempty(opts) && ~isempty(params)
            [y w] = sde_euler(f,g,tspan,y0,[],params);
        else
            [y w] = sde_euler(f,g,tspan,y0,opts,params);
        end
    end
    
    j = 0;
    while true  % timing loop
        j = j+1;
        
        % time an iteration
        if oneout{i}
            if isempty(opts) && isempty(params)
                tic
                y = sde_euler(f,g,tspan,y0);
                t(j) = toc;
            elseif ~isempty(opts) && isempty(params)
                tic
                y = sde_euler(f,g,tspan,y0,opts);
                t(j) = toc;
            elseif isempty(opts) && ~isempty(params)
                tic
                y = sde_euler(f,g,tspan,y0,[],params);
                t(j) = toc;
            else
                tic
                y = sde_euler(f,g,tspan,y0,opts,params);
                t(j) = toc;
            end
        else
            if isempty(opts) && isempty(params)
                tic
                [y w] = sde_euler(f,g,tspan,y0);
                t(j) = toc;
            elseif ~isempty(opts) && isempty(params)
                tic
                [y w] = sde_euler(f,g,tspan,y0,opts);
                t(j) = toc;
            elseif isempty(opts) && ~isempty(params)
                tic
                [y w] = sde_euler(f,g,tspan,y0,[],params);
                t(j) = toc;
            else
                tic
                [y w] = sde_euler(f,g,tspan,y0,opts,params);
                t(j) = toc;
            end
        end
        
        % check if times are sufficient
        if N < 20 && N >= 10
            if (j >= 10 && sum(t(1:j)) > 0.1 && std(t(1:j)) < Tol) || j > N
                mt(k) = mean(t(1:j));
                st(k) = std(t(1:j));
                break
            end
        elseif N < 10
            if j > N
                mt(k) = mean(t(1:j));
                st(k) = std(t(1:j));
                break2
            end
        else
            if (j >= 10 && sum(t(1:j)) > 0.1 && std(t(1:j)) < Tol) || ...
                    (j >= 20 && sum(t(1:j)) > 2) || j > N || sum(t(1:j)) > 5
                mt(k) = mean(t(1:j));
                st(k) = std(t(1:j));
                break
            end
        end
    end
    fprintf(1,['Test %s%d took an average of %1.3e sec. (%1.3e '...
               'sec./time-step) over %d runs.\n'],...
              sp(ones(1,ost-ceil(log10((i+1))))),i,mt(k),mt(k)/length(tspan),j);
end
fprintf(1,'Total mean time for all %d tests: %4.4f +/- %4.4f sec.\n',...
        lts,sum(mt),mean(st));

% reset path to prior state if we added toolbox
if PathAdded
    rmpath SDELab
end