function varargout=sde_euler_validate(dt,n,a,b,options)
%SDE_EULER_VALIDATE  Test SDE_EULER for performance and convergence order.
%   SDE_EULER_VALIDATE(DT,N)
%	YM = SDE_EULER_VALIDATE(DT,N)
%	YM = SDE_EULER_VALIDATE(DT,N,A,B)
%   YM = SDE_EULER_VALIDATE(DT,N,OPTIONS)
%	YM = SDE_EULER_VALIDATE(DT,N,A,B,OPTIONS)
%	[YM,YV] = SDE_EULER_VALIDATE(DT,N,...)
%   
%   See also: SDE_EULER, SDE_MILSTEIN_VALIDATE

%   For details of this validation method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9@case.edu, Created 11-1-10
%   Revision: 1.0, 4-11-12


close all
% Check inputs and outputs
if nargin < 2
    error('SDETools:sde_euler_validate:NotEnoughInputs',...
          'Not enough input arguments.');
else
    if nargin < 5
        if nargin == 3
            error('SDETools:sde_euler_validate:InvalidInputPattern',...
                 ['Not enough input arguments: both A and B must be '...
                  'specified.']);
        end
        options = [];
    elseif nargin == 5
        if isempty(options) && (ndims(options) ~= 2 || ...
                any(size(options) ~= 0) || ~(isstruct(options) || ...
                iscell(options) || isnumeric(options))) || ...
                ~isempty(options) && ~isstruct(options)
            error('SDETools:sde_euler_validate:InvalidSDESETStruct',...
                  'Invalid SDE options structure.  See SDESET.');
        end
    else
        error('SDETools:sde_euler_validate:TooManyInputs',...
              'Too many input arguments.');
    end
end
if nargout > 2
    error('SDETools:sde_euler_validate:TooManyOutputs',...
          'Too many output arguments.');
end

% Check DT
if ~isvector(dt) || ~isfloat(dt) || ~isreal(dt) || ~all(isfinite(dt))
    error(  'SDETools:sde_euler_validate:InvalidDT',...
            'DT must be a finite real floating-point vector.');
end
ldt = length(dt);
if length(dt) < 2
    error(  'SDETools:sde_euler_validate:BadInputSizeDT',...
            'Input vector DT must have length >= 2.');
end
dt = sort(dt);

% Check N
if ~isscalar(n) || ~isfloat(n) || ~isreal(n) || ~isfinite(n)
    error(  'SDETools:sde_euler_validate:InvalidN',...
            'N must be a finite real floating-point scalar.');
end
if isempty(n) || n < 1 || n-floor(n) ~= 0
    error(  'SDETools:sde_euler_validate:BadInputSizeN',...
            'Input N must be an integer >= 1.');
end

% Check A and B
if nargin >= 4
    if ~isscalar(a) || isempty(a) || ~isfloat(a) || ~isreal(a) || ~isfinite(a)
        error(  'SDETools:sde_euler_validate:InvalidA',...
                'A must be a finite real floating-point scalar.');
    end
    if ~isscalar(b) || isempty(b) || ~isfloat(b) || ~isreal(b) || ~isfinite(b)
        error(  'SDETools:sde_euler_validate:InvalidB',...
                'B must be a finite real floating-point scalar.');
    end
else
    a = 1;
    b = 1;
end

% Integration method is dependent on if SDE is Stratonovich or Ito form
Stratonovich = strcmp(sdeget(options,'SDEType','Stratonovich','flag'),...
        'Stratonovich');

% Create random number stream
RandSeed = sdeget(options,'RandSeed',[],'flag');
if ~isempty(RandSeed)
    if ~isscalar(RandSeed) || ~isnumeric(RandSeed) || ~isreal(RandSeed) || ...
            ~isfinite(RandSeed) || RandSeed >= 2^32 || RandSeed < 0
        error('SDETools:sde_euler_validate:InvalidRandSeed',...
              'RandSeed must be a non-negative integer value less than 2^32.');
    end
    % Create new stream based on seed value
	Stream = RandStream.create('mt19937ar','Seed',RandSeed);
else
    % Use default stream
    try
        Stream = RandStream.getGlobalStream;
    catch                                       %#ok<CTCH>
        Stream = RandStream.getDefaultStream;	%#ok<GETRS>
    end
end

% Override RandFun setting and specifiy stream via RANDN
options = sdeset(options,'RandFun',@(M,N)randn(Stream,M,N));

% Override non-diagonal noise, ConstFFUN, and ConstGFUN settings
options = sdeset(options,'DiagonalNoise','yes','ConstFFUN','no',...
                 'ConstGFUN','no');

t0 = 0;
tf = 20*dt(end);
y0 = ones(n,1);

f = @(t,y)a*y;
g = @(t,y)b*y;
Ym = zeros(ldt,1);
Yv = zeros(ldt,1);
if Stratonovich
    c = a;
    SDEType = 'Stratonovich';
else
    c = a-0.5*b^2;
    SDEType = 'Ito';
end

% Warm up for timing
[Y W] = sde_euler(f,g,[t0 t0+dt(1)],y0,options);    %#ok<NASGU,ASGLU>

% Loop through time-steps
ttotal = 0;
nsteps = 0;
for i=1:length(dt)
    t = t0:dt(i):tf;
    nsteps = nsteps+length(t);
    
    tic
    [Y W] = sde_euler(f,g,t,y0,options);
    ttotal = ttotal+toc;
    
    % Calculate error between analytic and simulated solutions
    Yerr = abs(y0(1)*exp(c*(t(end)-t0)+b*W(end,:))-Y(end,:));
    Ym(i) = mean(Yerr);
    Yv(i) = std(Yerr);
end

% Variable output
if nargout == 0
    disp(['Total simulation time: ' num2str(ttotal) ' seconds'])
    disp(['Mean of N simulations/time-step: ' num2str(ttotal/(nsteps)) ...
          ' seconds'])
else
    varargout{1} = Ym;
    if nargout == 2
        varargout{2} = Yv;
    end
end

% Plot results
figure
orders = [0.5 1.0 1.5 2.0]';
z = ones(length(orders),1);
xx = z*dt([1 end]);
logdt = log10(dt(end)/dt(1));
yy = Ym(1)*[z 10.^(orders*logdt)];
loglog(dt,Ym,'b.-',dt,Ym+Yv,'c',xx',yy','k')
text(xx(:,2)*10^(0.02*logdt),yy(:,2),cellstr(num2str(orders,'%1.1f')))
axis([dt(1) dt(end) Ym(1) yy(end,2)])
axis square
grid on
title(['SDE_EULER - ' SDEType ' - Convergence Order - ' num2str(n) ...
       ' simulations/time-step, A = ' num2str(a) ', B = ' num2str(b)],...
       'Interpreter','none')
xlabel('dt')
ylabel('Average Absolute Error')