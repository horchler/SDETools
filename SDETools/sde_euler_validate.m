function varargout=sde_euler_validate(dt,n,a,b,options)
%SDE_EULER_VALIDATE  Test SDE_EULER for performance and convergence order.
%   SDE_EULER_VALIDATE(DT,N)
%	YM = SDE_EULER_VALIDATE(DT,N)
%	YM = SDE_EULER_VALIDATE(DT,N,A,B)
%   YM = SDE_EULER_VALIDATE(DT,N,OPTIONS)
%	YM = SDE_EULER_VALIDATE(DT,N,A,B,OPTIONS)
%	[YM,YV] = SDE_EULER_VALIDATE(DT,N,...)
%
%   Example:
%       % Convergence order of Euler-Heun (Stratonovich) & Euler-Maruyama (Ito)
%       dt = logspace(-3,-1,3); n = 100; a = 1.5; b = 1;  % Try smaller b values
%       opts = sdeset('RandSeed',1);
%       sde_euler_validate(dt,n,a,b,opts);
%       opts = sdeset(opts,'SDEType','Ito');
%       sde_euler_validate(dt,n,a,b,opts);
%   
%   See also: SDE_EULER, SDE_GBM, SDE_MILSTEIN_VALIDATE, SDESET

%   For details of this validation method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, horchler @ gmail . com, Created 11-1-10
%   Revision: 1.2, 11-28-13


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
        if isempty(options) && (~sde_ismatrix(options) ...
                || any(size(options) ~= 0) || ~(isstruct(options) ...
                || iscell(options) || isnumeric(options))) ...
                || ~isempty(options) && ~isstruct(options)
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
    error('SDETools:sde_euler_validate:InvalidDT',...
          'DT must be a finite real floating-point vector.');
end
ldt = length(dt);
if length(dt) < 2
    error('SDETools:sde_euler_validate:BadInputSizeDT',...
          'Input vector DT must have length >= 2.');
end
dt = sort(dt);

% Check N
if ~isscalar(n) || ~isfloat(n) || ~isreal(n) || ~isfinite(n)
    error('SDETools:sde_euler_validate:InvalidN',...
          'N must be a finite real floating-point scalar.');
end
if isempty(n) || n < 1 || n ~= floor(n)
    error('SDETools:sde_euler_validate:BadInputSizeN',...
          'Input N must be an integer >= 1.');
end
% 20 batches of N trajectories
m = 20;
N = m*n;

% Check A and B
if nargin >= 4
    if ~isscalar(a) || isempty(a) || ~isfloat(a) || ~isreal(a) || ~isfinite(a)
        error('SDETools:sde_euler_validate:InvalidA',...
              'A must be a finite real floating-point scalar.');
    end
    if ~isscalar(b) || isempty(b) || ~isfloat(b) || ~isreal(b) || ~isfinite(b)
        error('SDETools:sde_euler_validate:InvalidB',...
              'B must be a finite real floating-point scalar.');
    end
else
    a = 1;
    b = 1;
end

% Check random number generation
if ~isempty(sdeget(options,'RandFUN',[],'flag'))
    error('SHCTools:sde_euler_validate:InvalidRandFUN',...
          'This function only supports the default random number stream.');
end
if strcmp(sdeget(options,'Antithetic','no','flag'),'yes')
    error('SHCTools:sde_euler_validate:Antithetic',...
          'This function does not support antithetic random variates.');
end

% Set random seed unless already specified
if isempty(sdeget(options,'RandSeed',[],'flag'))
    options = sdeset(options,'RandSeed',1);
end

% Override non-diagonal noise, ConstFFUN, and ConstGFUN settings
options = sdeset(options,'DiagonalNoise','yes','ConstFFUN','no',...
    'ConstGFUN','no');

% Get SDE type for plot
SDEType = sdeget(options,'SDEType','Stratonovich','flag');

t0 = 0;
tf = 20*max(dt);
y0 = ones(N,1);

f = @(t,y)a*y;
g = @(t,y)b*y;
Yerr(m,n) = 0;
Ym(ldt,1) = 0;
ehat(ldt,1) = 0;
Yv(ldt,1) = 0;
dehat(ldt,1) = 0;

% Warm up for timing
[Yeuler,W] = sde_euler(f,g,[t0 t0+dt(1)],y0,options);	%#ok<NASGU,ASGLU>

% Loop through time-steps
ttotal = 0;
nsteps = 0;
for i=1:ldt
    t = t0:dt(i):tf;
    nsteps = nsteps+length(t);
    
    tic
    [Yeuler,W] = sde_euler(f,g,t,y0,options);
    ttotal = ttotal+toc;
    
    %Ygbm = y0(1)*exp((a-b^2/2)*tf+b*W(end,:));
    Ygbm = sde_gbm(a,b,[t0 tf],y0,sdeset(options,'RandFUN',W([1 end],:)));
    
    % Calculate error between analytic and simulated solutions
    Yerr(:) = abs(Ygbm(end,:)-Yeuler(end,:));
    
    Ym(i) = mean(Yerr(:));
    
    ehatj = sum(Yerr,2)/n;
    ehat(i) = sum(ehatj)/m;
    Yv(i) = sum((ehatj-ehat(i)).^2)/(m-1);
    dehat(i) = tinv(0.95,m-1)*sqrt(Yv(i)/m);
end

% Variable output
if nargout == 0
    disp(['Total simulation time: ' num2str(ttotal) ' seconds']);
    disp(['Mean of ' int2str(N) ' simulations/time-step: ' ...
        num2str(ttotal/(nsteps)) ' seconds']);
else
    varargout{1} = Ym;
    if nargout == 2
        varargout{2} = Yv;
    end
end

% Plot results
hf=figure;
set(hf,'Color','w');
set(hf, 'Position', get(0,'Screensize'));
orders = [0.5 1.0 1.5 2.0]';
z = ones(length(orders),1);
xx = z*dt([1 end]);
logdt = log10(dt(end)/dt(1));
yy = Ym(1)*[z 10.^(orders*logdt)];
loglog(dt,Ym,'b.-',dt,ehat+dehat,'b--',dt,ehat-dehat,'b--',xx',yy','k')
text(xx(:,2)*10^(0.02*logdt),yy(:,2),cellstr(num2str(orders,'%1.1f')),...
    'FontSize',24)
axis([dt(1) dt(end) min([Ym(1) ehat(1)-dehat(1) yy(1)]) ...
      max([Ym(end) ehat(end)+dehat(end) yy(end,2)])])
set(gca,'FontSize',24)
axis square
grid on
title(['SDE\_EULER - ' SDEType ' - Convergence Order - ' int2str(N) ...
       ' simulations/time-step, \mu = ' num2str(a) ', \sigma = ' num2str(b)],...
       'Interpreter','Tex','FontSize',24)
xlabel('dt','FontSize',24)
ylabel('Average Absolute Error','FontSize',24)