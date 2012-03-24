function [Ym Yv] = sde_validate(sdefun,dt,n,a,b,SDEType,RandSeed)
%SDE_VALIDATE  Test SDE solver algorithms for performance and convergence order.
%   SDE_VALIDATE(SDEFUN,DT,N)
%	YM = SDE_VALIDATE(SDEFUN,DT,N)
%	YM = SDE_VALIDATE(SDEFUN,DT,N,A,B)
%	YM = SDE_VALIDATE(SDEFUN,DT,N,A,B,SDETYPE)
%	YM = SDE_VALIDATE(SDEFUN,DT,N,A,B,SDETYPE,RANDSEED)
%	[YM YV] = SDE_VALIDATE(SDEFUN,DT,N)

%   For details of this validation method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9@case.edu, 11-1-10
%   Revision: 1.0

close all
% Check inputs
if nargin < 3
    error(  'MATLAB:sde_validate:NotEnoughInputs',...
            'Not enough input arguments.  See SDE_VALIDATE.');
end

if ~isa(sdefun,'function_handle')
    error(  'MATLAB:sde_validate:NotAFunctionHandle',...
            'Input SDEFUN must be a function handle.  See SDE_VALIDATE.');
end

ldt = length(dt);
if length(dt) < 2
    error(  'MATLAB:sde_validate:BadInputSize',...
            'Input vector DT must have length >= 2.  See SDE_VALIDATE.');
end
dt = sort(dt);

if length(n) ~= 1 || n < 1
    error(  'MATLAB:sde_validate:BadInputSize',...
            'Input N must be a scalar >= 1.  See SDE_VALIDATE.');
end
n = floor(n);

t0 = 0;
tf = 5*dt(end);
y0 = ones(n,1);

if length(a) ~= 1
    error(  'MATLAB:sde_validate:BadInputSize',...
            'Input A must be a scalar.  See SDE_VALIDATE.');
end
if length(b) ~= 1
    error(  'MATLAB:sde_validate:BadInputSize',...
            'Input B must be a scalar.  See SDE_VALIDATE.');
end
amb = (a-0.5*b^2);
apb = (a+0.5*b^2);

if exist('SDEType','var')
    isIto = strncmpi(SDEType,'Ito',length(SDEType));
    if ~strncmpi(SDEType,'Stratonovich',length(SDEType)) && ~isIto
        error(	'MATLAB:sde_validate:ImproperValue',...
                'SDETYPE must be string: ''Stratonovich'' or ''Ito''.  See SDE_VALIDATE.');
    end
    if isIto
        SDEType = 'Ito';
    else
        SDEType = 'Stratonovich';
    end
else
    SDEType = 'Stratonovich';
    isIto = false;
end
options = sdeset('SDEType',SDEType,'DiagonalNoise','yes');

if exist('RandSeed','var')
    if RandSeed >= 2^32 || RandSeed < 0
        error(	'MATLAB:sde_validate:ImproperValue',...
                'RANDSEED must be positive and less than 2^32.  See SDE_VALIDATE.');
    end
else
    RandSeed = 1;
end

%Ye_Ito = y0*exp((a-0.5*b^2)*(t-t0) + b*(W-W0));
%Ye_Strat = y0*exp(a*(t-t0) + b*(W-W0));

tic
Stream = RandStream.create('mt19937ar','seed',RandSeed);
Ym = zeros(1,ldt);
Yv = zeros(1,ldt);
for i=1:length(dt)
    t = t0:dt(i):tf;
    %options = sdeset(options,'RandFun',@(M,N)randn(Stream,M,N));
    %options = sdeset(options,'RandFun',@(M,N)randn(Stream,M,N),'DGFUN',@(t,y)b*eye(n));
    options = sdeset(options,'RandSeed',RandSeed+i-1);
    %
    if isIto
        [Y W] = feval(sdefun,@(t,y)a*y,@(t,y)b*y,t,y0,options);
        %[Y W] = feval(sdefun,@(t,y)a*y,@(t,y)diag(b*y,0),t,y0,options);
        %Ye = y0(1)*exp(amb*(t(end)-t0)*ones(n,1) + b*sum(dW,2));
    else
        [Y W] = feval(sdefun,@(t,y)amb*y,@(t,y)b*y,t,y0,options);
        %[Y W] = feval(sdefun,@(t,y)amb*y,@(t,y)diag(b*y,0),t,y0,options);
        %[Y dW] = feval(sdefun,@(t,y)a*y,@(t,y)b*y,t,y0,options);
        %Ye = y0(1)*exp(a*(t(end)-t0)*ones(n,1) + b*sum(dW,2));
    end
    Ye = y0(1)*exp(amb*(t(end)-t0)*ones(n,1) + b*W(:,end));
    %
    %{
    if isIto
        [Y dW] = feval(sdefun,@(t,y)apb*y,@(t,y)b*y,t,y0,options);
    else
        [Y dW] = feval(sdefun,@(t,y)a*y,@(t,y)b*y,t,y0,options);
    end
    Ye = y0(1)*exp(a*(t(end)-t0)*ones(n,1) + b*sum(dW,2));
    %}
    YY = abs(Ye-Y(:,end));
    Ym(i) = mean(YY);
    Yv(i) = std(YY);
end
toc

% Plot results
figure
orders = [0.5 1.0 1.5 2.0]';
xx = ones(4,1)*dt([1 end]);
yy = Ym(1)*[ones(4,1) 10.^(orders*log10(dt(end)/dt(1)))];
loglog(dt,Ym,'b.-',dt,Ym+Yv,'c',xx',yy','k')
text(xx(:,2)*10^(0.02*log10(dt(end)/dt(1))),yy(:,2),cellstr(num2str(orders,'%1.1f')))
axis([dt(1) dt(end) Ym(1) yy(4,2)])
axis square
grid on
xlabel('dt')
ylabel('Average Absolute Error')

%{
Y1_Ito=sde_euler(@(t,y)f1(t,y,a),@(t,y)g1(t,y,b),t,y0,opts_Ito);
Y2_Ito=sde_milstein2(@(t,y)f1(t,y,a),@(t,y)g1(t,y,b),t,y0,opts_Ito);

tic
Y1_Ito=sde_euler(@(t,y)f1(t,y,a),@(t,y)g1(t,y,b),t,y0,opts_Ito);
toc

tic
Y2_Ito=sde_milstein2(@(t,y)f1(t,y,a),@(t,y)g1(t,y,b),t,y0,opts_Ito);
toc


tic
Y2_Ito=sde_milstein(@(t,y)f1(t,y,a),@(t,y)g1(t,y,b),t,y0,opts_Ito);
toc

tic
Y1_Strat=sde_euler(@(t,y)f1(t,y,a-0.5*b^2),@(t,y)g1(t,y,b),t,y0,opts_Strat);
toc

tic
Y2_Strat=sde_milstein(@(t,y)f1(t,y,a-0.5*b^2),@(t,y)g1(t,y,b),t,y0,opts_Strat);
toc

Ye_Ito=y0*exp((a-0.5*b^2)*(t-t0)+b*(W-W0));

figure
semilogy(t,Y1_Ito,'r',t,Y2_Ito,'g',t,Ye_Ito,'b')

figure
semilogy(t,Y1_Strat,'r',t,Y2_Strat,'g',t,Ye_Ito,'b')

figure
semilogy(t,abs(Y1_Ito-Y2_Ito),'r',t,abs(Y1_Strat-Y2_Strat),'g')
hold on
semilogy(t,abs(Ye_Ito-Y1_Ito),'m',t,abs(Ye_Ito-Y2_Ito),'c')
semilogy(t,abs(Ye_Ito-Y1_Strat),'k',t,abs(Ye_Ito-Y2_Strat),'b')


tic
Y1_Ito=sde_euler(@(t,y)f1(t,y,a+0.5*b^2),@(t,y)g1(t,y,b),t,y0,opts_Ito);
toc

tic
Y2_Ito=sde_milstein(@(t,y)f1(t,y,a+0.5*b^2),@(t,y)g1(t,y,b),t,y0,opts_Ito);
toc

tic
Y1_Strat=sde_euler(@(t,y)f1(t,y,a),@(t,y)g1(t,y,b),t,y0,opts_Strat);
toc

tic
Y2_Strat=sde_milstein(@(t,y)f1(t,y,a),@(t,y)g1(t,y,b),t,y0,opts_Strat);
toc

Ye_Strat=y0*exp(a*(t-t0)+b*(W-W0));

figure
semilogy(t,Y1_Ito,'r',t,Y2_Ito,'g',t,Ye_Strat,'b')

figure
semilogy(t,Y1_Strat,'r',t,Y2_Strat,'g',t,Ye_Strat,'b')

figure
semilogy(t,abs(Y1_Ito-Y2_Ito),'r',t,abs(Y1_Strat-Y2_Strat),'g')
hold on
semilogy(t,abs(Ye_Strat-Y1_Ito),'m',t,abs(Ye_Strat-Y2_Ito),'c')
semilogy(t,abs(Ye_Strat-Y1_Strat),'k',t,abs(Ye_Strat-Y2_Strat),'b')
%}

function y=f1(t,x,a)
y = a*x;

function y=g1(t,x,b)
y = b*x;