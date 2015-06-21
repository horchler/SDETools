function brusellator1_sde_euler
%BRUSELLATOR1_SDE_EULER  Demo of Brussellator model using sde_euler in SDETools.
%   Run a simulation of Brussellator model.
%
%   See also: SDE_EULER, SDE_SET

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-28-13
%   Revision: 1.0, 7-16-13


t = 0:1e-2:100;	% Time vector

n = 1;
N = n^2;
y0 = 2*rand(2*N,1);

a = 1;
b = 3;
sig = 5e-2;

% Initialize output function, specify Ito, set random seed
opts = sdeset('OutputFUN',@sdeplot,...
              'SDEType','Ito',...
              'RandSeed',2);

% Use Euler-Maruyama to integrate SDEs
sde_euler(@(t,y)f(t,y,a,b,N),@(t,y)g(t,y,sig,N),t,y0,opts);



function y=f(t,y,a,b,N)	%#ok<INUSL>
X = y(1:N);
Y = y(N+1:end);
y = [a-(b+1-Y.*X).*X;(b-Y.*X).*X];

function y=g(t,y,sig,N)	%#ok<INUSL>
y = sig*[-y(1:N);y(1:N)];