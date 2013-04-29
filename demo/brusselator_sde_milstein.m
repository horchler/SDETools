function brusselator_sde_milstein
%BRUSSELATOR_SDE_MILSTEIN  Noisy 2-D Brusselator demo using SDE_MILSTEIN.
%   Run a 2-D spatial simulation of noisy Brusselator equations.
%
%   See also: SDE_MILSTEIN, SDE_SET

%   Inspired by: http://en.wikipedia.org/wiki/File:Brusselator_space.gif

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-28-13
%   Revision: 1.0, 4-29-13


t = 0:0.05:800;     % Time vector
n = 64;             % Width and height of 
N = n^2;            % Number of Brusselators
y0 = 2*rand(2*N,1); % Randomize initial conditions

a = 1; b = 3;       % Reaction constants
sig = 5e-2;         % Noise magnitude
k = [0 0.25 0;
     0.25 -1 0.25;
     0 0.25 0];     % Laplace convolution kernel
dx = [0.2;0.02];    % X and Y diffusion coefficients

% Mask for borders
m = ones(n); m([1 n],1:n) = 0; m(2:n-1,[1 n]) = 0; m = [m(:);m(:)];

% Set derivative of diffusion, initialize output function, Ito, set random seed
opts = sdeset('DGFUN',dg(sig,m,N),...
              'OutputFUN',@(t,y)out(t,y,n,N),...
              'SDEType','Ito',...
              'RandSeed',1);

% Use Milstein method to integrate SDEs
sde_milstein(@(t,y)f(t,y,a,b,m,dx,k,n,N),@(t,y)g(t,y,sig,m,N),t,y0,opts);



function y=f(t,y,a,b,m,dx,k,n,N)   %#ok<INUSL>
% 2-D Brusselator drift function using convolution
X = y(1:N);
Y = y(N+1:end);
cx = conv2(reshape(X,[n n]),k,'same');
cy = conv2(reshape(Y,[n n]),k,'same');
y = m.*[a-(b+1-Y.*X).*X+dx(1)*cx(:);(b-Y.*X).*X+dx(2)*cy(:)];


function y=g(t,y,sig,m,N)  %#ok<INUSL>
% Brusselator diffusion function
y = sig*m.*[-y(1:N);y(1:N)];


function y=dg(sig,m,N)
% Derivative of Brusselator diffusion, constant
y = sig*m.*[-ones(N,1);ones(N,1)];


function out(t,y,n,N)
% Simple output function for 2-D Brusselator
persistent h z i
if isempty(h)
    i = 0;
    figure;
    z = zeros(n);
    z(:) = y(1:N);
    imagesc(z);
    axis square;
    title(['Brusselator, Time = ' sprintf('%.1f',t(1))])
    h = get(gca,'Children');
end
if mod(i,4) == 0
    z(:) = y(1:N);
    set(h,'CData',z);
    title(['Brusselator, Time = ' sprintf('%.1f',t(1))])
    drawnow('expose');
end
i = i+1;