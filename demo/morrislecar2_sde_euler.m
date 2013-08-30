function morrislecar2_sde_euler
%MORRISLECAR2_SDE_EULER  Noisy Morris-Lecar demo using SDE_EULER in SDETools.
%   Run a Langevin simulation of type II Morris-Lecar neurons at five different
%   numbers of K+ channels (noise levels).
%
%   See also: MORRISLECAR1_SDE_EULER, SDE_EULER, SDESET, SDEPLOT, ODEPLOT

%   John Rinzel and G. Bard Ermentrout, "Analysis of Neural Excitability and
%   Oscillations," in: "Methods in Neural Modeling: From Ions to Networks,"
%   2nd Ed., Eds. C. Koch, I. Segev, MIT Press, Cambridge, MA, 1998, Ch. 7,
%   251-292.  http://www.math.pitt.edu/~bard/classes/chapt7.ps

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-26-13
%   Revision: 1.2, 8-30-13


% Morris-Lecar Type I parameters (Rinzel & Ermentrout 1998)
V1 = -1.2; V2 = 18; V3 = 2; V4 = 30;
gCa = 4.4; gK = 8; gL = 2;
VK = -84; VL = -60; VCa = 120;
C = 20; phi = 0.04;

% Voltage-dependent funtions
lam = @(V)phi*cosh(0.5*(V-V3)/V4);
tv = @(V)tanh((V-V3)/V4);
winf = @(V)0.5*(1+tv(V));

t = 0:1e-2:100;         % Time vector
N = 5;                  % Number of neurons, different numbers of K+ channels
Ia = 93;                % Applied current
NK = logspace(2,4,N).';	% Vary number of K+ channels
V0 = -30+zeros(N,1);    % Initial Voltages
w0 = winf(V0);          % Initial fraction of open K+ channels

% Morris-Lecar Langevin drift and diffusion equations
V_drift = @(V,w)(Ia-0.5*(1+tanh((V-V1)/V2)).*(V-VCa)*gCa-w.*(V-VK)*gK-(V-VL)*gL)/C;
w_drift = @(V,w)lam(V).*(winf(V)-w);
w_diffusion = @(V,w)sqrt(lam(V).*(winf(V)-tv(V).*w)./NK);

% Combine V and w drift and diffusion equations, V = y(1:N), w = y(N+1:end)
f = @(t,y)[V_drift(y(1:N),y(N+1:end));w_drift(y(1:N),y(N+1:end))];
g = @(t,y)[zeros(N,1);w_diffusion(y(1:N),y(N+1:end))];

% Keep w >= 0, set random seed, Ito, set output function for V
opts = sdeset('NonNegative',N+1:2*N,...
              'RandSeed',1,...
              'SDEType','Ito',...
              'OutputFUN',@sdeplot,...
              'OutputYSelect',N+1:2*N);

% Use Euler-Maruyama to integrate SDEs
y = sde_euler(f,g,t,[V0(:);w0(:)],opts);

figure
plot(y(:,1:N),y(:,N+1:end))
title('Morris-Lecar Type II');
xlabel('V (mV)'); ylabel('w (Fraction open K+ channels)');