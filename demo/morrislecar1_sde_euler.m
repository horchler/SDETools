function morrislecar1_sde_euler
%MORRISLECAR1_SDE_EULER  Noisy Morris-Lecar demo using SDE_EULER in SDETools.
%   Run a Langevin simulation of type I Morris-Lecar neurons starting from a
%   range of different initial voltages.
%
%   See also: MORRISLECAR2_SDE_EULER, SDE_EULER, SDESET

%   John Rinzel and G. Bard Ermentrout, "Analysis of Neural Excitability and
%   Oscillations," in: "Methods in Neural Modeling: From Ions to Networks,"
%   2nd Ed., Eds. C. Koch, I. Segev, MIT Press, Cambridge, MA, 1998, Ch. 7,
%   251-292.  http://www.math.pitt.edu/~bard/classes/chapt7.ps

%   Andrew D. Horchler, horchler @ gmail . com, Created 9-11-12
%   Revision: 1.2, 8-30-13


% Morris-Lecar Type I parameters (Rinzel & Ermentrout 1998)
V1 = -1.2; V2 = 18; V3 = 12; V4 = 17.4;
gCa = 4; gK = 8; gL = 2;
VK = -84; VL = -60; VCa = 120;
C = 20; phi = 1/15;

% Voltage-dependent funtions
lam = @(V)phi*cosh(0.5*(V-V3)/V4);
tv = @(V)tanh((V-V3)/V4);
winf = @(V)0.5*(1+tv(V));

t = 0:1e-2:30;	% Time vector
Ia = 55;        % Applied current
V0 = -35:10:35; % Initial Voltages
w0 = winf(V0);  % Initial fraction of open K+ channels
N = numel(V0);  % Number of neurons
NK = 1e3;       % Number of K+ channels

% Morris-Lecar Langevin drift and diffusion equations
V_drift = @(V,w)(Ia-0.5*(1+tanh((V-V1)/V2)).*(V-VCa)*gCa-w.*(V-VK)*gK-(V-VL)*gL)/C;
w_drift = @(V,w)lam(V).*(winf(V)-w);
w_diffusion = @(V,w)sqrt(lam(V).*(winf(V)-tv(V).*w)/NK);

% Combine V and w drift and diffusion equations, V = y(1:N), w = y(N+1:end)
f = @(t,y)[V_drift(y(1:N),y(N+1:end));w_drift(y(1:N),y(N+1:end))];
g = @(t,y)[zeros(N,1);w_diffusion(y(1:N),y(N+1:end))];

% Keep w >= 0, set random seed, specify Ito
opts = sdeset('NonNegative',N+1:2*N,...
              'RandSeed',1,...
              'SDEType','Ito');

% Use Euler-Maruyama to integrate SDEs
y1 = sde_euler(f,g,t,[V0(:);w0(:)],opts);

figure
plot(y1(:,1:N),y1(:,N+1:end))
title('Morris-Lecar Type I');
xlabel('V (mV)'); ylabel('w (Fraction open K+ channels)');

y2 = sde_milstein(f,g,t,[V0(:);w0(:)],opts);

figure
plot(t,abs(y1-y2))