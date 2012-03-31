%SDELAB
%   Version 1.0, 3-31-12
%
%   Stochastic differential equation (SDE) solvers.
%       sde_euler       - Euler-Maruyama (Ito) and Euler-Heun (Stratonovich) 
%       sde_milstein    - Milstein method, derivative and derivative-free.
%
%   Stochastic interpolatolation and solver utilities.
%       sde_interp      - Brownian bridge interpolation.
%       sde_interpq     - Quick Brownian bridge interpolation.
%       sde_interpqn    - Quick linearly-spaced Brownian bridge interpolation.
%       sdeget          - Get SDE OPTIONS structure parameters.
%       sdeset          - Create/alter SDE OPTIONS structure.
%
%   Numerical validation and demos.
%       sde_validate    - Test SDE solver performance and convergence order.

%   Tested with Matlab 7.13.0.564 (R2011b)
%   Mac OS X 10.6.8 Build: 10K549, Java 1.6.0_29-b11-402-10M3527

%   Andrew D. Horchler, adh9@case.edu, Created 12-18-11
%   Revision: 1.0, 3-31-12