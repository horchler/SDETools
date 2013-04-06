%SDETools
%   Version 1.1, 4-5-13
%
%   Stochastic differential equation (SDE) solvers.
%       sde_euler       - Euler-Maruyama (Ito) and Euler-Heun (Stratonovich) 
%       sde_milstein    - Milstein method, derivative and derivative-free.
%
%   Stochastic differential equation (SDE) processes.
%       sde_bm          - Brownian motion process, analytic solution.
%       sde_gbm         - Geometric Brownian motion process, analytic solution.
%       sde_ou          - Ornstein-Uhlenbeck process, analytic solution.
%
%   Stochastic interpolatolation and SDE solver utilities.
%       sde_interp      - Brownian bridge interpolation.
%       sde_interpq     - Quick Brownian bridge interpolation.
%       sde_interpqn    - Quick linearly-spaced Brownian bridge interpolation.
%       sdeget          - Get SDE OPTIONS structure parameters.
%       sdeset          - Create/alter SDE OPTIONS structure.
%
%   Numerical validation and demos.
%       sde_euler_validate  - Test sde_euler performance and convergence order.

%   Tested with Matlab 8.0.0.783 (R2012b)
%   Mac OS X 10.8.2 (Build: 12D78), Java 1.6.0_43-b01-447-11M4203
%   Compatibility maintained back through Matlab 7.4 (R2007a)

%   Andrew D. Horchler, adh9 @ case . edu
%   Created: 12-18-11, Revision: 1.1, 4-5-13


%   Acknowledgment of support: This material is based upon work supported by the
%   National Science Foundation under Grant No. 1065489. Disclaimer: Any
%   opinions, findings, and conclusions or recommendations expressed in this
%   material are those of the author(s) and do not necessarily reflect the views
%   of the National Science Foundation.


%  Copyright © 2011-2013, Andrew D. Horchler
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%   * Redistributions of source code must retain the above copyright notice,
%     this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright notice,
%     this list of conditions and the following disclaimer in the documentation
%     and/or other materials provided with the distribution.
%   * Neither the name of Case Western Reserve University nor the names of its
%     contributors may be used to endorse or promote products derived from this
%     software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY
%  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.