%SDETools
%   Version 1.2, 9-21-13
%
%   Stochastic differential equation (SDE) solvers.
%       sde_euler       - Euler-Maruyama (Ito) and Euler-Heun (Stratonovich).
%       sde_milstein    - Milstein method, derivative and derivative-free.
%
%   SDE solver utilities.
%       sdeget          - Get SDE OPTIONS structure parameters.
%       sdephaseplot2   - 2-D phase plane SDE output function.
%       sdephaseplot3   - 3-D phase space SDE output function.
%       sdeplot         - Time series SDE output function.
%       sdeprint        - Command window printing SDE output function.
%       sdeset          - Create/alter SDE OPTIONS structure.
%
%   Stochastic differential equation (SDE) processes.
%       sde_bm          - Brownian motion process, analytic solution.
%       sde_gbm         - Geometric Brownian motion process, analytic solution.
%       sde_ou          - Ornstein-Uhlenbeck process, analytic solution.
%
%   Stochastic interpolatolation utilities.
%       sde_interp      - Brownian bridge interpolation.
%
%   Numerical validation.
%       sde_euler_validate  - Test sde_euler performance and convergence order.
%
%   Demos.
%       

%   This version tested with Matlab 9.0.0.341360 (R2016a)  
%   Mac OS X 10.11.4 (Build: 15E65), Java 1.7.0_75-b13  
%   Compatibility maintained back through Matlab 7.4 (R2007a)  

%   Andrew D. Horchler, horchler @ gmail . com
%   Created: 12-18-11, Revision: 1.2, 9-21-13


%   Acknowledgment of support: This material is based upon work supported by the
%   National Science Foundation under Grant No. 1065489. Disclaimer: Any
%   opinions, findings, and conclusions or recommendations expressed in this
%   material are those of the author(s) and do not necessarily reflect the views
%   of the National Science Foundation.


%  Copyright © 2011-2017, Andrew D. Horchler
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