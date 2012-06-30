SDETools
========
#####A Matlab Toolbox for the Numerical Solution of Stochastic Differential Equations (SDEs).#####
######Version 1.0, 6-30-12######
&nbsp;  

Stochastic differential equation solvers.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_euler```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_euler.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Euler-Maruyama (Ito) and Euler-Heun (Stratonovich)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_milstein```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_milstein.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Milstein method, derivative and derivative-free.

Stochastic differential equation processes.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_gbm```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_gbm.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Geometric Brownian motion process, analytic solution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_ou```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_ou.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Ornstein-Uhlenbeck process, analytic solution.

Stochastic interpolatolation and SDE solver utilities.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_interp```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_interp.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Brownian bridge interpolation.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_interpq```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_interpq.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Quick Brownian bridge interpolation.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_interpqn```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_interpqn.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Quick linearly-spaced Brownian bridge interpolation.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdeget```](https://github.com/horchler/SDETools/blob/master/SDETools/sdeget.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Get SDE OPTIONS structure parameters.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdeset```](https://github.com/horchler/SDETools/blob/master/SDETools/sdeset.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create/alter SDE OPTIONS structure.

Numerical validation and demos.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_euler_validate```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_euler_validate.m)&nbsp;- Test sde_euler performance and convergence order.  
&nbsp;  

--------

Andrew D. Horchler, *adh9 @ case . edu*, [biorobots.case.edu](http://biorobots.case.edu/)  
Created: 12-18-11, Revision: 1.0, 6-30-12  

This version tested with Matlab 7.14.0.739 (R2012a)  
Mac OS X 10.6.8 (Build: 10K549), Java 1.6.0_33-b03-424-10M3720  
Compatibility maintained back through Matlab 7.4 (R2007a)  
&nbsp;  

--------

Acknowledgment of support: This material is based upon work supported by the [National Science Foundation](http://www.nsf.gov/) under  
[Grant No.&nbsp;1065489](http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=1065489). Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.  
&nbsp;  

Copyright &copy; 2012, Andrew D. Horchler  
All rights reserved.  

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.