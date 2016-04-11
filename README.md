SDETools
========
#####A Matlab Toolbox for the Numerical Solution of Stochastic Differential Equations (SDEs).#####
######Version 1.2, 9-21-13######
#####Download Repository: [ZIP Archive](https://github.com/horchler/SDETools/archive/master.zip)#####

How to install (and uninstall) SDETools:  
 1. Download and expand the *[SDETools-master.zip](https://github.com/horchler/SDETools/archive/master.zip)* ZIP archive of the repository.  
 2. Move the resultant *SDETools-master* folder to the desired permanent location.  
 3. In Matlab, navigate to *SDETools-master/SDETools/* and run ```sde_install```. This adds the necessary files and folders to the Matlab search path. To uninstall SDETools, run ```sde_install('remove')```.  
 4. Minor edits and bug reports and fixes can be submitted by [filing an issue](https://github.com/horchler/SDETools/issues) or via email. To add new functionality or make propose major changes, please [fork the repository](https://help.github.com/articles/fork-a-repo). Any new features should be accompanied by some means of testing. Email or file an issue if you have any questions.  
&nbsp;  

--------

Stochastic differential equation (SDE) solvers.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_euler```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_euler.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Euler-Maruyama (Ito) and Euler-Heun (Stratonovich).  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_milstein```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_milstein.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Milstein method, derivative and derivative-free.

SDE solver utilities.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdeget```](https://github.com/horchler/SDETools/blob/master/SDETools/sdeget.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Get SDE OPTIONS structure parameters.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdephaseplot2```](https://github.com/horchler/SDETools/blob/master/SDETools/sdephaseplot2.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 2-D phase plane SDE output function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdephaseplot3```](https://github.com/horchler/SDETools/blob/master/SDETools/sdephaseplot3.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 3-D phase space SDE output function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdeplot```](https://github.com/horchler/SDETools/blob/master/SDETools/sdeplot.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Time series SDE output function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdeprint```](https://github.com/horchler/SDETools/blob/master/SDETools/sdeprint.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Command window printing SDE output function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sdeset```](https://github.com/horchler/SDETools/blob/master/SDETools/sdeset.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create/alter SDE OPTIONS structure.

Stochastic differential equation (SDE) processes.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_bm```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_bm.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Brownian motion process, analytic solution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_gbm```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_gbm.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Geometric Brownian motion process, analytic solution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_ou```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_ou.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Ornstein-Uhlenbeck process, analytic solution.

Stochastic interpolatolation utilities.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_interp```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_interp.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Brownian bridge interpolation.

Numerical validation.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```sde_euler_validate```](https://github.com/horchler/SDETools/blob/master/SDETools/sde_euler_validate.m)&nbsp;- Test sde_euler performance and convergence order.  
&nbsp;  

--------

Andrew D. Horchler, *adh9 @ case . edu*, [biorobots.case.edu](http://biorobots.case.edu/)  
Created: 12-18-11, Revision: 1.2, 9-21-13  

This version tested with Matlab 9.0.0.341360 (R2016a)  
Mac OS X 10.11.4 (Build: 15E65), Java 1.7.0_75-b13  
Compatibility maintained back through Matlab 7.4 (R2007a)  
&nbsp;  

--------

Acknowledgment of support: This material is based upon work supported by the [National Science Foundation](http://www.nsf.gov/) under [Grant No.&nbsp;1065489](http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=1065489). Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.  
&nbsp;  

Copyright &copy; 2011&ndash;2016, Andrew D. Horchler  
All rights reserved.  

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.