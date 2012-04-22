function sde_interpq_unittest
%SDE_INTERPQ_UNITTEST  
%

%   Andrew D. Horchler, adh9@case.edu, Created 3-3-12
%   Revision: 1.0, 4-21-12

% Scalar TI, vector X, M = 2
t=[1 2]';x=[3 4]';ti=-Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))	%#ok<*ISMAT>

t=[1 2]';x=[3 4]';ti=-1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2]';x=[3 4]';ti=0;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2]';x=[3 4]';ti=0.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2]';x=[3 4]';ti=1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 3)

t=[1 2]';x=[3 4]';ti=1.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 3.5)

t=[1 2]';x=[3 4]';ti=2;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 4)

t=[1 2]';x=[3 4]';ti=2.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2]';x=[3 4]';ti=Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2]';x=[3 4]';ti=NaN;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2]';x=[3 4]';ti=[];
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2]';x=[3 4]';ti=ones(0,1);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2]';x=[3 4]';ti=ones(1,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

% Scalar TI, matrix X, M = 2, N = 2
t=[1 2]';x=[3 4;5 6]';ti=-Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2]';x=[3 4;5 6]';ti=-1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2]';x=[3 4;5 6]';ti=0;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2]';x=[3 4;5 6]';ti=0.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2]';x=[3 4;5 6]';ti=1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [3 5]))

t=[1 2]';x=[3 4;5 6]';ti=1.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [3.5 5.5]))

t=[1 2]';x=[3 4;5 6]';ti=2;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [4 6]))

t=[1 2]';x=[3 4;5 6]';ti=2.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2]';x=[3 4;5 6]';ti=Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2]';x=[3 4;5 6]';ti=NaN;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2]';x=[3 4;5 6]';ti=[];
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2]';x=[3 4;5 6]';ti=ones(0,1);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2]';x=[3 4;5 6]';ti=ones(1,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

% Scalar TI, vector X, M = 3, linearly-spaced
t=[1 2 3]';x=[3 4 6]';ti=-Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 3]';x=[3 4 6]';ti=-1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 3]';x=[3 4 6]';ti=0;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 3]';x=[3 4 6]';ti=0.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 3]';x=[3 4 6]';ti=1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 3)

t=[1 2 3]';x=[3 4 6]';ti=1.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 3.5)

t=[1 2 3]';x=[3 4 6]';ti=2;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 4)

t=[1 2 3]';x=[3 4 6]';ti=2.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 5)

t=[1 2 3]';x=[3 4 6]';ti=3;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 6)

t=[1 2 3]';x=[3 4 6]';ti=3.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 3]';x=[3 4 6]';ti=Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 3]';x=[3 4 6]';ti=NaN;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 3]';x=[3 4 6]';ti=[];
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2 3]';x=[3 4 6]';ti=ones(0,1);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2 3]';x=[3 4 6]';ti=ones(1,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

% Scalar TI, vector X, M = 3, not linearly-spaced
t=[1 2 2.5]';x=[3 4 6]';ti=-Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 2.5]';x=[3 4 6]';ti=-1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 2.5]';x=[3 4 6]';ti=0;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 2.5]';x=[3 4 6]';ti=0.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 2.5]';x=[3 4 6]';ti=1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 3)

t=[1 2 2.5]';x=[3 4 6]';ti=1.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 3.5)

t=[1 2 2.5]';x=[3 4 6]';ti=2;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 4)

t=[1 2 2.5]';x=[3 4 6]';ti=2.25;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 5)

t=[1 2 2.5]';x=[3 4 6]';ti=2.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 6)

t=[1 2 2.5]';x=[3 4 6]';ti=3;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 2.5]';x=[3 4 6]';ti=Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 2.5]';x=[3 4 6]';ti=NaN;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && isnan(xi))

t=[1 2 2.5]';x=[3 4 6]';ti=[];
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2 2.5]';x=[3 4 6]';ti=ones(0,1);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2 2.5]';x=[3 4 6]';ti=ones(1,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

% Scalar TI, matrix X, M = 3, N = 2, linearly-spaced
t=[1 2 3]';x=[3 4 6;5 6 10]';ti=-Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=-1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=0;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=0.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [3 5]))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=1.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [3.5 5.5]))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=2;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [4 6]))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=2.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [5 8]))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=3;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [6 10]))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=3.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=NaN;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=[];
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=ones(0,1);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2 3]';x=[3 4 6;5 6 10]';ti=ones(1,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

% Scalar TI, matrix X, M = 3, N = 2, not linearly-spaced
t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=-Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=-1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=0;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=0.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=1;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [3 5]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=1.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [3.5 5.5]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=2;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [4 6]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=2.25;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [5 8]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=2.5;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(xi == [6 10]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=3;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=Inf;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=NaN;
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [1 2]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=[];
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=ones(0,1);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';ti=ones(1,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))


% Vector TI, vector X, M = 2, D = 2
t=[1 2]';x=[3 4]';ti=[-Inf -1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2]';x=[3 4]';ti=[-1 0]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2]';x=[3 4]';ti=[0 0.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2]';x=[3 4]';ti=[0.5 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2]';x=[3 4]';ti=[-Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2]';x=[3 4]';ti=[Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2]';x=[3 4]';ti=[NaN 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2]';x=[3 4]';ti=[1 1.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 3.5]'))

t=[1 2]';x=[3 4]';ti=[1 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 4]'))

t=[1 2]';x=[3 4]';ti=[1.5 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3.5 4]'))

t=[1 2]';x=[3 4]';ti=[1.5 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2]';x=[3 4]';ti=[1.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2]';x=[3 4]';ti=[1.5 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2]';x=[3 4]';ti=[2 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2]';x=[3 4]';ti=[2 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2]';x=[3 4]';ti=[2 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2]';x=[3 4]';ti=[2 2.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2]';x=[3 4]';ti=[2.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2]';x=[3 4]';ti=[Inf NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2]';x=[3 4]';ti=[NaN NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2]';x=[3 4]';ti=ones(0,2);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2]';x=[3 4]';ti=ones(2,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

% Vector TI, vector X, M = 3, D = 2, linearly-spaced
t=[1 2 3]';x=[3 4 6]';ti=[-Inf -1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6]';ti=[-1 0]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6]';ti=[0 0.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6]';ti=[0.5 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 3]';x=[3 4 6]';ti=[-Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 3]';x=[3 4 6]';ti=[Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 3]';x=[3 4 6]';ti=[NaN 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 3]';x=[3 4 6]';ti=[1 1.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 3.5]'))

t=[1 2 3]';x=[3 4 6]';ti=[1 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 4]'))

t=[1 2 3]';x=[3 4 6]';ti=[1.5 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3.5 4]'))

t=[1 2 3]';x=[3 4 6]';ti=[1.5 2.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3.5 5]'))

t=[1 2 3]';x=[3 4 6]';ti=[1.5 3]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3.5 6]'))

t=[1 2 3]';x=[3 4 6]';ti=[2 3]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [4 6]'))

t=[1 2 3]';x=[3 4 6]';ti=[1 3]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 6]'))

t=[1 2 3]';x=[3 4 6]';ti=[1.5 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2 3]';x=[3 4 6]';ti=[1.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2 3]';x=[3 4 6]';ti=[1.5 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2 3]';x=[3 4 6]';ti=[2 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 3]';x=[3 4 6]';ti=[2 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 3]';x=[3 4 6]';ti=[2 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 3]';x=[3 4 6]';ti=[2 3.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 3]';x=[3 4 6]';ti=[3.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6]';ti=[Inf NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6]';ti=[NaN NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 3]';x=[3 4 6]';ti=ones(0,2);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2 3]';x=[3 4 6]';ti=ones(2,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

% Vector TI, vector X, M = 3, D = 2, not linearly-spaced
t=[1 2 2.5]';x=[3 4 6]';ti=[-Inf -1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6]';ti=[-1 0]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6]';ti=[0 0.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6]';ti=[0.5 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 2.5]';x=[3 4 6]';ti=[-Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 2.5]';x=[3 4 6]';ti=[Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 2.5]';x=[3 4 6]';ti=[NaN 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && isnan(xi(1)) && xi(2) == 3)

t=[1 2 2.5]';x=[3 4 6]';ti=[1 1.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 3.5]'))

t=[1 2 2.5]';x=[3 4 6]';ti=[1 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 4]'))

t=[1 2 2.5]';x=[3 4 6]';ti=[1.5 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3.5 4]'))

t=[1 2 2.5]';x=[3 4 6]';ti=[1.5 2.25]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3.5 5]'))

t=[1 2 2.5]';x=[3 4 6]';ti=[1.5 2.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3.5 6]'))

t=[1 2 2.5]';x=[3 4 6]';ti=[2 2.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [4 6]'))

t=[1 2 2.5]';x=[3 4 6]';ti=[1 2.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(xi == [3 6]'))

t=[1 2 2.5]';x=[3 4 6]';ti=[1.5 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2 2.5]';x=[3 4 6]';ti=[1.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2 2.5]';x=[3 4 6]';ti=[1.5 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 3.5 && isnan(xi(2)))

t=[1 2 2.5]';x=[3 4 6]';ti=[2 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 2.5]';x=[3 4 6]';ti=[2 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 2.5]';x=[3 4 6]';ti=[2 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 2.5]';x=[3 4 6]';ti=[2 3]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && xi(1) == 4 && isnan(xi(2)))

t=[1 2 2.5]';x=[3 4 6]';ti=[3 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6]';ti=[Inf NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6]';ti=[NaN NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 1]) && all(isnan(xi)))

t=[1 2 2.5]';x=[3 4 6]';ti=ones(0,2);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

t=[1 2 2.5]';x=[3 4 6]';ti=ones(2,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 1]))

% Vector TI, matrix X, M = 2, N = 2, D = 2
t=[1 2]';x=[3 4;5 8]';ti=[-Inf -1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2]';x=[3 4;5 8]';ti=[-1 0]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2]';x=[3 4;5 8]';ti=[0 0.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2]';x=[3 4;5 8]';ti=[0.5 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2]';x=[3 4;5 8]';ti=[-Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2]';x=[3 4;5 8]';ti=[Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2]';x=[3 4;5 8]';ti=[NaN 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2]';x=[3 4;5 8]';ti=[1 1.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [3 3.5 5 6.5]'))

t=[1 2]';x=[3 4;5 8]';ti=[1 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [3 4 5 8]'))

t=[1 2]';x=[3 4;5 8]';ti=[1.5 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [3.5 4 6.5 8]'))

t=[1 2]';x=[3 4;5 8]';ti=[1.5 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [3.5 6.5]) && all(isnan(xi(2,:))))

t=[1 2]';x=[3 4;5 8]';ti=[1.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [3.5 6.5]) && all(isnan(xi(2,:))))

t=[1 2]';x=[3 4;5 8]';ti=[1.5 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [3.5 6.5]) && all(isnan(xi(2,:))))

t=[1 2]';x=[3 4;5 8]';ti=[2 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2]';x=[3 4;5 8]';ti=[2 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2]';x=[3 4;5 8]';ti=[2 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2]';x=[3 4;5 8]';ti=[2 2.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2]';x=[3 4;5 8]';ti=[2.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2]';x=[3 4;5 8]';ti=[Inf NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2]';x=[3 4;5 8]';ti=[NaN NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2]';x=[3 4;5 8]';ti=ones(0,2);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2]';x=[3 4;5 8]';ti=ones(2,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

% Vector TI, matrix X, M = 3, N = 2, D = 2, linearly-spaced
t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[-Inf -1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[-1 0]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[0 0.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[0.5 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[-Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[Inf 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[NaN 1]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(1,:))) && all(xi(2,:) == [3 5]))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[1 1.5]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [3 3.5 5 6.5]'))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[1 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [3 4 5 8]'))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[1.5 2]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [3.5 4 6.5 8]'))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[2 3]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [4 6 8 11]'))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[2.5 3]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [5 6 9.5 11]'))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[1 3]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == [3 6 5 11]'))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[1.5 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [3.5 6.5]) && all(isnan(xi(2,:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[1.5 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [3.5 6.5]) && all(isnan(xi(2,:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[1.5 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [3.5 6.5]) && all(isnan(xi(2,:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[2 -Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[2 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[2 NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[2 4]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(1,:) == [4 8]) && all(isnan(xi(2,:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[4 Inf]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[Inf NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=[NaN NaN]';
xi=sde_interpq(t,x,ti);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(isnan(xi(:))))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=ones(0,2);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

t=[1 2 3]';x=[3 4 6;5 8 11]';ti=ones(2,0);
xi=sde_interpq(t,x,ti);
assert(all(size(xi) == [0 2]))

disp('All tests passed.')