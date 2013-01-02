function sde_interpqn_unittest
%SDE_INTERPQN_UNITTEST  
%

%   Andrew D. Horchler, adh9@case.edu, Created 3-4-12
%   Revision: 1.0, 4-21-12


% set variance to zero for testing
eta=0;

% Vector X, M = 2
t=[1 2]';x=[3 4]';nt=-Inf;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))	%#ok<*ISMAT>

t=[1 2]';x=[3 4]';nt=-1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2]';x=[3 4]';nt=0;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2]';x=[3 4]';nt=1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [1 1]) && xi == 3.5)

t=[1 2]';x=[3 4]';nt=3;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [5 1]) && all(xi == [3 3.25 3.5 3.75 4]'))

t=[1 2]';x=[3 4]';nt=[];
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2]';x=[3 4]';nt=ones(0,1);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2]';x=[3 4]';nt=ones(1,0);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

% Matrix X, M = 2, N = 2
t=[1 2]';x=[3 4;5 6]';nt=-Inf;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2]';x=[3 4;5 6]';nt=-1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2]';x=[3 4;5 6]';nt=0;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [2 2]) && all(xi(:) == x(:)))

t=[1 2]';x=[3 4;5 6]';nt=1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [3 2]) && all(xi(:) == [3 3.5 4 5 5.5 6]'))

t=[1 2]';x=[3 4;5 6]';nt=3;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [5 2]) && all(xi(:) == [3 3.25 3.5 3.75 4 5 5.25 5.5 5.75 6]'))

t=[1 2]';x=[3 4;5 6]';nt=[];
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2]';x=[3 4;5 6]';nt=ones(0,1);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2]';x=[3 4;5 6]';nt=ones(1,0);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

% Vector X, M = 3, linearly-spaced
t=[1 2 3]';x=[3 4 6]';nt=-Inf;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 3]';x=[3 4 6]';nt=-1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 3]';x=[3 4 6]';nt=0;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [3 1]) && all(x == xi))

t=[1 2 3]';x=[3 4 6]';nt=1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [5 1]) && all(xi == [3 3.5 4 5 6]'))

t=[1 2 3]';x=[3 4 6]';nt=3;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [9 1]) && all(xi == [3 3.25 3.5 3.75 4 4.5 5 5.5 6]'))

t=[1 2 3]';x=[3 4 6]';nt=[];
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 3]';x=[3 4 6]';nt=ones(0,1);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 3]';x=[3 4 6]';nt=ones(1,0);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

% Vector X, M = 3, not linearly-spaced
t=[1 2 2.5]';x=[3 4 6]';nt=-Inf;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 2.5]';x=[3 4 6]';nt=-1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 2.5]';x=[3 4 6]';nt=0;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [3 1]) && all(x == xi))

t=[1 2 2.5]';x=[3 4 6]';nt=1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [5 1]) && all(xi == [3 3.5 4 5 6]'))

t=[1 2 2.5]';x=[3 4 6]';nt=3;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [9 1]) && all(xi == [3 3.25 3.5 3.75 4 4.5 5 5.5 6]'))

t=[1 2 2.5]';x=[3 4 6]';nt=[];
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 2.5]';x=[3 4 6]';nt=ones(0,1);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

t=[1 2 2.5]';x=[3 4 6]';nt=ones(1,0);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 1]))

% Matrix X, M = 3, N = 2, linearly-spaced
t=[1 2 3]';x=[3 4 6;5 6 10]';nt=-Inf;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 3]';x=[3 4 6;5 6 10]';nt=-1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 3]';x=[3 4 6;5 6 10]';nt=0;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [3 2]) && all(xi(:) == x(:)))

t=[1 2 3]';x=[3 4 6;5 6 10]';nt=1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [5 2]) && all(xi(:) == [3 3.5 4 5 6 5 5.5 6 8 10]'))

t=[1 2 3]';x=[3 4 6;5 6 10]';nt=3;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [9 2]) && all(xi(:) == [3 3.25 3.5 3.75 4 4.5 5 5.5 6 5 5.25 5.5 5.75 6 7 8 9 10]'))

t=[1 2 3]';x=[3 4 6;5 6 10]';nt=[];
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 3]';x=[3 4 6;5 6 10]';nt=ones(0,1);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 3]';x=[3 4 6;5 6 10]';nt=ones(1,0);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

% Matrix X, M = 3, N = 2, not linearly-spaced
t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=-Inf;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=-1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=0;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [3 2]) && all(xi(:) == x(:)))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=1;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [5 2]) && all(xi(:) == [3 3.5 4 5 6 5 5.5 6 8 10]'))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=3;
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [9 2]) && all(xi(:) == [3 3.25 3.5 3.75 4 4.5 5 5.5 6 5 5.25 5.5 5.75 6 7 8 9 10]'))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=[];
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=ones(0,1);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

t=[1 2 2.5]';x=[3 4 6;5 6 10]';nt=ones(1,0);
xi=sde_interpqn(t,x,nt,eta);
assert(ndims(xi) == 2 && all(size(xi) == [0 2]))

disp('All tests passed.')