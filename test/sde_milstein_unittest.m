function sde_milstein_unittest(tests)
%SDE_MILSTEIN_UNITTEST  Suite of unit tests for SDE_MILSTEIN solver function.
%   SDE_MILSTEIN_UNITTEST will perform all tests. Detailed feedback is provided
%   if any of the tests fail. See code (type 'edit SDE_MILSTEIN_UNITTEST' in the
%   command window) for test details or to add additional tests. The SDETools
%   Toolbox must be on the Matlab path or in the same directory as this
%   function.
%
%   SDE_MILSTEIN_UNITTEST(TESTS) where TESTS = [N1 N2 ... NM] is a vector of
%   specific test numbers (indices, 1 ... N) to run. This capability is useful
%   for rerunning a particular test after a failure. If multiple test numbers
%   are specified, they may be listed in any order, but they will be evaluated
%   in ascending order.
%
%   Not part of the the SDETools Toolbox; used only for development.
%   
%   See also: SDE_MILSTEIN, SDEARGUMENTS, SDE_MILSTEIN_BENCHMARK,
%       SDE_MILSTEIN_VALIDATE, SDE_EULER_UNITTEST, WAITTEXT

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-2-12
%   Revision: 1.2, 5-3-13


% Make sure toolbox on path, otherwise ensure we're in right location and add it
if ~isempty(strfind(path,'SDETools'))
    if exist('sde_milstein','file') ~= 2
        error('SDETools:sde_milstein_unittest:FunctionNotFound',...
             ['The SDETools Toolbox is appears to be on the Matlab path, '...
              'but the SDE_MILSTEIN solver function cannot be found.']);
    end
else
    if exist('../SDETools','dir') ~= 7
        error('SDETools:sde_milstein_unittest:ToolboxNotFound',...
             ['The SDETools Toolbox is not be on the Matlab path and the '...
              'root directory of the of the toolbox, SDETools, is not in '...
              'the same directory as the ''test'' directory.']);
    end
    addpath('../SDETools');
    Cleanup = onCleanup(@()rmpath('../SDETools'));  % Make sure path is reset
    if exist('sde_milstein','file') ~= 2
        error('SDETools:sde_milstein_unittest:FunctionNotFoundAfterAdd',...
             ['The SDETools Toolbox was added to the Matlab path, but the '...
              'SDE_MILSTEIN solver function cannot be found.']);
    end
end

% Validate input argument if it exists
if nargin == 1
    if isempty(tests) || ~isnumeric(tests) || ~all(isfinite(tests))
        error('SDETools:sde_milstein_unittest:InvalidArgument',...
              'Invalid argument.');
    end
    if any(tests < 1)
        error('SDETools:sde_milstein_unittest:NotAnIndex',...
              'Tests are numbered as indices, from 1 to N.');
    end
    runtests = true;
    tests = floor(tests);
else
    runtests = false;
end

lnum1 = cell(70,1);
cmd = cell(70,1);
msg = cell(70,1);

lnum = cell(234,1);
f = cell(234,1);
g = cell(234,1);
tspan = cell(234,1);
y0 = cell(234,1);
opts = cell(234,1);
twoout = cell(234,1);


% Error tests:

% Number of arguments:

% only one function
st = dbstack;
i = 1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sde_milstein:NotEnoughInputs';

% only one function with options
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,0:0.1:1,0,sdeset(''SDEType'',''Ito''));';
msg{i} = 'SDETools:sde_milstein:NotEnoughInputsOptions';

% no tspan
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0);';
msg{i} = 'SDETools:sde_milstein:NotEnoughInputs';

% no tspan with options
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0,sdeset(''SDEType'',''Ito''));';
msg{i} = 'SDETools:sde_milstein:NotEnoughInputsOptions';

% no ICs
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1);';
msg{i} = 'SDETools:sde_milstein:NotEnoughInputs';


% f:

% not a valid input, not []
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein('''',@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:InvalidFFUN';

% not a valid input, not float
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(uint8(0),@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:InvalidFFUN';

% not a valid input, not a matrix
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(ones(2,2,2),@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:InvalidFFUN';

% function doesn't exist
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@ff,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'MATLAB:UndefinedFunction';

% state argument not defined
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t)x+1,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNTooFewInputs';

% state output not specified
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@f1,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNNoOutput';

% state output not assigned in function
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@f2,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNUnassignedOutput';

% too many inputs required
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x,a)a*x+1,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNTooManyInputs';

% state output smaller
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x(1)+1,@(t,x)x+1,0:0.1:1,[0 0]);';
msg{i} = 'SDETools:sdearguments:FFUNDimensionMismatch';

% state output larger
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)[x;x]+1,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNDimensionMismatch';

% state output not a matrix
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)rand(1,1,2)+1,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNNotColumnVector';

% state output not a colum vector
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)[x,x]+1,@(t,x)x+1,0:0.1:1,[0 0]);';
msg{i} = 'SDETools:sdearguments:FFUNNotColumnVector';

% state output not non-empty
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)[],@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNNotColumnVector';

% state output not float
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)uint8(x)+1,@(t,x)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:FFUNNotColumnVector';


% g:

% not a function handle
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,'''',0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:InvalidGFUN';

% not a valid input, not float
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,uint8(0),0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:InvalidGFUN';

% not a valid input, not a matrix
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,ones(2,2,2),0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:InvalidGFUN';

% function doesn't exist
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@gg,0:0.1:1,0);';
msg{i} = 'MATLAB:UndefinedFunction';

% state argument not defined
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t)x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNTooFewInputs';

% state output not specified
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@g1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNNoOutput';

% state output not assigned in function
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@g2,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNUnassignedOutput';

% too many inputs required
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x,a)a*x+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNTooManyInputs';

% state output smaller
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x(1:2)+1,0:0.1:1,[0,0,0]);';
msg{i} = 'SDETools:sdearguments:GFUNDimensionMismatchDiagonal';

% state output larger
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)[x;x]+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNDimensionMismatchDiagonal';

% state output not a matrix
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)rand(1,1,2)+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNNotMatrix';

% state output number of rows not equal number of states for non-diagonal noise
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)rand(2)+1,0:0.1:1,0,sdeset(''Diagonal'',''no''));';
msg{i} = 'SDETools:sdearguments:GFUNDimensionMismatchNonDiagonal';

% state output not a column vector for diagonal noise
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)rand(2)+1,0:0.1:1,[0 0]);';
msg{i} = 'SDETools:sdearguments:GFUNDimensionMismatchDiagonal';

% state output not non-empty
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)[],0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNNotMatrix';

% state output not float
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)uint8(x)+1,0:0.1:1,0);';
msg{i} = 'SDETools:sdearguments:GFUNNotMatrix';


% tspan:

% scalar time
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0,0);';
msg{i} = 'SDETools:sdearguments:InvalidTSpanSize';

% const time
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,[0,0,0,0],0);';
msg{i} = 'SDETools:sdearguments:TspanNotMonotonic';

% non-monotonic time
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,sin(0:0.5*pi:2*pi),0);';
msg{i} = 'SDETools:sdearguments:TspanNotMonotonic';

% not float
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,uint8(0:9),0);';
msg{i} = 'SDETools:sdearguments:InvalidTSpanDataType';

% not non-empty
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,rand(0,9),0);';
msg{i} = 'SDETools:sdearguments:InvalidTSpanSize';

% not real
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,(0:0.1:1)+1i,0);';
msg{i} = 'SDETools:sdearguments:InvalidTSpanDataType';


% y0:

% not float
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,uint8(0));';
msg{i} = 'SDETools:sdearguments:Y0EmptyOrNotFloat';

% not non-empty
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,[]);';
msg{i} = 'SDETools:sdearguments:Y0EmptyOrNotFloat';


% options:

% invalid non-empty options
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,0);';
msg{i} = 'SDETools:sde_milstein:InvalidSDESETStruct';

% invalid empty options
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,'''');';
msg{i} = 'SDETools:sde_milstein:InvalidSDESETStruct';

% invalid empty options
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,zeros(0,1));';
msg{i} = 'SDETools:sde_milstein:InvalidSDESETStruct';

% invalid empty options
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,zeros(0,0,0));';
msg{i} = 'SDETools:sde_milstein:InvalidSDESETStruct';

% inavlid RandSeed, not a matrix
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandSeed'',ones(1,1,2)));';
msg{i} = 'SDETools:sderandfun:InvalidRandSeed';

% inavlid RandSeed, not a scalar
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandSeed'',[0,0]));';
msg{i} = 'SDETools:sderandfun:InvalidRandSeed';

% inavlid RandSeed, not finite
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandSeed'',NaN));';
msg{i} = 'SDETools:sderandfun:InvalidRandSeed';

% inavlid RandSeed, not real
st = dbstack;
i = i+1;
lnum{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandSeed'',1i));';
msg{i} = 'SDETools:sderandfun:InvalidRandSeed';

% inavlid RandSeed, not numeric
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandSeed'',logical(0)));';
msg{i} = 'SDETools:sderandfun:InvalidRandSeed';

% RandSeed too small
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandSeed'',-1));';
msg{i} = 'SDETools:sderandfun:InvalidRandSeed';

% RandSeed too large
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandSeed'',2^32));';
msg{i} = 'SDETools:sderandfun:InvalidRandSeed';


% RandFUN:

% undefined RandFUN
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandFUN'',@rr));';
msg{i} = 'MATLAB:UndefinedFunction';

% empty RandFUN output, D > N, nargout == 2
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = '[y,w]=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)[]));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray1';

% non-matrix RandFUN output, D > N, nargout == 2
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = '[y,w]=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)randn(p,q,2)));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray1';

% non-float RandFUN output, D > N, nargout == 2
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = '[y,w]=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)uint8(randn(p,q))));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray1';

% RandFUN output size mismatch, D > N, nargout == 2
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = '[y,w]=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)randn(1,q)));';
msg{i} = 'SDETools:sde_milstein:RandFUNDimensionMismatch1';

% empty RandFUN output, D > N, nargout == 1
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)[]));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray2';

% non-matrix RandFUN output, D > N, nargout == 1
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)randn(p,q,2)));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray2';

% non-float RandFUN output, D > N, nargout == 1
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)uint8(randn(p,q))));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray2';

% RandFUN output size mismatch, D > N, nargout == 1
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)[x,x]+1,0:0.1:1,0,sdeset(''Diagonal'',''no'',''RandFUN'',@(p,q)randn(p,1)));';
msg{i} = 'SDETools:sde_milstein:RandFUNDimensionMismatch2';

% empty RandFUN output, D <= N
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,[0,0],sdeset(''RandFUN'',@(p,q)[]));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray3';

% non-matrix RandFUN output, D <= N
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,[0,0],sdeset(''RandFUN'',@(p,q)randn(p,q,2)));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray3';

% non-float RandFUN output, D <= N
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,[0,0],sdeset(''RandFUN'',@(p,q)uint8(randn(p,q))));';
msg{i} = 'SDETools:sde_milstein:RandFUNNot2DArray3';

% RandFUN output size mismatch, D <= N
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,[0,0],sdeset(''RandFUN'',@(p,q)randn(1,q)));';
msg{i} = 'SDETools:sde_milstein:RandFUNDimensionMismatch3';

% RandFUN only has one input
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandFUN'',@(p)randn(p)));';
msg{i} = 'SDETools:sde_milstein:RandFUNTooFewInputs';

% RandFUN output not specified
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandFUN'',@r1));';
msg{i} = 'SDETools:sde_milstein:RandFUNNoOutput';

% RandFUN output not assigned
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandFUN'',@r2));';
msg{i} = 'SDETools:sde_milstein:RandFUNUnassignedOutput';

% RandFUN requires more than two inputs
st = dbstack;
i = i+1;
lnum1{i} = st.line;
cmd{i} = 'y=sde_milstein(@(t,x)x+1,@(t,x)x+1,0:0.1:1,0,sdeset(''RandFUN'',@(m,n,p)p*randn(m,n)));';
msg{i} = 'SDETools:sde_milstein:RandFUNTooManyInputs';



% normal unit tests:

% cases where no FOR loop is required:

% with numeric inputs

st = dbstack;
i = 1;
lnum{i} = st.line;
f{i} = '1';             % constant, autonomous, scalar
g{i} = '1';             % constant, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '1';             % constant, autonomous, scalar
g{i} = '1';             % constant, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1]';       	% constant, autonomous, vector
g{i} = '[1;1]';     	% constant, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% scalar ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1]';       	% constant, autonomous, vector
g{i} = '[1;1]';     	% constant, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';     	% scalar ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1]';       	% constant, autonomous, vector
g{i} = '[1,1;1,1]';   	% constant, autonomous, matrix, diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% scalar ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1]';       	% constant, autonomous, vector
g{i} = '[1,1;1,1]';    	% constant, autonomous, matrix, diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';     	% scalar ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1]';       	% constant, autonomous, vector
g{i} = '[1,1,1;1,1,1]';	% constant, autonomous, matrix, diagonal, D == 3
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% scalar ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [0 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1]';       	% constant, autonomous, vector
g{i} = '[1,1,1;1,1,1]';	% constant, autonomous, matrix, diagonal, D == 3
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';     	% scalar ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [0 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1;1]';      	% constant, autonomous, vector
g{i} = '[1,1;1,1;1,1]';	% constant, autonomous, matrix, diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% scalar ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [0 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '[1;1;1]';      	% constant, autonomous, vector
g{i} = '[1,1;1,1;1,1]';	% constant, autonomous, matrix, diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0,0]';     	% scalar ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [0 1];


% with function inputs

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';      	% constant, autonomous, scalar
g{i} = '@(t,x)1';     	% constant, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';     	% constant, autonomous, scalar
g{i} = '@(t,x)1';      	% constant, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';      	% constant, autonomous, scalar
g{i} = '@(t,x)[1,1]';  	% constant, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [0 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';      	% constant, autonomous, scalar
g{i} = '@(t,x)[1,1]';  	% constant, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [0 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';     	% constant, autonomous, scalar
g{i} = '@(t,x)1';      	% constant, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';      	% constant, autonomous, scalar
g{i} = '@(t,x)1';    	% constant, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1]';       % constant, autonomous, vector
g{i} = '@(t,x)ones(2,3)';	% constant, autonomous, scalar, non-diagonal, D == 3
tspan{i} = '0:0.5:10';   	% constant step-size, increasing
y0{i} = '[0,0]';          	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [0 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1]';       % constant, autonomous, vector
g{i} = '@(t,x)ones(2,3)';	% constant, autonomous, scalar, non-diagonal, D == 3
tspan{i} = '0:0.1:10';   	% non-constant step-size, increasing
y0{i} = '[0,0]';          	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [0 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1]';  	% constant, autonomous, vector
g{i} = '@(t,x)1';     	% constant, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1]';  	% constant, autonomous, vector
g{i} = '@(t,x)1';     	% constant, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';     	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1;1]';  	% constant, autonomous, vector
g{i} = '@(t,x)ones(3,2)';	% constant, autonomous, scalar, non-diagonal, D == 2
tspan{i} = '0:0.5:10';   	% constant step-size, increasing
y0{i} = '[0,0,0]';       	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1;1]';  	% constant, autonomous, vector
g{i} = '@(t,x)ones(3,2)';	% constant, autonomous, scalar, non-diagonal, D == 2
tspan{i} = '0:0.1:10';   	% non-constant step-size, increasing
y0{i} = '[0,0,0]';       	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1;1]';	% constant, autonomous, vector
g{i} = '@(t,x)ones(3)';	% constant, autonomous, scalar, non-diagonal, D == 3
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];
 
st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1;1]';	% constant, autonomous, vector
g{i} = '@(t,x)ones(3)';	% constant, autonomous, scalar, non-diagonal, D == 3
tspan{i} = '0:0.1:10'; 	% non-constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];


% default options:

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = '';
twoout{i} = [1 1];


% with options:

% non-diagonal, Stratonovich:

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.5:10';      % constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.1:10';      % non-constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0,0]';    	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';    	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'')';
twoout{i} = [1 1];


% diagonal, Ito:

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'')';
twoout{i} = [1 1];


% non-diagonal, Ito:

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.5:10';      % constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.1:10';      % non-constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0,0]';    	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';    	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'')';
twoout{i} = [1 1];


% diagonal, Stratonovich, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstFFUN'',''yes'')';
twoout{i} = [1 1];


% non-diagonal, Stratonovich, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.5:10';      % constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.1:10';      % non-constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0,0]';    	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';    	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];


% diagonal, Ito, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];


% non-diagonal, Ito, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';

twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.5:10';      % constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.1:10';      % non-constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0,0]';    	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';    	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstFFUN'',''yes'')';
twoout{i} = [1 1];


% diagonal, Stratonovich, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''ConstGFUN'',''yes'')';
twoout{i} = [1 1];


% non-diagonal, Stratonovich, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.5:10';      % constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.1:10';      % non-constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0,0]';    	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';    	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];


% diagonal, Ito, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1'; 	% non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0;0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];


% non-diagonal, Ito, constant FFUN

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1'; 	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.5:0';	% constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '10:-0.1:0';	% non-constant step-size, decreasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.5:5';	% constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '-5:0.1:5';	% non-constant step-size, increasing, negative to positive
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.5:-5';	% constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '5:-0.1:-5';	% non-constant step-size, decreasing, positive to negative
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';   	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';       	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.5:10';      % constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';         % non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 3
tspan{i} = '0:0.1:10';      % non-constant step-size, increasing
y0{i} = '[0,0]';            % vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0,0]';     	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)[x,x]+1';	% non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0,0]';    	% vector ICs, N == 3
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '[0,0]';     	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x(1)+1';	% non-homogeneous, autonomous, scalar, non-diagonal, D == 1
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '[0,0]';    	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, vector
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, matrix, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)x+1';    	% non-homogeneous, autonomous, scalar
g{i} = '@(t,x)x+1';     % non-homogeneous, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';            % scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''SDEType'',''Ito'',''ConstGFUN'',''yes'')';
twoout{i} = [1 1];


% cases that miss no FOR loop case, D > N, nargout == 1:

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';      	% constant, autonomous, scalar
g{i} = '@(t,x)[1,1]';  	% constant, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.5:10';	% constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 0];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)1';      	% constant, autonomous, scalar
g{i} = '@(t,x)[1,1]';  	% constant, autonomous, vector, non-diagonal, D == 2
tspan{i} = '0:0.1:10';	% non-constant step-size, increasing
y0{i} = '0';          	% scalar ICs, N == 1
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 0];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1]';       % constant, autonomous, vector
g{i} = '@(t,x)ones(2,3)';	% constant, autonomous, scalar, non-diagonal, D == 3
tspan{i} = '0:0.5:10';   	% constant step-size, increasing
y0{i} = '[0,0]';          	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 0];

st = dbstack;
i = i+1;
lnum{i} = st.line;
f{i} = '@(t,x)[1;1]';       % constant, autonomous, vector
g{i} = '@(t,x)ones(2,3)';	% constant, autonomous, scalar, non-diagonal, D == 3
tspan{i} = '0:0.1:10';   	% non-constant step-size, increasing
y0{i} = '[0,0]';          	% vector ICs, N == 2
opts{i} = ',sdeset(''Diagonal'',''no'',''ConstFFUN'',''yes'',''ConstGFUN'',''yes'')';
twoout{i} = [1 0];



nb =  char(160);    % Non-breaking space
bs = '\b';          % Backspace;
sp = ' ';           % Space

m = length(lnum1);
num = 0;

if runtests
    ts = sort(tests(tests <= m));
else
    ts = 1:m;
end
lts = length(ts);
lne = 2*length(num2str(lts));

if ~isempty(ts)
    waittext(0,'init');
end
for i = ts
    errv = true;
    try 
        eval(cmd{i});
    catch err
        if ~strcmp(err.identifier,msg{i})
            errmsg = ['Error test ' num2str(i) ' (of ' num2str(m) ') '...
                      'assertion failed. See line ' num2str(lnum1{i}-1) ...
                      ' of SDE_MILSTEIN_UNITTEST.m  for details.\n\n' nb nb ...
                      'Error ID:' nb nb nb msg{i} '\n' nb nb 'Expression:'...
                      nb cmd{i}];
            me = MException('SDETools:sde_milstein_unittest:AssertError',errmsg);
            rep = getReport(err,'extended');
            me2 = MException('SDETools:sde_milstein_unittest:AssertErrorCause',rep);
            throw(addCause(me,me2));
        end
        errv = false;
    end
    if errv
        errmsg = ['Error test ' num2str(i) ' (of ' num2str(m) ') not '...
                  'triggered. See line ' num2str(lnum1{i}-1) ' of '...
                  'SDE_MILSTEIN_UNITTEST.m for details.\n\n' nb nb 'Error ID:'...
                  nb nb nb msg{i} '\n' nb nb 'Expression:' nb cmd{i}];
        me = MException('SDETools:sde_milstein_unittest:EvaluationError',errmsg);
        rep = getReport(err,'extended');
        me2 = MException('SDETools:sde_milstein_unittest:EvaluationErrorCause',rep);
        throw(addCause(me,me2));
    end
    num = num+1;
    if mod(num,5) == 0 || num == 2 || num == lts && lts > 1
        wopts = struct('length',35,'prefix',[int2str(num) ' (of ' ...
            int2str(lts) ') error tests successfull... ']);
        waittext(num/lts,'waitbar',wopts);
    end
end
if num == 1
    waittext('1 error test performed successfully.');
elseif num > 1
    waittext(['All ' int2str(num) ' error tests performed successfully.']);
end


pm = m;
m = length(lnum);

if runtests
    ts = tests(tests > pm)-pm;
    ts = ts(ts <= m*14);
    lts = length(ts);
	mm = lts+1;
else
    mm = 0;
    for j = 1:m
        if twoout{j}(1) == 1
            mm = mm+6;
        end
        if twoout{j}(2) == 1
            mm = mm+8;
        end
    end
    lts = mm;
end
num = 0;
cnt = 0;

lnu = 2*length(num2str(lts));

if ~isempty(ts) || ~runtests
    waittext(0,'init');
    ms = 1:m;
else
    ms = [];
end
for j = ms
    n = 0;
    if twoout{j}(1) == 1
        n = n+6;
    end
    if twoout{j}(2) == 1
        n = n+8;
    end
    
    if runtests
        tts = ts(ts <= cnt+n);
        tts = mod(tts(tts > cnt)-1,n)+1;
    else
        tts = 1:n;
    end
    
    if ~isempty(tts)
        lnum2 = cell(n,1);
        cmd = cell(n,1);
        test1 = cell(n,1);
        test2 = cell(n,1);

        if twoout{j}(1) == 1
            % Output is not empty
            st = dbstack;
            i = 1;
            lnum2{i} = st.line;
            cmd{i} = ['y=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'isempty(y)';
            test2{i} = false;

            % Output is float
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['y=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'isfloat(y)';
            test2{i} = true;

            % Output is matrix
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['y=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'ndims(y)';
            test2{i} = 2;

            % Output is same dimension as Y0
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['y=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'size(y,2)';
            ic = eval(y0{j});
            test2{i} = size(ic(:),1);

            % Length of output same as TSPAN
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['y=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'size(y,1)';
            test2{i} = length(eval(tspan{j}));

            % First value is Y0
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['y=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'y(1,:)''';
            ic = eval(y0{j});
            test2{i} = ic(:);
        end
        if twoout{j}(2) == 1    % Output and check Wiener increments
            % Output is not empty
            st = dbstack;
            if twoout{j}(1) == 1
                i = i+1;
            else
                i = 1;
            end
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = '[isempty(y) isempty(w)]';
            test2{i} = false;

            % Output is float
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = '[isfloat(y) isfloat(w)]';
            test2{i} = true;

            % Output is matrix
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = '[ndims(y) ndims(w)]';
            test2{i} = 2;

            % Output of Y is same dimension as Y0
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'size(y,2)';
            ic = eval(y0{j});
            test2{i} = size(ic(:),1);

            % Output of W is same dimension as output of Y0 or GFUN
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'size(w,2)';
            t = eval(tspan{j});
            ic = eval(y0{j});
            eg = eval(g{j});
            if isa(eg,'function_handle')
                gout = feval(eg,t(1),ic(:));
            else
                gout = eg;
            end
            if isempty(strfind(opts{j},'''Diagonal'',''no'''))
                test2{i} = size(ic(:),1);
            else
                test2{i} = size(gout,2);
            end

            % Length of output same as TSPAN
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = '[size(y,1) size(w,1)]';
            test2{i} = length(eval(tspan{j}));

            % First value of Y is Y0
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'y(1,:)''';
            ic = eval(y0{j});
            test2{i} = ic(:);

            % First value of W is 0
            st = dbstack;
            i = i+1;
            lnum2{i} = st.line;
            cmd{i} = ['[y,w]=sde_milstein(' f{j} ',' g{j} ',' tspan{j} ',' y0{j} opts{j} ');'];
            test1{i} = 'w(1,:)';
            test2{i} = 0;
        end

        for i = tts
            try
                eval(cmd{i});
            catch err
                errmsg = ['Error during evaluation of unit test ' num2str(cnt+i) ...
                          '. See lines ' num2str(lnum{j}) ' and '...
                          num2str(lnum2{i}-1) ' of sde_milstein_UNITTEST.m for '...
                          'details.\n\n' nb nb 'Expression:' nb cmd{i}];
                me = MException('SDETools:sde_milstein_unittest:EvaluationError',errmsg);
                rep = getReport(err,'extended');
                me2 = MException('SDETools:sde_milstein_unittest:EvaluationErrorCause',rep);
                throw(addCause(me,me2));
            end
            errmsg = ['Unit test ' num2str(cnt+i) ' assertion failed. See lines '...
                      num2str(lnum{j}) ' and ' num2str(lnum2{i}-1) ' of '...
                      'SDE_MILSTEIN_UNITTEST.m for details.\n' nb nb 'Expression:'...
                      nb cmd{i} '\n\n' nb nb 'Assertion:' nb nb test1{i} nb '=='...
                      nb mat2str(test2{i})];
            assert(all(eval(test1{i}) == test2{i}),sprintf(errmsg));
            num = num+1;
            if mod(num,5) == 0 || num == 2 || num == lts && lts > 1
                wopts = struct('length',35,'prefix',[int2str(num) ' (of ' ...
                    int2str(lts) ') unit tests successfull... ']);
                waittext(num/lts,'waitbar',wopts);
            end
        end
    end
    cnt = cnt+n;
end
if num == 0
    if ~isempty(ms)
        waittext(0,'clear');
    end
elseif num == 1
    waittext('1 unit test performed successfully.');
else
    waittext(['All ' int2str(num) ' unit tests performed successfully.']);
end



function f1(t,x)        %#ok<*DEFNU,*INUSL>
y = x+1;                %#ok<*NASGU>


function y = f2(t,x)    %#ok<*STOUT>
x+1;                    %#ok<*VUNUS>


function g1(t,x)
y = x+1;


function y = g2(t,x)
x+1;


function r1(m,n)
y = randn(m,n);


function y = r2(m,n)
randn(m,n);