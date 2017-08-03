function Out = vcRobustMin(fcn,x0,MinParams)

% vcRobustMin
% 
% This routine builds on Chris Sims routine csminwel to make its results
% more robust.
%
% More in detail the robustness comes from running csminwel several times
% in succession, in which each guess will use the previous output, shock it
% and use the result as the new guess value. Furthermore the initial
% inverse hessian for each iteration is the inverse hessian from the
% previous one.
%
% Usage:
%   Out = vcRobustMin(fcn,x0)
%   Out = vcRobustMin(fcn,x0,MinParams)
%
% Inputs:
%   fcn
%     String with the name of the function to be minimized 
%   x0
%     Column vector with the guess values
%   MinParams (OPTIONAL): 
%     Can be set in one of two ways:
%     1) Cell array with input arguments for the function being minimized. 
%        No other minimization parameters are set.
%     2) Structure with new values for the parameters of the minimization 
%        and robustness check. The structure fields need to be among the 
%        following options, but do not have to include all of them - this 
%        way in a specific application only needs to set the parameters for 
%        which the default is not good.
%        2.1) parameters to be used in csminwel
%          MinParams.H0
%            Initial hessian. Default: eye(length(x0))
%          MinParams.grad
%            Name of gradient function. Default: []
%          MinParams.crit
%            Criterion for convergence. Default: 1e-7
%          MinParams.nit
%            Maximum number of iterations. Default: 1000
%          MinParams.varargin [cell array]
%            Additional arguments required to evaluate function. This is 
%            required if the function requires arguments. Default: {}
%        2.2) parameters needed for robustness part
%          MinParams.Rcrit 
%            Criterion for robustness convergence. Default: 1e-7
%          MinParams.Ritmax
%            Maximum number of robustness iterations. Default: 0 (implies
%            that no robustness is performed)
%          MinParams.Ritmin
%            Minimum number of robustness iterations. Default: 5
%          MinParams.Ritnoimp
%            Maximum number of iterations without improvement before 
%            robustness is declared. Default: 5
%          MinParams.RH0 
%            If set to 1 then the robust shock is based on the initial H0, 
%            if set to 0 then the robust shock is based on the latest H. 
%            Default: 0
%          MinParams.Rscaledown
%            Scaling of previous inverse hessian used for drawing candidate
%            guess value. Default: 0.1 (equivalent to 31.6% of SE)
%        2.3) Other parameters
%          MinParams.AltHessian
%            If set to 1, Hessian is computed using vcHessian,
%            if set to 0, Hessian is computed in csminwel.
%            Default: 0
% 
% Output:
%   The output of the minimization is introduced in a structure, as follows
%   1) Output typical of csminwel
%     Out.x
%       Argmin for this function
%     Out.f
%       Value of function at its minimum
%     Out.g
%       Gradient at minimum
%     Out.H
%       Hessian at minimium
%     Out.itct
%       Iteration count from last call of csminwel at minimum
%     Out.fc
%       Number of gradients computed in csminit at minimum
%     Out.rc
%       Return code from csminwel at minimum
%     Out.rcMsg
%       Message explaining return code from csminwel
%       rc = 0 => Normal solution
%       rc = 1 => Zero gradient
%       rc = 2 or rc = 4 => Back and forth on step length never finished
%       rc = 3 => Smallest step still improving too slow
%       rc = 5 => Largest step still improving too fast
%       rc = 6 => Smallest step still improving too slow, reversed gradient
%       rc = 7 => Warning: possible inaccuracy in H matrix
%       rc = 8 => Invalid initial guess
%   2) Output from robust analysis
%     Out.Ritct
%       Number of iterations of robustness
%     Out.Rrc
%       Return code for robustness
%     Out.RrcMsg
%       Message describing return code for robustness
%       Rrc = 0 => no robustness
%       Rrc = 1 => improvement in function below Rcrit for Ritnoimp
%                  consecutive iterations
%       Rrc = 2 => maximum number of iterations exceeded
%   3) Include the csminwel results for initial and robust iterations
%     Out.MinOutput
%       Structure array with 1+Ritct elements. The fields are the output of
%       csminwel for the corresponding minimization.
%   3) also list parameters used in minimization
%     Out.MinParams
%       Structure with all the parameters used, including guess value
%
%
% See also: 
% csminwel, csminit, numgrad, csolve, vcNumJacobian, vcHessian
%
% ..............................................................................
%
% Created: February 5, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2008-2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Default parameter values
if nargin==3
    if iscell(MinParams)
        MinParamsIn.varargin = MinParams;
    elseif isstruct(MinParams)
        MinParamsIn = MinParams; 
    else
        error('Wrong input: MinParams needs to be set as either a cell or structure array')
    end
    clear MinParams
end

% parameters needed for csminwel
MinParams.x0 = x0;
MinParams.H0 = eye(length(x0));
MinParams.grad = [];
MinParams.crit = 1e-7;
MinParams.nit = 1000;
MinParams.varargin = {};

% parameters needed for robustness part
MinParams.Rcrit = 1e-7;
MinParams.Ritmax = 0;
MinParams.Ritmin = 5;
MinParams.Ritnoimp = 5;
MinParams.RH0 = 0;
MinParams.Rscaledown = 0.1;

% other parameters
MinParams.AltHessian = 0;

%% Update parameter values
MinParamsFields = fieldnames(MinParams);
if nargin==1
    error('Need to specify at least function name and guess as separate arguments!')
elseif nargin==3
    MinParamsInFields = fieldnames(MinParamsIn);
    for j=1:length(MinParamsInFields)
        if ~ismember(MinParamsInFields{j},MinParamsFields)
            error('Wrong parameter field name: %s',MinParamsInFields{j})
        else
            MinParams.(MinParamsInFields{j}) = MinParamsIn.(MinParamsInFields{j});
        end
    end
elseif nargin>3
    error('Too many input arguments!')
end
Out.MinParams = MinParams;

%% Check consistency
% MinParams.Ritmin = min(MinParams.Ritmin,MinParams.Ritmax);
% MinParams.Ritnoimp = min(MinParams.Ritnoimp,MinParams.Ritmax);

%% Display parameters used:
disp(' ')
disp('Parameters used:')
disp(' ')
disp(MinParams)
disp(' ')
disp('Guess vector:')
x0
Out.x0 = x0;

%% -----------------------------------------------------------------------------

%% Check if it is good guess
f0 = feval(fcn,x0,MinParams.varargin{:})
Out.f0=f0;
if f0>1e50
    fprintf('\nWARNING: Invalid initial guess. Cannot proceed.')
    Out.x = x0;
    Out.f = f0;
    Out.g = [];
    Out.H = [];
    Out.itctr = 0;
    Out.fcr = 0;
    Out.rc = 8;
    Out.rcMsg = 'Invalid initial guess';
    Out.Ritct = 0;
    Out.Rrc = 0;
    Out.RrcMsg = 'No Robustness performed';
    Out.MinOutput = [];
    Out.MinParams = MinParams;
    return
end

%% Make first minimization
MinOutputFields = {'fh','xh','gh','Hh','itct','fc','rc'};
nMinOutputFields = length(MinOutputFields);
[fh,xh,gh,Hh,itct,fc,rc]=csminwel(fcn,x0,MinParams.H0,MinParams.grad,...
    MinParams.crit,MinParams.nit,MinParams.varargin{:});
if MinParams.AltHessian
    Hh = vcHessian(fcn,xh,MinParams.varargin{:});
end
for j=1:nMinOutputFields
    MinOutput(1).(MinOutputFields{j}) = eval(MinOutputFields{j});
end
MinOutput(1).rcMsg = rcErrorMsg(rc);

%% -----------------------------------------------------------------------------

%% Robustness

noimp = 0;
xhr = xh;
fhr = fh; fh0 = fh;
Hhr = Hh;
ghr = gh;
Hh_backup = MinParams.H0;
x0_backup = x0;
rcr = rc;
itctr = itct;
fcr = fc;

for rit=1:MinParams.Ritmax
    % check for semi positive definite inverse hessian
    [Ur,Dr] = eig((Hhr+Hhr')/2);
    Dr = diag(Dr);
    tol = eps(max(Dr)) * length(Dr);
    idxD = (abs(Dr) > tol);
    Dr = Dr(idxD);
    negeig = sum(Dr<0); % number of negative eigenvalues
    if negeig~=0
        disp('Warning: Previous Hessian is not positive semi definite!')
        disp('         Using earlier best alternative...')
        Hhr = Hh_backup;
    else
        Hh_backup = Hhr;
    end
    % choose H for robust shock
    if MinParams.RH0
        RShockVar = MinParams.H0;
    else
        RShockVar = Hhr;
    end
    x0r = mvnrnd(xhr,MinParams.Rscaledown*RShockVar)';
    
    % check whether it is an acceptable candidate
    GoodC = (feval(fcn,x0r,MinParams.varargin{:})<1e50);
    for j=1:1000
        if GoodC, break, end
        x0r = mvnrnd(xhr,MinParams.Rscaledown*RShockVar)';
        GoodC = (feval(fcn,x0r,MinParams.varargin{:})<1e50);
    end
    if ~GoodC
        disp('Warning: could not find suitable candidate after 1000 attempts.')
        disp('         Using last good candidate...')
        x0r = x0_backup;
    end
    
    [fh,xh,gh,Hh,itct,fc,rc]=csminwel(fcn,x0r,Hhr,MinParams.grad,...
        MinParams.crit,MinParams.nit,MinParams.varargin{:});
    if MinParams.AltHessian
        Hh = inv(vcHessian(fcn,xh,MinParams.varargin{:}));
    end
    for j=1:nMinOutputFields
        MinOutput(1+rit).(MinOutputFields{j}) = eval(MinOutputFields{j});
    end
    MinOutput(1+rit).rcMsg = rcErrorMsg(rc);

    disp(' ')
    disp(sprintf('Results from robustness iteration %.0f',rit))
    disp(sprintf('old function value: %.10f', fhr))
    disp(sprintf('new function value: %.10f', fh))
    dfh=fh-fhr;
    disp(sprintf('change in function value: %.10f', dfh))
    disp(' ')
    rf(rit) = fh;
    df(rit) = dfh;
    if -dfh>=MinParams.Rcrit
        xhr = xh;
        fhr = fh;
        Hhr = Hh;
        ghr = gh;
        rcr = rc;
        itctr = itct;
        fcr = fc;
        noimp=0;
    else
        noimp = noimp + 1;
        if (noimp>=MinParams.Ritnoimp) && (rit>=MinParams.Ritmin)
            Rrc = 1;
            RrcMsg = sprintf('Improvement in function below %g for %.0f consecutive iterations',...
                MinParams.Rcrit, MinParams.Ritnoimp);
            break
        end
    end
end

%% Checks
if MinParams.Ritmax==0
    Rrc = 0;
    RrcMsg = 'No robustness performed';
    rf = [];
    rit = 0;
elseif rit==MinParams.Ritmax
    Rrc = 2;
    RrcMsg = 'Maximum number of iterations exceeded';
end

%% print history of posterior values and improvements
% % fprintf('\nFunction value before robustness: %f\n', fh0)
% fprintf('Initial iterarion, function value: %.10f',fh0)
% fprintf('(Stopped at iteration %.0f)\n',MinOutput(1).itct)
% if Rrc~=0
%     fprintf('\nResults from robustness check:')
%     fprintf('\n------------------------------\n\n')
%     for rit=1:length(rf)
%         fprintf('Robustness iterarion %.0f, function value: %.10f, change %.10f',rit,rf(rit),df(rit))
%         fprintf('(Stopped at iteration %.0f)\n',MinOutput(1+rit).itct)
%     end
%     fprintf('\nFunction value after robustness: %f\n', fhr)
% end
% fprintf('\n%s\n',RrcMsg)

fprintf('\nRobustness analysis for minimization')
fprintf('\n------------------------------------\n')
fprintf('Iteration %3.0f: function value: %15.8f',0,MinOutput(1).fh)
% fprintf('         %15s','')
fprintf(' change: %15.8f',MinOutput(1).fh-f0)
fprintf(' stopped at iteration %.0f\n',MinOutput(1).itct)
itbest.fh = MinOutput(1).fh;
itbest.idx = 1;
for jr=2:length(MinOutput)
    fprintf('Iteration %3.0f: function value: %15.8f',jr-1,MinOutput(jr).fh)
    itchg = MinOutput(jr).fh-itbest.fh;
    fprintf(' change: %15.8f',itchg)
    fprintf(' stopped at iteration %.0f\n',MinOutput(jr).itct)
    if itchg<0
        itbest.fh = MinOutput(jr).fh;
        itbest.idx = jr;
    end
end
fprintf('\nBest iteration: %.0f',itbest.idx-1)
fprintf('\nBest iteration function value: %.8f\n',itbest.fh)
fprintf('\nRobustness change over initial min: %.8f',itbest.fh-MinOutput(1).fh)
fprintf('\nRobustness message: %s\n',RrcMsg)
fprintf('\nBest iteration message: %s\n\n',MinOutput(itbest.idx).rcMsg)

%% check for positive semi definite hessian
[U,D] = eig((Hhr+Hhr')/2);
D = diag(D);
tol = eps(max(D)) * length(D);
idxD = (abs(D) > tol);
D = D(idxD);
negeig = sum(D<0); % number of negative eigenvalues
if negeig~=0
    disp('Warning: Hessian is not positive semi definite!')
    disp('         Using last good one...')
    Hhr = Hh_backup;
    if all(all(Hh_backup==MinParams.H0))
        disp('Backup Hessian is initial guess!')
    end
end

%% -----------------------------------------------------------------------------

%% Prepare output

% Output typical of csminwel, for the optimal solution
Out.x = xhr;
Out.f = fhr;
Out.g = gh;
Out.H = Hhr;
Out.itctr = itctr;
Out.fcr = fcr;
Out.rc = rcr;
Out.rcMsg = rcErrorMsg(rcr);

% Output from robust analysis
if Rrc~=0
    Out.Ritct = rit;
else
    Out.Ritct = 0;
end
Out.Rrc = Rrc;
Out.RrcMsg = RrcMsg;

% for each minimization save the output:
Out.MinOutput = MinOutput;

% also list parameters used in minimization
Out.MinParams = MinParams;

%% -----------------------------------------------------------------------------

%% Subfunction: rcErrorMsg(rc)
function rcMsg = rcErrorMsg(rc)
if rc == 1
    rcMsg = 'Zero gradient';
elseif rc == 6
    rcMsg = 'Smallest step still improving too slow, reversed gradient';
elseif rc == 5
    rcMsg = 'Largest step still improving too fast';
elseif (rc == 4) || (rc==2)
    rcMsg = 'Back and forth on step length never finished';
elseif rc == 3
    rcMsg = 'Smallest step still improving too slow';
elseif rc == 7
    rcMsg = 'Warning: possible inaccuracy in H matrix';
else
    rcMsg = 'Normal solution';
end

%% -----------------------------------------------------------------------------



