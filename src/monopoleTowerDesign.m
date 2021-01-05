%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config_, node, elem] = monopoleTowerDesign (config, app)

if strcmp(config.units.LMT, 'metric')
    XLB = [0.1 2 0.5];
    XUB = [1 20 1];
else
    XLB = [0.1*3.28084 2 0.5];
    XUB = [1*3.28084 20 1];
end
X0 = XUB;

key = 'MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]';
tMax = 20;
tic;
problem.func = @(X)problem_function(X,config, app, tMax);
problem.o  = 1; % Number of objectives
problem.n  = 3; % Number of variables (in total)
problem.ni = 0; % Number of integer variables (0 <= ni <= n)
problem.m  = 2; % Number of constraints (in total)
problem.me = 0; % Number of equality constraints (0 <= me <= m)
problem.xl = XLB;
problem.xu = XUB;
problem.x  = X0;
option.maxeval   = 10000; % Maximum number of function evaluation (e.g. 1000000)
option.maxtime   = tMax;  % Maximum time limit in Seconds (e.g. 1 Day = 60*60*24)
option.printeval = 0;   % Print-Frequency for current best solution (e.g. 1000)
option.save2file = 0;   % Save SCREEN and SOLUTION to TXT-files [ 0=NO/ 1=YES]
option.param( 1) =  0;  % ACCURACY
option.param( 2) =  0;  % SEED
option.param( 3) =  0;  % FSTOP
option.param( 4) =  0;  % ALGOSTOP
option.param( 5) =  0;  % EVALSTOP
option.param( 6) =  0;  % FOCUS
option.param( 7) =  0;  % ANTS
option.param( 8) =  0;  % KERNEL
option.param( 9) =  0;  % ORACLE
option.param(10) =  0;  % PARETOMAX
option.param(11) =  0;  % EPSILON
option.param(12) =  0;  % BALANCE
option.param(13) =  0;  % CHARACTER
option.parallel = 1;
[ solution ] = midaco( problem, option, key);
[c, ~, ~, config_, node, elem] = evaluate (solution.x, config);
X = solution.x;
c

end

function [f, g] = problem_function( X, config0, app0, tMax0)

[c, ~, towerCost] = evaluate (X, config0);

% Objective functions F(X)
f = towerCost;

% Equality constraints G(X) = 0 MUST COME FIRST in g(1:me)
g = -c;

app0.Status.Text = ['Progress: ' num2str(min([(round(toc/tMax0*100)),100])) '%'];
drawnow();
end


function [c, ceq, towerCost, config0, node0, elem0] = evaluate (X, config0)

% build tower 
config0.IOfilename = num2str(randi(1e6));
config0.mast.nBeamElem = 10;
config0.mast.model = 'beam';
config0.mast.horBrace = false;
config0.mast.fwdBrace = false;
config0.mast.bwdBrace = false;
config0.mast.width = 0;
config0.mast.angle = 0;
config0.mast.nGuyLevels = 0;
config0.mast.Xsection = [];

config0.mast.taper = X(3);
config0.mast.leg.D = linspace(X(1), X(1)*X(3), config0.mast.nBeamElem);
config0.mast.leg.TR = X(2);

config0.mast.web.D = 0;

[config0, node0]        = createNodes (config0);
[config0, node0, elem0] = createElements (config0, node0);
[config0, node0, elem0] = assignProps (config0, node0, elem0);
[config0, node0, elem0] = assignLoads (config0, node0, elem0);

% solve
solve (config0, node0, elem0);
[config0, node0, elem0] = postProcess (config0, node0, elem0);

if ~config0.opt.fail
   
    % leg constraints
    cLegC  = max([elem0(:).stressLegCrit]) - 1;
    cLegTR = config0.mast.leg.TR/config0.mast.leg.TRMax - 1;
           
    % size foundation/anchors
    [config0, node0, elem0]  = sizeFoundation (config0, node0, elem0);
    [config0, node0, elem0]  = cost (config0, node0, elem0);
    
    % antenna displacement/guy attachment constraints
    cAntD = -1;
    if isfield(config0,'antenna')
        nIdx  = [config0.antenna(:).node];
        disp  = [node0(nIdx).disp];
        cAntD = max([cAntD,(max(abs(disp(4:6,:)))./[config0.antenna.maxDisp]-1)]);    
    end 
    
    ceq = [];
    c = [cLegC; cLegTR; cAntD];
    towerCost = config0.cost.total;
else
    ceq = [];
    c = [1;1;1]*100;
    config0.cost.total = 1e9;
    towerCost = 1e9;
end
end



