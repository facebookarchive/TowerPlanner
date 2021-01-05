%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config_, node, elem] = ssTowerDesign (config, app)

XLB = [0.5 1 1 1];
XUB = [1 10 20 20];
X0 = XUB;

key = 'MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]';
tMax = 20;
tic;
problem.func = @(X)problem_function(X,config, app, tMax);
problem.o  = 1; % Number of objectives
problem.n  = 4; % Number of variables (in total)
problem.ni = 3; % Number of integer variables (0 <= ni <= n)
problem.m  = 5; % Number of constraints (in total)
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

% build tower model
config0.IOfilename = num2str(randi(1e6));
config0.mast.nGuyLevels = 0;
config0.mast.model = 'lattice';

if strcmp(config0.units.LMT, 'metric')
pipeOD = [0.4050    0.5400    0.6750    0.8400    1.0500    1.3150    1.6600...
    1.9000    2.3750    2.8750    3.5000    4.0000    4.5000    5.5630...
    6.625 8.625 10.750 12.750 14 16]*0.0254;
else
pipeOD = [0.4050    0.5400    0.6750    0.8400    1.0500    1.3150    1.6600...
    1.9000    2.3750    2.8750    3.5000    4.0000    4.5000    5.5630...
    6.625 8.625 10.750 12.750 14 16]/12;
end

listMastWidth = linspace(config0.mast.height*0.05, config0.mast.height*0.2, 10);
config0.mast.taper = X(1);
config0.mast.width = listMastWidth(X(2));
config0.mast.leg.D = pipeOD(X(3));
config0.mast.web.D = pipeOD(X(4));

[config0, node0]        = createNodes (config0);
[config0, node0, elem0] = createElements (config0, node0);
[config0, node0, elem0] = assignProps (config0, node0, elem0);
[config0, node0, elem0] = assignLoads (config0, node0, elem0);

% solve
solve (config0, node0, elem0);
[config0, node0, elem0] = postProcess (config0, node0, elem0);

if ~config0.opt.fail
   
    % leg constraints
    if strcmp(config0.mast.model,'lattice')
        cLegC  = max([elem0(strcmp({elem0(:).type},'leg')).stressCrit]) - 1;
        cWebC  = max([elem0(strcmp({elem0(:).type},'web')).stressCrit]) - 1;
    elseif strcmp(config0.mast.model,'beam')
        cLegC  = max([elem0(idx).stressLegCrit]) - 1;
        cWebC  = max([elem0(idx).stressWebCrit]) - 1;
    end
    
    % SR constraints  
    dLeg = config0.mast.leg.D;
    tLeg = dLeg/config0.mast.leg.TR;
    rLeg   = sqrt(dLeg^2 + (dLeg-2*tLeg)^2)/4;    
    dWeb = config0.mast.web.D;
    rWeb = dWeb/4;
    
    legSR  = config0.mast.width*...
        tand(config0.mast.angle)/rLeg;
    webSR  = config0.mast.width/...
        cosd(config0.mast.angle)/rWeb;
    cLegSR = legSR/config0.mast.leg.SRMax - 1;   
    cWebSR = webSR/config0.mast.web.SRMax - 1;
            
    % size foundation/anchors
    [config0, node0, elem0]  = sizeFoundation (config0, node0, elem0);
    [config0, node0, elem0]  = cost (config0, node0, elem0);
    
    % antenna displacement constraints
    cAntD = -1;
    if isfield(config0,'antenna')
        nIdx  = [config0.antenna(:).node];
        disp  = [node0(nIdx).disp];
        cAntD = max([cAntD,(max(abs(disp(4:6,:)))./[config0.antenna.maxDisp]-1)]);   
    end 
    
    ceq = [];
    c = [cLegC;cWebC;cLegSR;cWebSR;cAntD];
    towerCost = config0.cost.total;
else
    ceq = [];
    c = [1;1;1;1;1;1;1]*100;
    config0.cost.total = 1e9;
    towerCost = 1e9;
end
end




