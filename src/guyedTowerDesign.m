%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config_, node, elem] = guyedTowerDesign (config, app)

XLB = [1 1 1];
XUB = [10 250 50];
X0 = XUB;

key = 'MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]';
tMax = 20;
tic;
problem.func = @(X)problem_function(X,config, app, tMax);
problem.o  = 1; % Number of objectives
problem.n  = 3; % Number of variables (in total)
problem.ni = 3; % Number of integer variables (0 <= ni <= n)
problem.m  = 8; % Number of constraints (in total)
problem.me = 0; % Number of equality constraints (0 <= me <= m)
problem.xl = XLB;
problem.xu = XUB;
problem.x  = X0;
option.maxeval   = 10000;    % Maximum number of function evaluation (e.g. 1000000)
option.maxtime   = tMax ; % Maximum time limit in Seconds (e.g. 1 Day = 60*60*24)
option.printeval = 0;  % Print-Frequency for current best solution (e.g. 1000)
option.save2file = 0;     % Save SCREEN and SOLUTION to TXT-files [ 0=NO/ 1=YES]
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

% build tower model

config0.IOfilename = num2str(randi(1e6));

if strcmp(config0.units.LMT, 'metric')
    pipeOD = [0.4050    0.5400    0.6750    0.8400    1.0500    1.3150    1.6600...
        1.9000    2.3750    2.8750    3.5000    4.0000    4.5000    5.5630]*0.0254;
    listMastWidth = linspace(0.25, 2, 5);
    listGuyDia = [3/16 1/4 5/16 3/8 7/16 1/2 9/16 5/8 3/4 7/8]*0.0254;
else
    pipeOD = [0.4050    0.5400    0.6750    0.8400    1.0500    1.3150    1.6600...
        1.9000    2.3750    2.8750    3.5000    4.0000    4.5000    5.5630]/12;
    listMastWidth = linspace(0.25, 2, 5)*3.28084;
    listGuyDia = [3/16 1/4 5/16 3/8 7/16 1/2 9/16 5/8 3/4 7/8]/12;
end

listN = 1:15;
listLegDia = pipeOD(5:end);
listWebDia = pipeOD(1:5);
elements = {listMastWidth, listLegDia, listWebDia};
combinations = cell(1, numel(elements));
[combinations{:}] = ndgrid(elements{:});
combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false);
listMastSet = [combinations{:}];


listGuyTension = linspace(0.10,0.20,5);
[A,B] = meshgrid(listGuyDia,listGuyTension);
c = cat(2,A',B');
listGuySet = reshape(c,[],2);

X = fix(X);

config0.mast.nGuyLevels = listN(X(1));
config0.mast.width      = listMastSet(X(2), 1);
config0.mast.leg.D      = listMastSet(X(2), 2);
config0.mast.web.D      = listMastSet(X(2), 3);

config0.guy.D   = listGuySet(X(3),1);
config0.guy.initialTension = listGuySet(X(3),2);

[config0, node0]        = createNodes (config0);
[config0, node0, elem0] = createElements (config0, node0);
[config0, node0, elem0] = assignProps (config0, node0, elem0);
[config0, node0, elem0] = assignLoads (config0, node0, elem0);

% solve

solve (config0, node0, elem0);

[config0, node0, elem0] = postProcess (config0, node0, elem0);

% objectives and constraints

if ~config0.opt.fail
   
    % leg constraints
    idx = strcmp({elem0(:).type},'leg');
    if strcmp(config0.mast.model,'lattice')
        cLegC  = max([elem0(idx).stressCrit]) - 1;
        cWebC  = cLegC;
    elseif strcmp(config0.mast.model,'beam')
        cLegC  = max([elem0(idx).stressLegCrit]) - 1;
        cWebC  = max([elem0(idx).stressWebCrit]) - 1;
    end
    
    % guy constraints
    idx = strcmp({elem0(:).type},'guy');
    cGuyS =  -min([elem0(idx).stressCrit]) - 0.0; %eps; - 0.20; -1; -0;
    cGuyT =  max([elem0(idx).stressCrit]) - 1;
    
    % mast constraints
    baseNode = find(strcmp({node0(:).l},'base'));
    idx = find([elem0(:).node1]==baseNode | [elem0(:).node2]==baseNode);
    A = sum([elem0(idx).A]);
    Lg = min([config0.guy.level]);
    P = 0;
    for i = 1:length(idx)
        P  = P + abs(elem0(idx(i)).intLoad.axMin);
    end
    sigY = config0.mat.steel.sigY.strut*config0.mat.steel.phiC;
    a = 1/7500;
    P = P - sum([elem0(idx).rho].*[elem0(idx).A]*Lg*config0.env.g);
    mastSR = abs(min(90,sqrt((sigY*A/P-1)/a)));
    
    dLeg = config0.mast.leg.D;
    tLeg = dLeg/config0.mast.leg.TR;
    rLeg   = sqrt(dLeg^2 + (dLeg-2*tLeg)^2)/4;
    
    dWeb = config0.mast.web.D;
    rWeb = dWeb/4;
    
    legSR  = config0.mast.width*...
        tand(config0.mast.angle)/rLeg;
    webSR = config0.mast.width/...
        cosd(config0.mast.angle)/rWeb;
    S  = Lg*sqrt(6)/mastSR;
    cMastSR = S/config0.mast.width - 1;
    cLegSR = legSR/mastSR - 1;   
    cWebSR = webSR/config0.mast.web.SRMax - 1;
    
    % antenna displacement/guy attachment constraints
    cAntD = -1;
    cAntG = -1;
    if isfield(config0,'antenna')
        nIdx  = [config0.antenna(:).node];
        disp  = [node0(nIdx).disp];
        cAntD = max([cAntD,(max(abs(disp(4:6,:)))./[config0.antenna.maxDisp]-1)]);      
        cAntG = max(-[config0.antenna(:).pairedWithGuy] + 0.5);
    end   

    % size foundation/anchors
    [config0, node0, elem0]  = sizeFoundation (config0, node0, elem0);
    [config0, node0, elem0]  = sizeAnchors (config0, node0, elem0);
    [config0, node0, elem0]  = cost (config0, node0, elem0);
    
    ceq = [];
    c = [cLegC;cWebC;cGuyT;cMastSR;cLegSR;cWebSR;cAntD;cAntG];
    towerCost = config0.cost.total;
else
    ceq = [];
    c = ones(9,1)*100;
    config0.cost.total = 1e9;
    towerCost = 1e9;
end
end




