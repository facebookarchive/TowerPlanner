%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function results = execute(app, mode)

if ispc
    copyfile 'thirdpartylibs/win/frame3dd.exe' 
elseif ismac
    copyfile 'thirdpartylibs/mac/frame3dd'
elseif isunix
    copyfile 'thirdpartylibs/unix/frame3dd' 
end

config = readData (app);

if strcmp(mode, 'analyze')
    [config, node, elem] = analyze(config);
elseif strcmp(mode, 'design')
     [config, node, elem] = design(config);       
end

% Structural Data
structData.elemType = {elem(:).type}';
structData.node1CoordX = [node([elem(:).node1]).x]';
structData.node1CoordY = [node([elem(:).node1]).y]';
structData.node1CoordZ = [node([elem(:).node1]).z]';

structData.node2CoordX = [node([elem(:).node2]).x]';
structData.node2CoordY = [node([elem(:).node2]).y]';
structData.node2CoordZ = [node([elem(:).node2]).z]';

structData.elemL = [elem(:).L]';
structData.elemE = [elem(:).E]';
structData.elemG = [elem(:).G]';
structData.elemA = [elem(:).A]';
structData.elemRho = [elem(:).rho]';
structData.elemEPA = [elem(:).EPA]';
structData.stressCrit = [elem(:).stressCrit]'*100;
structData.stressMode = {elem(:).stressCritMode}';
structData.stressMaxAllow = [elem(:).stressMaxAllow]';

% Cost Data
costData = config.cost.breakdown;


% Design Data
designData.LegDia = max(config.mast.leg.D);
designData.LegTaper = config.mast.taper;
designData.LegTRDesign = config.mast.leg.TR;
designData.WebDia =  config.mast.web.D;
designData.WebTRDesign = config.mast.web.TR;
designData.BracingAngle = config.mast.angle;
designData.MastWidth =  config.mast.width;
designData.nGuyLevels = config.mast.nGuyLevels;
designData.GuyDia = config.guy.D;
designData.GuyIT = config.guy.initialTension;
designData.SiteRadius = config.mast.siteRadius;

% csv File
T = table ({elem(:).type}', [node([elem(:).node1]).x]', [node([elem(:).node1]).y]',...
    [node([elem(:).node1]).z]', [node([elem(:).node2]).x]', [node([elem(:).node2]).y]',...
    [node([elem(:).node2]).z]', [elem(:).L]', [elem(:).E]', [elem(:).G]', [elem(:).A]', [elem(:).rho]',...
    [elem(:).EPA]', [elem(:).stressCrit]',  {elem(:).stressCritMode}', [elem(:).stressMaxAllow]',...
    'VariableNames',{'Member Type','X1','Y1','Z1','X2','Y2','Z2','Length','YoungsModulus','ShearModulus',...
    'Area', 'Density', 'EPA', 'FailureIndex','FailureMode','MaxAllowableStress'});
delete('results.csv');
delete('frame3dd*');
if ismac || isunix
delete('nul');
end
writetable(T,'results.csv');

results.structData = structData;
results.costData = costData;
results.designData = designData;

end

function config = readData (app)
            
            
            % general 
            config = [];
            switch app.TowerType
                case 'Monopole'
                     config.type = 'monopole';
                     app.Model = 'beam';
                     app.nGuyLevels = 0;
                case 'Self-Supported'
                    config.type = 'ss';
                    app.Model = 'lattice';
                    app.nGuyLevels = 0;
                case 'Guyed'
                    config.type = 'guyed';
                    app.Model = 'beam';
            end
            
            config.IOfilename = num2str(randi(1e6));
            
            % mast
            config.units.LMT = 'metric';
            config.mast.height = app.TowerHeight;
            config.mast.angle = app.WebAngle;
            config.mast.width = app.MastWidth;
            config.mast.taper = app.LegTaper;
            config.mast.Xsection = lower(app.CrossSection);
            config.mast.nGuyLevels = app.nGuyLevels;
            config.mast.siteRadius = app.SiteRadius;
            config.mast.model = lower(app.Model);
            
            config.mast.fwdBrace = false;
            config.mast.bwdBrace = false;       
            config.mast.horBrace = false;
            switch app.BracingType
                case 'X'
                    config.mast.fwdBrace = true;
                    config.mast.bwdBrace = true;
                case '/-'
                    config.mast.fwdBrace = true;
                    config.mast.horBrace = true;
                case '\-'
                    config.mast.bwdBrace = true;
                    config.mast.horBrace = true;
                case '-X-'
                    config.mast.fwdBrace = true;
                    config.mast.bwdBrace = true;
                    config.mast.horBrace = true;
            end
            
            if strcmp(config.type,'monopole')
                config.mast.nBeamElem = 10;
                config.mast.leg.D = linspace(app.LegDia, app.LegDia*config.mast.taper, config.mast.nBeamElem);
            else
                config.mast.leg.D = app.LegDia;
            end
            
            config.mast.leg.TR = app.LegTRDesign;
            config.mast.leg.TRMax = app.LegTRMax;
            config.mast.leg.SRMax = app.LegSRMax;
            config.mast.leg.taper = 1;             
            config.mast.web.D = app.WebDia;        
            config.mast.web.TR = app.WebTRDesign;
            config.mast.web.TRMax = app.WebTRMax;
            config.mast.web.SRMax = app.WebSRMax; 
          
            % guys
            config.anchor.nGuysPerAnchor = 5;
            config.guy.D = app.GuyDia;
            config.guy.initialTension = app.GuyIT;
            
            % antenna          
            if iscell(app.AntennaTable)
                AntennaTable = [];
                for k=1:length(app.AntennaTable)
                AntennaTable(k,:) = cell2mat(app.AntennaTable{k});
                end
                app.AntennaTable = AntennaTable;
            end
            
            idx = find((app.AntennaTable(:,1))>0);           
            for i = 1:length(idx)
            config.antenna(i).EPA    = (app.AntennaTable(idx(i),1));
            config.antenna(i).weight = (app.AntennaTable(idx(i),2));
            config.antenna(i).height = (app.AntennaTable(idx(i),3))/100*config.mast.height;            
            config.antenna(i).freq   = (app.AntennaTable(idx(i),4));
            config.antenna(i).dia    = (app.AntennaTable(idx(i),5));           
            config.antenna(i).power  = (app.AntennaTable(idx(i),6));
            
            if ~config.antenna(i).freq
                config.antenna(i).freq = 1.9;
            end          
            if ~config.antenna(i).dia
                config.antenna(i).dia = 0.1;
            end
            if ~config.antenna(i).height
                config.antenna(i).height = config.mast.height;
            end
            if ~config.antenna(i).power
                config.antenna(i).power = 10;
            end
            
            switch app.Units
                case 'Metric'
                    config.antenna(i).maxDisp = min(max(0.25, 16.2/(config.antenna(i).freq*config.antenna(i).dia)),5);
                case 'English'
                    config.antenna(i).maxDisp = min(max(0.25, 53.1/(config.antenna(i).freq*config.antenna(i).dia)),5);
            end
            config.antenna(i).pairedWithGuy = false;
            end
            
            if idx
            config.AntennaRFPropHeight = config.antenna(app.AntennaRFPropIdx).height;
            config.AntennaRFPropPower  = config.antenna(app.AntennaRFPropIdx).power;
            config.AntennaRFPropFreq   = config.antenna(app.AntennaRFPropIdx).freq;
            else
            config.AntennaRFPropHeight = 0;    
            end
                        
            config.TowerLat            = app.TowerLat;
            config.TowerLon            = app.TowerLon;
            
            % materials
            config.mat.steel.phiC = app.CompressionRF;
            config.mat.steel.phiT = app.TensionRF;
            config.mat.steel.sigY.strut = app.MastYield;  % grade 50
            config.mat.steel.sigY.guy = app.GuyYield;
            config.mat.steel.E = app.MastModulus;
            mu = 0.30;
            config.mat.steel.G = config.mat.steel.E/(2+mu);
            config.mat.steel.rho = app.MastDensity;
            config.mat.soil.bearingCapacity = app.SoilBearing;
            config.mat.soil.rho = app.SoilDensity;
            config.mat.soil.phi = app.SoilFriction;
            config.mat.soil.mu = 0.45;
            
            switch app.Units
                case 'Metric'
                    config.mat.soil.depth = 1;
                case 'English'
                    config.mat.soil.depth = 3.28084;
            end
           
            config.mat.concrete.rho = app.ConcreteDensity;
            
            % design
            config.designFactors.gustEff = app.GustEffectFactor;
            config.designFactors.windDirProb = app.WindDirectionProbability;
            
            switch app.StructureClass
                case 'I'
                     config.designFactors.IF = 0.87;
                case 'II'
                     config.designFactors.IF = 1;  
                case 'III'
                     config.designFactors.IF = 1.15;  
            end

            switch app.ExposureCategory
                case 'B'                   
            config.designFactors.Zg = 366;
            config.designFactors.alpha = 7;
                case 'C'                   
            config.designFactors.Zg = 274;
            config.designFactors.alpha = 9.5;            
                 case 'D'                   
            config.designFactors.Zg = 213;
            config.designFactors.alpha = 11.5;  
            end
            
            config.designFactors.windLoad = app.WindLoadMultiplier;
            config.designFactors.deadLoad = app.DeadLoadMultiplier;
            
            
            % environment 
            switch app.Units
                case 'Metric'
                    config.aero.rho = 1.225;
                    config.env.g    = 9.81;
                case 'English'
                    config.aero.rho =  0.0023769;
                    config.env.g    = 32.152231;
            end
            config.aero.windSpeed = app.WindSpeed;
            config.aero.windAngle = 0;                             
         
            % cost
            config.units.currency = 'USD';
            config.cost.mast = app.CostSteel;
            config.cost.guys = app.CostSteel;
            config.cost.concrete = app.CostConcrete;
            config.cost.erection = app.CostErection;
            switch app.Units
                case 'Metric'
                    config.cost.lightingHeight = 91;
                case 'English'
                    config.cost.lightingHeight = 300;
            end

            config.cost.lighting = app.CostLighting;
            config.cost.painting = app.CostPainting;
            config.cost.deadEnds = app.CostDeadends;
            config.cost.turnBuckles = app.CostTurnbuckles;
            config.cost.anchorRods = app.CostAnchorRods;
            config.cost.optionals = app.CostOptionals;
            config.cost.fencing = app.CostFencing;
            config.cost.preFabBuilding = app.CostPreFab;
            config.cost.siteExtras = app.CostSiteExtras;
            config.cost.logistics = app.CostLogistics;
            config.cost.installation = app.CostInstallation;
            config.cost.profit = app.CostProfit;
            config.cost.taxes  = app.CostTaxes;
            config.cost.contigency  = app.CostContigency;
            
end
      
function [config, node, elem] = analyze(config)
            
            disp('Creating nodes...');
            [config, node]       = createNodes (config);
            disp('Creating elements...');
            [config, node, elem] = createElements (config, node);
            disp('Assigning properites...');
            [config, node, elem] = assignProps (config, node, elem);
            disp('Assigning loads...');
            [config, node, elem] = assignLoads (config, node, elem);
            disp('Solving...');
            solve (config, node, elem);
            disp('Post-processing...');
            [config, node, elem] = postProcess (config, node, elem);
            [config, node, elem]  = sizeFoundation (config, node, elem);
            [config, node, elem]  = sizeAnchors (config, node, elem);
            [config, node, elem]  = cost (config, node, elem);
            
            
end

function [config, node, elem] = design(config)
            
            app.Status.Text = '';
            switch config.type
                case 'monopole'
                    config.mast.nGuyLevels = 0;
                    config.mast.horBrace = false;
                    config.mast.fwdBrace = false;
                    config.mast.bwdBrace = false;
                    [config, node, elem] = monopoleTowerDesign (config, app);

                case 'ss'
                    config.mast.angle = 45;
                    config.mast.nGuyLevels = 0;
                    config.mast.leg.TR = config.mast.leg.TRMax;
                    config.mast.web.TR = config.mast.web.TRMax;
                    [config, node, elem] = ssTowerDesign (config, app);

                case 'guyed'
                    config.mast.siteRadius = 0.8*config.mast.height;
                    config.mast.angle = 45;
                    config.mast.taper = 1.0;
                    config.mast.leg.TR = config.mast.leg.TRMax;
                    config.mast.web.TR = config.mast.web.TRMax;
                    [config, node, elem] = guyedTowerDesign (config, app);
            end
                      
end
       


