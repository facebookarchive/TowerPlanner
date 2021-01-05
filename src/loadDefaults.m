%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function app = loadDefaults ()

% tower
app.Units = 'Metric';
app.TowerType = 'Guyed';
app.CrossSection = 'triangle';
app.Model = 'Beam';
app.BracingType = 'X';
app.TowerHeight = 40;

% design
app.WebAngle = 45;
app.MastWidth = 1.5;
app.LegTaper = 1;
app.nGuyLevels = 5;
app.SiteRadius = 0.8*app.TowerHeight;
app.LegDia = 0.12;
app.WebDia = 0.04;
app.GuyDia = 0.02;
app.LegTRDesign = 12;
app.WebTRDesign = 2;
app.GuyIT = 0.1;

% materials
app.LegTRMax = 15;
app.LegSRMax = 150;
app.WebTRMax = 15;
app.WebSRMax = 200;
app.CompressionRF = 0.85;
app.TensionRF = 0.8;
app.MastYield = 5e8;
app.GuyYield = 1e9;
app.MastModulus = 2e11;
app.MastDensity = 7800;
app.SoilBearing = 1.4e5;
app.SoilDensity = 1730;
app.SoilFriction = 30;
app.ConcreteDensity = 2400;

% loading data
app.WindSpeed = 50;
app.AntennaTable(1,:)= [1,0,100,2,0,100]; %{'EPA'; 'Mass'; 'Height (%)'; 'Freq. (GHz)'; 'Diameter'; 'Power'};
app.AntennaRFPropIdx = 1;
app.TowerLat = 37.452961;
app.TowerLon = -122.181725;
app.GustEffectFactor = 0.85;
app.WindDirectionProbability = 0.85;
app.StructureClass = 'II';
app.ExposureCategory = 'C';
app.WindLoadMultiplier = 1.6;
app.DeadLoadMultiplier = 1.2;

% cost data
app.CostSteel = 5;
app.CostConcrete = 0.05;
app.CostErection = 200;
app.CostLighting = 200;
app.CostPainting = 25;
app.CostDeadends = 5000;
app.CostTurnbuckles = 1000;
app.CostAnchorRods = 5000;
app.CostOptionals = 10;
app.CostFencing = 25000;
app.CostPreFab = 75000;
app.CostLogistics = 10;
app.CostInstallation = 10;
app.CostProfit = 5;
app.CostSiteExtras = 20000;
app.CostTaxes = 10;
app.CostContigency = 5;

% environment
setenv('/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin');
 
end
