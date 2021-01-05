%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.


function [config, node, elem] = cost (config, node, elem)

%% Materials

% mast
idx = strcmp({elem(:).type},{'leg'}) | strcmp({elem(:).type},{'web'});
config.cost.breakdown.mast = sum([elem(idx).A].*[elem(idx).L].*...
    [elem(idx).rho]/config.designFactors.deadLoad)*config.cost.mast;

% guys
if config.mast.nGuyLevels>0
idx = strcmp({elem(:).type},{'guy'});
config.cost.breakdown.guy = sum([elem(idx).A].*[elem(idx).L].*...
    [elem(idx).rho])*config.cost.guys;
end

% foundation/anchors
concreteMass = config.anchor(end).mass;
if strcmp(config.type,'guyed')
concreteMass = concreteMass + sum([config.anchor(1:end-1).mass])*length(config.mast.Az);
end
config.cost.breakdown.foundation = concreteMass*config.cost.concrete;


%% Tower Accessories

if strcmp(config.type,'guyed')
config.cost.breakdown.deadends =  ...
    sum([config.guy.D]*config.cost.deadEnds.*...
    [config.guy.nPerAz]*length(config.mast.Az));
config.cost.breakdown.turnBuckles =  sum([config.guy.D]*...
    config.cost.turnBuckles.*[config.guy.nPerAz]*length(config.mast.Az));
config.cost.breakdown.anchorRods = (length(config.anchor)-1)*...
    length(config.mast.Az)*config.cost.anchorRods;
end

nLights = floor(config.mast.height/config.cost.lightingHeight);
config.cost.breakdown.lighting = nLights*config.cost.lighting;
config.cost.breakdown.painting = config.cost.painting*config.mast.height;
config.cost.breakdown.optionals = config.cost.optionals*config.mast.height;

%% Construction

config.cost.breakdown.erection = config.cost.erection*config.mast.height;
config.cost.breakdown.fencing = config.cost.fencing;
config.cost.breakdown.preFabBuilding = config.cost.preFabBuilding;
config.cost.breakdown.siteExtras = config.cost.siteExtras;

%% Contigencies

towerCost = sum(struct2array(config.cost.breakdown));
totalCost = towerCost*(1 + (config.cost.logistics + config.cost.installation + ...
    config.cost.profit + config.cost.contigency + config.cost.taxes)/100);
config.cost.breakdown.logistics = totalCost*config.cost.logistics/100;
config.cost.breakdown.installation = totalCost*config.cost.installation/100;
config.cost.breakdown.profit = totalCost*config.cost.profit/100;
config.cost.breakdown.taxes = totalCost*config.cost.taxes/100;
config.cost.breakdown.contigency = totalCost*config.cost.contigency/100;
%% Total

config.cost.total = totalCost;
