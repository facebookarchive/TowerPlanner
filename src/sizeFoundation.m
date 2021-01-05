%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config, node, elem] = sizeFoundation (config, node, elem)

baseNode = find(strcmp({node(:).l},'base'));
F =  [node(baseNode).reactLoad];
lateralLoad    = max(max(abs(F(1:2,:))));
columnLoad     = max(abs(F(3,:)));
overTurnMoment = max(abs(F(6,:)));

% first-order foundation sizing

% soil bearing capacity
padArea   = columnLoad/config.mat.soil.bearingCapacity;
padWidthB = sqrt(padArea);

% padwidth is atleast pole attachment size
if strcmp(config.mast.model,'lattice')
    padWidthP = config.mast.leg.D;
else
    padWidthP = max(config.mast.width);
end

% check for overturning
if overTurnMoment > 0
padWidthO = max(real(roots([1.5 -3*overTurnMoment/columnLoad ...
    -2*columnLoad/config.mat.soil.bearingCapacity])));
else
padWidthO = 0;    
end

% check for sliding
mu = config.mat.soil.mu;
phi = config.mat.soil.phi;
Kp = (tand(45+phi/2))^2;
gammaS = config.mat.soil.rho*config.env.g;
frictionResistance = mu*columnLoad;
F = lateralLoad - frictionResistance;
if F>0
padWidthS = ...
    max(real(roots([1 -2*padWidthP padWidthP^2 -8*(F)/(Kp*gammaS)])));
else
padWidthS = 0;
end

% check for shear
padWidth = max([padWidthB, padWidthP, padWidthO, padWidthS]);
padDepth = (padWidth - padWidthP)*0.5;
padArea  = padWidth^2;

k = length(config.anchor);
config.anchor(k+1).mass = padArea*padDepth*config.mat.concrete.rho*length(baseNode);
config.anchor(k+1).dimensions = [padWidth padWidth padDepth];
config.anchor(k+1).designLoad = [columnLoad lateralLoad overTurnMoment];
config.anchor(k+1).radius = 0;
config.anchor(k+1).nGuysPerAnchor = 0;
