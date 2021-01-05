%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.


function [config, node, elem] = assignProps (config, node, elem)

% leg
elemType = {elem(:).type};
idx = strcmp(elemType,'leg');
TR = config.mast.leg.TR;
P1 = [[node([elem(idx).node1]).x]' [node([elem(idx).node1]).y]' ...
    [node([elem(idx).node1]).z]'];
P2 = [[node([elem(idx).node2]).x]' [node([elem(idx).node2]).y]' ...
    [node([elem(idx).node2]).z]'];
L  = vecnorm(P1-P2,2,2);
zLoc = (P1(:,3) + P2(:,3))/2;

if strcmp(config.mast.model,'lattice')
    
    E = config.mat.steel.E;
    G = config.mat.steel.G;
    D  = (zLoc/config.mast.height*(config.mast.leg.taper - 1) + 1)*config.mast.leg.D;
    t  = D./TR;
    A  = pi*D.^2/4 - pi*(D - 2*t).^2/4;
    As = A*1e-3; 
    Jxx = pi/2*((D/2).^4 + (D/2 - t).^4); 
    Iyy = pi/4*((D/2).^4 + (D/2 - t).^4);
    Izz = pi/4*((D/2).^4 + (D/2 - t).^4); 
    rho = config.mat.steel.rho*config.designFactors.deadLoad;
    CD  = 1.2;
    EPA = CD*D;
    
elseif strcmp(config.mast.model,'beam')
    
    L                 = mean(L);
    a                 = config.mast.width;
    dLeg              = config.mast.leg.D;
    tLeg              = config.mast.leg.D/config.mast.leg.TR;
    dWeb              = config.mast.web.D;
    tWeb              = config.mast.web.D/config.mast.web.TR;
    aLeg              = pi*dLeg.^2/4 - pi*(dLeg - 2*tLeg).^2/4;
    aWeb              = pi*dWeb.^2/4 - pi*(dWeb - 2*tWeb).^2/4;
    rho               = config.mat.steel.rho*config.designFactors.deadLoad;
    theta             = 90 - config.mast.angle ;
    webDiagL          = a/cosd(config.mast.angle);
    nPanels           = L./(a*tand(config.mast.angle));
    E                 = config.mat.steel.E;
    G                 = config.mat.steel.G;
    CD                = 1.2;
   
    if ~strcmp(config.type, 'monopole')
        if (config.mast.bwdBrace || config.mast.fwdBrace) && config.mast.horBrace
            EA        = 3*E*aLeg;
            EI        = E*aLeg*a.^2/2;
            GJ        = 1/((4/(E*a.^2))*(1/(aLeg*(tand(theta))^2) + ...
                1/(aWeb*(sind(theta))^2*cosd(theta)) + tand(theta)/aWeb));
            GA        = 1/(2/(3*E*aWeb*(sind(theta))^2*cosd(theta)) + 2*tand(theta)/(3*E*aWeb));
            EPA       = CD*((3*sqrt(3)/4*tand(theta))*dWeb + ...
                ((cosd(theta))^2 + ...
                2/cosd(theta)*(1 - 0.25*(sind(theta))^2)^(3/2))*dWeb + 3*dLeg);
            nWebDiag  = nPanels; nWebHorz = nPanels;
            
            bracingMass   = nWebDiag*(aWeb*webDiagL)*rho  + nWebHorz*(aWeb.*a)*rho;
            legMass       = aLeg*rho;
            massPerLength = (legMass + bracingMass/L)*3;
            
        elseif config.mast.fwdBrace && config.mast.bwdBrace && ...
                ~config.mast.horBrace
            EA        = 3*E*aLeg;
            EI        = E*aLeg*a.^2/2;
            GJ        = E*a.^2*aWeb*(sind(theta))^2*cosd(theta)/2;
            GA        = 3*E*aWeb*(sind(theta))^2*cosd(theta);
            EPA       = CD*(2*((cosd(theta))^2 + ...
                2/cosd(theta)*(1 - 0.25*(sind(theta))^2)^(3/2))*dWeb + 3*dLeg);
            nWebDiag  = nPanels*2; nWebHorz = 0;
            bracingMass   = nWebDiag*(aWeb*webDiagL)*rho  + nWebHorz*(aWeb.*a)*rho;
            legMass       = aLeg*rho;
            massPerLength = (legMass + bracingMass/L)*3;
            
        elseif config.mast.fwdBrace && config.mast.bwdBrace &&...
                config.mast.horBrace
            EA        = 3*E*(aLeg + ...
                2*aWeb*aWeb*(cosd(theta)^3)/(aWeb + 2*aWeb*(sind(theta)^3)));
            EI        = E*a.^2/2*(aLeg + ...
                aWeb*aWeb/2*(cosd(theta))^3/(aWeb + 2*aWeb*(sind(theta)^3)));
            GA        = 3*E*aWeb*(sind(theta))^2*cosd(theta);
            GJ        = E*a.^2*aWeb*(sind(theta))^2*cosd(theta)/2;
            EPA       = CD*((3*sqrt(3)/4*tand(theta))*dWeb + ...
                2*((cosd(theta))^2 + ...
                2/cosd(theta)*(1 - 0.25*(sind(theta))^2)^(3/2))*dWeb + 3*dLeg);
            nWebDiag  = nPanels*2; nWebHorz = nPanels;
            bracingMass   = nWebDiag*(aWeb*webDiagL)*rho  + nWebHorz*(aWeb.*a)*rho;
            legMass       = aLeg*rho;
            massPerLength = (legMass + bracingMass/L)*3;
        end
        
    else
        EA        = E*aLeg;
        EI        = E*pi/4*((dLeg/2).^4 + (dLeg/2 - tLeg).^4);
        GJ        = G*pi/2*((dLeg/2).^4 + (dLeg/2 - tLeg).^4);
        EPA       = CD*dLeg;
        massPerLength = aLeg*rho;
    end
        
    A             = EA/E;
    As            = GA/G;
    D             = nan;
    Jxx           = GJ/G;
    Izz           = EI/E;
    Iyy           = EI/E;
    rho           = massPerLength./A;
    
end

L = num2cell(L);
E = num2cell(E);
G = num2cell(G);
A = num2cell(A);
As = num2cell(As);
D = num2cell(D);
Jxx = num2cell(Jxx);
Iyy = num2cell(Iyy);
Izz = num2cell(Izz);
rho = num2cell(rho);
EPA = num2cell(EPA);


[elem(idx).L]    = deal(L{:});
[elem(idx).E]    = deal(E{:});
[elem(idx).G]    = deal(G{:});
[elem(idx).A]    = deal(A{:});
[elem(idx).As]   = deal(As{:});
[elem(idx).D]    = deal(D{:});
[elem(idx).Jxx]  = deal(Jxx{:});
[elem(idx).Iyy]  = deal(Iyy{:});
[elem(idx).Izz]  = deal(Izz{:});
[elem(idx).rho]  = deal(rho{:});
[elem(idx).EPA]  = deal(EPA{:});


% web
if strcmp(config.mast.model,'lattice')
    elemType = {elem(:).type};
    idx = strcmp(elemType,'web');
    D  = config.mast.web.D*ones(length(find(idx>0)),1);
    TR = config.mast.web.TR;
    E = config.mat.steel.E;
    G = config.mat.steel.G;
    rho = config.mat.steel.rho*config.designFactors.deadLoad;
    P1 = [[node([elem(idx).node1]).x]' [node([elem(idx).node1]).y]'...
        [node([elem(idx).node1]).z]'];
    P2 = [[node([elem(idx).node2]).x]' [node([elem(idx).node2]).y]' ...
        [node([elem(idx).node2]).z]'];
    L  = vecnorm(P1-P2,2,2);
    t  = D/TR;
    
    D((P1(:,3) - P2(:,3)==0)) = D((P1(:,3) - P2(:,3)==0));          
    
    A  = pi*D.^2/4 - pi*(D - 2*t).^2/4;
    As = A*1e-3; 
    Jxx = pi/2*((D/2).^4 + (D/2 - t).^4); 
    Iyy = pi/4*((D/2).^4 + (D/2 - t).^4);
    Izz = pi/4*((D/2).^4 + (D/2 - t).^4); 
    CD  = 1.2;
    EPA = CD*D;
    
    L = num2cell(L);
    A = num2cell(A);
    As = num2cell(As);   
    D = num2cell(D);
    Jxx = num2cell(Jxx);
    Iyy = num2cell(Iyy);
    Izz = num2cell(Izz);
    EPA = num2cell(EPA);
    
    [elem(idx).L]    = deal(L{:});
    [elem(idx).D]    = deal(D{:});
    [elem(idx).E]    = deal(E);
    [elem(idx).G]    = deal(G);
    [elem(idx).A]    = deal(A{:});
    [elem(idx).As]   = deal(As{:});
    [elem(idx).Jxx]  = deal(Jxx{:});
    [elem(idx).Iyy]  = deal(Iyy{:});
    [elem(idx).Izz]  = deal(Izz{:});
    [elem(idx).rho]  = deal(rho);
    [elem(idx).EPA]  = deal(EPA{:});
end

% guy
if config.mast.nGuyLevels > 0
    elemType = {elem(:).type};
    idx = strcmp(elemType,'guy');
    D  = config.guy(1).D;
    Em = config.mat.steel.E;
    G = 1e-6;
    rho = config.mat.steel.rho*config.designFactors.deadLoad;
    P1 = [[node([elem(idx).node1]).x]' [node([elem(idx).node1]).y]' ...
        [node([elem(idx).node1]).z]'];
    P2 = [[node([elem(idx).node2]).x]' [node([elem(idx).node2]).y]' ...
        [node([elem(idx).node2]).z]'];
    b  = max(P1(:,3),P2(:,3));
    c  = vecnorm(P1-P2,2,2);
    a  = sqrt(c.^2 - b.^2);
    for i = 1:length(b)
        [~, gIdx(i)] = min(abs([config.guy(:).level] - b(i)));
    end
    N  = [config.guy(gIdx).nPerAz];   
    A  = pi*D^2/4;
    As = A*1e-3;
    w  = A*rho*config.env.g/config.designFactors.deadLoad;
    sigY = config.mat.steel.sigY.guy;
    IT = [config.guy(gIdx).initialTension]*A*sigY;
    EAeq  = (Em*A./(1 + (w^2*a.^2*A*Em./(12*(IT(:).^3)))));   % ignore initial tension effect on stiffness, since it's effect is ~10^3 smaller.
    Jxx = 1e-5;
    Iyy = 1e-5;
    Izz = 1e-5;
    CD  = 1.2;
    EPA = CD*D*N;
    
    L   = num2cell(c);
    Lv  = num2cell(b);
    Lh  = num2cell(a);
    IT  = num2cell(IT(:).*N(:));
    E   = num2cell(EAeq(:)/A);
    A   = num2cell(A*N);
    N   = num2cell(N);
    As  = num2cell(As);  
    EPA  = num2cell(EPA); 
    
    [elem(idx).L]    = deal(L{:});
    [elem(idx).E]    = deal(E{:});
    [elem(idx).G]    = deal(G);
    [elem(idx).A]    = deal(A{:});
    [elem(idx).As]   = deal(As{:});
    [elem(idx).D]    = deal(D);
    [elem(idx).Jxx]  = deal(Jxx);
    [elem(idx).Iyy]  = deal(Iyy);
    [elem(idx).Izz]  = deal(Izz);
    [elem(idx).rho]  = deal(rho);
    [elem(idx).EPA]  = deal(EPA{:});
    [elem(idx).IT]   = deal(IT{:});
    [elem(idx).Lv]   = deal(Lv{:});
    [elem(idx).Lh]   = deal(Lh{:});
    [elem(idx).Lh]   = deal(Lh{:});
    [elem(idx).N]    = deal(N{:});
end

