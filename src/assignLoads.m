%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config, node, elem] = assignLoads (config, node, elem)

% loading
Gh  = config.designFactors.gustEff;
Kd  = config.designFactors.windDirProb;
I   = config.designFactors.IF;
airDensity = config.aero.rho;
q   = 0.5*airDensity*(config.aero.windSpeed)^2*Gh*Kd*I*...
    config.designFactors.windLoad; 

% wind loading
for i = 1:length(elem)
    
    P1 = [[node([elem(i).node1]).x]' [node([elem(i).node1]).y]' ...
        [node([elem(i).node1]).z]'];
    P2 = [[node([elem(i).node2]).x]' [node([elem(i).node2]).y]' ...
        [node([elem(i).node2]).z]'];
    C  = (P2 - P1)./norm(P2-P1);
    EPA = elem(i).EPA;
    Cx = C(:,1);Cy = C(:,2); Cz = C(:,3);
    
    % coordinate transformation: global to local (copied from frame3dd)
    t1=0;t2=0;t3=0;t4=0;t5=0;t6=0;t7=0;t8=0;t9=0;
    
    if abs(Cz) == 1.0
        t3 =  Cz;
        t5 =  1;
        t7 = -Cz;
    else
        den = sqrt ( 1.0 - Cz^2);
        t1 = Cx;
        t2 = Cy;
        t3 = Cz;
        t4 = -Cy/den;
        t5 = Cx/den;
        t7 = -Cx*Cz/den;
        t8 = -Cy*Cz/den;
        t9 = den;
    end
    
    Q = [[t1;t2;t3],[t4;t5;t6],[t7;t8;t9]];    
    h = (node([elem(i).node1]).z + node([elem(i).node2]).z)/2;
    Kz  = 2.01*(h/config.designFactors.Zg)^(2/config.designFactors.alpha);
    F = q*Kz.*EPA;
    alpha = config.aero.windAngle; 
    F1 = F*cosd(alpha);
    F2 = F*sind(alpha);
    F3 = 0;
    QF = Q'*[F1; F2; F3];
    theta = acosd(QF(1)/rms([QF(2),QF(3)]));
    elem(i).appLoad(1,:)   = 0;
    elem(i).appLoad(2:3,:) = [QF(2)*sind(theta)^2; QF(3)*sind(theta)^2];
end

% antenna loading 
if isfield(config,'antenna')
    for i = 1:length(config.antenna)
        nIdx = config.antenna(i).node;
        h = config.antenna(i).height;
        Kz  = 2.01*(h/config.designFactors.Zg)^(2/config.designFactors.alpha);
        F = q*Kz*config.antenna(i).EPA;
        alpha = config.aero.windAngle;
        node(nIdx).appLoad(1,:) = F*cosd(alpha);
        node(nIdx).appLoad(2,:) = F*sind(alpha);
        node(nIdx).appLoad(3,:) = ...
            -config.antenna(i).weight*(config.env.g)*ones(1,length(alpha))*...
            config.designFactors.deadLoad;
    end
end

% guy initial tension loading 
if strcmp(config.type, 'guyed')
    elemType = {elem(:).type};
    idx = find(strcmp(elemType,'guy')>0);
    for i = 1:length(idx)
        node1 = elem(idx(i)).node1;
        node2 = elem(idx(i)).node2;
        if node(node2).z > node(node1).z
            nIdx = node2;
        else
            nIdx = node1;
        end
        alpha = config.aero.windAngle; 
        Tw    = elem(idx(i)).IT*elem(idx(i)).Lv/elem(idx(i)).L;
        node(nIdx).appLoad(3,:) = node(nIdx).appLoad(3,:) + ...
            -Tw*ones(1,length(alpha));
    end    
end
