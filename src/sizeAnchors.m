%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config, node, elem] = sizeAnchors(config, node, elem)

if config.mast.nGuyLevels > 0
    % first-order anchor sizing method
    for i = 1:length(config.anchor)-1
        % anchor reactions
        tol = 1e-3;
        anchorNode = find(strcmp({node(:).l},'anchor') & ...
            abs([node(:).r]-config.anchor(i).radius)<tol);
        F =  [node(anchorNode).reactLoad];
        Fx = (abs(F(1,:)));
        Fy = (abs(F(2,:)));
        Fz = (abs(F(3,:)));
        
        Tv = Fz;
        Th = rms([Fx, Fy]);
        
        % add initial tension
        gTv = []; gTh = [];
        for j = 1:length(anchorNode)
            guyElem = ismember([elem(:).node1], anchorNode(j)) |...
                ismember([elem(:).node2], anchorNode(j));
            T = [elem(guyElem).IT];
            zGuy = max([node([elem(guyElem).node1]).z]',...
                [node([elem(guyElem).node2]).z]');
            rGuy =  max([node([elem(guyElem).node1]).r]',...
                [node([elem(guyElem).node2]).r]');
            theta = atan(zGuy./rGuy);
            gTv(j) =  T*sin(theta);
            gTh(j) =  T*cos(theta);
        end
        
        Tv = max(Tv + gTv);
        Th = max(Th + gTh);
        
        phi = config.mat.soil.phi;
        hSoil = config.mat.soil.depth;
        gammaS = config.mat.soil.rho*config.env.g;
        gammaC = config.mat.concrete.rho*config.env.g;
        Kp = (tand(45+phi/2))^2;
        df =  2*hSoil*gammaS/gammaC;
        lf = 3*Th/(Kp*df^2*(gammaS + gammaC));
        bf = Tv/(df*lf*gammaC);
        
        config.anchor(i).mass = df*lf*bf*config.mat.concrete.rho;
        config.anchor(i).dimensions = [lf bf df]';
        config.anchor(i).designLoad = [Tv Th 0];
        
    end
end
