%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config, node] = createNodes (config)

% cross-section
if strcmp(config.mast.Xsection,'triangle')
    config.mast.Az = [0 120 240];
    config.mast.R  = config.mast.width/2/cosd(30);
elseif strcmp(config.mast.Xsection,'square')
    config.mast.Az = [45 135 225 315];
    config.mast.R  = config.mast.width/2/cosd(45);    
end

% associate guy heights to specific anchor. Anchor radius fixed at mean of
% guy heights, i.e., at 45 deg. 
if config.mast.nGuyLevels>0
    guyLevels = linspace(config.mast.height, 0, config.mast.nGuyLevels+1);
    guyLevels(end) = [];
    guyCnt = 0;
    ancIdx = 1;
    
    for i = 1:length(guyLevels)
        
        % determine no. of guys per Az depending on proximity to antenna
        config.guy(i).level = guyLevels(i);
        config.guy(i).nPerAz = 1;
        
        if isfield(config,'antenna')
            [val, idx] = min(abs(config.guy(i).level - [config.antenna.height]));
            if val < 0.05 && ~config.antenna(idx).pairedWithGuy
                config.guy(i).nPerAz = 2;
                config.antenna(idx).pairedWithGuy = true;
            end
        end
        
        guyCnt = guyCnt + config.guy(i).nPerAz;
        config.guy(i).anchorIdx = ancIdx;
        config.guy(i).D = config.guy(1).D;
        config.guy(i).initialTension = config.guy(1).initialTension;
        
        if guyCnt == config.anchor(1).nGuysPerAnchor
            guyIdx = find([config.guy(:).anchorIdx] == ancIdx);
            config.anchor(ancIdx).radius = mean([config.guy(guyIdx).level]);
            guyCnt = 0;
            ancIdx = ancIdx + 1;
        elseif guyCnt > config.anchor(1).nGuysPerAnchor
            guyIdx = find([config.guy(1:end-1).anchorIdx] == ancIdx);
            config.anchor(ancIdx).radius = mean([config.guy(guyIdx).level]);
            guyCnt = config.guy(i).nPerAz;
            ancIdx = ancIdx + 1;                                        
        elseif i == length(guyLevels) && ...
                length(guyLevels)-i < config.anchor(1).nGuysPerAnchor  
            guyIdx = find([config.guy(:).anchorIdx] == ancIdx);
            config.anchor(ancIdx).radius = mean([config.guy(guyIdx).level]);
        end
        
    end
    config.anchor(1).radius = config.mast.siteRadius;
end

% mast
X = []; Y = []; Z = [];R = [];T = []; L = [];
if strcmp(config.mast.model,'lattice')
    
    %Rm = config.mast.R;  % unused for now
    Az = config.mast.Az;
    
    for i = 1:length(Az)
        
        panelHeight = mean(config.mast.width)*tand(config.mast.angle);
        z  = config.mast.height:-panelHeight:0;
        z  = z(z>0);
        Rm = linspace(config.mast.R*config.mast.taper,...
            config.mast.R, length(z));
        x = Rm*cosd(Az(i));
        y = Rm*sind(Az(i));
        r = Rm;
        t = Az(i)*ones(1,length(z));
        l = cell(1, length(z)); l(:) = {'mast'};
        X = [X, x];
        Y = [Y, y];
        Z = [Z, z];
        R = [R, r];
        T = [T, t];
        L = [L, l];
    end
    
    if config.mast.nGuyLevels>0
        X = [X, 0]; Y = [Y, 0]; Z = [Z, 0];
        R = [R, 0]; T = [T, 0]; L = [L, {'base'}];
    else
        X  = [X, config.mast.R*cosd(Az)];
        Y  = [Y, config.mast.R*sind(Az)];
        Z  = [Z, zeros(1, length(Az))];
        R  = [R, config.mast.R*ones(1,length(Az))];
        T  = [T, Az];
        l  = cell(1, length(Az)); l(:) = {'base'};
        L  = [L, l];
    end
    
elseif strcmp(config.mast.model,'beam')
    
    if strcmp(config.type,'guyed')
        nBeamElem             = config.mast.nGuyLevels*10;
        if config.mast.nGuyLevels > 1   
            dec          = ...
                1/2^(max(round(log(nBeamElem/config.mast.nGuyLevels)/log(2)),0));
            nodeCoords   = ...
                interp1(1:config.mast.nGuyLevels+1,...
                [0 sort([config.guy(:).level])], 1:dec:config.mast.nGuyLevels+1);
        else
            nodeCoords   = ...
                linspace (0, config.mast.height,  nBeamElem + 1);  
        end
    else 
        nodeCoords    = ...
            linspace (0, config.mast.height,  config.mast.nBeamElem + 1);
    end
    
    z = nodeCoords;
    x = zeros(1,length(z));
    y = zeros(1,length(z));
    r = zeros(1,length(z));
    t = zeros(1,length(z));
    l = cell(1, length(z)); l(:) = {'mast'};l(1) = {'base'};
    X = [X, x];
    Y = [Y, y];
    Z = [Z, z];
    R = [R, r];
    T = [T, t];
    L = [L, l];
end
    
% anchor nodes
if config.mast.nGuyLevels>0
    
    Az = config.mast.Az;
    for i = 1:length(Az)       
        x = [config.anchor.radius]*cosd(Az(i));
        y = [config.anchor.radius]*sind(Az(i));
        z = [config.anchor.radius]*0;
        r = [config.anchor.radius];
        t = Az(i)*ones(1,length(z));
        l = cell(1, length(z)); l(:) = {'anchor'};
        X = [X, x];
        Y = [Y, y];
        Z = [Z, z];
        R = [R, r];
        T = [T, t];
        L = [L, l];
    end
    
end

node = struct('x',num2cell(X),'y',num2cell(Y),'z',...
    num2cell(Z),'r',num2cell(R),'t',num2cell(T),'l',L);

% assign no loads
[node(:).appLoad]    = deal([0;0;0;0;0;0]);

% associate antenna nodes
if isfield(config,'antenna')
    for i = 1:length(config.antenna)
        h = config.antenna(i).height;
        [~, hIdx] = min(abs([node(:).z] - h));
        config.antenna(i).node = hIdx;
        config.antenna(i).offset = node(config.antenna(i).node).r;
    end
end
