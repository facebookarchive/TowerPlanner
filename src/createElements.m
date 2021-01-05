%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.


function [config, node, elem] = createElements (config, node)

elem = [];

% leg elements
if strcmp(config.mast.model,'lattice')
    
    Az = config.mast.Az;
    tol = 1e-3;
    
    for i = 1:length(Az)
        
        lIdx = find(strcmp({node(:).l},'mast') & abs([node.t] - Az(i))<tol);
        [~, sIdx] = sort([node(lIdx).z]);
        nodeIDs = lIdx(sIdx);
        nodePairs = [nodeIDs(1:end-1)' nodeIDs(2:end)'];
        elemType  = cell(1,size(nodePairs,1));
        elemType(:) = {'leg'};
        elem = [elem; struct('node1',num2cell(nodePairs(:,1)),...
            'node2',num2cell(nodePairs(:,2)),'type', elemType')];
    end
    
elseif strcmp(config.mast.model,'beam')
    
    lIdx = find(strcmp({node(:).l},'mast'));
    [~, sIdx] = sort([node(lIdx).z]);
    nodeIDs = lIdx(sIdx);
    nodePairs = [nodeIDs(1:end-1)' nodeIDs(2:end)'];
    elemType  = cell(1,size(nodePairs,1));
    elemType(:) = {'leg'};
    elem = [elem; struct('node1',num2cell(nodePairs(:,1)),...
        'node2',num2cell(nodePairs(:,2)),'type', elemType')];
    
end

% base elements
if config.mast.nGuyLevels>0 && strcmp(config.mast.model,'lattice')
    
    tol = 1e-3;
    mNodes = node(strcmp({node(:).l},'mast'));
    minMastHeight = min([mNodes(:).z]);
    aIdx = find(abs([node(:).z] - minMastHeight)<tol);
    bIdx = find(strcmp({node(:).l},'base'));    
    nodePairs = [aIdx' repmat(bIdx,length(aIdx),1)];
    elemType  = cell(1,size(nodePairs,1));
    elemType(:) = {'leg'};
    elem = [elem; struct('node1',num2cell(nodePairs(:,1)),...
        'node2',num2cell(nodePairs(:,2)),'type', elemType')];
else
    
    tol = 1e-3;
    mNodes = node(strcmp({node(:).l},'mast'));
    minMastHeight = min([mNodes(:).z]);
    aIdx = find(abs([node(:).z] - minMastHeight)<tol);
    [~,sIdx] = sort([node(aIdx).t]);
    aIdx     = aIdx(sIdx);
    bIdx     = find(strcmp({node(:).l},'base'));
    [~,sIdx] = sort([node(bIdx).t]);
    bIdx     = bIdx(sIdx);   
    
    nodePairs = [aIdx' bIdx'];
    elemType  = cell(1,size(nodePairs,1));
    elemType(:) = {'leg'};
    elem = [elem; struct('node1',num2cell(nodePairs(:,1)),...
        'node2',num2cell(nodePairs(:,2)),'type', elemType')];
        
end

% guy elements
if config.mast.nGuyLevels>0
    Az = config.mast.Az;
    tol = 1e-3;
    for i = 1:length(Az)
        for j = 1:length([config.guy.level])
            
            if strcmp(config.mast.model,'lattice')
                mastNodeAz = find(strcmp({node(:).l},'mast') &...
                    abs([node.t] - Az(i))<tol);
            else
                mastNodeAz = find(strcmp({node(:).l},'mast'));
            end
            [~, idx]   = min(abs([node(mastNodeAz).z] - ...
                config.guy(j).level));
            mastNode   = mastNodeAz(idx);
                       
            anchorRadius = config.anchor(config.guy(j).anchorIdx).radius;
            anchorNode = find([node.z] <tol &...
                abs([node.r] - anchorRadius)<tol &  abs([node.t] - ...
                Az(i))<tol);
            
            elem = [elem; struct('node1',mastNode,...
                'node2',anchorNode,'type', 'guy')];
            
        end
    end
end

% bracing elements
if  strcmp(config.mast.model,'lattice')
    Rm = config.mast.R;     % seems like unused
    Az = [config.mast.Az' [config.mast.Az(2:end) config.mast.Az(1)]'];
    tol = 1e-3;
    
    for i = 1:length(Az)
        
        if strcmp(config.type, 'guyed')
            aIdx = find((strcmp({node(:).l},'mast')) &...
                abs([node.t] - Az(i,1))<tol);
            [~, sIdx] = sort([node(aIdx).z]);
            aIDs = aIdx(sIdx);
            
            bIdx = find((strcmp({node(:).l},'mast')) &...
                abs([node.t] - Az(i,2))<tol);
            [~, sIdx] = sort([node(bIdx).z]);
            bIDs = bIdx(sIdx);
        else
            aIdx = find((strcmp({node(:).l},'mast') | strcmp({node(:).l},'base')) &...
                abs([node.t] - Az(i,1))<tol);
            [~, sIdx] = sort([node(aIdx).z]);
            aIDs = aIdx(sIdx);
            
            bIdx = find((strcmp({node(:).l},'mast') | strcmp({node(:).l},'base')) &...
                abs([node.t] - Az(i,2))<tol);
            [~, sIdx] = sort([node(bIdx).z]);
            bIDs = bIdx(sIdx);
        end
        
        if config.mast.horBrace
            nodePairs = [aIDs' bIDs'];
            elemType  = cell(1,size(nodePairs,1));
            elemType(:) = {'web'};
            elem = [elem; struct('node1',num2cell(nodePairs(:,1)),...
                'node2',num2cell(nodePairs(:,2)),'type', elemType')];
        end
        
        if config.mast.fwdBrace
            nodePairs = [aIDs(1:end-1)' bIDs(2:end)'];
            elemType  = cell(1,size(nodePairs,1));
            elemType(:) = {'web'};
            elem = [elem; struct('node1',num2cell(nodePairs(:,1)),...
                'node2',num2cell(nodePairs(:,2)),'type', elemType')];
        end
        
        if config.mast.bwdBrace
            nodePairs = [aIDs(2:end)' bIDs(1:end-1)'];
            elemType  = cell(1,size(nodePairs,1));
            elemType(:) = {'web'};
            elem = [elem; struct('node1',num2cell(nodePairs(:,1)),...
                'node2',num2cell(nodePairs(:,2)),'type', elemType')];
        end
        
    end
end
