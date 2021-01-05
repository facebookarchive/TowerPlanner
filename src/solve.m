%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function solve (config, node, elem)
    
% nodes
nNodes = length(node);
XYZ = [[node(:).x];[node(:).y];[node(:).z];zeros(1,nNodes)];

% element connectivity
ELT = [[elem(:).node1];[elem(:).node2]];

% reactions
RCT = zeros(6, nNodes);	
if strcmp(config.type,'guyed')
idx = find(strcmp({node(:).l},{'base'}));    
RCT(1:6,idx) = repmat([1;1;1;0;0;0],1,length(idx));
idx = find(strcmp({node(:).l},{'anchor'}));    
RCT(1:6,idx) = repmat([1;1;1;1;1;1],1,length(idx));
else
    idx = find(strcmp([{node(:).l}],{'base'}));  
    RCT(1:6,idx) = repmat([1;1;1;1;1;1],1,length(idx));
end


% element properties 
nElem       = length(elem);
EAIJ(1,:)   = [elem(:).A];
EAIJ(2,:)   = [elem(:).As];
EAIJ(3,:)   = [elem(:).As];
EAIJ(4,:)   = [elem(:).Jxx];
EAIJ(5,:)   = [elem(:).Iyy];
EAIJ(6,:)   = [elem(:).Izz];
EAIJ(7,:)   = [elem(:).E];
EAIJ(8,:)   = [elem(:).G];
EAIJ(9,:)   = zeros(1,nElem);
EAIJ(10,:)  = [elem(:).rho];

% point loads
P        = zeros(6, nNodes);
P(:,:)   = [node(:).appLoad];

% distributed loads
U        = [elem(:).appLoad];

% prescribed displacements
D     = zeros(6, nNodes);	

% run frame 3dd
frame_3dd(XYZ,ELT,RCT,EAIJ,P,U,D, config.IOfilename);


end
