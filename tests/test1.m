%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

%% Example run
% !! Run this from the root folder !!  
function results = test1()
% Default Values 
app = loadDefaults();

% Modify Defaults
app.TowerType = 'Self-Supported';
app.Model = 'lattice';
app.WindSpeed = 30;
app.TowerHeight = 40;
app.AntennaTable(1,:)= [1,10,100,2,0,100];

% Execute

mode = 'analyze';
p = mfilename('fullpath');
pp = pwd;
filepath = fileparts(p);
cd(filepath);cd('../src');
results = execute(app, mode);
cd(pp);
end

