%Copyright (c) Facebook, Inc. and its affiliates.

%This source code is licensed under the MIT license found in the
%LICENSE file in the root directory of this source tree.

function [config, node, elem] = postProcess (config, node, elem)

% read data
fid = fopen(['IO' config.IOfilename '.out']);
lCnt = 0;
tline = fgetl(fid);
while lCnt<1e3
    if strcmp(tline,'N O D E   D I S P L A C E M E N T S  					(global)')
        break;
    end
    tline = fgetl(fid);
    lCnt = lCnt + 1;
end

C2  = textscan(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f','HeaderLines',1);
C3  = textscan(fid,'%f\t%f\t%f%c\t%f\t%f\t%f\t%f\t%f','HeaderLines',2);
C4  = textscan(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f','HeaderLines',2);
C5  = textscan(fid,'%f\t%c%c%c\t%f\t%f\t%f\t%f\t%f\t%f','HeaderLines',4);
fclose(fid);

if length(C2{1}) < 2 | any(isnan(cell2mat(C4)))
    config.opt.fail = true;
    return;
else
    config.opt.fail = false;
end

% nodal displacements
C2 = cell2mat(C2);
for i = 1:size(C2,1)
    n = C2(i,1);
    node(n).disp = C2(i,2:7)';
end

% nodal reaction loads
C4 = cell2mat(C4);
for i = 1:size(C4,1)
    n = C4(i,1);
    node(n).reactLoad = C4(i,2:7)';
end


% element internal loads
axMax = [C5{5}(1:2:end)];
axMin = [C5{5}(2:2:end)];
axMode = [C3{4}(1:2:end)];
Vy   = [C5{6}(1:2:end)];
Vz   = [C5{7}(1:2:end)];
Txx  = [C5{8}(1:2:end)];
Myy  = [C5{9}(1:2:end)];
Mzz  = [C5{10}(1:2:end)];

shear    = abs([[C5{6}(1:2:end)] [C5{6}(2:2:end)] [C5{7}(1:2:end)] [C5{7}(2:2:end)]]);
torsion    = abs([[C5{8}(1:2:end)] [C5{8}(2:2:end)]]);
moment    = abs([[C5{9}(1:2:end)] [C5{9}(2:2:end)] [C5{10}(1:2:end)] [C5{10}(2:2:end)]]);

for i = 1:length(elem)
    elem(i).intLoad.axMax   = axMax(i);
    elem(i).intLoad.axMin   = axMin(i);
    elem(i).intLoad.shear   = max(shear(i,:));
    elem(i).intLoad.torsion = max(torsion(i,:));
    elem(i).intLoad.moment  = max(moment(i,:));
    elem(i).intLoad.axMode  = axMode(i);
end

% evaluate failure
for i = 1:length(elem)
    
    if strcmp(config.mast.model,'lattice') &&...
            ~strcmp(elem(i).type,'guy')
        
        D = elem(i).D;
        L = elem(i).L;
        TR = config.mast.([elem(i).type]).TR;
        t  = D/TR;
        A  = pi*D^2/4 - pi*(D - 2*t)^2/4;
        
        [cMax, tMax] = ...
            memberStrength (D, t, config.mat.steel.phiC, config.mat.steel.phiT, ...
            config.mat.steel.E, L, config.mat.steel.sigY.strut);
        axMode       = elem(i).intLoad.axMode;
        [~, idx] = max([elem(i).intLoad.axMax/tMax, -elem(i).intLoad.axMin/cMax]);
        
        if axMode == 't'
            elem(i).stressCrit = abs(elem(i).intLoad.axMax)/tMax;
            elem(i).stressMaxAllow = tMax/A;
            elem(i).stressCritMode = 'Tension';
        else   %if axMode == 'c'
            elem(i).stressCrit = abs(elem(i).intLoad.axMin)/cMax;
            elem(i).stressMaxAllow = cMax/A;
            elem(i).stressCritMode = 'Compression';
        end
        
        
    elseif strcmp(config.mast.model,'beam') &&...
            ~strcmp(elem(i).type,'guy')
        
        a    = config.mast.width;
        dLeg = config.mast.leg.D;
        if length(dLeg) > 1, dLeg = dLeg(i); end
        tLeg = dLeg/config.mast.leg.TR;
        aLeg = pi*dLeg^2/4 - pi*(dLeg - 2*tLeg)^2/4;
        lLeg = a*tand(config.mast.angle);
        dWeb = config.mast.web.D;
        tWeb = config.mast.web.D/config.mast.web.TR;
        aWeb = pi*dWeb^2/4 - pi*(dWeb - 2*tWeb)^2/4;
        lWeb = a/cosd(config.mast.angle);
        if config.mast.fwdBrace && config.mast.bwdBrace
            lWeb = lWeb*0.5;
        end
        
        cMaxLeg = ...
            memberStrength (dLeg, tLeg, config.mat.steel.phiC, config.mat.steel.phiT, ...
            config.mat.steel.E, lLeg, config.mat.steel.sigY.strut);
        
        cMaxWeb = ...
            memberStrength (dWeb, tWeb, config.mat.steel.phiC, config.mat.steel.phiT, ...
            config.mat.steel.E, lWeb, config.mat.steel.sigY.strut);
        
        
       F     = elem(i).intLoad;
        if strcmp(config.type, 'guyed')
            
            % leg stress
            axMax = (abs(F.axMin)/3 + abs(F.moment)*2/(sqrt(3)*a));
            lfIdx = axMax/cMaxLeg;
            
            % web stress           
            Fmax = 0;
            for alpha = 0 %[0 90 270]
                Fz        = F.shear*sind(alpha);
                Fy        = -F.shear*cosd(alpha);
                F1        = 2*F.torsion/(3^0.5*a) + 2/3*Fz;
                F2        = 2*F.torsion/(3^0.5*a) - 1/3*Fz - Fy/(2*sind(60));
                F3        = 2*F.torsion/(3^0.5*a) - 1/3*Fz + Fy/(2*sind(60));
                Fmax      = max([Fmax abs([F1 F2 F3])]);
            end
            
            rLeg          = sqrt(dLeg^2 + (dLeg - 2*tLeg)^2)/4;
            legSR         = lLeg/rLeg;
            Fs            = axMax;
            Ps            = (1.5 + (legSR-60)/60)*Fs/100;
            Ps            = max(min(Ps, 2.5*Fs/100), 1.5*Fs/100);
            
            if config.mast.fwdBrace && config.mast.bwdBrace
                Ps = Ps/2;
                Fmax = Fmax/2;
            end
            
            wfIdx  = max([Fmax, Ps])/cMaxWeb/cosd(config.mast.angle);
            elem(i).stressLegCrit = lfIdx;
            elem(i).stressLegMaxAllow = cMaxLeg/aLeg;
            elem(i).stressWebCrit = wfIdx;
            elem(i).stressWebMaxAllow = cMaxWeb/aWeb;
            
            if lfIdx > wfIdx
                elem(i).stressCrit = lfIdx;
                elem(i).stressMaxAllow = cMaxLeg/aLeg;
                elem(i).stressCritMode = 'Leg-Compression';
            else
                elem(i).stressCrit = wfIdx;
                elem(i).stressMaxAllow = cMaxWeb/aWeb;
                elem(i).stressCritMode = 'Web-Compression';
            end
            
        elseif strcmp(config.type, 'monopole') || strcmp(config.type, 'ss')
            
            S = pi/32*(dLeg^4 - (dLeg - 2*tLeg)^4)/dLeg;
            axMax = abs(F.axMin) + abs(F.moment)/S;
            lfIdx = axMax/cMaxLeg;
            elem(i).stressLegCrit = lfIdx;
            elem(i).stressLegMaxAllow = cMaxLeg/aLeg;
            elem(i).stressCrit = lfIdx;
            elem(i).stressMaxAllow = cMaxLeg/aLeg;
            elem(i).stressCritMode = 'Compresssion';
            
        end
        
    end
        
    if strcmp(elem(i).type,'guy')
        sigY = config.mat.steel.sigY.guy;
        phiTG = config.mat.steel.phiT*0.75;
        A    = elem(i).A;
        IT   = elem(i).IT;
        axMode = elem(i).intLoad.axMode;
        
        if axMode == 't' 
            T = IT + elem(i).intLoad.axMax;
        elseif axMode == 'c' 
            T = IT + elem(i).intLoad.axMin;
        end
        
        elem(i).stressCrit     = T/(sigY*phiTG*A);
        elem(i).stressMaxAllow = sigY*phiTG*A;  
        if T > 0
        elem(i).stressCritMode = 'Tension';
        else
        elem(i).stressCritMode = 'Slack';    
        end
    end
    
end

% clear data files
delete(['IO' config.IOfilename '*']);

end

function [compStrength, tensStrength] = ...
    memberStrength (D, t, phiT, phiC, E, L, sigY)

r    = sqrt(D^2 + (D-2*t)^2)/4;
A    = pi*D^2/4 - pi*(D - 2*t)^2/4;

% compressive strength
if D/t ~= 2
    if D/t < 0.114*E/sigY
        Fy = sigY;
    elseif D/t > 0.114*E/sigY && D/t < 0.448*E/sigY
        Fy = (0.0379*E/sigY/(D/t) + 2/3)*sigY;
    else
        Fy = 0.337*E/(D/t);
    end
else
    Fy = sigY;
end

lamda  = L/(r*pi)*sqrt(Fy/E);                     
if lamda <= 1.5
    Fcr = (0.658^(lamda^2))*Fy;
else
    Fcr = (0.877/(lamda^2))*Fy;
end
compStrength = A*Fcr*phiC;

% tensile strength
tensStrength = A*sigY*phiT;

end
