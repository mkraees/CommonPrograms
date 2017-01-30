% make colour patch
% Vinay Shirhatti, 31 May 2016
% 
% Input arguments for the patch:
% pass a structure variable stim which is a gaborStimulus structure 
% OR optionally with the following fields: 
% (a,e) => centre of the patch 
% r => radius (2 dimensional if both inner and outer radii are passed)
% (h,s,v) : (hue, saturation, value) for the patch
%
%**************************************************************************

function [colourPatch] = makeColourPatchHSV(hplot,stim,passhsv,showStim)

if ~exist('showStim','var')
    showStim = 0;
end

if ~exist('passhsv','var')
    passhsv = 0;
end

if passhsv
    if isfield(stim,'a')
        azi = stim.a;
    else
        azi = 0;
    end
    if isfield(stim,'e')
        ele = stim.e;
    else
        ele = 0;
    end
    ori = stim.h;
    sf = stim.s;
    c = stim.v;
    
    if length(stim.r)==1    % Gabor
        radMax = stim.r;
        radMin = 0;
    else
        radMax = stim.r(2);
        radMin = stim.r(1);
    end
else
    azi = stim.azimuthDeg;
    ele = stim.elevationDeg;
    sf  = stim.spatialFreqCPD;
    ori = stim.orientationDeg;
    c   = stim.contrastPC;
    
    if length(stim.radiusDeg)==1    % Gabor
        radMax = stim.radiusDeg;
        radMin = 0;
    else
        radMax = stim.radiusDeg(2);
        radMin = stim.radiusDeg(1);
    end
end

if radMax==radMin
    radMax=0; radMin=0;
end

numTheta = 60;
theta = 0:(2*pi/numTheta):(2*pi-2*pi/numTheta);
dtheta = (theta(2)-theta(1))/2;
theta = theta-repmat(dtheta,1,length(theta)); % so that theta is at the centre of the sector

innercircleX = radMin.*cos(theta)+azi;
innercircleY = radMin.*sin(theta)+ele;

outercircleX = radMax.*cos(theta)+azi;
outercircleY = radMax.*sin(theta)+ele;

x = [innercircleX innercircleX(1) outercircleX outercircleX(1) innercircleX(1)];
y = [innercircleY innercircleY(1) outercircleY outercircleY(1) innercircleY(1)];

hue = ori/360; sat = sf; val = c/100; % hsvcolor = [hue sat val]
rgbcolor = hsv2rgb([hue sat val]);

colourPatch = patch('XData',x,'YData',y,'FaceColor',rgbcolor,...
            'linestyle','none','parent',hplot);
        
if showStim
    figure;
    patch('XData',x,'YData',y,'FaceColor',rgbcolor,...
    'linestyle','none','parent',hplot);
end


end