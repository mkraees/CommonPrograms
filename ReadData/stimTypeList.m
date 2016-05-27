function [timeStartFromBaseLineList,deltaTList] = stimTypeList

timeStartFromBaseLineList(1) = -0.55; deltaTList(1) = 1.024;  % stimOn = 0.2; stimOff = 0.3; Used for RF mapping protocols. 
timeStartFromBaseLineList(2) = -1.148; deltaTList(2) = 2.048; % stimOn = 0.4; stimOff = 0.6;
timeStartFromBaseLineList(3) = -1.5; deltaTList(3) = 4.096;   % stimOn = 1.5, stimOff = 1.5;
timeStartFromBaseLineList(4) = -0.848; deltaTList(4) = 2.048; % stimOn = 0.8; stimOff = 0.7;

timeStartFromBaseLineList(5) = -1; deltaTList(5) = 3.2768; % stimOn = 0.8; stimOff = 0.7;  For BP data with Fs = 2500: Compatible with matching pursuit
timeStartFromBaseLineList(6) = -0.5; deltaTList(6) = 2.048; % stimOn = 0.8; stimOff = 0.7; For BR data (Fs = 2000) and EG Data (Fs = 1000): Compatible with matching pursuit