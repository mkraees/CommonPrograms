% Fit a curve to a set of data points using fminsearch
% Vinay Shirhatti, 27 April 2017
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [x,yfit,fval,exitflag,output] = fminsearchCurveFit(fun,startPt,xData,yData,evalnfun,opts)

if ~exist('opts','var')
    opts = optimset('TolX',1e-6,...
                    'TolFun',1e-6,...
                    'MaxIter',5000,...
                    'Display','off',...
                    'LargeScale','off',...
                    'MaxFunEvals',5000);
end

if ~exist('evalnfun','var'); evalnfun = 'mse'; end

switch evalnfun
    case 'mse'
        [x,fval,exitflag,output] = fminsearch(@(x) mseFunction2(x,yData,xData,fun),startPt,opts);
    case 'abs'
        [x,fval,exitflag,output] = fminsearch(@(x) absFunction(x,yData,xData,fun),startPt,opts);
    otherwise
        [x,fval,exitflag,output] = fminsearch(@(x) mseFunction2(x,yData,xData,fun),startPt,opts);
end

switch fun
    case 'linear'
        yfit = linearFunction(x,xData);
    case 'power'
        yfit = powerFunction(x,xData);
    case 'quad'
        yfit = parabolaFunction(x,xData);
    otherwise
        yfit = linearFunction(x,xData);
end
end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Evaluation function (mean square error)
function mse = mseFunction2(x,yData,xData,fun)    
switch fun
    case 'linear'
        yHat = linearFunction(x,xData);
    case 'power'
        yHat = powerFunction(x,xData);
    case 'quad'
        yHat = parabolaFunction(x,xData);
    otherwise
        yHat = linearFunction(x,xData);
end
mse = sum((yData-yHat).^2);
end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Evaluation function (mean square error)
function numpoints = absFunction(x,yData,xData,fun)    
switch fun
    case 'linear'
        yHat = linearFunction(x,xData);
    case 'power'
        yHat = powerFunction(x,xData);
    case 'quad'
        yHat = parabolaFunction(x,xData);
    otherwise
        yHat = linearFunction(x,xData);
end
% intersectPoints = length(intersect(yHat,yData));
% numpoints = length(xData)-intersectPoints;

numpoints = sum(abs(yData-yHat));

end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Curve fitting options
function y=powerFunction(x,xData)
    A=x(1);X=x(2);C=x(3);
    y=A*(xData).^(X) + max(C,0);
end

function y=linearFunction(x,xData)
    A=x(1);C=x(2);
    y=A*(xData) + C;
end

function y=parabolaFunction(x,xData)
    A=x(1);X=x(2);B=x(3); C=x(3);
    y=-A*(xData-X).^B + C;
end