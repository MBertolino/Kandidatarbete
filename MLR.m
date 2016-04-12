clear;

% Load Data
load('KexJobbData.mat')
depMarket = 1;                 % Dependent Market
indepMarket = 2:40;            % Possible independent markets
nPred = 3;                     % How many predictors    
predTime = 40;                 % How many days to predict
lag = 1;                       % How many days ago we look at the indep markets


% Specify Time Period
startTrain = '01-May-2004';          % Start
endTrain = '01-Jul-2008';            % End

% Use this time period
[dates, closingPrice] = removeNaN(dates, closingPrice);
indS = find(dates > datenum(startTrain, 'dd-mmm-yyyy'), 1);
indE = find(dates > datenum(endTrain, 'dd-mmm-yyyy'), 1);
datesTrain = dates(indS:indE);
datesPred = dates(indE+1:indE+predTime);
clPr = closingPrice(indS:indE+predTime, [depMarket indepMarket]);

%% Regression
yTrain = clPr(1+lag:end-predTime,1);
xTrain = clPr(1:end-predTime-lag,2:end);
%yTrain = diff(yTrain);
%xTrain = diff(xTrain);

% Standardize data
xTrain = zscore(xTrain);

% Normal regress
method{1} = 'Regress';
[b1, yHat(:,1), pred1] = NormalRegress(yTrain, xTrain, nPred);

% Lsq Ridge
method{2} = 'Identity Ridge';
[b2, yHat(:,2), pred2] = RidgeRegress(yTrain, xTrain, nPred);

% Quadratic solution
method{3} = 'Generalized Ridge';
[b3, yHat(:,3), pred3] = GeneralizedRegress(yTrain, xTrain, nPred);



%% Prediction
yVal = clPr(end-predTime+1:end,1);
xVal = clPr(end-predTime+1-lag:end-lag,pred1);
%yVal = diff(yVal);
%xVal = diff(xVal);

% Standardize data
xVal = zscore(xVal);

% Prediction
XVal = [ones(size(xVal(:,1))) xVal];
yPred(:,1) = XVal*b1;
yPred(:,2) = XVal*b2;
yPred(:,3) = XVal*b3;



%% Plots
figure()
for ip = 1:3
    p(ip) = subplot(3,1,ip);
    hold on;
    plot(datesTrain(1+lag:end), [yTrain yHat(:,ip)])
    ylabel('$$$');
    xlabel('Time [Days]');
    title(['Multiple Linear Regression Fit, ' method(ip)]);
    legend('Real data','Prediction', 'Location', 'Northeast');
    linkaxes(p,'x','y');
    datetick('x')
    
end
figure()
for ip = 1:3
    p2(ip) = subplot(3,1,ip);
    hold on;
    plot(datesPred(1:end), [yVal yPred(:,ip)])
    ylabel('$$$');
    xlabel('Time [Days]');
    title(['Multiple Linear Regression Fit, ' method(ip)]);
    legend('Real data','Prediction', 'Location', 'Northeast');
    linkaxes(p2,'x','y');
    datetick('x')
end


% Clear variables
clear closingPrice begT endT indS indE nIndep p p2

%}