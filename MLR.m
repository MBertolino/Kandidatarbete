clear;

% Load Data
load('KexJobbData.mat')
depMarket = 1;                 % Dependent Market
indepMarket = 1:40;            % Possible independent markets
predTime = 20;                 % How many days to predict
lag = 1;                       % How many days ago we look at the indep markets


% Specify Time Period
startTrain = '01-Jun-2006';
endTrain = '01-Jul-2008';

% Use this time period
[dates, closingPrice] = removeNaN(dates, closingPrice);
indS = find(dates > datenum(startTrain, 'dd-mmm-yyyy'), 1);
indE = find(dates > datenum(endTrain, 'dd-mmm-yyyy'), 1);
datesTrain = dates(indS:indE);
datesPred = dates(indE+1:indE+predTime);
clPr = closingPrice(indS:indE+predTime, [depMarket indepMarket]);

%% Regression
yTrain = clPr(1+lag:end-predTime,1);
xTrain = clPr(1:end-predTime-lag,:);
%yTrain = diff(yTrain);
%xTrain = diff(xTrain);
XTrain = [ones(size(xTrain(:,1))) xTrain];

% Standardize data
[xTrain, mu, sigma] = zscore(xTrain);

% Normal regress
method{1} = 'Regress';
[b1, yHat(:,1)] = NormalRegress(yTrain, XTrain);

% Lsq Ridge
method{2} = 'Identity Ridge';
[b2, yHat(:,2)] = RidgeRegress(yTrain, XTrain);


%% Prediction
yVal = clPr(end-predTime+1:end,1);
xVal = clPr(end-predTime+1-lag:end-lag,:);
%yVal = diff(yVal);
%xVal = diff(xVal);

% Standardize data
%xVal = zscore(xVal); % ANVÄNDA GAMLA MU OCH SIGMA
%xVal = (xVal - mu')./sigma';

% Prediction
XVal = [ones(size(xVal(:,1))) xVal];
yPred(:,1) = XVal*b1;
yPred(:,2) = XVal*b2;


%% Plots
figure()
for ip = 1:2
    p(ip) = subplot(2,1,ip);
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
for ip = 1:2
    p2(ip) = subplot(2,1,ip);
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