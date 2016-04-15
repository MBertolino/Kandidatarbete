clear;

% Load Data
load('KexJobbData.mat')
depMarket = 1;                 % Dependent Market
indepMarket = 1:40;               % Possible independent markets
predTime = 1;                  % How many days to predict
lag = 1;                       % How many days ago we look at the indep markets


% Use this time period
[dates, clPr] = removeNaN(dates, closingPrice);

%% Regression

% Initialize
% yTrain = zeros(length(lag)+predTime,1);
% yHat = yTrain;
% yVal = [yTrain; zeros(predTime,1)];
% yPred = yVal;


trend = zeros(size(dates));
ret = trend;

%for i = 1+lag(end)+predTime:length(dates)-predTime
for i = 1+lag(end)+predTime:1+lag(end)+predTime+500
    yTrain(i-lag(end)-predTime,:) = clPr(i,depMarket) - clPr(i-predTime,depMarket);
    xTrain(i-lag(end)-predTime,:) = clPr(i-predTime,indepMarket) - clPr(i-lag-predTime,indepMarket);
    
    
    %     %% Prediction
    yVal(i+predTime) = clPr(i+predTime,depMarket) - clPr(i,depMarket);
    xVal(i+predTime) = clPr(i,indepMarket) - clPr(i-lag,indepMarket);
    %
    %
    %     % Positioning every predTime'th day
    %     if ~mod(i,predTime)
    %         trend(i) = sign(yPred(i+predTime));
    %     end
    
    
end


% Standardize data and add intercept
[xTrain, mu, sigma] = zscore(xTrain);
XTrain = [ones(size(xTrain(:,1))) xTrain];
xVal = (xVal - mu')./sigma';
XVal = [ones(size(xVal(:,1))) xVal];

% Normal regress
method{1} = 'Regress';
[b1, yHat1] = NormalRegress(yTrain, XTrain);

% % Lsq Ridge
% method{2} = 'Identity Ridge';
% [b2, yHat2] = RidgeRegress(yTrain, XTrain);


% Prediction
yPred = XVal*b1;
% yPred(:,2) = XVal*b2;

% % Calculate profit
% gamma = zeros(size(trend));
% gamma(trend > 0) = 1;
% gamma(trend < 0) = -1;
% ret = yVal.*gamma;
% profit = cumsum(ret);

% % Make it same length as dates
% yTrain = [yTrain; NaN*ones(predTime,1)];
% yHat = [yHat; NaN*ones(predTime,1)];

%% Plots
figure()
for ip = 1:1
    %p(ip) = subplot(2,1,ip);
    hold on;
    plot(dates, [yTrain yHat(:,ip)])
    ylabel('$$$');
    xlabel('Time [Days]');
    title(['Multiple Linear Regression Fit, ' method(ip)]);
    legend('Real data','Prediction', 'Location', 'Northeast');
    %linkaxes(p,'x','y');
    datetick('x')
    
end


%{
figure()
for ip = 1:1
    %p2(ip) = subplot(2,1,ip);
    hold on;
    plot(dates, [yVal yPred(:,ip)])
    ylabel('$$$');
    xlabel('Time [Days]');
    title(['Multiple Linear Regression Fit, ' method(ip)]);
    legend('Real data','Prediction', 'Location', 'Northeast');
    %linkaxes(p2,'x','y');
    datetick('x')
end


% Clear variables
clear closingPrice begT endT indS indE nIndep p p2

%}