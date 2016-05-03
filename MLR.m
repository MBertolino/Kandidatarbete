clear;
%{
y = b1x1 + b2x2 + ... + bnxn + e

y - Dependent assets in which to choose a long or short position
    Price change from today to predTime days ahead

x - Assumed independent assets used to make a prediction
    Price change from day t-1 to t in the range of lag(end) days back
%}
tic;


%% Setup
% Load Data
load('KexJobbData.mat')

% Assets
class = 1;
name = assetClassNames(class);
depAsset = find(assetClass == class); % 1:13 14:22 23:35 36:40
indepAsset = depAsset;
Ld = length(depAsset);
Li = length(indepAsset);

% Prediction Param
lag = 1:21;                       % How many days ago we look at the indep assets
predTime = 21;                    % How many days to predict
trainTime = 650;
lambda = [0 1e2 1e3 1e4];
Ll = length(lambda);

% Investment Param
bankStart = 10000;
risk = 0.05;

% Remove NaN's
% Start at 02-Jan-2009
% End at 28-Jan-2016
[dates, clPr] = removeNaN(dates(7447 - trainTime - 2*predTime + 1:end), ...
    closingPrice(7447 - trainTime - 2*predTime + 1:end, :));
tradePeriods = floor((length(dates) - trainTime - predTime)/predTime) - 1;
diffClPr = diff(clPr);

% Pre-allocating
yTrain = zeros(trainTime - lag(end) - predTime, Ld);
xTrain = zeros(trainTime - lag(end) - predTime, lag(end)*Li);
yVal = zeros(tradePeriods, Ld);
yPred = zeros(tradePeriods, Ll*Ld);
sigmay = yVal;
holding = zeros(tradePeriods, Ll*Ld);
holdingTot = zeros(tradePeriods, Ll);
datez = holdingTot;
holding(1,:) = bankStart;


%% Regression
% Initial step
% Data relating assets to be predicted, y, and assets to use as
% predictors, x. For each day in a tradePeriods sized window y and x are
% related.
% yTrain - is the training set of the to-be-predicted y
% xTrain - is the training set of the predictors x
for i = 1:trainTime - predTime - lag(end)
    yTrain(i,:) = clPr(i + lag(end) + predTime, depAsset) ...
        - clPr(i + lag(end), depAsset);
    xTemp = diffClPr(i + lag, indepAsset)';
    xTrain(i,:) = reshape(xTemp.', 1, []);
end

% Sliding window
% After a prediction the training window is moved so only the latest data
% points are used as training data as they are assumed to be more accurate
for j = 1:tradePeriods
    yTrain(1:end - predTime, :) = yTrain(predTime + 1:end, :);
    xTrain(1:end - predTime, :) = xTrain(predTime + 1:end, :);
    for i = predTime:trainTime - lag(end)
        yTrain(i,:) = clPr(i + lag(end) + j*predTime, depAsset) ...
            - clPr(i + lag(end) + (j-1)*predTime, depAsset);
        xTemp = diffClPr(i + lag, indepAsset)';
        xTrain(i,:) = reshape(xTemp.', 1, []);
    end
    
    % Standardize data and add intercept
    [yTrain, muy, sigmay] = zscore(yTrain);
    [xTrain, mux, sigmax] = zscore(xTrain);
    
    % For every invested asset, calculate the regression coefficients
    % using both OLS and Ridge
    b(:,1:Ll*Ld) = RidgeRegress(yTrain, xTrain, lambda);
    
    
    %% Prediction & Validation
    % Prediction of the change in price of each asset
    % XVal - are the predictors
    % yPred - is the predicted change in price of each asset
    xTemp = diffClPr(i + lag + j*predTime, indepAsset)';
    xVal = reshape(xTemp.', 1, []);
    xVal = (xVal - mux)./sigmax;
    yPred(j,:) = xVal*b;
    
    % Validation
    % yVal - is the actual price change measured at the end of the
    % prediction time
    yVal(j,:) = clPr(i + lag(end) + (j+1)*predTime, depAsset) ...
        - clPr(i + lag(end) + j*predTime, depAsset);
    yVal(j,:) = (yVal(j,:) - muy)./sigmay;
    
    % Dates adjustment
    % At each predicted day, the date is extracted
    datez(j,:) = dates(i + lag(end) + (j+1)*predTime);
end


%% Strategy
% gamma - is the position (1 or -1) for each asset
% ret - is the risk adjusted return for taking a position
% infoRet - is the information quotient for a strategy
for i = 1:Ll
    gamma(:,1+(i-1)*Ld:i*Ld) = sign(yPred(:,1+(i-1)*Ld:i*Ld));
    ret(:,1+(i-1)*Ld:i*Ld) = bsxfun (@rdivide, yVal.*gamma(:,1+(i-1)*Ld:i*Ld), std(yVal));
    infoRet(i) = mean(ret(:,1+(i-1)*Ld:i*Ld))/std(ret(:,1+(i-1)*Ld:i*Ld)) ...
        * sqrt(250); % Annualized
end
infoRet

% Calculate the evolution of a holding for each asset
for ii = 2:length(ret)
    holding(ii,:) = holding(ii - 1, :).*(1 + risk*ret(ii - 1, :));
end

% Calculate the evolution of the total holding
for il = 1:Ll
    holdingTot(:,il) = sum(holding(:, (il-1)*Ld + 1: il*Ld), 2)/Ld;
end

%% Plots
%Plot expected returns for each asset
% figure()
% hold on;
% plot(datez(:,1), yVal)
% hold on;
% plot(datez(:,1), yPred)
% title('Expected return')
% ylabel('Risk-adjusted return')
% xlabel('Time [years]')
% legend('Real data','Prediction', 'Location', 'Northeast');
% datetick('x')

% Plot the evolution of the total holding
figure()
plot(datez, holdingTot)
ylabel('Holding [$]') % ;)
xlabel('Time [Days]')
title(['Holding using MLR on ' name])
str = cellstr(num2str(lambda', 'lambda = %d'));
legend(str, 'Location', 'NorthWest');
datetick('x')

hold off;
toc;