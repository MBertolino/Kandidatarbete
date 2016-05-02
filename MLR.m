clear;
%{
y = b1x1 + b2x2 + ... + bnxn + e

y - Dependent assets in which to choose a long or short position
    Price change from today to predTime days ahead

x - Assumed independent assets used to make a prediction
    Price change from day t-1 to t in the range of lag(end) days back
%}
tic;

% Load Data
load('KexJobbData.mat')

% Assets
depAsset = 36:40;
indepAsset = 36:40;
Ld = length(depAsset);
Li = length(indepAsset);

% Prediction Param
lag = 1:21;                       % How many days ago we look at the indep assets
predTime = 21;                    % How many days to predict
trainTime = 300;
lambda = 2e2;

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
yPred = zeros(tradePeriods, 2*Ld);
sigmay = yVal;
b = zeros(lag(end)*Li + 1, 2*Ld);
profitTot = zeros(tradePeriods, 2);
profit = zeros(tradePeriods, 2*Ld);
datez = profitTot;
profit(1,:) = bankStart;


%% Regression
% Initial step
for i = 1:trainTime - predTime - lag(end)
    yTrain(i,:) = clPr(i + lag(end) + predTime, depAsset) ...
        - clPr(i + lag(end), depAsset);
    xTemp = diffClPr(i + lag, indepAsset)';
    xTrain(i,:) = reshape(xTemp.', 1, []);
end

% Sliding window
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
    XTrain = [ones(size(xTrain(:,1))) xTrain];
    
    % For every invested market, calculate the regression coefficients
    % using both OLS and Ridge
    for m = 1:Ld
        % OLS Regression
        method{1} = 'OLS Regression';
        b(:,m) = regress(yTrain(:,m), XTrain);
        
        % Ridge Regression
        method{2} = 'Ridge Regression';
%         b(:,Ld+m) = RidgeRegress(yTrain(:,m), XTrain, lambda);
        b(:,Ld+m) = ridge(yTrain(:,m), XTrain, lambda);
    end
    
    
    %% Prediction & Validation
    % Prediction
    xTemp = diffClPr(i + lag + j*predTime, indepAsset)';
    xVal = reshape(xTemp.', 1, []);
    xVal = (xVal - mux)./sigmax;
    XVal = [ones(size(xVal(:,1))) xVal];
    yPred(j,:) = XVal*b;
    
    % Validation
    yVal(j,:) = clPr(i + lag(end) + (j+1)*predTime, depAsset) ...
        - clPr(i + lag(end) + j*predTime, depAsset);
    yVal(j,:) = (yVal(j,:) - muy)./sigmay;
    
    % Dates adjustment
    datez(j,:) = dates(i + lag(end) + (j+1)*predTime);
end


%% Strategy
for i = 1:2
    gamma(:,1+(i-1)*Ld:i*Ld) = sign(yPred(:,1+(i-1)*Ld:i*Ld));
    ret(:,1+(i-1)*Ld:i*Ld) = bsxfun (@rdivide, yVal.*gamma(:,1+(i-1)*Ld:i*Ld), std(yVal));
    infoRet(i) = mean(ret(:,1+(i-1)*Ld:i*Ld))/std(ret(:,1+(i-1)*Ld:i*Ld)) ...
        * sqrt(250); % Annualized
end

% Calculate the profit for using the stategy
for ii = 2:length(ret)
    profit(ii,:) = profit(ii - 1, :).*(1 + risk*ret(ii - 1, :));
end
profitTot(:,1) = sum(profit(:, 1:Ld), 2)/Ld;
profitTot(:,2) = sum(profit(:, Ld+1:2*Ld), 2)/Ld;


%% Plots
% Plot accumulated profit for each market
% figure()
% hold on;
% plot(datez(:,1), profit)
% title('Accumulated Profit')
% ylabel('Risk-adjusted profit')
% xlabel('Time [years]')
% legend('OLS', 'Ridge')
% datetick('x')
%
% % Plot accumulated total profit
% figure()
% plot(datez(:,1), profitTot)
% title('Accumulated Total Profit')
% ylabel('Risk-adjusted profit')
% xlabel('Time [years]')
% legend('OLS', 'Ridge')
% datetick('x')
%
% %Plot expected returns for each market
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

% Plot the investment
figure()
plot(datez, profitTot)
ylabel('Profit [$$$]') % ;)
xlabel('Time [Days]')
title('Profit using MLR in dollars')
legend('OLS',['Ridge, lambda = ' num2str(lambda)])
datetick('x')

toc;