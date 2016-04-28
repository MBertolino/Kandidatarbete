clear;
tic;

% Parameters
depMarket = 36:37;
indepMarket = 36:40;
Ld = length(depMarket);
Li = length(indepMarket);
lag = 1:20;                       % How many days ago we look at the indep markets
lambda = 2e2;

predTime = 21;                   % How many days to predict
trainTime = 300;

% Load Data
load('KexJobbData.mat')

% Remove NaN's
% Start at 02-Jan-2009
% End at 28-Jan-2016
[dates, clPr] = removeNaN(dates(7447-trainTime-2*predTime+1:end), ...
    closingPrice(7447-trainTime-2*predTime+1:end,:));
tradePeriods = floor((length(dates) - trainTime - predTime)/predTime)-1;


%% Regression
% Pre-allocating
yTrain = zeros(trainTime-lag(end)-predTime,Ld);
xTrain = zeros(trainTime-lag(end)-predTime, lag(end)*Li);
yVal = zeros(tradePeriods, Ld);
yPred = zeros(tradePeriods, 2*Ld);
sigmay = yVal;
datez = yVal;
b1 = zeros(lag(end)*Li+1,Ld);
b2 = b1;

% Initial step
for i = 1:trainTime-predTime-lag(end)
    yTrain(i,:) = clPr(i+lag(end)+2*predTime,depMarket) - clPr(i+lag(end)+1*predTime,depMarket);
    xTemp = repmat(clPr(i+lag(end)+1*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+1*predTime,indepMarket);
    xTrain(i,:) = reshape(xTemp.',1,[]);
end

% Sliding window
for j = 1:tradePeriods
    yTrain(1:end-predTime,:) = yTrain(predTime+1:end,:);
    xTrain(1:end-predTime,:) = xTrain(predTime+1:end,:);
    for i = predTime:trainTime-lag(end)
        yTrain(i,:) = clPr(i+lag(end)+(j+1)*predTime,depMarket) - clPr(i+lag(end)+j*predTime,depMarket);
        xTemp = repmat(clPr(i+lag(end)+j*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+j*predTime,indepMarket);
        xTrain(i,:) = reshape(xTemp.',1,[]);
    end
    
    % Standardize data and add intercept
    [xTrain, mux, sigmax] = zscore(xTrain);
    XTrain = [ones(size(xTrain(:,1))) xTrain];
    [yTrain, muy, sigmay] = zscore(yTrain);
    
    for m = 1:Ld
        % OLS Regression
        method{1} = 'OLS Regression';
        b1(:,m) = regress(yTrain(:,m), XTrain);
        
        % Ridge Regression
        method{2} = 'Ridge Regression';
        b2(:,m) = RidgeRegress(yTrain(:,m), XTrain, lambda);
        %         b2(:,m) = ridge(yTrain(:,m), XTrain, lambda);
    end
    
    
    %% Prediction & Validation
    % Prediction
    xTemp = repmat(clPr(i+lag(end)+(j+1)*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+(j+1)*predTime,indepMarket);
    xVal = reshape(xTemp.',1,[]);
    xVal = (xVal - mux)./sigmax;
    XVal = [ones(size(xVal(:,1))) xVal];
    
    yPred(j,1:Ld) = XVal*b1;                            % Regress
    yPred(j,Ld+1:2*Ld) = XVal*b2;                       % Ridge
    
    % Validation
    yVal(j,:) = clPr(i+lag(end)+(j+2)*predTime,depMarket) - clPr(i+lag(end)+(j+1)*predTime,depMarket);
    yVal(j,:) = (yVal(j,:) - muy)./sigmay;
    
    % Dates adjustment
    datez(j) = dates(i+lag(end)+(j+1)*predTime);
end

% Strategy
for i = 1:Ld
    gamma(:,1+(i-1)*Ld:i*Ld) = sign(yPred(:,1+(i-1)*Ld:i*Ld));
    ret(:,1+(i-1)*Ld:i*Ld) = yVal.*gamma(:,1+(i-1)*Ld:i*Ld);
    profit(:,1+(i-1)*Ld:i*Ld) = cumsum(ret(:,1+(i-1)*Ld:i*Ld));
    profitTot(:,i) = sum(profit(:,1+(i-1)*Ld:i*Ld),2);
    infoRet(i) = mean(ret(:,1+(i-1)*Ld:i*Ld))/std(ret(:,1+(i-1)*Ld:i*Ld)) ...
        * sqrt(250); % Annualized
end

%% Plots

% Plot accumulated profit for each market
figure()
hold on;
plot(datez(:,1), profit)
title('Accumulated Profit')
ylabel('Risk-adjusted profit')
xlabel('Time [years]')
legend('OLS', 'Ridge')
datetick('x')


% Plot accumulated total profit
figure()
plot(datez(:,1), profitTot)
title('Accumulated Total Profit')
ylabel('Risk-adjusted profit')
xlabel('Time [years]')
legend('OLS', 'Ridge')
datetick('x')


%Plot expected returns for each market
figure()
hold on;
plot(datez(:,1),yVal)
hold on;
plot(datez(:,1),yPred)
title('Expected return')
ylabel('Risk-adjusted return')
xlabel('Time [years]')
legend('Real data','Prediction', 'Location', 'Northeast');
datetick('x')
%}
toc;