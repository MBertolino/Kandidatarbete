clear;

% Load Data
load('KexJobbData.mat')
[dates, clPr] = removeNaN(dates, closingPrice);

% Parameters
depMarket = 1;
indepMarket = 1:40;
lag = 1:4;                       % How many days ago we look at the indep markets

trainTime = 400;
predTime = 21;                   % How many days to predict
tradePeriods = 160;


%% Regression

% Pre-allocating
yTrain = zeros(trainTime-lag(end)-predTime,1);
xTrain = zeros(trainTime-lag(end)-predTime, lag(end)*length(indepMarket));
yVal = zeros(tradePeriods,1);
yPred1 = yVal;
yPred2 = yVal;
datez = yVal;
b = zeros(lag(end)*length(indepMarket)+1,2);

for m = 1:length(depMarket)
    
    % Sliding window
    for j = 1:tradePeriods
        for i = 1:trainTime-predTime-lag(end)
            yTrain(i,:) = clPr(i+lag(end)+(j+1)*predTime,depMarket(m)) - clPr(i+lag(end)+j*predTime,depMarket(m));
            xTemp = repmat(clPr(i+lag(end)+j*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+j*predTime,indepMarket);
            xTrain(i,:) = reshape(xTemp.',1,[]);
        end
          
        % Standardize data and add intercept
        [xTrain, mux, sigmax] = zscore(xTrain);
        XTrain = [ones(size(xTrain(:,1))) xTrain];
        [yTrain, muy, sigmay] = zscore(yTrain);
        
        % OLS Regression
        method{1} = 'OLS Regression';
        b(:,1) = regress(yTrain,XTrain);
        
        % Ridge Regression
        method{2} = 'Ridge Regression';
        b(:,2) = RidgeRegress(yTrain, XTrain);
        
        
        %% Prediction & Validation
        % Prediction
        xTemp = repmat(clPr(i+lag(end)+(j+1)*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+(j+1)*predTime,indepMarket);
        xVal = reshape(xTemp.',1,[]);
        xVal = (xVal - mux)./sigmax;
        XVal = [ones(size(xVal(:,1))) xVal];
        
        yPred1(j) = XVal*b(:,1);             % Regress
        yPred2(j) = XVal*b(:,2);             % Ridge
        
        % Validation
        yVal(j) = clPr(i+lag(end)+(j+2)*predTime,depMarket(m)) - clPr(i+lag(end)+(j+1)*predTime,depMarket(m));
        yVal(j) = (yVal(j) - muy)/sigmay;
        
        % Dates adjustment
        datez(j) = dates(i+lag(end)+(j+1)*predTime);
    end
    
    % Calculate profit
    yPred = yPred2;
    gamma = sign(yPred);
    ret = yVal.*gamma;
    profit(:,m) = cumsum(ret);
end

profit = sum(profit,2);

%% Plots

figure()
plot(datez, profit)
ylabel('$$$')
xlabel('Time [days]')
datetick('x')

figure()
hold on;
plot(datez,yVal)
hold on;
plot(datez,yPred)
legend('Real data','Prediction', 'Location', 'Northeast');
datetick('x')