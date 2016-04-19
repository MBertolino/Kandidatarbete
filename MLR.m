clear;

% Load Data
load('KexJobbData.mat')
depMarket = 1;                   % Dependent Market
indepMarket = 1:40;              % Possible independent markets
lag = 1:4;                       % How many days ago we look at the indep markets

% Use this time period
[dates, clPr] = removeNaN(dates, closingPrice);

trainTime = 400;
predTime = 21;                   % How many days to predict
tradePeriods = 150;




%% Regression

% Pre-allocating
yTrain = zeros(trainTime-lag(end)-predTime,1);
xTrain = zeros(trainTime-lag(end)-predTime, lag(end)*length(indepMarket));
yVal = zeros(tradePeriods,1);
yPred1 = yVal;
yPred2 = yVal;
datez = yVal;
b = zeros(lag(end)*length(indepMarket)+1,2);

% Sliding window
for j = 1:tradePeriods
    for i = 1:trainTime-predTime-lag(end)
        yTrain(i,:) = clPr(i+lag(end)+(j+1)*predTime,depMarket) - clPr(i+lag(end)+j*predTime,depMarket);
        xTemp = repmat(clPr(i+lag(end)+j*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+j*predTime,indepMarket);
        xTrain(i,:) = reshape(xTemp.',1,[]);
    end
    
    
    % Standardize data and add intercept
    [xTrain, mux, sigmax] = zscore(xTrain);
    XTrain = [ones(size(xTrain(:,1))) xTrain];
    [yTrain2, muy, sigmay] = zscore(yTrain);
    
    
    % Normal regress
    method{1} = 'Regress';
    b(:,1) = regress(yTrain,XTrain);
    
    % Lsq Ridge
    method{2} = 'Identity Ridge';
    b(:,2) = RidgeRegress(yTrain, XTrain);
    
    
    %% Prediction
    
    % Prediction
    xTemp = repmat(clPr(i+lag(end)+(j+1)*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+(j+1)*predTime,indepMarket);
    xVal = reshape(xTemp.',1,[]);
    xVal = xVal - repmat(mux,length(xVal(j)),1);
    xVal = xVal./sigmax;
    XVal = [ones(size(xVal(:,1))) xVal];
    
    yPred1(j,1) = XVal*b(:,1); % Regress
    yPred1(j,2) = yPred1(j,1) + sigmay; % + 1 std
    yPred1(j,3) = yPred1(j,1) - sigmay; % - 1 std
    yPred2(j) = XVal*b(:,2); % Ridge
    
    
    % Validation
    yVal(j) = clPr(i+lag(end)+(j+2)*predTime,depMarket) - clPr(i+lag(end)+(j+1)*predTime,depMarket);
    
    
    % Dates adjustment
    datez(j) = dates(i+lag(end)+(j+1)*predTime);
end

% Calculate profit
yPred = yPred1;
gamma = sign(yPred);
for k = 1:3
    ret(:,k) = yVal.*gamma(:,k);
    profit(:,k) = cumsum(ret(:,k));
end

%% Plots

figure()
plot(datez, profit)
ylabel('$$$')
xlabel('Time [days]')
legend('OMX')
datetick('x')

figure()
hold on;
plot(datez,yPred)
hold on;
% plot(datez,yVal)
legend('Real data','Prediction', 'Location', 'Northeast');
datetick('x')