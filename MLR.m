clear;

% Load Data
load('KexJobbData.mat')
[dates, clPr] = removeNaN(dates, closingPrice);

% Parameters
depMarket = 1:2;
indepMarket = 1:40;
lag = 1:3;                       % How many days ago we look at the indep markets

trainTime = 300;
predTime = 21;                   % How many days to predict
tradePeriods = 160;


%% Regression

% Pre-allocating
yTrain = zeros(trainTime-lag(end)-predTime,length(depMarket));
xTrain = zeros(trainTime-lag(end)-predTime, lag(end)*length(indepMarket));
yVal = zeros(tradePeriods,length(depMarket));
yPred1 = yVal;
yPred2 = yVal;
sigmay = yVal;
datez = yVal;
b1 = zeros(lag(end)*length(indepMarket)+1,length(depMarket));
b2 = b1;
    
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
        [yTrain, muy, sigmay] = zscore(yTrain);
        
        for m = 1:length(depMarket)
            % OLS Regression
            method{1} = 'OLS Regression';
%            b1(:,m) = regress(yTrain(:,m),XTrain);
            
            % Ridge Regression
            method{2} = 'Ridge Regression';
            b2(:,m) = RidgeRegress(yTrain(:,m), XTrain);
        end
        
        
        %% Prediction & Validation
        % Prediction
        xTemp = repmat(clPr(i+lag(end)+(j+1)*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+(j+1)*predTime,indepMarket);
        xVal = reshape(xTemp.',1,[]);
        xVal = (xVal - mux)./sigmax;
        XVal = [ones(size(xVal(:,1))) xVal];
        
%         yPred1(j,:) = XVal*b1;             % Regress
        yPred2(j,:) = XVal*b2;             % Ridge
        
               
        % Validation
        yVal(j,:) = clPr(i+lag(end)+(j+2)*predTime,depMarket) - clPr(i+lag(end)+(j+1)*predTime,depMarket);
        yVal(j,:) = (yVal(j,:) - muy)./sigmay;  
        
        % Dates adjustment
        datez(j) = dates(i+lag(end)+(j+1)*predTime);
     end
    
    % Calculate profit
    yPred = yPred2;
    gamma = sign(yPred).*repmat(sigmay,tradePeriods,1);
    gamma = gamma./sum(sigmay,2);
    ret = yVal.*gamma;
    profit = cumsum(ret);

profitTot = sum(profit,2);

%% Plots

% Plot accumulated profit for each market
figure()
plot(datez(:,1), profit)
title('Accumulated Profit')
ylabel('$$$')
xlabel('Time [days]')
datetick('x')

% Plot accumulated total profit
figure()
plot(datez(:,1), profitTot)
title('Accumulated Total Profit')
ylabel('$$$')
xlabel('Time [days]')
datetick('x')

% Plot expected returns for each market
figure()
hold on;
plot(datez(:,1),yVal)
hold on;
plot(datez(:,1),yPred)
title('Expected return')
legend('Real data','Prediction', 'Location', 'Northeast');
datetick('x')