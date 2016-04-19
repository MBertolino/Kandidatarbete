clear;

% Load Data
load('KexJobbData.mat')
depMarket = 1;                   % Dependent Market
indepMarket = 1:40;              % Possible independent markets
lag = 1:4;                       % How many days ago we look at the indep markets

% Use this time period
[dates, clPr] = removeNaN(dates, closingPrice);

trainTime = 200;
predTime = 21;                   % How many days to predict
tradePeriods = 150;


%% Regression

% Sliding window
for j = 1:tradePeriods
    for i = 1:trainTime-predTime-lag(end)
        yTrain(i,:) = clPr(i+lag(end)+(j+1)*predTime,depMarket) - clPr(i+lag(end)+j*predTime,depMarket);
        xTemp = repmat(clPr(i+lag(end)+j*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+j*predTime,indepMarket);
        % xTemp = repmat(clPr(i+lag(end)-predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag-predTime,indepMarket);
        xTrain(i,:) = reshape(xTemp.',1,[]);
    end
    
    % Standardize data and add intercept
    [xTrain, mux, sigmax] = zscore(xTrain);
    XTrain = [ones(size(xTrain(:,1))) xTrain];
    [yTrain2, muy, sigmay] = zscore(yTrain);
    
    
    % Normal regress
    method{1} = 'Regress';
    [b1, yHat] = NormalRegress(yTrain, XTrain);
    
    % Lsq Ridge
    method{2} = 'Identity Ridge';
    [b2, yHat2] = RidgeRegress(yTrain, XTrain);
    
    
    %% Prediction
    % Validation
    yVal(j) = clPr(i+lag(end)+(j+2)*predTime,depMarket) - clPr(i+lag(end)+(j+1)*predTime,depMarket);
    
    xTemp = repmat(clPr(i+lag(end)+(j+1)*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+(j+1)*predTime,indepMarket);
    xVal = reshape(xTemp.',1,[]);
    [nRow, nCol] = size(xVal);
    xVal = (xVal - repmat(mux,length(yVal(j)),1));
    xVal = xVal./(sigmax'*ones(1,nRow))';
    XVal = [ones(size(xVal(:,1))) xVal];
    
    
    % Prediction
    yPred1(j) = XVal*b1; % Regress
    yPred2(j) = XVal*b2; % Ridge
    yPred = yPred2;
    
    % Dates adjustment
    datez(j) = dates(i+lag(end)+(j+1)*predTime);
end

% Calculate profit
gamma = sign(yPred);
ret = yVal.*gamma;
profit = cumsum(ret);

%% Plots

figure()
plot(datez, profit)
ylabel('$$$')
xlabel('Time [days]')
legend('OMX')
datetick('x')

figure()
hold on;
plot(datez,yPred,'-*')
hold on;
plot(datez,yVal)
datetick('x')

%{
figure()
for ip = 1:1
    %p(ip) = subplot(2,1,ip);
    hold on;
    plot(dates(2+lag(end)+predTime:end), [yTrain yHat(:,ip)])
    ylabel('$$$');
    xlabel('Time [Days]');
    title(['Multiple Linear Regression Fit, ' method(ip)]);
    legend('Real data','Prediction', 'Location', 'Northeast');
    %linkaxes(p,'x','y');
    datetick('x')
    
end



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

%}