clear;

% Load Data
load('KexJobbData.mat')
depMarket = 10;                   % Dependent Market
indepMarket = 1:40;              % Possible independent markets
lag = 1:5;                       % How many days ago we look at the indep markets


trainTime = 300;
predTime = 10;                   % How many days to predict
tradePeriods = 200;

% Use this time period
[dates, clPr] = removeNaN(dates, closingPrice);

%% Regression

% Sliding window
for j = 1:tradePeriods
    for i = 1:trainTime
        yTrain(i,:) = clPr(i+lag(end)+j*predTime,depMarket) - clPr(i+lag(end)+(j-1)*predTime,depMarket);
        xTemp = repmat(clPr(i+lag(end)+(j-1)*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+(j-1)*predTime,indepMarket);
        % xTemp = repmat(clPr(i+lag(end)-predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag-predTime,indepMarket);
        xTrain(i,:) = reshape(xTemp.',1,[]);
    end
    
    % Standardize data and add intercept
    [xTrain, mu, sigma] = zscore(xTrain);
    XTrain = [ones(size(xTrain(:,1))) xTrain];
    
    
    % Normal regress
    method{1} = 'Regress';
    [b1, yHat] = NormalRegress(yTrain, XTrain);
    
    %     % Lsq Ridge
    %     method{2} = 'Identity Ridge';
    %     [b2, yHat2] = RidgeRegress(yTrain, XTrain);
    
    %% Prediction
    % Validation
    yVal(j) = clPr(i+lag(end)+(j+1)*predTime,depMarket) - clPr(i+lag(end)+j*predTime,depMarket);
    xTemp = repmat(clPr(i+lag(end)+j*predTime,indepMarket),lag(end),1) - clPr(i+lag(end)-lag+j*predTime,indepMarket);
    xVal = reshape(xTemp.',1,[]);
    [nr, nc] = size(xVal);
    xVal = (xVal - repmat(mu,length(yVal(j)),1));
    xVal = xVal./(sigma'*ones(1,nr))';
    XVal = [ones(size(xVal(:,1))) xVal];
    
    % Prediction
    yPred(j) = XVal*b1;
    %yPred2(j) = XVal*b2; yVal]
    
    % Dates adjustment
    datez(j) = dates(i+lag(end)+(j+1)*predTime);
end

% Calculate profit
trend = sign(yPred);
gamma = zeros(size(trend));
gamma(trend > 0) = 1;
gamma(trend < 0) = -1;
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


% Clear variables
clear closingPrice begT endT indS indE nIndep p p2

%}