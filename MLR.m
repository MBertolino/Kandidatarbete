clear;
%{
y = b1x1 + b2x2 + ... + bnxn + e

y - Dependent assets in which to choose a long or short position
    Price change from today to predTime days ahead

x - Assumed independent assets used to make a prediction
    Price change from day t-1 to t in the range of lag days back
%}
tic;


%% Setup
% Load Data
load('KexJobbData.mat')

% Assets
class = 4;
name = {'Equities', 'Fixed Income'};
% depAsset = find(assetClass == class); % 1:13 14:22 23:35 36:40
% indepAsset = depAsset;
depAsset = [1; 13]';
indepAsset = [1; 13]';
Ld = length(depAsset(:,1));
Li = length(indepAsset(:,1));

% Prediction Param
lag = [1 2];
predTime = 1;                    % How many days to predict
trainTime = 100;
lambda = [0];
Ll = length(lambda);

% Investment Param
bankStart = 10000;
risk = 0.05;

% Remove NaN's
% Start at 02-Jan-2009
% End at 28-Jan-2016 7447
[dates, clPr] = removeNaN(dates(7447 - trainTime - 2*predTime + 1:end), ...
    closingPrice(7447 - trainTime - 2*predTime + 1:end, :));
tradePeriods = floor((length(dates) - trainTime - predTime)/predTime) - 1;
diffClPr = diff(clPr);

<<<<<<< HEAD
for l = 1:length(lag)
    % Pre-allocating
    b = zeros(lag(l)*Li, 1);
    yTrain = zeros(trainTime - lag(l) - predTime, Ld);
    xTrain = zeros(trainTime - lag(l) - predTime, lag(l)*Li);
    yVal = zeros(tradePeriods, Ld);
    yPred = zeros(tradePeriods, Ll*Ld);
    sigmay = zeros(1, Ld);
    holding = zeros(tradePeriods, Ll*Ld);
    holdingTot = zeros(tradePeriods, Ll);
    datez = holdingTot;
    holding(1,:) = bankStart;
    gamma = zeros(tradePeriods, Ll*Ld);
    
    %% Regression
    % Data relating assets to be predicted, y, and assets to use as
    % predictors, x. For each day in a tradePeriods sized window y and x are
    % related.
    % After a prediction the training window is moved so only the latest data
    % points are used as training data as they are assumed to be more accurate
    
    h = waitbar(0,'Total time...');
    for iDep = 1:length(depAsset(1,:))
        % Sliding window
        for j = 1:tradePeriods
            %     yTrain(1:end - predTime, :) = yTrain(predTime + 1:end, :);
            %     xTrain(1:end - predTime, :) = xTrain(predTime + 1:end, :);
            for i = 1 + lag(l):trainTime - predTime
                yTrain(i-lag(l), :) = clPr(i + j*predTime, depAsset(:,iDep)) ...
                    - clPr(i + (j-1)*predTime, depAsset(:,iDep));
                xTemp = diffClPr(i - lag(l) + (j-1)*predTime : ...
                    i - 1 + (j-1)*predTime, indepAsset(:,iDep))';
                xTrain(i-lag(l), :) = reshape(xTemp.', 1, []);
            end
            
            % Standardize data and add intercept
            [yTrain, muy, sigmay] = zscore(yTrain);
            [xTrain, mux, sigmax] = zscore(xTrain);
            
            % For every invested asset, calculate the regression coefficients
            % using both OLS and Ridge
            b(:,1:Ll*Ld) = RidgeRegress(yTrain, xTrain, lambda);
            
            
            %% Prediction & Validation
            % Prediction of the change  in price of each asset
            % XVal - are the predictors
            % yPred - is the predicted change in price of each asset
            xTemp = diffClPr(i - lag(l) + j*predTime : ...
                i - 1 + j*predTime, indepAsset(:,iDep))';
            xVal = reshape(xTemp, 1, []);
            xVal = (xVal - mux)./sigmax;
            yPred(j,:) = xVal*b;
            
            % Smart positioning (optional)
            if j > 99
                gamma(j,:) = yPred(j,:)./mean(abs(yPred(j-99:j,:)));
            else
                gamma(j,:) = yPred(j,:)./mean(abs(yPred(1:j,:)));
            end
            
            % Validation
            % yVal - is the actual price change measured at the end of the
            % prediction time
            yVal(j,:) = clPr(i + (j+1)*predTime, depAsset(:,iDep)) ...
                - clPr(i + j*predTime, depAsset(:,iDep));
            yVal(j,:) = (yVal(j,:) - muy)./sigmay;
            
            % Dates adjustment
            % At each predicted day, the date is extracted
            datez(j,:) = dates(i + (j+1)*predTime);
            
            waitbar(j/tradePeriods)
        end
        
        
        %% Strategy
        % gamma - is the position to take for each asset
        % ret - is the risk adjusted return for taking a position
        % sharpe - is the information quotient for a strategy
        
        % Position and return
        gamma(abs(gamma) > 1) = sign(gamma(abs(gamma) > 1)); % Smart
        %     gamma = sign(yPred);                           % +/- 1
        ret = bsxfun (@rdivide, repmat(yVal,1,Ll).*gamma, repmat(std(yVal),1,Ll));
        
        % Calculate the evolution of a holding for each asset
        for ih = 2:length(ret)
            holding(ih,:) = holding(ih - 1, :).*(1 + risk*ret(ih - 1, :));
        end
        
        % Calculate the evolution of the total holding
        for il = 1:Ll
            sharpe(il,iDep) = mean(ret(:, 1 + (il-1)*Ld:il*Ld)) ...
                /std(ret(:, 1 + (il-1)*Ld:il*Ld)) ...
                * sqrt(250); % Annualized
            holdingTot(:,il) = holdingTot(:,il) + sum(holding(:, (il-1)*Ld + 1:il*Ld), 2)/Ld;
        end
        
        
        %% Plots
        % Plot the evolution of the total holding
        figure()
        plot(datez, holdingTot)
        ylabel('Holding [$]') % ;)
        xlabel('Time [Days]')
        title(['Holding using MLR on ' name(iDep) ' with lag ' lag(l)])
        str = cellstr(num2str(lambda', 'lambda = %d'));
        legend(str, 'Location', 'NorthWest');
        datetick('x')
        
    end
    holdingTot = holdingTot/iDep;
    disp(['Sharpe ratio for lag ' sharpe])
    close(h);
=======
% Pre-allocating
yTrain = zeros(trainTime - lag - predTime, Ld);
xTrain = zeros(trainTime - lag - predTime, lag*Li);
yVal = zeros(tradePeriods, Ld);
yPred = zeros(tradePeriods, Ll*Ld);
sigmay = zeros(1, Ld);
holding = zeros(tradePeriods, Ll*Ld);
holdingTot = zeros(tradePeriods, Ll);
datez = holdingTot;
holding(1,:) = bankStart;


%% Regression
% Data relating assets to be predicted, y, and assets to use as
% predictors, x. For each day in a tradePeriods sized window y and x are
% related.
% After a prediction the training window is moved so only the latest data
% points are used as training data as they are assumed to be more accurate

for iDep = 1:length(depAsset(1,:))
    % Sliding window
    for j = 1:tradePeriods
        %     yTrain(1:end - predTime, :) = yTrain(predTime + 1:end, :);
        %     xTrain(1:end - predTime, :) = xTrain(predTime + 1:end, :);
        for i = 1 + lag:trainTime - predTime
            yTrain(i-lag, :) = clPr(i + j*predTime, depAsset(:,iDep)) ...
                - clPr(i + (j-1)*predTime, depAsset(:,iDep));
            xTemp = diffClPr(i - lag + (j-1)*predTime : ...
                i - 1 + (j-1)*predTime, indepAsset(:,iDep))';
            xTrain(i-lag, :) = reshape(xTemp.', 1, []);
        end
        
        % Standardize data and add intercept
        [yTrain, muy, sigmay] = zscore(yTrain);
        [xTrain, mux, sigmax] = zscore(xTrain);
        
        % For every invested asset, calculate the regression coefficients
        % using both OLS and Ridge
        b(:,1:Ll*Ld) = RidgeRegress(yTrain, xTrain, lambda);
        
        
        %% Prediction & Validation
        % Prediction of the change  in price of each asset
        % XVal - are the predictors
        % yPred - is the predicted change in price of each asset
        xTemp = diffClPr(i - lag + j*predTime : ...
            i - 1 + j*predTime, indepAsset(:,iDep))';
        xVal = reshape(xTemp, 1, []);
        xVal = (xVal - mux)./sigmax;
        yPred(j,:) = xVal*b;
        
        % Validation
        % yVal - is the actual price change measured at the end of the
        % prediction time
        yVal(j,:) = clPr(i + (j+1)*predTime, depAsset(:,iDep)) ...
            - clPr(i + j*predTime, depAsset(:,iDep));
        yVal(j,:) = (yVal(j,:) - muy)./sigmay;
        
        % Dates adjustment
        % At each predicted day, the date is extracted
        datez(j,:) = dates(i + lag + (j+1)*predTime);
    end
    
    
    %% Strategy
    % gamma - is the position (1 or -1) for each asset
    % ret - is the risk adjusted return for taking a position
    % infoRet - is the information quotient for a strategy
    for i = 1:Ll
        gamma(:,1+(i-1)*Ld:i*Ld) = sign(yPred(:,1+(i-1)*Ld:i*Ld));
        ret(:,1+(i-1)*Ld:i*Ld) = bsxfun (@rdivide, yVal.*gamma(:,1+(i-1)*Ld:i*Ld), std(yVal));
        infoRet(i,iDep) = mean(ret(:,1+(i-1)*Ld:i*Ld))/std(ret(:,1+(i-1)*Ld:i*Ld)) ...
            * sqrt(250); % Annualized
    end
    infoRet
    
    % Calculate the evolution of a holding for each asset
    for ii = 2:length(ret)
        holding(ii,:) = holding(ii - 1, :).*(1 + risk*ret(ii - 1, :));
    end
    
    % Calculate the evolution of the total holding
    for il = 1:Ll
        holdingTot(:,il) = holdingTot(:,il) + sum(holding(:, (il-1)*Ld + 1: il*Ld), 2)/Ld;
    end
>>>>>>> master
end
holdingTot = holdingTot/iDep;

<<<<<<< HEAD
=======
%% Plots
% Plot expected returns for each asset
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

>>>>>>> master
hold off;
toc;
