% try

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

% Prediction Param
trainTime = 2200;
predTime = 21;                 % How many days to predict
timeFrame = 7448;              % Frame to remove NaN
lag = [100 200 300];
lambda = [0 1e0 1e1 1e2 1e3 1e4 1e5];
Ll = length(lambda);

% Investment Param
assetIndex = 1:7;
bankStart = 10000;
risk = 0.05;
smart = 1;                     % 0 or 1, 1 uses smart positioning, 0 don't

% Assets
name = {'Equities 1', 'Equities 2', 'Fixed Income 1', 'Fixed Income 2', ...
    'Commodities 1', 'Commodities 2', 'Foreign Exchange'};
% Calculate all asset classes at ones
for asset = assetIndex
    switch(asset)
        case(1)
            depAsset = [1:6]';
            indepAsset = [1:6]';
        case(2)
            depAsset = [7:13]';
            indepAsset = [7:13]';
        case(3)
            depAsset = [14:18]';
            indepAsset = [14:18]';
        case(4)
            depAsset = [19:22]';
            indepAsset = [19:22]';
        case(5)
            depAsset = [23:29]';
            indepAsset = [23:29]';
        case(6)
            depAsset = [30:35]';
            indepAsset = [30:35]';
        case(7)
            depAsset = [36:40]';
            indepAsset = [36:40]';
    end
    Ld = length(depAsset(:,1));
    Li = length(indepAsset(:,1));
    
    
    for l = 1:length(lag)
        % Remove NaN's
        % Start at 02-Jan-2009
        % End at 06-Jan-2016
        [datesNoNaN, clPr] = removeNaN(dates(timeFrame-predTime -trainTime:end), ...
            closingPrice(timeFrame-predTime-trainTime:end, :));
        tradePeriods = floor((length(datesNoNaN) - trainTime)/predTime);
        diffClPr = diff(clPr);
        
        % Pre-allocating lag-dependent variables
        b = zeros(lag(l)*Li, 1);
        yTrain = zeros(trainTime - lag(l) - predTime, Ld);
        xTrain = zeros(trainTime - lag(l) - predTime, lag(l)*Li);
        yVal = zeros(tradePeriods, Ld);
        yPred = zeros(tradePeriods, Ll*Ld);
        sigmay = zeros(1, Ld);
        holding = zeros(tradePeriods, Ll*Ld);
        holdingTot = zeros(tradePeriods, Ll);
        holding(1,:) = bankStart;
        holdingTot(1,:) = bankStart;
        datesAdjusted = holdingTot;
        gamma = zeros(tradePeriods, Ll*Ld);
        
        % Not an attractive solution but much faster to keep ridgeEye
        % outside RidgeRegression.m since it is redone unneccesarily
        ridgeEye = diag(repelem(lambda, 1, lag(l)*Ld)*eye(Ll*lag(l)*Ld));
        
        
        %% Regression
        % Data relating assets to be predicted, y, and assets to use as
        % predictors, x. For each day in a tradePeriods sized window y and x are
        % related.
        % After a prediction the training window is moved so only the latest data
        % points are used as training data as they are assumed to be more accurate
        
        % Create a waitbar to show calculation time
        h = waitbar(0,['Lag: ' num2str(l) '/' num2str(length(lag)) ...
            ', Asset class: ' num2str(asset) '/' num2str(assetIndex(end))]);
        
        % Sliding window
        for j = 1:tradePeriods
            % Speed Up - reuse data
            if j > 1
                yTrain(1:end - predTime, :) = yTrain(predTime + 1:end, :);
                xTrain(1:end - predTime, :) = xTrain(predTime + 1:end, :);
                start = trainTime - 2*predTime;
            else
                start = 1 + lag(l);
            end
            
            % Train the model
            for i = start:trainTime - predTime
                yTrain(i-lag(l), :) = clPr(i + j*predTime, depAsset) ...
                    - clPr(i + (j-1)*predTime, depAsset);
                xTemp = diffClPr(i - lag(l) + (j-1)*predTime : ...
                    i - 1 + (j-1)*predTime, indepAsset)';
                xTrain(i-lag(l), :) = reshape(xTemp.', 1, []);
            end
            
            % Standardize data
            [yTrainStd, muy, sigmay] = zscore(yTrain);
            [xTrainStd, mux, sigmax] = zscore(xTrain);
            
            % For every invested asset, calculate the regression coefficients
            % using both OLS and Ridge
            b(:,1:Ll*Ld) = RidgeRegress(yTrainStd, xTrainStd, lambda, ridgeEye);
            
            
            %% Prediction & Validation
            % Prediction of the change in price of each asset
            % xPred - are the predictors
            % yPred - is the predicted change in price of each asset
            xTemp = diffClPr(i - lag(l) + j*predTime : ...
                i - 1 + j*predTime, indepAsset)';
            xPred = reshape(xTemp, 1, []);
            xPred = (xPred - mean([xTrain; xPred]))./std([xTrain; xPred]);
            yPred(j,:) = xPred*b;
            
            % Smart positioning (optional)
            if smart > 0
                if j > 99
                    gamma(j,:) = yPred(j,:)./mean(abs(yPred(j-99:j,:)));
                else
                    gamma(j,:) = yPred(j,:)./mean(abs(yPred(1:j,:)));
                end
            end
            
            % Validation
            % yVal - is the actual price change measured at the end of the
            % prediction time
            yVal(j,:) = clPr(i + (j+1)*predTime, depAsset) ...
                - clPr(i + j*predTime, depAsset);
            yVal(j,:) = (yVal(j,:) - muy)./sigmay;
            
            % Dates adjustment
            % At each predicted day, the date is extracted
            datesAdjusted(j,:) = datesNoNaN(i + (j+1)*predTime);
            
            waitbar(j/tradePeriods);
        end
        
        
        %% Strategy
        % gamma - is the position to take for each asset
        % ret - is the risk adjusted return for each asset
        % retTot - is the total r.a return for each ridge tuning param
        % sharpe - is the information quotient for a strategy
        % holding - is the evolution of a holding for each asset
        % holdingTot - is the evolution of a holding for each asset
        % group
        
        % Positioning
        if smart > 0
            gamma(abs(gamma) > 1) = sign(gamma(abs(gamma) > 1));
        else
            gamma = sign(yPred);                           % +/- 1
        end
        
        % Returns and Sharpe
        ret = repelem(yVal,1,Ll).*gamma;
        retTot = cell2mat(arrayfun(@(x) sum(ret(:, x:Ll:end), 2), 1:Ll, 'uni', 0));
        sharpe = mean(retTot)./std(retTot)*sqrt(250/predTime);
        
        % Calculate the development of the total holding for each lambda
        for ih = 2:length(ret(:,1))
            holdingTot(ih,:) = holdingTot(ih - 1, :).*(1 + risk*retTot(ih - 1, :));
        end
        
        
        %% Plots
        % Plot the evolution of the total holding
        figure()
        plot(datesAdjusted, holdingTot)
        ylabel('Holding [$]') % ;)
        xlabel('Time [Days]')
        title(['Holding using MLR on ' name(asset) ' with lag ' lag(l)])
        str = cellstr(num2str(lambda', 'lambda = %d'));
        legend(str, 'Location', 'NorthWest');
        datetick('x')
        disp(['Sharpe ratio for lag ' num2str(lag(l)) ', and asset ' num2str(asset) ': ' num2str(sharpe)])
        close(h);
    end
end
hold off;
toc;
load('handel')
sound(y,Fs)

% catch
%     setpref('Internet','E_mail','mattias.bertolino@gmail.com');
%     sendmail('mattias.bertolino@gmail.com','Calculation failed.')
% end