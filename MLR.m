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
trainTime = 1640; %1640
predTime = 21;                    % How many days to predict
timeFrame = [5787];               % Frame to remove NaN

% Setup Param
lag = [100];
assetIndex = 1:7;
lambda = [0 1e8];
Ll = length(lambda);

% Investment Param
bankStart = 10000;
risk = 0.01;

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
        [datesNoNaN, clPr] = removeNaN(dates(5787:end), closingPrice(5787:end, :));
        tradePeriods = floor((length(datesNoNaN) - trainTime)/predTime);
        diffClPr = diff(clPr);
        
        % Pre-allocating
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
        
        
        %% Regression
        % Data relating assets to be predicted, y, and assets to use as
        % predictors, x. For each day in a tradePeriods sized window y and x are
        % related.
        % After a prediction the training window is moved so only the latest data
        % points are used as training data as they are assumed to be more accurate
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
            b(:,1:Ll*Ld) = RidgeRegress(yTrain, xTrain, lambda);
            
            %% Prediction & Validation
            % Prediction of the change  in price of each asset
            % xPred - are the predictors
            % yPred - is the predicted change in price of each asset
            xTemp = diffClPr(i - lag(l) + j*predTime : ...
                i - 1 + j*predTime, indepAsset)';
            xPred = reshape(xTemp, 1, []);
            %                         xPred = (xPred - mux)./sigmax;
            xPred = (xPred - mean([xTrain(end-99:end, :); xPred]))./std([xTrain(end-99:end, :); xPred]);
            yPred(j,:) = xPred*b;
            
            % Smart positioning (optional)
            if j > 99
                gamma(j,:) = yPred(j,:)./mean(abs(yPred(j-99:j,:)));
            else
                gamma(j,:) = yPred(j,:)./mean(abs(yPred(1:j,:)));
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
        % ret - is the risk adjusted return for taking a position
        % sharpe - is the information quotient for a strategy
        % holding - is the evolution of a holding for each asset
        % holdingTot - is the evolution of a holding for each asset
        % group
        
        % Position and return
        gamma(abs(gamma) > 1) = sign(gamma(abs(gamma) > 1)); % Smart
        %         gamma = sign(yPred);                           % +/- 1
%                 ret = bsxfun (@rdivide, repelem(yVal,1,Ll).*gamma, repelem(std(yVal),1,Ll));
        ret = repelem(yVal,1,Ll).*gamma;
        retTot = cell2mat(arrayfun(@(x) sum(ret(:, x:Ll:end), 2), 1:Ll, 'uni', 0));
        sharpe = mean(retTot)./std(retTot)*sqrt(250/predTime);
        
        % Calculate the evolution of a holding for each asset
        for ih = 2:length(ret(:,1))
            holdingTot(ih,:) = holdingTot(ih - 1, :).*(1 + risk*retTot(ih - 1, :));
        end
        
        %         % Calculate the evolution of the total holding for every lambda
        %         for il = 1:Ll
        %             holdingTot(:,il) = sum(holding(:, il:Ll:Ll*Ld), 2)/Ld;
        %         end
        
        
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