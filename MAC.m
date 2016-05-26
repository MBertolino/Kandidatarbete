% Moving Average Crossover - main program
% This code applies moving average crossover on 40 assets for three
% different positioning strategies. The filter lengths are varied between
% 50 and 2000 respectively 25 and 500.
% The assets were divided into seven groups in order to compare the results
% to the results from MLR.
%
% 2016 Iliam Barkino, Mattias Bertolino


clear
% Load Data
load('KexJobbData.mat')


%% Process data to ajust for NaNs
clPr = closingPrice;
[dates2, clPr] = removeNaN(dates, clPr);
dates2 = dates2(2062:end); % Starting from 2009
clPr = clPr(2062:end,:);


%% Parameters
longVector = (50:50:2000); % Different lengths of long filter
shortVector = (25:25:500); % Different lengths of short filter

stdevDays = 21; % Number of days used to calculate the standard deviation

iWeights = 3; % Number of positioning strategies
iMarkets = 7; % Number of asset groups
iGamma = 1; % Binary switch to determine what positioning strategy to use

% Matrix containing all obtained maximal values of information ratio
maxIR = zeros(iWeights, iMarkets, 1);

% Matrix showing what combinations of long and short filters that yielded
% the highest information ratios
longShort = zeros(iWeights, iMarkets, 2);

% Matrix containing all developements of invested capital for different
% combinations of long and short filters
holdingsMatrix = zeros(length(dates2), length(longVector),...
    length(shortVector));
risk = 0.05; % risk aversion

% Matrix containing all information ratios for different combinations of
% long and short filters
infoRatios = zeros(length(longVector), length(shortVector));


%% Starting looping for all possible combinations
for iw = 1:iWeights
    for im = 1:iMarkets
        
        % Bar that shows status of the run
        h = waitbar(0,['Weight type: ' num2str(iw) '/' num2str(iWeights) ...
            ', Market cathegory: ' num2str(im) '/' num2str(iMarkets)]);
        waitbar(im/(iMarkets*iWeights) + (iw - 1)/iWeights)
        
        switch(im) % Deciding which group of assets to evaluate
            case(1)
                clPr2 = clPr(:,1:6);
            case(2)
                clPr2 = clPr(:,7:13);
            case(3)
                clPr2 = clPr(:,14:18);
            case(4)
                clPr2 = clPr(:,19:22);
            case(5)
                clPr2 = clPr(:,23:29);
            case(6)
                clPr2 = clPr(:,30:35);
            case(7)
                clPr2 = clPr(:,36:40);
        end
        
        for l = 1:length(longVector)
            long = longVector(l);
            for s = 1:length(shortVector)
                short = shortVector(s);
                if short < long % long filter must be larger than short filter
                    
                    
                    %% Weights
                    if iw == 1
                        % Normal Weights
                        wLong = ones(1, long)/long;
                        wShort = ones(1, short)/short;
                    elseif iw == 2
                        % Exponential Weights
                        alphaLong = 2/(long + 1); % Smoothing Param
                        wLong = repmat(1-alphaLong, 1, long).^(1:long);
                        wLong = wLong/sum(wLong);
                        alphaShort = 2/(short + 1); % Smoothing Param
                        wShort = repmat(1-alphaShort, 1, short).^(1:short);
                        wShort = wShort/sum(wShort);
                    else
                        % Linear Weights (sums of digits)
                        wLong = 1/((long+1)*(long/2)) * flipud((1:long)');
                        wShort = 1/((short+1)*(short/2)) * flipud((1:short)');
                    end
                    
                    % Filters
                    avgClL = filter(wLong, 1, clPr2);
                    avgClS = filter(wShort, 1, clPr2);
                    
                    
                    %% Positioning
                    % To be able to redo the matrixes later
                    [row, col] = size(clPr2);
                    trend = avgClS - avgClL;
                    
                    % Choose positioning method (1 = smart, 0 = simple)
                    if iGamma > 0
                        absTrend = abs(trend);
                        absMeans = zeros(row-short:col);
                        for mm = 1:row-short
                            absMeans(mm,:) = mean(absTrend(mm:mm + short, :));
                        end
                        trendMean = [ones(short,col); absMeans];
                        gamma = trend./trendMean;
                        gamma(abs(gamma) > 1) = sign(trend(abs(gamma) > 1));
                    else
                        gamma = sign(trend); % gamma(i,j) = +/- 1
                        
                        % Calculating position changes, if
                        % long average = short average, the position i kept
                        [row2, col2] = find(gamma == 0);
                        for i = 1:length(row2)
                            gamma(row2(i), col2(i)) = ...
                                gamma(row2(i) - 1, col2(i));
                        end
                    end
                    
                    
                    %% Returns
                    % One day price difference
                    deltaP = diff(clPr2); % daily return
                    % row of zeros for later calculations
                    deltaP = [deltaP; zeros(1, col)];
                    
                    % Standard deviation for returns
                    stdev1 = zeros(length(clPr2) - stdevDays, col);
                    for i = 1:length(clPr2) - stdevDays
                        A = deltaP(i:i + stdevDays - 1, :);
                        stdev1(i,:) = std(A);
                    end
                    % Dimension fit
                    stdev1 = [ones(stdevDays, col); stdev1];
                    
                    % Adjust for standard deviation dimension
                    deltaP(1:stdevDays, :) = zeros(stdevDays, col);
                    
                    % Return of each asset on each day
                    ret1 = deltaP.*gamma./stdev1;
                    % Mean return of all assets on each day
                    retTot1 = sum(ret1, 2)/col;
                    
                    
                    %% Investing
                    % Developement of invested capital
                    holdings = zeros(length(dates2), col);
                    holdings(1,:) = 10000; % Starting att $10,000
                    
                    for ii = 2:length(ret1)
                        % Calculating developement
                        holdings(ii,:) = holdings(ii - 1,:).*...
                            (1 + risk*ret1(ii,:));
                    end
                    
                    % Storing the mean developement for all assets in a
                    % specific group
                    holdingsMatrix(:,l,s) = mean(holdings, 2);
                    
                    %Calculating mean return
                    meanProffit = mean((ret1(22:end, :)));
                    %Calculating standard deviation of returns
                    proffitStd = std((ret1(22:end, :)));
                    
                    %Information ratio on daily bais
                    meanInfoRet = meanProffit/proffitStd;
                    
                    %Annualizing the information ratio
                    infoRatios(l,s) = meanInfoRet*sqrt(250);
                end
            end
            %{
            holdingsMatrix (holdingsMatrix == 0) = nan;
            Holdings = zeros(length(dates2), 1);
            
            for y = 1:length(dates2)
                Holdings(y) = nanmean(nanmean(holdingsMatrix(y,:,:)));
            end
            %}
        end
        
        % Finding maximal Information Ratio between aall combinations of
        % long and short filter length
        IR = infoRatios;
        IR(IR == 0) = nan;
        maxIR(iw,im) = nanmax(nanmax(nanmax(IR)));
        
        % Finding what combination of long and short filter length that
        % yielded maximal Information Ratio
        [iLong, iShort] = find (IR == maxIR(iw,im));
        
        if length(iLong) > 1
            iLong = floor(mean(iLong));
            iShort = floor(mean(iShort));
        end
        % Storing filter combinations with highest Information Ratio
        longShort(iw,im,1) = longVector(iLong);
        longShort(iw,im,2) = shortVector(iShort);
        
        close(h);
    end
end
