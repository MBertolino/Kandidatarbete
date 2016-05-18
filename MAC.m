clear
% Load Data
tic;
load('KexJobbData.mat')

%% Process data to ajust for NaNs
clPr = closingPrice;
[dates2, clPr] = removeNaN(dates, clPr);
dates2 = dates2(2062:end);
clPr = clPr(2062:end,:);


%% Moving Averages
% Parameters
% longVector = (50:50:2000);
% shortVector = (25:25:500);
longVector = 500;
shortVector = 200;
stdevDays = 21;

iWeights = 1;
iMarkets = 1;

MAXSharpes = zeros(iWeights, iMarkets, 1);
longShort = zeros(iWeights, iMarkets, 2);

% Investment
HOLDINGS = zeros(length(dates2), length(longVector), length(shortVector));
risk = 0.05;
MEANPROFIT = zeros(length(longVector), length(shortVector));
MEANSTD = MEANPROFIT;
MEANSHARPE = MEANPROFIT;

meanProfit = 0;
meanStd = 0;
meanSharpe = 0;

for iw = 1:iWeights
    disp(iw);
    
    % Calculate all asset classes at ones
    for im = 1:iMarkets
        disp(im)
        h = waitbar(0,['Weight type: ' num2str(iw) '/' num2str(iWeights) ...
            ', Market cathegory: ' num2str(im) '/' num2str(iMarkets)]);
        waitbar(im/(iMarkets*iWeights) + (iw - 1)/iWeights)
        
        switch(im)
            case(1)
                ClPR = clPr(:,1:6);
            case(2)
                ClPR = clPr(:,7:13);
            case(3)
                ClPR = clPr(:,14:18);
            case(4)
                ClPR = clPr(:,19:22);
            case(5)
                ClPR = clPr(:,23:29);
            case(6)
                ClPR = clPr(:,30:35);
            case(7)
                ClPR = clPr(:,36:40);
            case(8)
                ClPR = clPr(:,:);
        end
        
        for l = 1:length(longVector)
            long = longVector(l);
            for s = 1:length(shortVector)
                short = shortVector(s);
                if short < long
                    
                    
                    %% Weights
                    if iw == 1
                        % Normal weights
                        wLong = ones(1, long)/long;
                        wShort = ones(1, short)/short;
                    
                    elseif iw == 2       
                        % Exponential weights
                        alphaLong = 2/(long + 1); % Smoothing Param
                        wLong = repmat(1-alphaLong, 1, long).^(1:long);
                        wLong = wLong/sum(wLong);
                        alphaShort = 2/(short + 1); % Smoothing Param
                        wShort = repmat(1-alphaShort, 1, short).^(1:short);
                        wShort = wShort/sum(wShort);
                        
                    else
                        % Linear weights (sums of digits)
                        wLong = 1/((long+1)*(long/2)) * flipud((1:long)');
                        wShort = 1/((short+1)*(short/2)) * flipud((1:short)');          
                    end
                    
                    % Calculate filter
                    avgClL = filter(wLong, 1, ClPR);
                    avgClS = filter(wShort, 1, ClPR);      
                    
                    
                    %% Positioning
                    [row, col] = size(ClPR); % To be able to redo the matrices later
                    trend = avgClS - avgClL;
                    
                    %{
                    absTrend = abs(trend);
                    absMeans = zeros(row - short:col);
                    for mm = 1:row - short
                        absMeans(mm,:) = mean(absTrend(mm:mm + short, :));
                    end
                    trendMean = [ones(short,col); absMeans];
                    gamma = trend./trendMean;
                    gamma(abs(gamma) > 1) = sign(trend(abs(gamma) > 1));
                    %}
                    
                    gamma = sign(trend); % gamma(i,j) = +/- 1
                    
                    %Calculating position changes, if long average = short average, the
                    %position is kept
                    [row2, col2] = find(gamma == 0);
                    for i = 1:length(row2)
                        gamma(row2(i), col2(i)) = gamma(row2(i) - 1,  col2(i));
                    end
                    
                    
                    %% Returns
                    % One day price difference
                    deltaP = diff(ClPR); % daily return
                    deltaP = [deltaP; zeros(1, col)];
                    
                    %Standard deviation for returns
                    stdev1 = zeros(length(ClPR) - stdevDays, col);
                    for i = 1:length(ClPR) - stdevDays
                        A = deltaP(i:i + stdevDays - 1, :);
                        stdev1(i,:) = std(A);
                        
                    end
                    stdev1 = [ones(stdevDays, col); stdev1];
                         
                    %Return
                    deltaP(1:stdevDays, :) = zeros(stdevDays, col);
                    
                    ret1 = deltaP.*gamma./stdev1;
                    retTot1 = sum(ret1, 2)/col;
                    profit1 = cumsum(ret1);
                    profitTot1 = cumsum(retTot1);
                        
                    holdings = zeros(length(dates2), col);
                    holdings(1,:) = 1000;
                    
                    
                    %% Investing
                    for ii = 2:length(ret1)
                        holdings(ii,:) = holdings(ii-1,:).*(1 + risk*ret1(ii,:));
                    end
                    holdings = mean(holdings, 2);
                    HOLDINGS(:,l,s) = holdings;
                    
                    meanProfit = mean((ret1(200:end,:)));
                    meanStd = std((ret1(200:end,:)));
                    meanSharpe = meanProfit/meanStd;
                    
                    MEANSHARPE(l,s) = meanSharpe*sqrt(250);
                end
            end
            
            HOLDINGS (HOLDINGS == 0) = nan;
            Holdings = zeros(length(dates2), 1);
            
            for y = 1:length(dates2)
                Holdings(y) = nanmean(nanmean(HOLDINGS(y,:,:)));
            end
%             %{
            figure()
            plot(dates2, Holdings);
            ylabel 'Holdings [$]'
            xlabel 'Year'
            title 'Value of holdings over time, using MAC'
            datetick;
%             %}
            
        end
        g = MEANSHARPE;
        g(g == 0) = nan;
        maxSharpe = nanmax(nanmax(nanmax(g)));
        MAXSharpes(iw, im) = maxSharpe;
        [iLong, iShort] = find (g == maxSharpe);
        if length(iLong) > 1
            iLong = floor(mean(iLong));
            iShort = floor(mean(iShort));
        end
        longShort(iw, im, 1) = longVector(iLong);
        longShort(iw, im, 2) = shortVector(iShort);
        close(h);
    end
end

% surf(longVector, shortVector, g);
% xlabel 'Length of long mean'
% ylabel 'Length of short mean'
% zlabel 'Sharpe Ratio'

% stem(reshape(longShort(:,:,1), 24, 1));

t = toc;
