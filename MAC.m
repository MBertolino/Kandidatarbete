clear

% Load Data
tic;
load('KexJobbData2.mat')
MAXSharpes = zeros(40);
MINSharpes = MAXSharpes;

markets = 1;
<<<<<<< HEAD

=======
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
ClPR = ClPr(:,markets);

%Positionsbyten, kostnad --> 1/20-dels standardavvikelse.
%Portföljvikt genom std

<<<<<<< HEAD
=======

>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
%% Process data to ajust for NaNs
%[dates2, ClPr] = removeNaN(dates, ClPr);
%dates2 = dates2(2862:end-275);
%ClPr = ClPr(2862:end-275,:);
dates2 = dates2(2062:end);
ClPR = ClPR(2062:end,:);
<<<<<<< HEAD
=======


>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
%% Moving Averages
% Parameters
longVector = [127 252 378 504 756];
long = 200;
shortVector = [63 126 188 251];
short = 50;
stdevDays = 21;

HOLDINGS = zeros(length(dates2), length(longVector), length(shortVector));
% Investment
risk = 0.05;

MEANPROFFIT = zeros(length(longVector), length(shortVector));
MEANSTD = MEANPROFFIT;
MEANSHARPE = MEANPROFFIT;

meanProffit = 0;
meanStd = 0;
meanSharpe = 0;

%accProfitsNormal = zeros(length(dates2),length(longVector), length(shortVector));
%dayProfitsNormal = accProfitsNormal;

%{
endAccProfits = zeros(length(longVector), length(shortVector));
slopesProfits = endAccProfits;
stdevProfits = endAccProfits;
slopesOverStdev = endAccProfits;
%}
for l = 1:length(longVector)
    long = longVector(l);
    for s= 1:length(shortVector)
        short = shortVector(s);
        if short<long
<<<<<<< HEAD
            
            
            %% Weights   
            % Normal Weights
            wLong = ones(1, long)/long;
            wShort = ones(1, short)/short;
            %{
% Exponential Weights
alphaLong = 2/(long + 1); % Smoothing Param
wLong = repmat(1-alphaLong, 1, long).^(1:long);
wLong = wLong/sum(wLong);
alphaShort = 2/(short + 1); % Smoothing Param
wShort = repmat(1-alphaShort, 1, short).^(1:short);
wShort = wShort/sum(wShort);
            %}
            avgClL = filter(wLong, 1, ClPR);
            avgClS = filter(wShort, 1, ClPR);           
            
            % Linear Weights (sums of digits)
            %{
wLong = 1/((long+1)*(long/2)) * flipud((1:long)');
wShort = 1/((short+1)*(short/2)) * flipud((1:short)');

avgClL = filter(wLong, 1, ClPr);
avgClS = filter(wShort, 1, ClPr);
            %}
            
                      
=======
            
            
            %% Weights
            % Normal Weights
            wLong = ones(1, long)/long;
            wShort = ones(1, short)/short;
            
            %{
                % Exponential Weights
                alphaLong = 2/(long + 1); % Smoothing Param
                wLong = repmat(1-alphaLong, 1, long).^(1:long);
                wLong = wLong/sum(wLong);
                alphaShort = 2/(short + 1); % Smoothing Param
                wShort = repmat(1-alphaShort, 1, short).^(1:short);
                wShort = wShort/sum(wShort);
            
                % Linear Weights (sums of digits)
                wLong = 1/((long+1)*(long/2)) * flipud((1:long)');
                wShort = 1/((short+1)*(short/2)) * flipud((1:short)');

                avgClL = filter(wLong, 1, ClPr);
                avgClS = filter(wShort, 1, ClPr);
            %}
            
            avgClL = filter(wLong, 1, ClPR);
            avgClS = filter(wShort, 1, ClPR);
            
            
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
            %% Positioning
            [row, col] = size(ClPR); % To be able to redo the matrixes later
            trend = avgClS - avgClL;
            gamma = sign(trend); % gamma(i,j) = +/- 1
            
<<<<<<< HEAD
            %Calculating position changes, if long average = short average, the
            %position i kept
            [row2, col2] = find(gamma == 0);
            for i = 1:length(row2)
                gamma(row2(i),col2(i)) = gamma(row2(i)-1, col2(i));
=======
            % Calculating position changes, if long average = short average,
            % the position is kept
            [row2, col2] = find(gamma == 0);
            for i = 1:length(row2)
                gamma(row2(i), col2(i)) = gamma(row2(i) - 1, col2(i));
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
            end
            
            
            %% Transaction costs
            %Vector with number of position changes for each day
            posChangeDay = sum((abs(diff(gamma))/2), 2);
            posChangeTot = sum(posChangeDay);
            
            
            %% Returns
<<<<<<< HEAD
            %One day price difference
            deltaP = diff(ClPR); % daily return
            deltaP = [deltaP; zeros(1,col)];        
=======
            % One day price difference
            deltaP = diff(ClPR); % daily return
            deltaP = [deltaP; zeros(1, col)];
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
            
            % Standard deviation for returns
            stdev1 = zeros(length(ClPR) - stdevDays, col);
            for i = 1:length(ClPR) - stdevDays
                A = deltaP(i:i + stdevDays - 1, :);
<<<<<<< HEAD
                stdev1(i,:) = std(A);               
            end
            stdev1 = [ones(stdevDays,col); stdev1];

            deltaP(1:stdevDays,:) = zeros(stdevDays, col);
            
            %  Return
            ret1 = deltaP .*gamma./stdev1;
            retTot1 = sum(ret1,2)/col;
=======
                stdev1(i,:) = std(A);
            end
            stdev1 = [ones(stdevDays, col); stdev1];
            deltaP(1:stdevDays, :) = zeros(stdevDays, col);
            
            % Return
            ret1 = deltaP.*gamma./stdev1;
            retTot1 = sum(ret1, 2)/col;
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
            profit1 = cumsum(ret1);
            profitTot1 = cumsum(retTot1);
            
            % dayProfitsNormal(:,l,s) = retTot1;
            % accProfitsNormal(:,l,s) = profitTot1;
            
            holdings = zeros(length(dates2), col);
            holdings(1,:) = 10000;
            
            % Investing
            for ii = 2:length(ret1)
<<<<<<< HEAD
                holdings(ii,:) = holdings(ii-1,:).*(1 + risk*ret1(ii,:));
            end
            holdings = mean(holdings,2);
            %}
=======
                holdings(ii,:) = holdings(ii-1,:).*(1+risk*ret1(ii,:));
            end
            holdings = mean(holdings, 2);
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
            HOLDINGS(:,l,s) = holdings;
            
            %{
% Returns without weighting
figure(1);
a(1) = subplot(3,1,1);
plot(dates2, profitTot1)
datetick('x')
ylabel('Return')
xlabel('Time [Days]');
title('Accumulated Profit')

a(2) = subplot(3,1,2);
plot(dates2, retTot1')
datetick('x')
ylabel('Return')
xlabel('Time [Days]')
title('Daily return')

% Plot the investment
a(3) = subplot(3,1,3);
plot(dates2, holdings)
datetick('x')
ylabel('Holdings [$]') % ;)
xlabel('Time [Days]')
title('Profit using MAC in dollars')
linkaxes([a(1) a(2) a(3)],'x');
            %}
<<<<<<< HEAD
        
=======
            
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
            meanProffit = mean(diff(holdings));
            meanStd = std(diff(holdings));
            meanSharpe = meanProffit/meanStd;
            
<<<<<<< HEAD
            MEANSHARPE(1,l,s) = meanSharpe*sqrt(250);                
=======
            MEANSHARPE(1,l,s) = meanSharpe*sqrt(250);
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
        end
    end
    
    HOLDINGS (HOLDINGS == 0) = nan;
<<<<<<< HEAD
    
    Holdings = zeros(length(dates2),1);
=======
    Holdings = zeros(length(dates2), 1);
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
    for y = 1:length(dates2)
        Holdings(y) = nanmean(nanmean(HOLDINGS(y,:,:)));
    end
    
    figure(2)
    plot(dates2, Holdings);
    ylabel 'Holdings [$]'
    xlabel 'Year'
    title 'Value of holdings over time, using MAC'
    datetick;
    g = MEANSHARPE;
    g(g == 0) = nan;
    minSharp = min(min(min(g)));
    maxSharp = max(max(max(g)));
    
    MAXSharpe(y) = maxSharp;
    MINSharpe(y) = minSharp;
end

<<<<<<< HEAD
t=toc;
=======
t = toc;
>>>>>>> 3167cf9dfdeceaf033f37d1a19cc00a82c2a5fd7
