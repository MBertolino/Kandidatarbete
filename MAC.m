%clear
% Load Data
load('KexJobbData.mat')
clPr = closingPrice(:,1);

%% Process data to ajust for NaNs

[dates2, clPr] = removeNaN(dates, clPr);

%% Dependence

Cov = cov(clPr);
Corr = corr(clPr);

%% Moving Averages
% Parameters
long = 200;
short = 50;

% % Normal No Weights (SMA)
% coeffL = ones(1, long)/long;
% coeffS = ones(1, short)/short;

% Weighted (WMA)
% coeffL = (long + 1 - (1:long))/sum(1:long);
% coeffS = (short + 1 - (1:short))/sum(1:short);

% % Exponential Weights (EWMA)
% lambdaL = 2/(long + 1); % Smoothing Param
% coeffL = repmat(1-lambdaL, 1, long).^(1:long);
% coeffL = coeffL/sum(coeffL);
% lambdaS = 2/(short + 1); % Smoothing Param
% coeffS = repmat(1-lambdaS, 1, short).^(1:short);
% coeffS = coeffS/sum(coeffS);

avgClL = filter(coeffL, 1, clPr);
avgClS = filter(coeffS, 1, clPr);

% Positioning
trend = avgClS - avgClL;
gamma = ones(size(trend(2:end,:)))*(-1);
gamma(trend(2:end,:) > 0) = 1;
deltaP = clPr(2:end,:)-clPr(1:end-1,:);
ret = deltaP.*gamma;
profit = cumsum(ret);

%% Plot 
% Closing Price and Moving Averages
figure();
hold on;
plot(dates2,clPr)
plot(dates2, [avgClL avgClS]);
legend(['Closing Prices'], ['Moving Average ' num2str(long)], ...
    ['Moving Average ' num2str(short)], 'location','best');
ylabel('$$$');
xlabel('Time [Days]');
title('OMX Closing Price');
datetick 'x'

% Returns
figure();
a(1) = subplot(2,1,1);
plot(dates2(2:end), profit)
datetick('x')
ylabel('$$$')
xlabel('Time [Days]');
title('Accumulated Profit')

a(2) = subplot(2,1,2);
plot(dates2(2:end), ret)
datetick('x')
ylabel('$$$')
xlabel('Time [Days]')
title('Daily return')
linkaxes([a(1) a(2)],'x');
%}