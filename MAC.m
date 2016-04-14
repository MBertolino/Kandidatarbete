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
long = 30;
short = 20;

% % Normal No Weights (SMA)
% wLong = ones(1, long)/long;
% wShort = ones(1, short)/short;


%Linear Weights (sums of digits)
wLong = 1/((long+1)*(long/2)) * flipud((1:long)');
wShort = 1/((short+1)*(short/2)) * flipud((1:short)');

% % Exponential Weights (EWMA)
% alphaL = 2/(long + 1); % Smoothing Param
% wLong = repmat(1-alphaL, 1, long).^(1:long);
% wLong = wLong/sum(wLong);
% alphaS = 2/(short + 1); % Smoothing Param
% wShort = repmat(1-alphaS, 1, short).^(1:short);
% wShort = wShort/sum(wShort);

avgClL = filter(wLong, 1, clPr);
avgClS = filter(wShort, 1, clPr);

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
