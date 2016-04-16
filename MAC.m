%clear
% Load Data
load('KexJobbData.mat')
ClPr = closingPrice(:,1:40);

%% Process data to ajust for NaNs

[dates2, ClPr] = removeNaN(dates, ClPr);
%ClPr = ClPr(:,1);

%% Dependence

Cov = cov(ClPr);
Corr = corr(ClPr);

%% Moving Averages
% Parameters
long = 3;
short = 2;


% Normal Weights
% wLong = ones(1, long)/long;
% wShort = ones(1, short)/short;

% Exponential Weights
alphaLong = 2/(long + 1); % Smoothing Param
wLong = repmat(1-alphaLong, 1, long).^(1:long);
wLong = wLong/sum(wLong);
alphaShort = 2/(short + 1); % Smoothing Param
wShort = repmat(1-alphaShort, 1, short).^(1:short);
wShort = wShort/sum(wShort);

avgClL = filter(wLong, 1, ClPr);
avgClS = filter(wShort, 1, ClPr);


%Linear Weights (sums of digits)
wLong = 1/((long+1)*(long/2)) * flipud((1:long)');
wShort = 1/((short+1)*(short/2)) * flipud((1:short)');


% Positioning
[row, col] = size(ClPr); % To be able to redo the matrixes later
trend = avgClS - avgClL;
gamma = sign(trend); % gamme(i,j) = +/- 1

%One day price difference
deltaP = diff(ClPr); % daily return
deltaP = [zeros(1,col);deltaP];

%Two days price difference
%deltaP = ClPr(3:end,:)-ClPr(1:end-2,:);
%deltaP = [zeros(2,col);deltaP];

STDEV = zeros(length(ClPr)-21,col);

%{
%If w use absolute return
for i = 1:length(ClPr)-21
    A = deltaP(i:i+21,:);
    STDEV(i,:) = std(A);
    
end
%}

%{
%If w use relative return
change = ClPr(2:end,:)./ClPr(1:end-1,:);
change = [zeros(1,col);change];
for i = 1:length(ClPr)-21
    A = change(i:i+21,:);
    STDEV(i,:) = std(A);
    
end
%}


%If w use static data points
for i = 1:length(ClPr)-21
    A = ClPr(i:i+21,:);
    STDEV(i,:) = std(A);
    
end
%}

STDEV = [ones(21,col);STDEV];
deltaP(1:21,:) = zeros(21,col);
ret = deltaP .*gamma./STDEV;
profit = cumsum(ret);
PROFIT = sum(profit');
PROFIT = PROFIT';

%% Plot 
% Closing Price and Moving Averages
figure(1);
%hold on;
%plot(dates2,ClPr)
ClPr2 = ClPr./STDEV;
plot(dates2, [ClPr2 avgClL avgClS ]);
legend(['Closing Prices'], ['Moving Average ' num2str(long)], ...
    ['Moving Average ' num2str(short)], 'location','best');
ylabel('$$$');
xlabel('Time [Days]');
title('Market Closing Prices');
datetick 'x'
%hold off;

% Returns
figure(2);
a(1) = subplot(2,1,1);
plot(dates2, profit)
datetick('x')
ylabel('$$$')
xlabel('Time [Days]');
title('Accumulated Profit')

a(2) = subplot(2,1,2);
plot(dates2, ret)
datetick('x')
ylabel('$$$')
xlabel('Time [Days]')
title('Daily return')
linkaxes([a(1) a(2)],'x');
%}
