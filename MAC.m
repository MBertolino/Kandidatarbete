%clear
% Load Data
load('KexJobbData.mat')
ClPr = closingPrice(:,1:40);

%% Process data to ajust for NaNs

[dates2, ClPr] = removeNaN(dates, ClPr);


%% Moving Averages

% Parameters
PosChanges = zeros(50,50);
retTot = PosChanges;
long = 200;
short = 50;
stdevDays = 21;
portfolioWeights = ones(size(ClPr));


%% Weights

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
%{
wLong = 1/((long+1)*(long/2)) * flipud((1:long)');
wShort = 1/((short+1)*(short/2)) * flipud((1:short)');

avgClL = filter(wLong, 1, ClPr);
avgClS = filter(wShort, 1, ClPr);
%}


%% Positioning
[row, col] = size(ClPr); % To be able to redo the matrixes later
trend = avgClS - avgClL;
gamma = sign(trend); % gamma(i,j) = +/- 1

%Calculating position changes, if long average = short average, the
%position i kept
[row2, col2] = find(gamma == 0);
for i = 1:length(row2)
    gamma(row2(i),col2(i)) = gamma(row2(i)-1, col2(i));
end


%% Transaction costs

%Vector with number of position changes for each day
posChangeDay = sum((abs(diff(gamma))/2)');
posChangeTot = sum(posChangeDay);


%% Returns
%One day price difference
deltaP = diff(ClPr); % daily return
deltaP = [deltaP;zeros(1,col)];


%Standard deviation
stdev = zeros(length(ClPr)-stdevDays,col);
%If w use absolute return
for i = 1:length(ClPr)-stdevDays
    A = deltaP(i:i+stdevDays-1,:);
    stdev(i,:) = std(A);
    
end

stdev = [ones(stdevDays,col);stdev];

deltaP(1:stdevDays,:) = zeros(stdevDays,col);
ret = deltaP .*gamma./stdev;
retSum=sum(ret);

profit = cumsum(ret);
profitTot = cumsum(retTot);



%{
%% Plot

% Closing Price and Moving Averages
figure(1);
%hold on;
%plot(dates2,ClPr)
ClPr2 = deltaP./stdev;
plot(dates2, [ClPr avgClL avgClS ]);
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
plot(dates2, profitTot)
datetick('x')
ylabel('$$$')
xlabel('Time [Days]');
title('Accumulated Profit')

a(2) = subplot(2,1,2);
plot(dates2, retTot')
datetick('x')
ylabel('$$$')
xlabel('Time [Days]')
title('Daily return')
linkaxes([a(1) a(2)],'x');
%}
