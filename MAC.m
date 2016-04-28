%clear
% Load Data
load('KexJobbData.mat')
ClPr = closingPrice(:,1:40);

%% Process data to ajust for NaNs

[dates2, ClPr] = removeNaN(dates, ClPr);


%% Moving Averages

% Parameters
long = 200;
short = 50;
stdevDays = 21;

% Investment
risk = 1;
bank(1) = 10000;


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
posChangeDay = sum((abs(diff(gamma))/2),2);
posChangeTot = sum(posChangeDay);


% Returns
%One day price difference
deltaP = diff(ClPr); % daily return
deltaP = [deltaP;zeros(1,col)];


%Standard deviation for returns
stdev1 = zeros(length(ClPr)-stdevDays,col);
for i = 1:length(ClPr)-stdevDays
    A = deltaP(i:i+stdevDays-1,:);
    stdev1(i,:) = std(A);
    
end
stdev1 = [ones(stdevDays,col);stdev1];


%Standard deviation for portfolio weights
stdev2 = zeros(length(ClPr)-stdevDays,col);
changeFactor = [ClPr(2:end,:)./ClPr(1:end-1,:);zeros(1,col)];
for i = 1:length(ClPr)-stdevDays
    B = changeFactor(i:i+stdevDays-1,:);
    stdev2(i,:) = std(B);
end
stdev2 = [ones(stdevDays,col);stdev2];
stdev2Sum = sum(stdev2,2);

wPortfolio = repmat(stdev2Sum,1,col)./stdev2;
wPortfolio = wPortfolio./repmat(sum(wPortfolio,2),1,col);


deltaP(1:stdevDays,:) = zeros(stdevDays,col);

%Return with equal portfolio weighting
ret1 = deltaP .*gamma./stdev1;
retTot1 = sum(ret1,2)/col;
profit1 = cumsum(ret1);
profitTot1 = cumsum(retTot1);
infoRet1 = mean(ret1)/std(ret1) ...
        * sqrt(250);

%Return with risk-based portfolio weighting
ret2 = col*ret1.*wPortfolio;
retTot2 = sum(ret2,2)/col;
profit2 = cumsum(ret2);
profitTot2 = cumsum(retTot2);

% Investing
for ii = 2:length(ret)
    bank(ii) = bank(ii-1)*(1+risk*ret(ii,1)); % Om vi delar ret med 10 så får vi rimliga siffror!
end


%% Plot

% Returns without weighting
figure();
a(1) = subplot(2,1,1);
plot(dates2, profitTot1)
datetick('x')
ylabel('Return')
xlabel('Time [Days]');
title('Accumulated Profit with equal Portfolio Weights')

a(2) = subplot(2,1,2);
plot(dates2, retTot1')
datetick('x')
ylabel('Return')
xlabel('Time [Days]')
title('Daily return with equal Portfolio Weights')
linkaxes([a(1) a(2)],'x');

% Returns with weighting
figure();
a(1) = subplot(2,1,1);
plot(dates2, profitTot2)
datetick('x')
ylabel('Return')
xlabel('Time [Days]');
title('Accumulated Profit with stdev Portfolio Weights')

a(2) = subplot(2,1,2);
plot(dates2, retTot2')
datetick('x')
ylabel('Return')
xlabel('Time [Days]')
title('Daily return with stdev Portfolio Weights')
linkaxes([a(1) a(2)],'x');

% Plot the investment
figure()
plot(bank)
ylabel('Profit [$$$]') % ;)
xlabel('Time [Days]')
title('Profit using MAC in dollars')
%}
