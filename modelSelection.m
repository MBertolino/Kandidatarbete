clear;

% Load Data
load('KexJobbData.mat')

timeFrame = 7448;              % Frame to remove NaN 7448

% Prediction Param
trainTime = 300;              % 1000 1640 2040
predTime = 21;                 % How many days to predict
timeFrame = 7448;              % Frame to remove NaN 7448
lag = 10;
Li = 5;
dep = [1];
nPred = 5;

% Remove NaN's
% Start at 02-Jan-2009
% End at 06-Jan-2016
[datesNoNaN, clPr] = removeNaN(dates(timeFrame-predTime-trainTime:end), ...
    closingPrice(timeFrame-predTime-trainTime:end, :));
diffClPr = diff(clPr);

yTrain = zeros(trainTime - lag - predTime, 1);
RSS = zeros(size(clPr, 2), 1);
indepOut = [];

for iDep = dep
    indep = [];
    while length(indep) < nPred
        h = waitbar(0,['']);
        for k = 1:size(clPr, 2)
            
            ridgeEye = zeros(lag*(1 + length(indep)) + 1);
            xTrain = zeros(trainTime -lag - predTime, lag*(1 + length(indep)) + 1);
            
            if sum(indep == k) == 0
                for i = 1 + lag:trainTime - predTime
                    if k < 2
                        yTrain(i - lag, :) = clPr(i + predTime, iDep) ...
                            - clPr(i, iDep);
                    end
                    
                    xTemp = diffClPr(i - lag : i - 1, [k indep])';
                    xTrain(i - lag, :) = [1 reshape(xTemp.', 1, [])];
                end
                
                b = RidgeRegress(yTrain, xTrain, 0, ridgeEye);
                RSS(k) = mean((yTrain - xTrain*b).^2);
            end
            
            waitbar(k/size(clPr,2));
        end
        close(h);
        next = find(RSS == min(RSS));
        indep = [indep next];
    end
    indepOut = [indepOut indep'];
end
