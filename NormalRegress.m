function [b ,yHat, pred] = NormalRegress(yTrain, xTrain, nPred)

% Choose the nPred best predictors
[xTrain, pred] = choosePred(xTrain, yTrain, nPred); % Välj nPred stycken bästa marknader med crossval (function)
XTrain = [ones(size(xTrain(:,1))) xTrain];

% Solve with regress
b = regress(yTrain,XTrain);
yHat = XTrain*b;
end
