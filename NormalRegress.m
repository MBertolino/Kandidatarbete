function [b ,yHat, pred] = NormalRegress(yTrain, xTrain, nPred)

% Choose the nPred best predictors
[pred] = choosePred(xTrain, yTrain, nPred); % Välj nPred stycken bästa marknader med crossval (function)
xTrain = xTrain(:, pred);
XTrain = [ones(size(xTrain(:,1))) xTrain];

% Solve with regress
b = regress(yTrain,XTrain);
yHat = XTrain*b;

end
