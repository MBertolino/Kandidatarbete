function [b ,yHat] = NormalRegress(yTrain, XTrain)

% Solve with regress
b = regress(yTrain,XTrain);
yHat = XTrain*b;

end
