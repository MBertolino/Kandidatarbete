function [b, yHat] = RidgeRegress(yTrain, XTrain)

[row, col] = size(XTrain);

% Penalty factor
lambda = 1e5;
ridgeEye = eye(col);
ridgeEye(1,1) = 0;

% Ridge data
% XRidge = [XTrain; (lambda*ridgeEye)];
% yRidge = [yTrain; zeros(col,1)];

% Solving with lsqr
b = lsqr(XRidge,yRidge);

% Fit to training data
yHat = XTrain*b;

end