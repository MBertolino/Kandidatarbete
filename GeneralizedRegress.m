function [b, yHat, pred] = GeneralizedRegress(yTrain, xTrain, nPred)

% Choose the nPred best predictors
pred = choosePred(xTrain, yTrain, nPred, 'Ridge');
xTrain = xTrain(:,pred);
XTrain = [ones(size(xTrain(:,1))) xTrain];

% Implement Ridge, don't penalize the intercept
lambda = 1e-3;
ridgeEye = eye(nPred+1);
ridgeEye(1,1) = 0;
XRidge = [XTrain; (lambda*ridgeEye)];
yRidge = [yTrain; zeros(nPred+1,1)];


% Solve the system using quadprog
H = XRidge'*XRidge;
f = -XRidge'*yRidge;
% lb = zeros(length(pred)+1,1);          % Possible lower bound
% ub = Inf*ones(length(pred)+1,1);       % Possible upper bound
b = quadprog(H,f);
yHat = XTrain*b;

end