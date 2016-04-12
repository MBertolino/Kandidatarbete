function [b, yHat, pred] = GeneralizedRegress(yTrain, xTrain, nPred)

% Choose the nPred best predictors
pred = choosePred(xTrain, yTrain, nPred);
xTrain = xTrain(:,pred);
XTrain = [ones(size(xTrain(:,1))) xTrain];  

% Implement Ridge
lambda = 1e-3;
XRidge = [XTrain; (lambda*eye(nPred+1))];
yRidge = [yTrain; zeros(nPred+1,1)];


% Solve the system using quadprog
H = XRidge'*XRidge;
f = -XRidge'*yRidge;
% lb = zeros(length(pred)+1,1);          % Possible lower bound
% ub = Inf*ones(length(pred)+1,1);       % Possible upper bound
b = quadprog(H,f);
yHat = XTrain*b;

end