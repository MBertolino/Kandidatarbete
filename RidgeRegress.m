function [b, yHat, pred] = RidgeRegress(yTrain, xTrain, nPred)

% Choose the nPred best predictors
pred = choosePred(xTrain, yTrain, nPred, 'Ridge');
xTrain = xTrain(:,pred);
XTrain = [ones(size(xTrain(:,1))) xTrain];

% Find the best lambda with cross validation
lambda = linspace(0,1e-1,1e3);
%lambda = 1e-3;
cvMSE = zeros(size(lambda));
ridgeEye = eye(nPred+1);
ridgeEye(1,1) = 0;
for i = 1:length(lambda)
    xRidge = [XTrain; (lambda(i)*ridgeEye)];
    yRidge = [yTrain; zeros(nPred+1,1)];
    XRidge = [ones(size(xRidge(:,1)))  xRidge];
    b(:,i) = lsqr(XRidge,yRidge);
    
    yFit = @(XRidge,yRidge,XTest)(XTest*b(:,i)); % Anropa crossval
    cvMSE(i) = crossval('MSE',XRidge,yRidge,'predfun',yFit);
end
[junk, idx] = min(cvMSE);

% Solve with lsqr
XRidge = [XTrain; (lambda(idx)*ridgeEye)];
yRidge = [yTrain; zeros(nPred+1,1)];
b = lsqr(XRidge,yRidge);

yHat = XTrain*b;

end