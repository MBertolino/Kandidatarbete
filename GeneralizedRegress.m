function [b, yHat, pred] = GeneralizedRegress(yTrain, xTrain, nPred)

% Choose the nPred best predictors
pred = choosePred(xTrain, yTrain, nPred, 'Ridge');
xTrain = xTrain(:,pred);
XTrain = [ones(size(xTrain(:,1))) xTrain];

yRidge = [yTrain; zeros(nPred+1,1)];



% Find the best lambda with cross validation
expon = -25:2;
lambda = [0 10.^expon];
%lambda = 1e-3;

cvMSE = zeros(size(lambda));
ridgeEye = eye(nPred+1);
ridgeEye(1,1) = 0;
for i = 1:length(lambda)
    xRidge = [XTrain; (lambda(i)*ridgeEye)];
    XRidge = [ones(size(xRidge(:,1))) xRidge];
    
    % Use Crossvalidation to determine optimal lambda
    K = 10;
    C = cvpartition(length(XRidge),'KFold',K);
    err = zeros(C.NumTestSets,1);
    for icv = 1:C.NumTestSets
        trIdx = C.training(icv);
        teIdx = C.test(icv);
        
        % Solve with quadprog
        Hicv = XRidge'*XRidge;
        ficv = -XRidge'*yRidge;
        bicv = quadprog(Hicv,ficv);
        ytest = XRidge(teIdx,:)*bicv;
        err(icv) = mean(abs(ytest-yRidge(teIdx)).^2);
    end
    cvMSE(i) = mean(err);
        
end
[junk, idx] = min(cvMSE);               % Optimal lambda index


% Solve with quadprog
XRidge = [XTrain; (lambda(idx)*ridgeEye)];
H = XRidge'*XRidge;
f = -XRidge'*yRidge;
% lb = zeros(length(pred)+1,1);          % Possible lower bound
% ub = Inf*ones(length(pred)+1,1);       % Possible upper bound
b = quadprog(H,f);
yHat = XTrain*b;

end