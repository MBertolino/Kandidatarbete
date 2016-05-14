function bOut = RidgeRegress(yTrain, xTrain, lambda)

[rowX, colX] = size(xTrain);
[rowy, coly] = size(yTrain);
lambdaLength = length(lambda);

[U, D, V] = svd(xTrain);

% Resize to match number of tuning parameters lambda'
ridgeEye = [];
for i = 1:lambdaLength
    ridgeEye = blkdiag(ridgeEye, lambda(i)*eye(colX));
end

diagD = diag(repmat(diag(D'*D),lambdaLength,1));

% Regress for all tuning parameters and dependent variables at once
b = (diagD + ridgeEye)\repmat(D'*U'*yTrain,lambdaLength,1);
b = reshape(b, colX, coly*lambdaLength);
bOut = V*b;

end