function bOut = RidgeRegress(yTrain, xTrain, lambda)

[rowX, colX] = size(xTrain);
[rowy, coly] = size(yTrain);
lambdaLength = length(lambda);

ridgeEye = [];
for i = 1:lambdaLength
ridgeEye = blkdiag(ridgeEye, lambda(i)*eye(colX));
end

[U, D, V] = svd(xTrain);
diagD = diag(repmat(diag(D'*D),lambdaLength,1));

b = reshape((diagD + ridgeEye)\repmat(D'*U'*yTrain,lambdaLength,1), colX, coly*lambdaLength);
bOut = V*b; 

end