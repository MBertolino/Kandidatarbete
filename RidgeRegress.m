function bOut = RidgeRegress(yTrain, XTrain, lambda)

[rowX, colX] = size(XTrain);
[rowy, coly] = size(yTrain);
lambdaLength = length(lambda);

% Penalty factor
ridgeEye = eye(colX);

% Solve with svd
[U, D, V] = svd(XTrain);

bOut = [];
b = zeros(colX, coly);
for j = 1:lambdaLength
    for m = 1:coly
        b(:,m) = V*((D'*D + lambda(j)*ridgeEye)\D'*U'*yTrain(:,m));
    end
    bOut = [bOut b];
end

end