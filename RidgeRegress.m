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

% Solve with lsqr
% ridgeEye(1,1) = 0;
% Ridge data
% XRidge = [XTrain; (lambda*ridgeEye)];
% yRidge = [yTrain; zeros(col,1)];
% b = lsqr(XRidge,yRidge);

% H = V*(D'*D + lambda*ridgeEye)*V';
% H = (H+H')/2;           % Why necessary
% f = -V*D'*U'*yTrain;
% b = quadprog(H,f);

end