function b = RidgeRegress(yTrain, XTrain)

[row, col] = size(XTrain);

% Penalty factor
lambda = 1e1;
ridgeEye = eye(col);

% Solve with svd
[U, D, V] = svd(XTrain);
b = V*((D'*D + lambda*ridgeEye)\D'*U'*yTrain);



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