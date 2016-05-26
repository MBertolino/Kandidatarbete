function bOut = RidgeRegress(yTrain, xTrain, lambda, ridgeEye)
% RidgeRegress calculates regression coefficients.
% bOut = RidgeRegress(yTrain, xTrain, lambda, ridgeEye) calculates
% regression coefficients parallelly for multiple dependent variables
% with shared design matrix using Singular Value Decomposition.
%
% Input:
% yTrain (matrix)   - is a matrix with each column representing a dependent 
%                       variable and each row representing an observation.
% xTrain (matrix)   - is the shared design matrix with each column 
%                       representing a covariate and each row representing
%                       an observation.
% lambda (vector)   - is a vector of each tuning parameter
% ridgeEye (matrix) - is diagonal with size length(lambda) x col(yTrain)
%                       with each element of lambda repeated length(lambda)
%                       times along the diagonal
%
% Output:
% bOut (matrix)     - is the resulting regression coefficent vector for each
%                       dependent variable repeated for each lambda
%
% 2016 Iliam Barkino, Mattias Bertolino

[rowX, colX] = size(xTrain);
[rowy, coly] = size(yTrain);
lambdaLength = length(lambda);

[U, D, V] = svd(xTrain);

% Resize to match number of tuning parameters lambda
diagD = diag(repmat(diag(D'*D), lambdaLength, 1));

% Regress for all tuning parameters and dependent variables at once
b = (diagD + ridgeEye)\repmat(D'*U'*yTrain, lambdaLength, 1);
b = reshape(b, colX, coly*lambdaLength);
bOut = V*b;

end