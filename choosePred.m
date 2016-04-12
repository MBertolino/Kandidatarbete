function [bestPred] = choosePred(x, y, nPred)

[rowx, colx] = size(x);

% Find subset of predictors: bestPred
cvMSE = zeros(length(colx),1);
bestPred = [];
for ip = 1:nPred
    for in = 1:colx
        if ~sum(ismember(bestPred,in)) % Don't reuse predictor
            X = [ones(size(x(:,1))) x(:,[bestPred in])];
            
            b = regress(y,X);          % Want to use different methods
            yFit = @(X,y,XTest)(XTest*b);
            cvMSE(in) = crossval('MSE',X,y,'predfun',yFit);
        end
    end
    [junk, nextPred] = min(cvMSE);
    bestPred = [bestPred nextPred];
end

end