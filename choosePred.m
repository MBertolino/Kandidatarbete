function [bestPred] = choosePred(x, y, nPred, variation)

[rowx, colx] = size(x);

% Find subset of predictors: bestPred
cvMSE = zeros(length(colx),1);
bestPred = [];
for ip = 1:nPred
    for in = 1:colx
        if ~sum(ismember(bestPred,in)) % Don't reuse predictor
            X = [ones(size(x(:,1))) x(:,[bestPred in])];
            
            % Which variation to use
            if strcmpi(variation,'Normal')
                b = regress(y,X);
                
            elseif strcmpi(variation,'Ridge')
                lambda = 1e-3;
                [rowPred, colPred] = size([bestPred in]);
                XRidge = [X; (lambda*eye(colPred+1))];
                yRidge = [y; zeros(colPred+1,1)];
                b = lsqr(XRidge,yRidge);
                
            else
                % If not Normal or Ridge - Wrong input
                ME = MException('MyComponent:noSuchVariable', ...
                    'Input %s is not valid',variation);
                throw(ME)
            end
            
            yFit = @(X,y,XTest)(XTest*b);
            cvMSE(in) = crossval('MSE',X,y,'predfun',yFit);
        end
    end
    [junk, nextPred] = min(cvMSE);
    bestPred = [bestPred nextPred];
end

end