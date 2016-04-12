function [bestPred] = choosePred(x, y, nPred)

[rowx, colx] = size(x);

cvMSE = zeros(length(colx),1);     % Crossval Mean Squared Error
bestPred = [];
for ip = 1:nPred                   % Vi vill ha nPred stycken predictors
    for in = 1:colx                % Vi har colx stycken att välja på
        if ~sum(ismember(bestPred,in))     % Se till att vi inte åter%använder en predictor
            S = sqrt(sum(x(:,in).^2)/rowx);  % Normera%
            x(:,in) = x(:,in)./S;
            X = [ones(size(x(:,1))) x(:,[bestPred in])];
            
            b = regress(y,X);          % koeff, (OBS regress kanske är fel)
            yFit = @(X,y,XTest)(XTest*b); % Anropa crossval
            cvMSE(in) = crossval('MSE',X,y,'predfun',yFit);
        end
    end
    [junk, nextPred] = min(cvMSE);  % Bäst predictordelmängd är den med lägst
    bestPred = [bestPred nextPred]; % Lägg till den bästa
end
bestPred = bestPred;            % Vi vill inte ha dependent market med

end