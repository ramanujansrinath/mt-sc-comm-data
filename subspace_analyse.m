function out = subspace_analyse_online(targ,src)
    % RRR
    [op,y,cvLoss,cvLossFull] = regressAndPlot(targ',src');
    out.optDim = op;
    out.perf = y';
    out.optPerf = y(op);
    out.cvLoss = cvLoss;
    out.cvLossFull = cvLossFull;
    
    % targ fa
    [outEigVal,outDim] = getFA_eig_dim(targ');
    out.outEigVal = outEigVal;
    out.outDim = outDim;
    
    % src fa
    [outEigVal,outDim] = getFA_eig_dim(src');
    out.srcEigVal = outEigVal;
    out.srcDim = outDim;
    
    % rsc
    [rsc1,rsc2,rscS] = getSharedRsc(src',targ');
    out.rsc = rsc2;
    out.srcRsc = rsc1;
    out.sharedRsc = rscS;
    
    % rate
    out.targRate = nanmean(targ,2);
    out.srcRate = nanmean(src,2);
end

function [op,y,cvLoss,cvLossFull] = regressAndPlot(Y,X)
    [cvLoss,op] = mtscRegress(Y,X);
    cvLossFull = mtscRegress_full(Y,X);
    y = 1-cvLoss(1,:);
    cvLoss = cvLoss(2,:);
end

function [cvLoss,optDim,B] = mtscRegress(Y,X) % mostly from semedo
    %%% Cross-validate Reduced Rank Regression

    % Vector containing the interaction dimensionalities to use when fitting
    % RRR. 0 predictive dimensions results in using the mean for prediction.
    numDimsUsedForPrediction = 1:min([10 size(Y,2)]);

    % Number of cross validation folds.
    cvNumFolds = min([10 size(Y,2)]);

    % Initialize default options for cross-validation.
    cvOptions = statset('crossval');

    % If the MATLAB parallel toolbox is available, uncomment this line to
    % enable parallel cross-validation.
    % cvOptions.UseParallel = true;

    % Regression method to be used.
    regressMethod = @ReducedRankRegress;

    % Auxiliary function to be used within the cross-validation routine (type
    % 'help crossval' for more information). Briefly, it takes as input the
    % the train and test sets, fits the model to the train set and uses it to
    % predict the test set, reporting the model's test performance. Here we
    % use NSE (Normalized Squared Error) as the performance metric. MSE (Mean
    % Squared Error) is also available.
    
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
        numDimsUsedForPrediction, 'LossMeasure', 'NSE');
    
    % do ridge
    % dMaxShrink = .5:.01:1;
    % lambda = GetRidgeLambda(dMaxShrink, X);
    % cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
    %     (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
    %     'LossMeasure', 'NSE');

    % Cross-validation routine.
    cvl = crossval(cvFun, Y, X, ...
          'KFold', cvNumFolds, ...
        'Options', cvOptions);

    % Stores cross-validation results: mean loss and standard error of the
    % mean across folds.
    cvLoss = [ mean(cvl); std(cvl)/sqrt(cvNumFolds) ];

    % To compute the optimal dimensionality for the regression model, call
    % ModelSelect:
    optDim = ModelSelect...
        (cvLoss, numDimsUsedForPrediction);
end

function cvLoss = mtscRegress_full(Y,X)
    cvNumFolds = min([10 size(Y,2)]);
    cvOptions = statset('crossval');
    % do ridge
    dMaxShrink = .5:.01:1;
    lambda = GetRidgeLambda(dMaxShrink, X);
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
        'LossMeasure', 'NSE');

    % Cross-validation routine.
    cvl = crossval(cvFun, Y, X, ...
          'KFold', cvNumFolds, ...
        'Options', cvOptions);

    % Stores cross-validation results: mean loss and standard error of the
    % mean across folds.
    cvLoss = [ mean(cvl); std(cvl)/sqrt(cvNumFolds) ];

    lambSel = ModelSelect(cvLoss, lambda);
    cvLoss = cvLoss(:,lambda==lambSel);
end

function [eigVals,outDim] = getFA_eig_dim(Y)
    nFolds = 2;
    nDims = size(Y,2); % min(5,size(Y,2));
    dim = crossvalidate_fa(Y','zDimList',nDims,'showPlots',false,'numFolds',nFolds);
    L=dim.estParams.L;
    LL=L*L';
    [~,D]=eig(LL);
    la=diag(D);
    la=sort(la,'descend');
    eigVals=la(1:nDims);
    [~,outDim] = getDim(eigVals);
end

function [rsc1,rsc2,rscS] = getSharedRsc(Y1,Y2)
    idx = [ones(size(Y1,2)); 2*ones(size(Y2,2),size(Y1,2))];
    idx = [idx [3*ones(size(Y1,2),size(Y2,2)); 4*ones(size(Y2,2))]];

    Y = [Y1 Y2];
    rsc = corrcoef(Y,'rows','pairwise');
    rsc = tril(rsc,-1);
    
    c1 = rsc(idx == 1);
    c2 = rsc(idx == 4);
    cs = rsc(idx == 2);
    
    rsc1 = (c1(c1~=0));
    rsc2 = (c2(c2~=0));
    rscS = (cs(cs~=0));
end

function [vEx,outDim] = getDim(eigVals)
    vEx = eigVals./repmat(sum(eigVals),size(eigVals,1),1);
    
    y = cumsum(vEx);
    ind2 = find(y>0.95,1,'first');
    if ind2 == 1
        outDim = 1;
    else
        ind1 = ind2-1;
        c = (y(ind1)*ind2 - ind1*y(ind2));
        m = y(ind2)-y(ind1);

        outDim = (0.95-c)/m;
    end
end
