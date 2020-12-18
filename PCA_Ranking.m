function  [mse, Pt, St] = PCA_Ranking(X, alpha, W, attriNum)

    [mse, Pt, St, optFlag] = PCAlgorithm(X, alpha, W, attriNum);
    while optFlag == 0
        disp('re-sorting')
        [mse, Pt, St, optFlag] = PCAlgorithm(X, alpha, W, attriNum);
    end

