function [mse, Pt, St,optFlag] = PCAlgorithm(X,alpha, W, attriNum)

    optFlag = 1;
    P00 = 0.5*(1 - alpha);
    P03 = 0.5*(1 + alpha);  
    Pt = [P00, rand(attriNum,1), rand(attriNum,1), P03];
    M = [1 -3 3 -1; 0, 3, -6, 3; 0, 0, 3, -3; 0, 0, 0, 1];
    [lastmse, ~, St] = nearestPointsOnCurve(X, Pt, M,W);
    deltaM = 10000;
    
    while deltaM > 0.00001        
        Zt = [ones(1, size(X,2)); St; St.^2; St.^3];
        Pt = getControlPoints(Zt, M, X, P00, P03, W,attriNum);
        [mse, ~, St] = nearestPointsOnCurve(X, Pt, M, W);       
        deltaM =  lastmse - mse;
        if deltaM < 0
            optFlag = 0;
            break;
        end
        lastmse = mse; 
    end    
