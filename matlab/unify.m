%% Testing

coeff = -0.7;
n = 20;
W = [0.5;0.5];
AttriNum = 2;
M = [1 -3 3 -1; 0, 3, -6, 3; 0, 0, 3, -3; 0, 0, 0, 1];

[u, ~] = random2dAttributeValues(coeff, n);
X = u';
size(X)

alpha = [1, 1]';
if coeff < 0
    alpha(2) = -1;
end

tic
[~, controlPoints, score] = PCA_Ranking(X,alpha, W, AttriNum);
toc

[~, pc, s] = nearestPointsOnCurve(X, controlPoints, M, W);
s

function [u, curcoeff] = random2dAttributeValues(coeff, n)

    cur_diff = 0.1;
    while cur_diff > 0.0001
        u = copularnd('Gaussian',coeff,n);
        xx = corrcoef(u);
        curcoeff = xx(1,2);
%         curcoeff = corr(u(:,1), u(:,2), 'Type', 'Spearman');
        cur_diff = abs(curcoeff - coeff);
    end
end

function  [mse, Pt, St] = PCA_Ranking(X, alpha, W, attriNum)

    [mse, Pt, St, optFlag] = PCAlgorithm(X, alpha, W, attriNum);
    while optFlag == 0
        disp('re-sorting')
        [mse, Pt, St, optFlag] = PCAlgorithm(X, alpha, W, attriNum);
    end
end

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
end

function [msedist, pc, s] = nearestPointsOnCurve(X, P, M, W)
    t = 0:0.0001:1;
    Z = [ones(1,length(t)); t; t.^2; t.^3];
    bzc = P * M * Z;
    distsq = zeros(length(t),size(X,2));
    for i = 1 : size(X,1)
        idelta = bsxfun(@minus,bzc(i,:)',X(i,:));
%         size(idelta)
        distsq = distsq + W(i) * idelta.^2;
    end
    size(distsq)
    [dists, k] = min(distsq,[],1);
    pc = bzc(:,k)';
    s = t(k);
    msedist = sum(dists);
end

function p = getControlPoints(Zt, M, X, P00, P03, W, AttriNum)
    cvx_begin quiet
        variable p(4,AttriNum) 
        minimize( sum((Zt'*M'*p -X').^2*W)) 
        subject to
            0 <= p <= 1
            [1,0,0,0]*p == P00'
            [0,0,0,1]*p == P03'
    cvx_end
    p = p';
end




