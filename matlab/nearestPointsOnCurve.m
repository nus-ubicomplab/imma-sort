function [msedist, pc, s] = nearestPointsOnCurve(X, P, M, W)
    t = 0:0.0001:1;
    Z = [ones(1,length(t)); t; t.^2; t.^3];
    bzc = P * M * Z;
    distsq = zeros(length(t),size(X,2));
    for i = 1 : size(X,1)
        idelta = bsxfun(@minus,bzc(i,:)',X(i,:));
        distsq = distsq + W(i) * idelta.^2;
    end
    [dists, k] = min(distsq,[],1);
    pc = bzc(:,k)';
    s = t(k);
    msedist = sum(dists);
