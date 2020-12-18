%% Testing

coeff = -0.7;
n = 20;
W = [0.5;0.5];
AttriNum = 2;
M = [1 -3 3 -1; 0, 3, -6, 3; 0, 0, 3, -3; 0, 0, 0, 1];

[u, ~] = random2dAttributeValues(coeff, n);
X = u';

alpha = [1, 1]';
if coeff < 0
    alpha(2) = -1;
end

tic
[~, controlPoints, score] = PCA_Ranking(X,alpha, W, AttriNum);
toc

[~, pc, s] = nearestPointsOnCurve(X, controlPoints, M, W);