function [u, curcoeff] = randomMoreThan2AttributeValues(coeff, n, attriNum)

% Rho = ones(attriNum,attriNum) * coeff *(-1);
% Rho(1:(attriNum+1):end) = 1;
% Rho(1,2) = coeff;
% Rho(2,1) = coeff;
% top2 = -1;
% 
% while(top2 < (coeff-0.05) || top2 >(coeff+0.05))
% 
% Z = mvnrnd(zeros(1, attriNum), Rho, n);
% u = normcdf(Z,0,1);
% curcoeff = corr(u, 'Type', 'Spearman');
% top2 = curcoeff(1,2);
% 
% end

Rho = ones(attriNum,attriNum) * coeff;
Rho(1:(attriNum+1):end) = 1;
% Rho(1,2) = coeff;
% Rho(2,1) = coeff;
top2 = -1;

while(top2 < (coeff-0.05) || top2 >(coeff+0.05))

Z = mvnrnd(zeros(1, attriNum), Rho, n);
u = normcdf(Z,0,1);
curcoeff = corr(u, 'Type', 'Spearman');
top2 = curcoeff(1,2);

end