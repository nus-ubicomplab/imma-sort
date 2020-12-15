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