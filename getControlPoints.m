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