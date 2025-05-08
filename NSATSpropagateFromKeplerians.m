function XSats = NSATSpropagateFromKeplerians(Ki,mu,etR,numsats)
%NSATSPROPAGATEFROMKEPLERIANS 
    [~, nT] = size(etR);
    XSats(:,:,numsats) = zeros(6,nT);
    Kns = KnsfromK(Ki,numsats,mu);

    for n = 1:numsats
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end
end

