function [Xs] = propagateFromKeplerians(Ki,mu,etR)
    %PROPAGATEFROMKEPLERIANS Returns Cartesian Trajectories from initial
    %orbital elements. All orbital elements are assumed constant except for the
    %true anomaly.

    [~, nT] = size(etR);
    dt = (etR(nT)-etR(1))/(nT-1);
    [a,e,i,Om,w,f0i] = unpackKeplerian(Ki);

    Xs(:,1) = Keplerian2Cartesian(a,e,i,Om,w,f0i,mu);
    Xs(:,nT) = Xs(:,1); % preallocating
    
    f0last = f0i;
    for n = 2:nT
        f0new = updateTrueAnomaly(a, e, i, Om, w, f0last, mu, dt);
        Xs(:,n) = Keplerian2Cartesian(a,e,i,Om,w,f0new,mu);
        f0last = f0new;
    end
end

