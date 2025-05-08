function Kns = KnsfromK_Inverted(Ki,numsats,mu)
%KNSFROMK 
    [a,e,i,Om,w,f0i] = unpackKeplerian(Ki);
    Kns(:,numsats) = zeros(6,1);
    % Assume orbital symmetry around the Sun for now
    for n = 1:numsats
        Kn = Ki;
        spacing = 2*pi/numsats;
        if e > 1e-4 % Orbit non-circular, change argument of periapsis w
            Kn(5) = w + spacing*n;
        end

        % Regardless, stagger mean anomaly linearly out in time using
        % updateTrueAnomaly
        T = 2*pi*sqrt(a^3/mu); % Orbital Period
        dt = T-T*n/numsats;
        Kn(6) = updateTrueAnomaly(a,e,i,Om,w,f0i,mu,dt);

        Kns(:,n) = Kn;
    end
end

