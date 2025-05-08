function Kns = Circular_Kns(Ki, numsats, mu)
%Similar to KnsfromK, but explicitly for circular orbits
    [~,~,~,~,~,f0i] = unpackKeplerian(Ki);
    Kns(:,numsats) = zeros(6,1);

    for n = 1:numsats
        Kn = Ki;
        spacing = 2*pi/numsats;

        % Orbit is circular, let's stagger the true anomaly (= mean anomaly
        % in this case)
        Kn(6) = f0i + (n-1)*spacing;
        Kns(:,n) = Kn;
    end
end