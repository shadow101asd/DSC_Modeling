function f0new = updateTrueAnomaly(a, e, i, Om, w, f0, mu, dt)
%UPDATETRUEANOMALY Summary of this function goes here
%   Detailed explanation goes here

% Compute current eccentric anomaly

Ecur = 2*atan(sqrt((1-e)/(1+e))*tan(f0/2));

% Compute Current Mean Anomaly from Kepler's Formula

Mecur = Ecur - e*sin(Ecur);

% Update Mean Anomaly

n = sqrt(mu/a^3); % Mean motion
Meupdated = Mecur + n*dt;

% Compute Updated Eccentric Anomaly from Meupdated
% Numerical method:

fun = @(E) E - e*sin(E) - Meupdated;
Eupdated = fzero(fun, Ecur); % Ecur used as initial guess

% Compute Updated True Anomaly from Eupdated

f0new = atan2((sin(Eupdated)*sqrt(1-e^2)),(cos(Eupdated)-e));
end

