function [a, enorm, i, Om, w, f0] = Cartesian2Keplerian(X, mu)
%CARTESIAN2KEPLERIAN: For now, only return eccentricity norm!!
    % Unpack X
    r = X(1:3,:);
    v = X(4:6,:);

    n = size(r(1,:));
    n = n(1,2);
    % Vallado algorithm
    h = cross(r, v);
    H = (h(1,:).^2 + h(2,:).^2 + h(3,:).^2).^0.5;
    K = [zeros(2,n); ones(1,n)];
    N = cross(K, h);
    Nnorm = (N(1,:).^2 + N(2,:).^2 + N(3,:).^2).^0.5;
    R = (r(1,:).^2 + r(2,:).^2 + r(3,:).^2).^0.5;
    e = 1/mu * cross(v, h) - r./R; % Eccentricity
    enorm = (e(1,:).^2 + e(2,:).^2 + e(3,:).^2).^0.5;

    V = (v(1,:).^2 + v(2,:).^2 + v(3,:).^2).^0.5;
    E = 0.5*V.^2 - mu*R.^-1;

    if norm(enorm-1) > 1e-6 * sqrt(n)
        a = -mu./(2*E);
        p = a.*(ones(1,n)-e.^2);
    else
        a = inf;
        p = h.^2 ./mu;
    end

    i = acos(h(3,:)./H);
    Om = acos(N(1,:)./Nnorm);
    w = acos(dot(N,e)./(Nnorm.*enorm));
    f0 = acos(dot(e,r)./(enorm.*R));

    % Check conditions for angle Keplerian elements

    for j = 1:n
        if dot(r(:,j),v(:,j)) < 0
            f0(1,j) = 2*pi - f0(1,j);
        end
        if N(2,j) < 0
            Om(1,j) = 2*pi - Om(1,j);
        end
        if e(3,j) < 0
            w(1,j) = 2*pi - w(1,j);
        end
    end
    % Special Cases
    if norm(enorm) < 1e-6 * sqrt(n) % if circular
        if norm(i) < 1e-6 * sqrt(n) % if equatorial
            f0 = acos(r(1, :)./R); % lambda_true
        else
            f0 = acos(dot(N,r)./(Nnorm.*R)); % u
        end
    else
        if norm(i) < 1e-6 * sqrt(n) % if equatorial
            f0 = acos(r(1, :)./R); % lambda_true
        end
    end
end

