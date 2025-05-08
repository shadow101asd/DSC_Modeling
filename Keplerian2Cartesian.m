function y = Keplerian2Cartesian(a, e, i, Om, w, f0, mu)
    % Convert Keplerian Orbital Elements to Cartesian Coordinates
    if size(e) > 1 % Does this matter?
        enorm = (e(1,:).^2 + e(2,:).^2 + e(3,:).^2).^0.5;
    else
        enorm = e;
    end

    p = a.*(1-enorm.^2);
    [~, n] = size(a);
    % Special Cases
    if norm(enorm) < 1e-6 * sqrt(n) % if circular
        w = 0*w;
        if norm(i) < 1e-6 * sqrt(n) % if equatorial
            %f0 = acos(r(1, :)./R); % lambda_true
            Om = Om * 0;
        else
            %f0 = acos(dot(N,r)./(Nnorm.*R)); % u
        end
    else
        if norm(i) < 1e-6 * sqrt(n) % if equatorial
            %f0 = acos(r(1, :)./R); % lambda_true
            Om = Om * 0;
        end
    end
    
    % Storing Variables for Computational Efficiency
    cosnu = cos(f0);
    sinnu = sin(f0);
    rootmup = sqrt(mu.*p.^-1);
    rPQW = [p.*cosnu./(1+ enorm.*cosnu); p.*sinnu./(1+ enorm.*cosnu); zeros(1,n)];
    vPQW = [-rootmup.*sinnu; rootmup.*(enorm + cosnu); zeros(1,n)];
    
    % Computing IJK/PQW
    IJKOverPQW1 = cos(Om).*cos(w) - sin(Om).*sin(w).*cos(i);
    IJKOverPQW2 = sin(Om).*cos(w) + cos(Om).*sin(w).*cos(i);
    IJKOverPQW3 = sin(w).*sin(i);
    IJKOverPQW4 = -cos(Om).*sin(w) - sin(Om).*cos(w).*cos(i);
    IJKOverPQW5 = -sin(Om).*sin(w) + cos(Om).*cos(w).*cos(i);
    IJKOverPQW6 = cos(w).*sin(i);
    IJKOverPQW7 = sin(Om).*sin(i);
    IJKOverPQW8 = -cos(Om).*sin(i);
    IJKOverPQW9 = cos(i);
    
    for j = 1:n
        IJKOverPQW = [IJKOverPQW1(1, j), IJKOverPQW4(1, j), IJKOverPQW7(1, j);
                        IJKOverPQW2(1, j), IJKOverPQW5(1, j), IJKOverPQW8(1, j); 
                        IJKOverPQW3(1, j), IJKOverPQW6(1, j), IJKOverPQW9(1, j)];
        rPQW(:,j) = IJKOverPQW * rPQW(:,j);
        vPQW(:,j) = IJKOverPQW * vPQW(:,j);
    end

    y = [rPQW;vPQW];
end

