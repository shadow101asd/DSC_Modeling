function dy = twobodyode(t, y, GM)
% Two body problem with one mass much larger than the other.
r = sqrt(y(1)^2 + y(2)^2+ y(3)^2);
dy = [y(4:6); -GM*y(1:3)./r^3];
end