function Ds = distanceBetweenXs(X1,X2)
%DISTANCEBETWEENXS Distances between points of the 6-vector format. If one
%of the Xs has three dimensions, enter it as X2.
    Ds = sqrt((X1(1,:)-X2(1,:,:)).^2 + (X1(2,:)-X2(2,:,:)).^2 + (X1(3,:)-X2(3,:,:)).^2);
end

