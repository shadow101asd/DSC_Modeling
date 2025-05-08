function meanDistance2 = computeLinkBudget(X1,X2)
    %COMPUTELINKBUDGET Summary of this function goes here
    %   Detailed explanation goes here
    meanDistance2 = mean((X1(1,:)-X2(1,:)).^2 + (X1(2,:)-X2(2,:)).^2 + (X1(3,:)-X2(3,:)).^2);
end

