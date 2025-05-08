function Xs = collateXs(X1,X2,XSats)
%COLLATEXS 
    [~,nT,numsats] = size(XSats);
    Xs = zeros(6,nT,numsats+2);
    
    Xs(:,:,1:numsats) = XSats;
    Xs(:,:,numsats+1) = X1;
    Xs(:,:,numsats+2) = X2;
end

