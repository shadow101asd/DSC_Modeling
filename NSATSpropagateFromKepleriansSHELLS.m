function XSats = NSATSpropagateFromKepleriansSHELLS(Ki,mu,etR,shells,es,cartwheels)
%NSATSPROPAGATEFROMKEPLERIANS 
    [~, nT] = size(etR);
    numsats_shell = sum(shells);
    XSats(:,:,numsats_shell) = zeros(6,nT);
    Kns = [];
    
    for i = 1:length(shells)
        shellnum = shells(i);
        K = [Ki(1); es(i); Ki(3:6)];
        Kni = KnsfromK_Inverted(K,shellnum,mu);
        Kns = [Kns, Kni];
    end
    
    if cartwheels > 1
        for k = 2:cartwheels % Add cartwheel symmetricals from the Earth-centric cartwheel
            for n = 1:numsats_shell
                Kns2 = [Kns(1:4,n);Kns(5,n)+2*pi*(k-1)/(cartwheels);Kns(6,n)];
                Kns = [Kns,Kns2];
            end
        end
    end

    for n = 1:(numsats_shell*cartwheels)
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end

    
end

