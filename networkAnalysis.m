function [graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis(X1,X2,name1,name2,XSats)
    % Network Analysis
    
    % Create Adjacency Matrices
    Xs = collateXs(X1,X2,XSats);
    D12 = distanceBetweenXs(X1,X2);
    [~,nT,numsats] = size(XSats);
    
    A = createAdjacencyMatrix(Xs,@distanceBetweenXs,D12);
    
    % Create graphs from adjacency matrices
    
    names(numsats+1) = name1;
    names(numsats+2) = name2;
    for j=1:numsats
        names(j) = strcat("Sat", num2str(j));
    end
    graphslist{nT} = graph(A(:,:,nT),names);
    
    numedges = zeros(nT,1);
    
    for t = 1:nT
        g = graph(A(:,:,t),names);
        graphslist{t} = g;
        numedges(t) = g.numedges;
    end
    
    % Compute bottleneck shortest path for each graph
    paths{nT} = {};
    Dpaths = zeros(nT,1);
    Dequivs = zeros(nT,1);
    for t = 1:nT
        g = graphslist{t};
        [path, Dpath, edgepath] = bottleneckshortestpath2(g, name1, name2);
        paths{t} = path;
        Dpaths(t) = Dpath;
        Dbottleneck = bottleneckDistance(g,path);
        Dequivs(t) = Dbottleneck;
    end
end

