function Dmax = bottleneckDistance(graph,path)
%BOTTLENECKDISTANCE Summary of this function goes here
%   Detailed explanation goes here
    [~,n] = size(path);
    n=n-1;
    Dmax = 0.0;
    for i = 1:n
        edgelength = graph.Edges.Weight(findedge(graph,path(i),path(i+1)));
        Dmax = max(Dmax,edgelength);
    end
end

