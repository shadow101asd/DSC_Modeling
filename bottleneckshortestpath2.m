function [path, D, edgepath] = bottleneckshortestpath2(g, nodename1, nodename2)
%BOTTLENECKSHORTESTPATH
    T = minspantree(g, Method="dense", Root=nodename1);
    [path1, ~, ~] = shortestpath(T,nodename1,nodename2);

    % Return and check for more direct, roughly equivalent bottleneck paths
    margin = 1.05;
    D_Bottleneck = bottleneckDistance(g,path1);
    DMax = margin*D_Bottleneck;
    
    A2 = adjacency(g, "weighted");
    A2(A2>DMax) = 0.0; % Prune graph
    graph2 = graph(A2, g.Nodes);

    [path,D,edgepath] = shortestpath(graph2,nodename1,nodename2);
end

