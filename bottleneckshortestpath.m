function [path, D, edgepath] = bottleneckshortestpath(graph, nodename1, nodename2)
%BOTTLENECKSHORTESTPATH
    T = minspantree(graph, Method="dense", Root=nodename1);
    [path, D, edgepath] = shortestpath(T,nodename1,nodename2);
end

