function A = createAdjacencyMatrix(Xs,pairwiseFunc,cutoff)
%Create a time dependent adjacency matrix for a graph where all nodes are 
% connected, and edge weights are computed with a passed in pairwiseFunc

    [~,nT,N] = size(Xs);
    A = zeros(N,N,nT);
    for a = 1:N
        for b = 1:N
            if a == b
                A(a,b) = 0; % No self-loops
            else 
                if A(b,a,1) ~= 0.0
                    A(a,b,:) = A(b,a,:); % Avoids recomputing since A is symmetric
                else
                    A(a,b,:) = pairwiseFunc(Xs(:,:,a),Xs(:,:,b));
                end
            end
        end
    end
    % Remove edges at every time step with weight larger than cutoff of t
    for t = 1:nT
        At = A(:,:,t);
        At(At>cutoff(t)) = 0.0;
        A(:,:,t) = At;
    end
end

