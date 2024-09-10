function [EBC_pos,EBC_neg,BC_pos,BC_neg] = edge_betweenness_wei_sign(G)
% Edited version of Rubinov's function (01/25/16 - Spielberg)
%       Now computed for positive and negative weights separately
%
%
%EDGE_BETWEENNESS_WEI_SIGN    Edge betweenness centrality
%
%   EBC_pos = edge_betweenness_wei_sign(G);
%   [EBC_pos,EBC_neg,BC_pos,BC_neg] = edge_betweenness_wei_sign(G);
%
%   Edge betweenness centrality is the fraction of all shortest paths in 
%   the network that contain a given edge. Edges with high values of 
%   betweenness centrality participate in a large number of shortest paths.
%
%   Input:      G,      Directed/undirected connection-length matrix.
%
%   Output:     EBC_pos = edge betweenness centrality matrix for pos weights
%               EBC_neg = edge betweenness centrality matrix for neg weights
%               BC_pos  = nodal betweenness centrality vector for pos weights
%               BC_neg  = nodal betweenness centrality vector for neg weights
%
%   Notes:
%       The input matrix must be a connection-length matrix, typically
%   obtained via a mapping from weight to length. For instance, in a
%   weighted correlation network higher correlations are more naturally
%   interpreted as shorter distances and the input matrix should
%   consequently be some inverse of the connectivity matrix. 
%       Betweenness centrality may be normalised to the range [0,1] as
%   BC/[(N-1)(N-2)], where N is the number of nodes in the network.
%
%   Reference: Brandes (2001) J Math Sociol 25:163-177.
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2012


n=length(G);
% E=find(G); G(E)=1./G(E);        %invert weights
BC_pos           = zeros(n,1);                 %vertex betweenness
BC_neg           = zeros(n,1);                 %vertex betweenness
EBC_pos          = zeros(n);                   %edge betweenness
EBC_neg          = zeros(n);                   %edge betweenness
G_pos            = threshold_absolute(G,0);
G_pos(1:n+1:end) = 0;
G_neg            = threshold_absolute(G*-1,0);
G_neg(1:n+1:end) = 0;

% For positive weights
for u = 1:n
    D_pos     = inf(1,n);
    D_pos(u)  = 0;          %distance from u
    NP_pos    = zeros(1,n);
    NP_pos(u) = 1;          %number of paths from u
    S_pos     = true(1,n);  %distance permanence (true is temporary)
    P_pos     = false(n);   %predecessors
    Q_pos     = zeros(1,n);
    q_pos     = n;          %order of non-increasing distance
    G1_pos    = G_pos;
    V_pos     = u;
    while 1
        S_pos(V_pos)    = 0; %distance u->V is now permanent
        G1_pos(:,V_pos) = 0; %no in-edges as already shortest
        for v = V_pos
            Q_pos(q_pos) = v;
            q_pos        = q_pos-1;
            W_pos        = find(G1_pos(v,:)); %neighbours of v
            for w = W_pos
                Duw_pos = D_pos(v)+G1_pos(v,w);       %path length to be tested
                if Duw_pos<D_pos(w)                   %if new u->w shorter than old
                    D_pos(w)   = Duw_pos;
                    NP_pos(w)  = NP_pos(v);           %NP(u->w) = NP of new path
                    P_pos(w,:) = 0;
                    P_pos(w,v) = 1;                   %v is the only predecessor
                elseif Duw_pos==D_pos(w)              %if new u->w equal to old
                    NP_pos(w)  = NP_pos(w)+NP_pos(v); %NP(u->w) sum of old and new
                    P_pos(w,v) = 1;                   %v is also a predecessor
                end
            end
        end

        minD_pos = min(D_pos(S_pos));
        if isempty(minD_pos)
            break %all nodes reached, or
        elseif isinf(minD_pos) %...some cannot be reached:
            Q_pos(1:q_pos) = find(isinf(D_pos));
            break %...these are first-in-line
        end
        V_pos = find(D_pos==minD_pos);
    end

    DP_pos = zeros(n,1);                          %dependency
    for w = Q_pos(1:n-1)
        BC_pos(w) = BC_pos(w)+DP_pos(w);
        for v = find(P_pos(w,:))
            DPvw_pos     = (1+DP_pos(w)).*NP_pos(v)./NP_pos(w);
            DP_pos(v)    = DP_pos(v)+DPvw_pos;
            EBC_pos(v,w) = EBC_pos(v,w)+DPvw_pos;
        end
    end
end

% For negative weights
for u = 1:n
    D_neg     = inf(1,n);
    D_neg(u)  = 0;          %distance from u
    NP_neg    = zeros(1,n);
    NP_neg(u) = 1;          %number of paths from u
    S_neg     = true(1,n);  %distance permanence (true is temporary)
    P_neg     = false(n);   %predecessors
    Q_neg     = zeros(1,n);
    q_neg     = n;          %order of non-increasing distance
    G1_neg    = G_neg;
    V_neg     = u;
    while 1
        S_neg(V_neg)    = 0; %distance u->V is now permanent
        G1_neg(:,V_neg) = 0; %no in-edges as already shortest
        for v = V_neg
            Q_neg(q_neg) = v;
            q_neg        = q_neg-1;
            W_neg        = find(G1_neg(v,:)); %neighbours of v
            for w = W_neg
                Duw_neg = D_neg(v)+G1_neg(v,w);       %path length to be tested
                if Duw_neg<D_neg(w)                   %if new u->w shorter than old
                    D_neg(w)   = Duw_neg;
                    NP_neg(w)  = NP_neg(v);           %NP(u->w) = NP of new path
                    P_neg(w,:) = 0;
                    P_neg(w,v) = 1;                   %v is the only predecessor
                elseif Duw_neg==D_neg(w)              %if new u->w equal to old
                    NP_neg(w)  = NP_neg(w)+NP_neg(v); %NP(u->w) sum of old and new
                    P_neg(w,v) = 1;                   %v is also a predecessor
                end
            end
        end

        minD_neg = min(D_neg(S_neg));
        if isempty(minD_neg)
            break %all nodes reached, or
        elseif isinf(minD_neg) %...some cannot be reached:
            Q_neg(1:q_neg) = find(isinf(D_neg));
            break %...these are first-in-line
        end
        V_neg = find(D_neg==minD_neg);
    end

    DP_neg = zeros(n,1);                          %dependency
    for w = Q_neg(1:n-1)
        BC_neg(w) = BC_neg(w)+DP_neg(w);
        for v = find(P_neg(w,:))
            DPvw_neg     = (1+DP_neg(w)).*NP_neg(v)./NP_neg(w);
            DP_neg(v)    = DP_neg(v)+DPvw_neg;
            EBC_neg(v,w) = EBC_neg(v,w)+DPvw_neg;
        end
    end
end