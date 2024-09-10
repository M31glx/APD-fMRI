function [lambda_pos,lambda_neg,efficiency_pos,efficiency_neg,ecc_pos,ecc_neg,radius_pos,radius_neg,diameter_pos,diameter_neg] = charpath_sign(L,diagonal_dist,infinite_dist)
%CHARPATH_SIGN  Characteristic path length, global efficiency and related statistics
%
%   [lambda_pos,lambda_neg]                                = charpath_sign(L);
%   [lambda_pos,lambda_neg,efficiency_pos,efficiency_neg]  = charpath_sign(L);
%   [lambda_pos,lambda_neg,efficiency_pos,efficiency_neg,ecc_pos,ecc_neg,radius_pos,radius_neg,diameter_pos,diameter_neg] = charpath_sign(L,diagonal_dist,infinite_dist);
%
%   The network characteristic path length is the average shortest path
%   length between all pairs of nodes in the network. The global efficiency
%   is the average inverse shortest path length in the network. The nodal
%   eccentricity is the maximal path length between a node and any other
%   node in the network. The radius is the minimal eccentricity, and the
%   diameter is the maximal eccentricity.
%
%   Input:      L,              length matrix (NOTE: different from the 
%                               original function, which needed a distance 
%                               matrix)
%               diagonal_dist   optional argument
%                               include distances on the main diagonal
%                                   (default: diagonal_dist=0)
%               infinite_dist   optional argument
%                               include infinite distances in calculation
%                                   (default: infinite_dist=1)
%
%   Outputs:    lambda,         network characteristic path length
%               efficiency,     network global efficiency
%               ecc,            nodal eccentricity
%               radius,         network radius
%               diameter,       network diameter
%
%   Notes:
%       The input distance matrix may be obtained with any of the distance
%   functions, e.g. distance_bin, distance_wei.
%       Characteristic path length is defined here as the mean shortest
%   path length between all pairs of nodes, for consistency with common
%   usage. Note that characteristic path length is also defined as the
%   median of the mean shortest path length from each node to all other
%   nodes.
%       Infinitely long paths (i.e. paths between disconnected nodes) are
%   included in computations by default. This behavior may be modified with
%   via the infinite_dist argument.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008
%   Mika Rubinov, U Cambridge, 2010/2015

%   Modification history
%   2002: original (OS)
%   2010: incorporation of global efficiency (MR)
%   2015: exclusion of diagonal weights by default (MR)
%   2016: inclusion of infinite distances by default (MR)
%   2016: computed for positive and negative weights separately

D_pos = distance_wei(threshold_absolute(L,0));
D_neg = distance_wei(threshold_absolute((L*-1),0));

n = size(L,1);
if any(any(isnan(D_pos))) || any(any(isnan(D_neg)))
    error('The distance matrix must not contain NaN values');
end
if ~exist('diagonal_dist','var') || ~diagonal_dist || isempty(diagonal_dist)
    D_pos(1:n+1:end) = NaN;     % set diagonal distance to NaN
    D_neg(1:n+1:end) = NaN;     % set diagonal distance to NaN
end
if  exist('infinite_dist','var') && ~infinite_dist
    D_pos(isinf(D_pos)) = NaN;  % ignore infinite path lengths
    D_neg(isinf(D_neg)) = NaN;  % ignore infinite path lengths
end

Dv_pos = D_pos(~isnan(D_pos));  % get non-NaN indices of D
Dv_neg = D_neg(~isnan(D_neg));  % get non-NaN indices of D

% Mean of entries of D(G)
lambda_pos = mean(Dv_pos);
lambda_neg = mean(Dv_neg);

% Efficiency: mean of inverse entries of D(G)
efficiency_pos = mean(1./Dv_pos);
efficiency_neg = mean(1./Dv_neg);

% Eccentricity for each vertex
ecc_pos        = nanmax(D_pos,[],2);
ecc_neg        = nanmax(D_neg,[],2);

% Radius of graph
radius_pos     = min(ecc_pos);
radius_neg     = min(ecc_neg);

% Diameter of graph
diameter_pos   = max(ecc_pos);
diameter_neg   = max(ecc_neg);
