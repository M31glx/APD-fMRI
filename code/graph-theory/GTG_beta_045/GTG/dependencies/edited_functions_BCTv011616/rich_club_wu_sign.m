function [Rw_pos,Rw_neg] = rich_club_wu_sign(CIJ,varargin)
%%RICH_CLUB_WU_SIGN Rich club coefficients curve (weighted undirected graph)
%
%   [Rw_pos,Rw_neg] = rich_club_wu_sign(CIJ,varargin) % rich club curve for weighted graph
%
%   The weighted rich club coefficient, Rw, at level k is the fraction of
%   edge weights that connect nodes of degree k or higher out of the
%   maximum edge weights that such nodes might share.
%
%   Inputs:
%       CIJ:        weighted directed connection matrix
%
%       k-level:    (optional) max level of RC(k).
%                   (by default k-level quals the maximal degree of CIJ)
%                
%   Output:
%       Rw:         rich-club curve
%
%
%   References:     
%       T Opsahl et al. Phys Rev Lett, 2008, 101(16)
%       M van den Heuvel, O Sporns, J Neurosci 2011 31(44)
%
%   Martijn van den Heuvel, University Medical Center Utrecht, 2011

%   Modification History:
%   2011: Original
%   2015: Expanded documentation (Mika Rubinov)
%   2016: Computed for pos and neg weights separately (JMS)


NofNodes = size(CIJ,2);     %number of nodes
CIJ(1:NofNodes+1:end) = 0;
CIJ_pos = CIJ.*(CIJ>0);
CIJ_neg = -CIJ.*(CIJ<0);

NodeDegree_pos = sum((CIJ_pos~=0)); %define degree of each node
NodeDegree_neg = sum((CIJ_neg~=0)); %define degree of each node

%define to which level rc should be computed
if size(varargin,2)==1
    klevel_pos = varargin{1};
    klevel_neg = varargin{1};
elseif isempty(varargin)
    klevel_pos = max(NodeDegree_pos);
    klevel_neg = max(NodeDegree_neg);
else
    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
end

% For positive weights
%wrank contains the ranked weights of the network, with strongest connections on top
wrank_pos = sort(CIJ_pos(:), 'descend');

%loop over all possible k-levels
for kk = 1:klevel_pos
    
    SmallNodes_pos = find(NodeDegree_pos<kk);
    
    if isempty(SmallNodes_pos);
        Rw_pos(kk)=NaN;             %#ok<*AGROW>
        continue
    end
    
    %remove small nodes with NodeDegree<kk
    CutoutCIJ_pos                   = CIJ_pos;
    CutoutCIJ_pos(SmallNodes_pos,:) = [];
    CutoutCIJ_pos(:,SmallNodes_pos) = [];
    
    %total weight of connections in subset E>r
    Wr_pos = sum(CutoutCIJ_pos(:));
    
    %total number of connections in subset E>r
    Er_pos = length(find(CutoutCIJ_pos~=0));
    
    %E>r number of connections with max weight in network
    wrank_r_pos = wrank_pos(1:1:Er_pos);
    
    %weighted rich-club coefficient
    Rw_pos(kk) = Wr_pos/sum(wrank_r_pos); 
end

% For negative weights
%wrank contains the ranked weights of the network, with strongest connections on top
wrank_neg = sort(CIJ_neg(:), 'descend');

%loop over all possible k-levels
for kk = 1:klevel_neg
    
    SmallNodes_neg = find(NodeDegree_neg<kk);
    
    if isempty(SmallNodes_neg);
        Rw_neg(kk)=NaN;             %#ok<*AGROW>
        continue
    end
    
    %remove small nodes with NodeDegree<kk
    CutoutCIJ_neg                   = CIJ_neg;
    CutoutCIJ_neg(SmallNodes_neg,:) = [];
    CutoutCIJ_neg(:,SmallNodes_neg) = [];
    
    %total weight of connections in subset E>r
    Wr_neg = sum(CutoutCIJ_neg(:));
    
    %total number of connections in subset E>r
    Er_neg = length(find(CutoutCIJ_neg~=0));
    
    %E>r number of connections with max weight in network
    wrank_r_neg = wrank_neg(1:1:Er_neg);
    
    %weighted rich-club coefficient
    Rw_neg(kk) = Wr_neg/sum(wrank_r_neg); 
end
