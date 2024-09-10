function [r_pos,r_neg] = assortativity_wei_sign(CIJ,flag)
% Edited version of Sporns et al.'s function (01/25/16 - Spielberg)
%       Now computed for positive and negative weights separately
%       (directionality calpability removed). Also, forces diagonal to be
%       0, consistent with the help file. 
%
%%ASSORTATIVITY_WEI_SIGN      Assortativity coefficient
%
%   [r_pos,r_neg] = assortativity_wei_sign(CIJ,flag);
%
%   The assortativity coefficient is a correlation coefficient between the
%   strengths (weighted degrees) of all nodes on two opposite ends of a link.
%   A positive assortativity coefficient indicates that nodes tend to link to
%   other nodes with the same or similar strength.
%
%   Inputs:           CIJ,    weighted directed/undirected connection matrix
%                    flag,    0, undirected graph: strength/strength correlation
%                             1, directed graph: out-strength/in-strength correlation
%                             2, directed graph: in-strength/out-strength correlation
%                             3, directed graph: out-strength/out-strength correlation
%                             4, directed graph: in-strength/in-strength correlation
%
%   Outputs: r_pos, r_neg,     assortativity coefficient
%
%   Notes: The main diagonal should be empty. For flag 1 the function computes 
%   the directed assortativity described in Rubinov and Sporns (2010) NeuroImage.
%
%   Reference:  Newman (2002) Phys Rev Lett 89:208701
%               Foster et al. (2010) PNAS 107:10815-10820
%
%   Olaf Sporns, Indiana University, 2007/2008
%   Vassilis Tsiaras, University of Crete, 2009
%   Murray Shanahan, Imperial College London, 2012
%   Mika Rubinov, University of Cambridge, 2012

if ~exist('flag','var')
    flag = 0;
end

n = length(CIJ);
CIJ_pos = threshold_absolute(CIJ,0);
CIJ_pos(1:n+1:end) = 0;
CIJ_neg = threshold_absolute(CIJ*-1,0);
CIJ_neg(1:n+1:end) = 0;

if (flag==0)                        % undirected version
    str_pos       = strengths_und(CIJ_pos);
    str_neg       = strengths_und(CIJ_neg);
    
    [i_pos,j_pos] = find(triu(CIJ_pos,1)>0);
    [i_neg,j_neg] = find(triu(CIJ_neg,1)>0);
    
    K_pos         = length(i_pos);
    K_neg         = length(i_neg);
    
    stri_pos      = str_pos(i_pos);
    stri_neg      = str_neg(i_neg);
    
    strj_pos      = str_pos(j_pos);
    strj_neg      = str_neg(j_neg);
else                                % directed versions
    [is_pos,os_pos] = strengths_dir(CIJ_pos);
    [is_neg,os_neg] = strengths_dir(CIJ_neg);
    
    [i_pos,j_pos]   = find(CIJ_pos>0);
    [i_neg,j_neg]   = find(CIJ_neg>0);
    
    K_pos           = length(i_pos);
    K_neg           = length(i_neg);

    switch flag
        case 1
            stri_pos = os_pos(i_pos);
            stri_neg = os_neg(i_neg);
            
            strj_pos = is_pos(j_pos);
            strj_neg = is_neg(j_neg);
        case 2
            stri_pos = is_pos(i_pos);
            stri_neg = is_neg(i_neg);
            
            strj_pos = os_pos(j_pos);
            strj_neg = os_neg(j_neg);
        case 3
            stri_pos = os_pos(i_pos);
            stri_neg = os_neg(i_neg);
            
            strj_pos = os_pos(j_pos);
            strj_neg = os_neg(j_neg);
        case 4
            stri_pos = is_pos(i_pos);
            stri_neg = is_neg(i_neg);
            
            strj_pos = is_pos(j_pos);
            strj_neg = is_neg(j_neg);
    end
end

% compute assortativity
r_pos = (sum(stri_pos.*strj_pos)/K_pos-(sum(0.5*(stri_pos+strj_pos))/ ...
    K_pos)^2)/(sum(0.5*(stri_pos.^2+strj_pos.^2))/K_pos-(sum(0.5* ...
    (stri_pos+strj_pos))/K_pos)^2);

r_neg = (sum(stri_neg.*strj_neg)/K_neg-(sum(0.5*(stri_neg+strj_neg))/ ...
    K_neg)^2)/(sum(0.5*(stri_neg.^2+strj_neg.^2))/K_neg-(sum(0.5* ...
    (stri_neg+strj_neg))/K_neg)^2);

