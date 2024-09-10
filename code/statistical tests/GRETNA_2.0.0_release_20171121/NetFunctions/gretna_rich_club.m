function phi = gretna_rich_club(A)
%==========================================================================
% This function is used to calculate rich club metric for a binary graph or
% network G.
%
%
% Syntax:  function [rc] = gretna_rich_club(A)
%
% Input:
%        A:
%                The adjencent matrix of G.
%
% Outputs:
%        phi:
%                Rich club coefficient of G.
%                1 by N-1 Array: if k>k_max or k<k_min. rc(i)=NaN
%
% References:
% 1. Martijn P. van den Heuvel and Olaf Sporns (2011) Rich-Club Organization
%    of the Human Connectome. The Journal of Neuroscience, 31:15775C15786.
%
% Jinhui WANG, CCBD, HZNU, HangZhou, 2013/04/25, Jinhui.Wang.1982@gmail.com
%
% Revised by Xindi Wang 20151231
%==========================================================================

A=A - diag(diag(A));
A=abs(A);
N=size(A, 1);
phi=nan(1, N-1);

K=sum(A);
kmin = max([1 min(K)]); kmax = max(K);

k = kmin:kmax-1;

for i = 1:length(k)
    ind = find(K <= k(i));
    
    net = A;
    net(ind,:) = [];
    net(:,ind) = [];
    
    if sum(net(:)) == 0, break, end
    
    %rc.deg(i,1) = k(i);
    phi(1, k(i)) = sum(net(:))/length(net)/(length(net)-1);
end