%input the residuals and calculate the between ROIs correlation by
%partialling out the effect of other ROIs
function [RR,PP]=PartialOut(resid)
m = size(resid,2); 
RR = eye(m); %create identity matrix of X
% set up another matrix to store p values
PP = RR;
for i = 1:m 
    for j = 1:i-1 
       k = setdiff(1:m,[i j]); %select columns to be controlled for
       [RR(i,j) PP(i,j)] = partialcorr(resid(:,i),resid(:,j),resid(:,k));
       RR(j,i) = RR(i,j);
       PP(j,i) = PP(i,j);
    end 
end
