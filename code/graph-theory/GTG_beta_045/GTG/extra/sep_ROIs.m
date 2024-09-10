function sep = sep_ROIs(comb)

% Create a separate binary matrix for each ROI in the input matrix. 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 03.07.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

labels = unique(comb(:));
labels(labels==0) = [];
sep = zeros([size(comb),length(labels)]);

for ROI = 1:length(labels)
    temp = zeros(size(comb));
    temp(comb==labels(ROI)) = 1;
    sep(:,:,:,ROI) = temp;
end
