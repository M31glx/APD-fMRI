function [out,NaN_ROIs] = rem_nan_ROIs_ts(out)

% Remove ROIs with NaN entries from GTG output structure. 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 11.25.14
%
% WARNING: This is a beta version. There no known bugs, but only limited
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights


% Find NaNs
NaN_ROIs = [];
for sub = 1:length(out.subs)
    NaN_inds = find(sum(isnan(out.ts{sub}),2)>0);
    if ~isempty(NaN_inds)
        NaN_ROIs = [NaN_ROIs,NaN_inds(~ismember(NaN_inds,NaN_ROIs))'];
    end
end
NaN_ROIs = sort(NaN_ROIs);

% Remove NaN ROIs
for sub = 1:length(out.subs)
    out.ts{sub}(NaN_ROIs,:) = [];
end
out.ROI_labels(NaN_ROIs) = [];
out.nROI = length(out.ROI_labels);
