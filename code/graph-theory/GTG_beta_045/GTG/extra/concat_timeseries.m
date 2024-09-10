function [out_tot] = concat_timeseries(out_1,out_2,manip_flag,varargin)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 02.29.16
%
% manip_flag = what type of data manipulation should occur before
%              concatenation. 
%              Specify 0 to do nothing
%              Specify 1 to mean center each timeseries individually before
%              concatenation (to avoid mean signal differences between runs
%              driving association measures)
%              Specify 2 to zscore each timeseries individually before
%              concatenation (to additionally avoid higher variance
%              timeseries having a greater impact on association measures)
%
% We recommend against specifying 0, as this may have a strong (and 
% differential by participant) impact on association measures.  
%
% Note that the output data will be sorted according to participant IDs.
% This is due to the fact that, if the number of participants differs for 
% the two input files, the script will match up participants.  
%
%
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if ~isempty(varargin)
    num_cens_vols_1 = varargin{1};
    num_cens_vols_2 = varargin{2};
    max_num_cens_vols = varargin{3};
    max_vols = varargin{4};
    
    num_cens_vols_1(isnan(num_cens_vols_1)) = max_vols;
    num_cens_vols_2(isnan(num_cens_vols_2)) = max_vols;
    num_cens_vols_tot = num_cens_vols_1 + num_cens_vols_2;
    
    out_1.ts(logical(num_cens_vols_tot > max_num_cens_vols)) = {NaN};
    out_2.ts(logical(num_cens_vols_tot > max_num_cens_vols)) = {NaN};
end

if ~issorted(out_1.subs)
    [out_1.subs,idx_1] = sort(out_1.subs);
    out_1.ts = out_1.ts(idx_1);
end

if ~issorted(out_2.subs)
    [out_2.subs,idx_2] = sort(out_2.subs);
    out_2.ts = out_2.ts(idx_2);
end

subs_2ts = intersect(out_1.subs,out_2.subs);
subs_1ts = setdiff(out_1.subs,out_2.subs);

nROIs = size(out_1.ts{1},1);

if nROIs ~= size(out_2.ts{1},1)
    error('Timeseries must have the same number of ROIs')
end

if isfield(out_1,'ROI_labels')
    out_tot.ROI_labels = out_1.ROI_labels;
end

for sub = 1:length(subs_2ts)
    curr_ts1            = out_1.ts{ismember(out_1.subs,subs_2ts{sub})};
    curr_ts2            = out_2.ts{ismember(out_2.subs,subs_2ts{sub})};
    out_tot.subs{sub,1} = subs_2ts{sub};
    
    switch manip_flag
        case 1
            try
                curr_ts1 = curr_ts1 - nanmean(curr_ts1(:));
                curr_ts2 = curr_ts2 - nanmean(curr_ts2(:));
            catch %#ok<*CTCH>
                curr_ts1 = curr_ts1 - mean(curr_ts1(:));
                curr_ts2 = curr_ts2 - mean(curr_ts2(:));
            end
                
        case 2
            try
                curr_ts1 = (curr_ts1 - nanmean(curr_ts1(:)))./nanstd(curr_ts1(:));
                curr_ts2 = (curr_ts2 - nanmean(curr_ts2(:)))./nanstd(curr_ts2(:));
            catch
                curr_ts1 = (curr_ts1 - mean(curr_ts1(:)))./std(curr_ts1(:));
                curr_ts2 = (curr_ts2 - mean(curr_ts2(:)))./std(curr_ts2(:));
            end
    end
    
    if size(curr_ts1,1) == 1 && size(curr_ts2,1) == 1
        out_tot.ts{sub,1} = NaN;
        out_tot.bounddiff_1val(sub,1) = NaN;
        out_tot.bounddiff_5val(sub,1) = NaN;
    elseif size(curr_ts1,1) == 1
        out_tot.ts{sub,1} = curr_ts2;
        out_tot.bounddiff_1val(sub,1) = NaN;
        out_tot.bounddiff_5val(sub,1) = NaN;
    elseif size(curr_ts2,1) == 1
        out_tot.ts{sub,1} = curr_ts1;
        out_tot.bounddiff_1val(sub,1) = NaN;
        out_tot.bounddiff_5val(sub,1) = NaN;
    else
        out_tot.ts{sub,1} = [curr_ts1,curr_ts2];
        
        try
            out_tot.bounddiff_1val(sub,1) = nanstd(curr_ts1(:,end) - curr_ts2(:,end))./nanmean([nanstd(curr_ts1,0,2);nanstd(curr_ts2,0,2)]);
            out_tot.bounddiff_5val(sub,1) = nanstd(nanmean(curr_ts1(:,(end-4):end),2)-nanmean(curr_ts2(:,(end-4):end),2))./nanmean([nanstd(curr_ts1,0,2);nanstd(curr_ts2,0,2)]);
        catch
            out_tot.bounddiff_1val(sub,1) = std(curr_ts1(:,end) - curr_ts2(:,end))./mean([std(curr_ts1,0,2);std(curr_ts2,0,2)]);
            out_tot.bounddiff_5val(sub,1) = std(mean(curr_ts1(:,(end-4):end),2) - mean(curr_ts2(:,(end-4):end),2))./mean([std(curr_ts1,0,2);std(curr_ts2,0,2)]);
        end
    end
end

for sub = 1:length(subs_1ts)
    if any(ismember(out_1.subs,subs_1ts{sub}))
        curr_ts = out_1.ts{ismember(out_1.subs,subs_2ts{sub})};
    elseif any(ismember(out_2.subs,subs_1ts{sub}))
        curr_ts = out_2.ts{ismember(out_2.subs,subs_2ts{sub})};
    end
    subind                 = sub+length(subs_2ts);
    out_tot.subs{subind,1} = subs_1ts{sub};
    
    switch manip_flag
        case 1
            try
                curr_ts = curr_ts - nanmean(curr_ts(:));
            catch %#ok<*CTCH>
                curr_ts = curr_ts - mean(curr_ts(:));
            end
                
        case 2
            try
                curr_ts = (curr_ts - nanmean(curr_ts(:)))./nanstd(curr_ts(:));
            catch
                curr_ts = (curr_ts - mean(curr_ts(:)))./std(curr_ts(:));
            end
    end
    
    if size(curr_ts,1) == 1
        out_tot.ts{subind,1} = NaN;
    else
        out_tot.ts{subind,1} = curr_ts;
    end
    out_tot.bounddiff_1val(subind,1) = NaN;
    out_tot.bounddiff_5val(subind,1) = NaN;
end

if ~issorted(out_tot.subs)
    [out_tot.subs,idx] = sort(out_tot.subs);
    out_tot.ts = out_tot.ts(idx);
    out_tot.bounddiff_1val = out_tot.bounddiff_1val(idx);
    out_tot.bounddiff_5val = out_tot.bounddiff_5val(idx);
end
