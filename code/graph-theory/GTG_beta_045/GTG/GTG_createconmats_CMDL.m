function GTG_createconmats_CMDL(config_filename)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: Beta 0.45 (03.30.16)
%
% Usage:    GTG_createconmats_CMDL('config_filename')
% 
% Input:    config_filename: Name (with path if not in current directory) 
%                            of the .mat file saved by the user. Note that
%                            this file must have all the desired options
%                            selected
% 
% History:
% 03.30.16 - Beta 0.45 - initial release of command line version% 
%
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but modifed, 
% redistributed versions must have the same rights

load(config_filename);

% Check whether inputs have been specified
if ~exist('out','var')
    msgbox('Enter output from previous stage (structure with timeseries)','Error','error')
    return
elseif ~isfield(out,'ts') && ~isfield(out,'deconv_ts')%#ok<*NODEF>
    error('No timeseries found in input structure')
elseif ~isfield(out,'subs')
    error('Participant IDs not found in input structure')
end
if ~isfield(out,'div_by_cond')
    out.div_by_cond = 'No';
elseif (strcmp(out.div_by_cond,'Yes - Block Design') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use canonical HRF') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use node-specific mean HRF')) && ~isfield(out,'timing_info')
    msgbox('Specify timing info','Error','error')
    return
end
if strcmp(out.div_by_cond,'Yes - Block Design') && ~isfield(out,'block_detr')
    out.block_detr = '0 (mean centering)';
end
if strcmp(out.div_by_cond,'Yes - Block Design') && ~isfield(out,'partial_des')
    out.partial_des = 'Yes';
end
if ~isfield(out,'TR')
    out.TR = 2;
end
if ~isfield(out,'connect_type')
    msgbox('Specify a type of connectivity to calculate','Error','error')
    return
end
if ~isfield(out,'min_TP')
    out.min_TP = 50;
else
    if out.min_TP<2
        out.min_TP = 2;
    end
end
if ~isfield(out,'nROI')
    out.nROI = length(out.ROI_labels);
end
if strcmp(out.div_by_cond,'Yes - Block Design') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use node-specific mean HRF') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use canonical HRF')
    out.data_type = 'Original timeseries';
elseif ~isfield(out,'data_type')
    if isfield(out,'ts')
        out.data_type = 'Original timeseries';
    elseif isfield(out,'deconv_ts')
        out.data_type = 'Deconvolved timeseries';
    end
elseif strcmpi(out.data_type,'Compute Conmats From:') && isfield(out,'ts')
    out.data_type = 'Original timeseries';
elseif strcmpi(out.data_type,'Compute Conmats From:') && ~isfield(out,'ts')
    out.data_type = 'Deconvolved timeseries';
elseif (strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both')) && ~isfield(out,'ts')
    msgbox('Connectivity was specified to be computed from the original timeseries, but they were not found','Error','error')
    return
elseif (strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both')) && ~isfield(out,'deconv_ts')
    msgbox('Connectivity was specified to be computed from the deconvolved timeseries, but they were not found','Error','error')
    return
end

if (strcmp(out.div_by_cond,'Yes - Block Design') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use canonical HRF') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use node-specific mean HRF'))
    if size(out.timing_info,2)>1
        if size(out.timing_info,2)>length(out.subs)
            msgbox('Too many sets of timing information were specified (i.e., > #participants)','Error','error')
            return
        elseif size(out.timing_info,2)<length(out.subs)
            msgbox('Too few sets of timing information were specified (i.e., < #participants)','Error','error')
            return
        end
        out.ind_timing_info = 1;
    else
        out.ind_timing_info = 0;
    end
end

if isempty(strfind(out.outname,'/'))
    out.outname = [pwd,'/',out.outname];
elseif out.outname(end)=='/'
    out.outname = [out.outname,'out'];
end
if ~strcmpi(out.outname(end-3:end),'.mat')
    out.outname = [out.outname,'.mat'];
end

if isfield(out,'num_conditions')
    for cond = 1:out.num_conditions
        for sub = 1:size(out.timing_info,2)
            out.timing_info{cond,sub} = out.timing_info{cond,sub}/out.TR;
        end
    end
end

use_parfor      = out.use_parfor;
num_par_workers = out.num_par_workers;

if use_parfor
    if isempty(gcp('nocreate'))
        if num_par_workers>feature('numCores')
            num_par_workers = feature('numCores');
        end
        if num_par_workers>1
            try
                parpool('local',num_par_workers);
            catch
                matlabpool('open',num_par_workers); %#ok<DPOOL>
            end
        else
            use_parfor = false;
        end
    end
end

nROI       = out.nROI;
min_TP     = out.min_TP;
ROI_labels = out.ROI_labels;

if isfield(out,'deconv_ts')
    totels  = 0;
    totnans = 0;
    for cl = 1:size(handles.out.deconv_ts,2)
        for cs = 1:size(handles.out.deconv_ts,1)
            totels  = totels+numel(handles.out.deconv_ts{cs,cl});
            totnans = totnans+sumn(isnan(handles.out.deconv_ts{cs,cl}));
        end
    end
    if totels==totnans
        out = rmfield(out,'deconv_ts');
    end
end

if strcmp(out.div_by_cond,'Yes - Block Design')
    if ~isfield(out,'deconv_ts')
        if strcmp(out.partial_des,'Yes')
            hrf_task    = spm_hrf(out.TR);
            task_desmat = [ones(size(out.ts{1},2),1),zeros(size(out.ts{1},2),out.num_conditions)];
            if out.ind_timing_info==0
                for cond = 1:out.num_conditions
                    for block = 1:size(out.timing_info{cond,1},2)
                        curr_onset = round(out.timing_info{cond,1}(1,block));
                        if curr_onset<1
                            curr_onset = 1;
                        end
                        curr_offset = round(out.timing_info{cond,1}(1,block)+out.timing_info{cond,1}(2,block))-1;
                        if curr_offset<curr_onset
                            curr_offset = curr_onset;
                        end
                        task_desmat(curr_onset:curr_offset,(cond+1)) = 1;
                    end
                    conved                  = conv(task_desmat(:,(cond+1)),hrf_task);
                    task_desmat(:,(cond+1)) = conved(1:size(task_desmat(:,(cond+1)),1));
                end
                for currsub = 1:length(out.subs)
                    for currROI = 1:nROI
                        part_vals                  = regstats(out.ts{currsub}(currROI,:),task_desmat,eye(size(task_desmat,2)),{'r'});
                        out.ts{currsub}(currROI,:) = part_vals.r;
                    end
                end
            elseif out.ind_timing_info==1
                for currsub = 1:length(out.subs)
                    for cond = 1:out.num_conditions
                        for block = 1:size(out.timing_info{cond,currsub},2)
                            curr_onset = round(out.timing_info{cond,currsub}(1,block));
                            if curr_onset<1
                                curr_onset = 1;
                            end
                            curr_offset = round(out.timing_info{cond,currsub}(1,block)+out.timing_info{cond,currsub}(2,block))-1;
                            if curr_offset<curr_onset
                                curr_offset = curr_onset;
                            end
                            task_desmat(curr_onset:curr_offset,(cond+1)) = 1;
                        end
                        conved                  = conv(task_desmat(:,(cond+1)),hrf_task);
                        task_desmat(:,(cond+1)) = conved(1:size(task_desmat(:,(cond+1)),1));
                    end
                    for currROI = 1:nROI
                        part_vals                  = regstats(out.ts{currsub}(currROI,:),task_desmat,eye(size(task_desmat,2)),{'r'});
                        out.ts{currsub}(currROI,:) = part_vals.r;
                    end
                end
            end
        end
        dt          = out.TR/16;
        NT          = out.TR/dt;
        hrf         = spm_hrf(dt);
        out.ts_orig = out.ts;
        out.ts      = cell(length(out.subs),out.num_conditions);
        nTP         = size(out.ts_orig{1},2);
        progressbar('Progress Deconvolving')
        for currsub = 1:length(out.subs)
            currsub_ts = out.ts_orig{currsub};
            if size(currsub_ts,2)>1
                decon_ts = zeros(nROI,round(nTP.*NT));
                if use_parfor
                    parfor currROI = 1:nROI
                        P      = {};
                        ROI_ts = currsub_ts(currROI,:)';
                        xb     = spm_dctmtx(nTP*NT+128,nTP);
                        P{1}.X = zeros(nTP,nTP);
                        for i = 1:nTP
                            Hx          = conv(xb(:,i),hrf);
                            P{1}.X(:,i) = Hx((1:NT:nTP*NT)+128);
                        end
                        xb                  = xb(129:end,:);
                        P{1}.C              = speye(nTP,nTP)/4;
                        P{2}.X              = sparse(nTP,1);
                        P{2}.C              = speye(nTP,nTP)*nTP/trace(P{1}.X'*P{1}.X);
                        C                   = spm_PEB(ROI_ts,P);
                        decon_ts(currROI,:) = xb*C{2}.E(1:nTP);
                    end
                else
                    for currROI = 1:nROI
                        P      = {};
                        ROI_ts = currsub_ts(currROI,:)';
                        xb     = spm_dctmtx(nTP*NT+128,nTP);
                        P{1}.X = zeros(nTP,nTP);
                        for i = 1:nTP
                            Hx          = conv(xb(:,i),hrf);
                            P{1}.X(:,i) = Hx((1:NT:nTP*NT)+128);
                        end
                        xb                  = xb(129:end,:);
                        P{1}.C              = speye(nTP,nTP)/4;
                        P{2}.X              = sparse(nTP,1);
                        P{2}.C              = speye(nTP,nTP)*nTP/trace(P{1}.X'*P{1}.X);
                        C                   = spm_PEB(ROI_ts,P);
                        decon_ts(currROI,:) = xb*C{2}.E(1:nTP);
                    end
                end
            end
            for cond = 1:out.num_conditions
                if out.ind_timing_info==0
                    block_on  = round(out.timing_info{cond,1}(1,:)*NT)+1;
                    block_off = round(((out.timing_info{cond,1}(1,:)+out.timing_info{cond,1}(2,:))-1)*NT)+1;
                elseif out.ind_timing_info==1
                    block_on  = round(out.timing_info{cond,currsub}(1,:)*NT)+1;
                    block_off = round(((out.timing_info{cond,currsub}(1,:)+out.timing_info{cond,currsub}(2,:))-1)*NT)+1;
                end
                if block_on<1
                    block_on = 1;
                end
                if block_off<block_on
                    block_off = block_on;
                end
                for block = 1:length(block_on)
                    if ~strcmp(out.block_detr,'Polynomial Detrending Order') && ~strcmp(out.block_detr,'No detrending') && ~strcmp(out.block_detr,'0 (mean centering)')
                        decon_ts(:,block_on(block):block_off(block)) = spm_detrend(decon_ts(:,block_on(block):block_off(block)),str2double(out.block_detr));
                    elseif strcmp(out.block_detr,'0 (mean centering)')
                        decon_ts(:,block_on(block):block_off(block)) = decon_ts(:,block_on(block):block_off(block))-mean(decon_ts(:,block_on(block):block_off(block)));
                    end
                    out.ts{currsub,cond} = [out.ts{currsub,cond},decon_ts(:,block_on(block):block_off(block))];
                end
            end
            prog = currsub/length(out.subs);
            progressbar(prog)
        end
    else
        out.deconv_ts_orig = out.deconv_ts;
        out.ts_orig        = out.ts;
        out.ts             = cell(length(out.subs),out.num_conditions);
        for currsub = 1:length(out.subs)
            decon_ts = out.deconv_ts_orig{currsub};
            for cond = 1:out.num_conditions
                if out.ind_timing_info==0
                    block_on  = round(out.timing_info{cond,1}(1,:))+1;
                    block_off = round(((out.timing_info{cond,1}(1,:)+out.timing_info{cond,1}(2,:))-1))+1;
                elseif out.ind_timing_info==1
                    block_on  = round(out.timing_info{cond,currsub}(1,:))+1;
                    block_off = round(((out.timing_info{cond,currsub}(1,:)+out.timing_info{cond,currsub}(2,:))-1))+1;
                end
                if block_on<1
                    block_on = 1;
                end
                if block_off<block_on
                    block_off = block_on;
                end
                for block = 1:length(block_on)
                    if ~strcmp(out.block_detr,'Polynomial Detrending Order') && ~strcmp(out.block_detr,'No detrending') && ~strcmp(out.block_detr,'0 (mean centering)')
                        decon_ts(:,block_on(block):block_off(block)) = spm_detrend(decon_ts(:,block_on(block):block_off(block)),str2double(out.block_detr));
                    elseif strcmp(out.block_detr,'0 (mean centering)')
                        decon_ts(:,block_on(block):block_off(block)) = decon_ts(:,block_on(block):block_off(block))-mean(decon_ts(:,block_on(block):block_off(block)));
                    end
                    out.ts{currsub,cond} = [out.ts{currsub,cond},decon_ts(:,block_on(block):block_off(block))];
                end
            end
            prog = currsub/length(out.subs);
            progressbar(prog)
        end
    end
elseif strcmp(out.div_by_cond,'Yes - Event-Related Design - Use canonical HRF')
    progressbar('Progress Calculating Beta Series')
    hrf = spm_hrf(out.TR);
    nTP = size(out.ts{1},2);
    if out.ind_timing_info==0
        beta_desmats = cell(out.num_conditions,1);
        for cond = 1:out.num_conditions
            for trial = 1:size(out.timing_info{cond,1}(1,:),2)
                pred_ct = zeros(nTP,1);
                pred_ot = zeros(nTP,1);
                for ccond = 1:out.num_conditions
                    for ctrial = 1:size(out.timing_info{cond,1}(1,:),2)
                        trial_start = round(out.timing_info{ccond,1}(1,ctrial));
                        if trial_start<1
                            trial_start = 1;
                        end
                        trial_end = round(out.timing_info{ccond,1}(1,ctrial)+out.timing_info{ccond,1}(2,ctrial))-1;
                        if trial_end<trial_start
                            trial_end = trial_start;
                        end
                        if ccond==cond && ctrial==trial
                            pred_ct(trial_start:trial_end) = 1;
                        else
                            pred_ot(trial_start:trial_end) = 1;
                        end
                    end
                end
                beta_desmats{cond,trial}(1:nTP,1) = ones(nTP,1);
                convpred                          = conv(pred_ct,hrf);
                beta_desmats{cond,trial}(1:nTP,2) = convpred(1:nTP);
                convpred                          = conv(pred_ot,hrf);
                beta_desmats{cond,trial}(1:nTP,3) = convpred(1:nTP);
            end
        end
        out.beta_ts = cell(length(out.subs),out.num_conditions);
        for currsub = 1:length(out.subs)
            for cond = 1:out.num_conditions
                for currROI = 1:nROI
                    for trial = 1:size(out.timing_info{cond,1}(1,:),2)
                        B                                        = (beta_desmats{cond,trial}'*beta_desmats{cond,trial})\beta_desmats{cond,trial}'*out.ts{currsub}(currROI,:)';
                        out.beta_ts{currsub,cond}(currROI,trial) = B(2);
                    end
                end
            end
            prog = currsub/length(out.subs);
            progressbar(prog)
        end
    elseif out.ind_timing_info==1
        out.beta_ts = cell(length(out.subs),out.num_conditions);
        for currsub = 1:length(out.subs)
            beta_desmats = cell(out.num_conditions,1);
            for cond = 1:out.num_conditions
                for trial = 1:size(out.timing_info{cond,currsub}(1,:),2)
                    pred_ct = zeros(nTP,1);
                    pred_ot = zeros(nTP,1);
                    for ccond = 1:out.num_conditions
                        for ctrial = 1:size(out.timing_info{cond,currsub}(1,:),2)
                            trial_start = round(out.timing_info{ccond,currsub}(1,ctrial));
                            if trial_start<1
                                trial_start = 1;
                            end
                            trial_end = round(out.timing_info{ccond,currsub}(1,ctrial)+out.timing_info{ccond,currsub}(2,ctrial))-1;
                            if trial_end<trial_start
                                trial_end = trial_start;
                            end
                            if ccond==cond && ctrial==trial
                                pred_ct(trial_start:trial_end) = 1;
                            else
                                pred_ot(trial_start:trial_end) = 1;
                            end
                        end
                    end
                    beta_desmats{cond,trial}(1:nTP,1) = ones(nTP,1);
                    convpred                          = conv(pred_ct,hrf);
                    beta_desmats{cond,trial}(1:nTP,2) = convpred(1:nTP);
                    convpred                          = conv(pred_ot,hrf);
                    beta_desmats{cond,trial}(1:nTP,3) = convpred(1:nTP);
                end
            end
            for cond = 1:out.num_conditions
                for currROI = 1:nROI
                    for trial = 1:size(out.timing_info{cond,currsub}(1,:),2)
                        B                                        = (beta_desmats{cond,trial}'*beta_desmats{cond,trial})\beta_desmats{cond,trial}'*out.ts{currsub}(currROI,:)';
                        out.beta_ts{currsub,cond}(currROI,trial) = B(2);
                    end
                end
            end
            prog = currsub/length(out.subs);
            progressbar(prog)
        end
    end
elseif strcmp(out.div_by_cond,'Yes - Event-Related Design - Use node-specific mean HRF')
    progressbar('Progress Calculating Beta Series')
    nTP = size(out.ts{1},2);
    
    if out.ind_timing_info==0
        hrf_length   = length(spm_hrf(out.TR));
        hrf_length   = hrf_length+round(out.timing_info{1,1}(2,1))-1;
        FIR_preds    = zeros(nTP,hrf_length);
        trial_starts = [];
        for cond = 1:out.num_conditions
            trial_starts = [trial_starts round(out.timing_info{cond,1}(1,:))]; %#ok<AGROW>
        end
        for cimp = 1:hrf_length
            FIR_preds(trial_starts,cimp) = 1;
            trial_starts                 = trial_starts+1;
        end
        out.spec_hrfs = zeros(length(out.subs),nROI,hrf_length);
        for currsub = 1:length(out.subs)
            for currROI = 1:nROI
                out.spec_hrfs(currsub,currROI,:) = regress(out.ts{currsub}(currROI,:)',FIR_preds);
            end
        end
        out.beta_ts = cell(length(out.subs),out.num_conditions);
        for currsub = 1:length(out.subs)
            for currROI = 1:nROI
                beta_desmats = cell(out.num_conditions,1);
                for cond = 1:out.num_conditions
                    for trial = 1:size(out.timing_info{cond,1}(1,:),2)
                        pred_ct = zeros(nTP,1);
                        pred_ot = zeros(nTP,1);
                        for ccond = 1:out.num_conditions
                            for ctrial = 1:size(out.timing_info{cond,1}(1,:),2)
                                trial_start = round(out.timing_info{ccond,1}(1,ctrial));
                                if trial_start<1
                                    trial_start = 1;
                                end
                                trial_end = round(out.timing_info{ccond,1}(1,ctrial)+out.timing_info{ccond,1}(2,ctrial))-1;
                                if trial_end<trial_start
                                    trial_end = trial_start;
                                end
                                if ccond==cond && ctrial==trial
                                    pred_ct(trial_start:trial_end) = 1;
                                else
                                    pred_ot(trial_start:trial_end) = 1;
                                end
                            end
                        end
                        beta_desmats{cond,trial}(1:nTP,1) = ones(nTP,1);
                        convpred                          = conv(pred_ct,squeeze(out.spec_hrfs(currsub,currROI,:)));
                        beta_desmats{cond,trial}(1:nTP,2) = convpred(1:nTP);
                        convpred                          = conv(pred_ot,squeeze(out.spec_hrfs(currsub,currROI,:)));
                        beta_desmats{cond,trial}(1:nTP,3) = convpred(1:nTP);
                    end
                end
                for cond = 1:out.num_conditions
                    for trial = 1:size(out.timing_info{cond,1}(1,:),2)
                        B                                        = (beta_desmats{cond,trial}'*beta_desmats{cond,trial})\beta_desmats{cond,trial}'*out.ts{currsub}(currROI,:)';
                        out.beta_ts{currsub,cond}(currROI,trial) = B(2);
                    end
                end
            end
            prog = currsub/length(out.subs);
            progressbar(prog)
        end
    elseif out.ind_timing_info==1
        out.spec_hrfs = zeros(length(out.subs),nROI,hrf_length);
        out.beta_ts   = cell(length(out.subs),out.num_conditions);
        for currsub = 1:length(out.subs)
            hrf_length   = length(spm_hrf(out.TR));
            hrf_length   = hrf_length+round(out.timing_info{1,currsub}(2,1))-1;
            FIR_preds    = zeros(nTP,hrf_length);
            trial_starts = [];
            for cond = 1:out.num_conditions
                trial_starts = [trial_starts round(out.timing_info{cond,currsub}(1,:))]; %#ok<AGROW>
            end
            for cimp = 1:hrf_length
                FIR_preds(trial_starts,cimp) = 1;
                trial_starts                 = trial_starts+1;
            end
            for currROI = 1:nROI
                out.spec_hrfs(currsub,currROI,:) = regress(out.ts{currsub}(currROI,:)',FIR_preds);
            end
            
            for currROI = 1:nROI
                beta_desmats = cell(out.num_conditions,1);
                for cond = 1:out.num_conditions
                    for trial = 1:size(out.timing_info{cond,currsub}(1,:),2)
                        pred_ct = zeros(nTP,1);
                        pred_ot = zeros(nTP,1);
                        for ccond = 1:out.num_conditions
                            for ctrial = 1:size(out.timing_info{cond,currsub}(1,:),2)
                                trial_start = round(out.timing_info{ccond,currsub}(1,ctrial));
                                if trial_start<1
                                    trial_start = 1;
                                end
                                trial_end = round(out.timing_info{ccond,currsub}(1,ctrial)+out.timing_info{ccond,currsub}(2,ctrial))-1;
                                if trial_end<trial_start
                                    trial_end = trial_start;
                                end
                                if ccond==cond && ctrial==trial
                                    pred_ct(trial_start:trial_end) = 1;
                                else
                                    pred_ot(trial_start:trial_end) = 1;
                                end
                            end
                        end
                        beta_desmats{cond,trial}(1:nTP,1) = ones(nTP,1);
                        convpred                          = conv(pred_ct,squeeze(out.spec_hrfs(currsub,currROI,:)));
                        beta_desmats{cond,trial}(1:nTP,2) = convpred(1:nTP);
                        convpred                          = conv(pred_ot,squeeze(out.spec_hrfs(currsub,currROI,:)));
                        beta_desmats{cond,trial}(1:nTP,3) = convpred(1:nTP);
                    end
                end
                for cond = 1:out.num_conditions
                    for trial = 1:size(out.timing_info{cond,currsub}(1,:),2)
                        B                                        = (beta_desmats{cond,trial}'*beta_desmats{cond,trial})\beta_desmats{cond,trial}'*out.ts{currsub}(currROI,:)';
                        out.beta_ts{currsub,cond}(currROI,trial) = B(2);
                    end
                end
            end
            prog = currsub/length(out.subs);
            progressbar(prog)
        end
    end
end

switch out.connect_type
    case 'Pearson Correlation'
        mat_outname        = strrep(out.outname,'.mat','_fullcorr.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_fullcorr.mat');
        logfile_outname    = strrep(out.outname,'.mat','_fullcorr_logfile.txt');
    case 'Partial Correlation'
        mat_outname        = strrep(out.outname,'.mat','_partialcorr.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_partialcorr.mat');
        logfile_outname    = strrep(out.outname,'.mat','_partialcorr_logfile.txt');
    case 'Regularized Partial Correlation'
        mat_outname        = strrep(out.outname,'.mat','_regpartialcorr.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_regpartialcorr.mat');
        logfile_outname    = strrep(out.outname,'.mat','_regpartialcorr_logfile.txt');
    case 'Mutual Information'
        mat_outname        = strrep(out.outname,'.mat','_mutualinfo.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_mutualinfo.mat');
        logfile_outname    = strrep(out.outname,'.mat','_mutualinfo_logfile.txt');
    case 'Robust Correlation'
        mat_outname        = strrep(out.outname,'.mat','_robustcorr.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_robustcorr.mat');
        logfile_outname    = strrep(out.outname,'.mat','_robustcorr_logfile.txt');
    case 'Kendall''s Tau'
        mat_outname        = strrep(out.outname,'.mat','_kendalltau.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_kendalltau.mat');
        logfile_outname    = strrep(out.outname,'.mat','_kendalltau_logfile.txt');
    case 'Spearman''s Rho'
        mat_outname        = strrep(out.outname,'.mat','_spearmanrho.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_spearmanrho.mat');
        logfile_outname    = strrep(out.outname,'.mat','_spearmanrho_logfile.txt');
    case 'Gaussian Copula'
        mat_outname        = strrep(out.outname,'.mat','_gausscopula.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_gausscopula.mat');
        logfile_outname    = strrep(out.outname,'.mat','_gausscopula_logfile.txt');
    case 't Copula'
        mat_outname        = strrep(out.outname,'.mat','_tcopula.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_tcopula.mat');
        logfile_outname    = strrep(out.outname,'.mat','_tcopula_logfile.txt');
    case 'Frank Copula'
        mat_outname        = strrep(out.outname,'.mat','_frankcopula.mat');
        deconv_mat_outname = strrep(out.outname,'.mat','_deconv_frankcopula.mat');
        logfile_outname    = strrep(out.outname,'.mat','_frankcopula_logfile.txt');
end

if size(out.ts,1)~=length(out.subs) && size(out.ts,2)==length(out.subs)
    out.ts = out.ts';
elseif size(out.ts,1)~=length(out.subs) && size(out.ts,2)~=length(out.subs)
    msgbox('The number of participant timeseries differs from the number of participant IDs','Error','error')
    return
end

if isfield(out,'num_conditions')
    out.num_rep_levs = out.num_conditions;
elseif ~isfield(out,'beta_ts')
    out.num_rep_levs = size(out.ts,2);
else
    out.num_rep_levs = size(out.beta_ts,2);
end

if isfield(out,'ts') && (strcmpi(out.data_type,'Compute Conmats From:') || strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both'))
    outmats = zeros(nROI,nROI,length(out.subs),out.num_rep_levs);
end
if isfield(out,'deconv_ts') && (strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both'))
    outmats_deconv = zeros(nROI,nROI,length(out.subs),out.num_rep_levs);
end
progressbar('Progress Calculating Connectivity Matrices')

if ~strcmp(out.div_by_cond,'Yes - Event-Related Design - Use node-specific mean HRF') && ~strcmp(out.div_by_cond,'Yes - Event-Related Design - Use canonical HRF')
    switch out.connect_type
        case {'Pearson Correlation','Mutual Information','Robust Correlation','Kendall''s Tau','Spearman''s Rho','Gaussian Copula','t Copula','Frank Copula'}
            calc_cells = logical(triu(ones(nROI),1));
            for currsub = 1:length(out.subs)
                if iscell(out.subs)
                    sub = out.subs{currsub};
                else
                    sub = num2str(out.subs(currsub));
                end
                for rep_lev = 1:out.num_rep_levs
                    if isfield(out,'ts') && (strcmpi(out.data_type,'Compute Conmats From:') || strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both'))
                        currsub_ts = squeeze(out.ts{currsub,rep_lev});
                        if size(currsub_ts,2)>1
                            for rowvar = 1:nROI
                                rowts = currsub_ts(rowvar,:)';
                                if use_parfor
                                    parfor colvar = 1:nROI
                                        if calc_cells(rowvar,colvar)==1
                                            colts   = currsub_ts(colvar,:)';
                                            usevols = logical((~isnan(rowts)).*(~isnan(colts)));
                                            if sum(usevols)>=min_TP
                                                switch out.connect_type
                                                    case 'Pearson Correlation'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts(usevols),colts(usevols));
                                                    case 'Mutual Information'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = kernelmi(rowts(usevols)',colts(usevols)');
                                                    case 'Robust Correlation'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = bendcorr(rowts(usevols),colts(usevols),0);
                                                    case 'Kendall''s Tau'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts(usevols),colts(usevols),'type','Kendall');
                                                    case 'Spearman''s Rho'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts(usevols),colts(usevols),'type','Spearman');
                                                    case 'Gaussian Copula'
                                                        rowts_kdens                            = ksdensity(rowts(usevols),rowts(usevols),'function','cdf');
                                                        colts_kdens                            = ksdensity(colts(usevols),colts(usevols),'function','cdf');
                                                        tempcop                                = copulafit('Gaussian',[rowts_kdens,colts_kdens]);
                                                        outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 't Copula'
                                                        rowts_kdens                            = ksdensity(rowts(usevols),rowts(usevols),'function','cdf');
                                                        colts_kdens                            = ksdensity(colts(usevols),colts(usevols),'function','cdf');
                                                        tempcop                                = copulafit('t',[rowts_kdens,colts_kdens]);
                                                        outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 'Frank Copula'
                                                        rowts_kdens                            = ksdensity(rowts(usevols),rowts(usevols),'function','cdf');
                                                        colts_kdens                            = ksdensity(colts(usevols),colts(usevols),'function','cdf');
                                                        outmats(rowvar,colvar,currsub,rep_lev) = copulafit('Frank',[rowts_kdens,colts_kdens]);
                                                end
                                            else
                                                logfile_fid                            = fopen(logfile_outname,'a');
                                                outmats(rowvar,colvar,currsub,rep_lev) = NaN;
                                                fprintf(logfile_fid,'Less than %i usable timepoints for computing the correlation between %s and %s for %s\n',min_TP,ROI_labels{rowvar},ROI_labels{colvar},sub); %#ok<*PFBNS>
                                                fclose(logfile_fid);
                                            end
                                        end
                                    end
                                else
                                    logfile_fid = fopen(logfile_outname,'a');
                                    for colvar = 1:nROI
                                        if calc_cells(rowvar,colvar)==1
                                            colts   = currsub_ts(colvar,:)';
                                            usevols = logical((~isnan(rowts)).*(~isnan(colts)));
                                            if sum(usevols)>=min_TP
                                                switch out.connect_type
                                                    case 'Pearson Correlation'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts(usevols),colts(usevols));
                                                    case 'Mutual Information'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = kernelmi(rowts(usevols)',colts(usevols)');
                                                    case 'Robust Correlation'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = bendcorr(rowts(usevols),colts(usevols),0);
                                                    case 'Kendall''s Tau'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts(usevols),colts(usevols),'type','Kendall');
                                                    case 'Spearman''s Rho'
                                                        outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts(usevols),colts(usevols),'type','Spearman');
                                                    case 'Gaussian Copula'
                                                        rowts_kdens                            = ksdensity(rowts(usevols),rowts(usevols),'function','cdf');
                                                        colts_kdens                            = ksdensity(colts(usevols),colts(usevols),'function','cdf');
                                                        tempcop                                = copulafit('Gaussian',[rowts_kdens,colts_kdens]);
                                                        outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 't Copula'
                                                        rowts_kdens                            = ksdensity(rowts(usevols),rowts(usevols),'function','cdf');
                                                        colts_kdens                            = ksdensity(colts(usevols),colts(usevols),'function','cdf');
                                                        tempcop                                = copulafit('t',[rowts_kdens,colts_kdens]);
                                                        outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 'Frank Copula'
                                                        rowts_kdens                            = ksdensity(rowts(usevols),rowts(usevols),'function','cdf');
                                                        colts_kdens                            = ksdensity(colts(usevols),colts(usevols),'function','cdf');
                                                        outmats(rowvar,colvar,currsub,rep_lev) = copulafit('Frank',[rowts_kdens,colts_kdens]);
                                                end
                                            else
                                                outmats(rowvar,colvar,currsub,rep_lev) = NaN;
                                                fprintf(logfile_fid,'Less than %i usable timepoints for computing the correlation between %s and %s for %s\n',min_TP,out.ROI_labels{rowvar},out.ROI_labels{colvar},sub);
                                            end
                                        end
                                    end
                                    fclose(logfile_fid);
                                end
                            end
                            outmats(:,:,currsub,rep_lev) = outmats(:,:,currsub,rep_lev)+outmats(:,:,currsub,rep_lev)';
                        else
                            outmats(:,:,currsub,rep_lev) = NaN;
                        end
                    end
                    if isfield(out,'deconv_ts') && (strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both'))
                        currsub_deconv_ts = squeeze(out.deconv_ts{currsub,rep_lev});
                        if size(currsub_deconv_ts,2)>1
                            for rowvar = 1:nROI
                                row_deconv_ts = currsub_deconv_ts(rowvar,:)';
                                if use_parfor
                                    connect_type = out.connect_type;
                                    parfor colvar = 1:nROI
                                        if calc_cells(rowvar,colvar)==1
                                            col_deconv_ts  = currsub_deconv_ts(colvar,:)';
                                            deconv_usevols = logical((~isnan(row_deconv_ts)).*(~isnan(col_deconv_ts)));
                                            if sum(deconv_usevols)>=min_TP
                                                switch connect_type
                                                    case 'Pearson Correlation'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = corr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols));
                                                    case 'Mutual Information'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = kernelmi(row_deconv_ts(deconv_usevols)',col_deconv_ts(deconv_usevols)');
                                                    case 'Robust Correlation'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = bendcorr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),0);
                                                    case 'Kendall''s Tau'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = corr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'type','Kendall');
                                                    case 'Spearman''s Rho'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = corr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'type','Spearman');
                                                    case 'Gaussian Copula'
                                                        rowts_kdens                                   = ksdensity(row_deconv_ts(deconv_usevols),row_deconv_ts(deconv_usevols),'function','cdf');
                                                        colts_kdens                                   = ksdensity(col_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'function','cdf');
                                                        tempcop                                       = copulafit('Gaussian',[rowts_kdens,colts_kdens]);
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 't Copula'
                                                        rowts_kdens                                   = ksdensity(row_deconv_ts(deconv_usevols),row_deconv_ts(deconv_usevols),'function','cdf');
                                                        colts_kdens                                   = ksdensity(col_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'function','cdf');
                                                        tempcop                                       = copulafit('t',[rowts_kdens,colts_kdens]);
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 'Frank Copula'
                                                        rowts_kdens                                   = ksdensity(row_deconv_ts(deconv_usevols),row_deconv_ts(deconv_usevols),'function','cdf');
                                                        colts_kdens                                   = ksdensity(col_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'function','cdf');
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = copulafit('Frank',[rowts_kdens,colts_kdens]);
                                                end
                                            else
                                                logfile_fid                                   = fopen(logfile_outname,'a');
                                                outmats_deconv(rowvar,colvar,currsub,rep_lev) = NaN;
                                                fprintf(logfile_fid,'Less than %i usable deconvovlved timepoints for computing the correlation between %s and %s for %s\n',min_TP,ROI_labels{rowvar},ROI_labels{colvar},sub); %#ok<*PFBNS>
                                                fclose(logfile_fid);
                                            end
                                        end
                                    end
                                else
                                    logfile_fid = fopen(logfile_outname,'a');
                                    for colvar = 1:nROI
                                        if calc_cells(rowvar,colvar)==1
                                            col_deconv_ts  = currsub_deconv_ts(colvar,:)';
                                            deconv_usevols = logical((~isnan(row_deconv_ts)).*(~isnan(col_deconv_ts)));
                                            if sum(deconv_usevols)>=min_TP
                                                switch out.connect_type
                                                    case 'Pearson Correlation'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = corr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols));
                                                    case 'Mutual Information'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = kernelmi(row_deconv_ts(deconv_usevols)',col_deconv_ts(deconv_usevols)');
                                                    case 'Robust Correlation'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = bendcorr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),0);
                                                    case 'Kendall''s Tau'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = corr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'type','Kendall');
                                                    case 'Spearman''s Rho'
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = corr(row_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'type','Spearman');
                                                    case 'Gaussian Copula'
                                                        rowts_kdens                                   = ksdensity(row_deconv_ts(deconv_usevols),row_deconv_ts(deconv_usevols),'function','cdf');
                                                        colts_kdens                                   = ksdensity(col_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'function','cdf');
                                                        tempcop                                       = copulafit('Gaussian',[rowts_kdens,colts_kdens]);
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 't Copula'
                                                        rowts_kdens                                   = ksdensity(row_deconv_ts(deconv_usevols),row_deconv_ts(deconv_usevols),'function','cdf');
                                                        colts_kdens                                   = ksdensity(col_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'function','cdf');
                                                        tempcop                                       = copulafit('t',[rowts_kdens,colts_kdens]);
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                                    case 'Frank Copula'
                                                        rowts_kdens                                   = ksdensity(row_deconv_ts(deconv_usevols),row_deconv_ts(deconv_usevols),'function','cdf');
                                                        colts_kdens                                   = ksdensity(col_deconv_ts(deconv_usevols),col_deconv_ts(deconv_usevols),'function','cdf');
                                                        outmats_deconv(rowvar,colvar,currsub,rep_lev) = copulafit('Frank',[rowts_kdens,colts_kdens]);
                                                end
                                            else
                                                outmats_deconv(rowvar,colvar,currsub,rep_lev) = NaN;
                                                fprintf(logfile_fid,'Less than %i usable deconvolved timepoints for computing the correlation between %s and %s for %s\n',min_TP,out.ROI_labels{rowvar},out.ROI_labels{colvar},sub);
                                            end
                                        end
                                    end
                                    fclose(logfile_fid);
                                end
                            end
                            outmats_deconv(:,:,currsub,rep_lev) = outmats_deconv(:,:,currsub,rep_lev)+outmats_deconv(:,:,currsub,rep_lev)';
                        else
                            outmats_deconv(:,:,currsub,rep_lev) = NaN;
                        end
                    end
                    prog = currsub/length(out.subs);
                    progressbar(prog)
                end
            end
        case 'Partial Correlation'
            if out.min_TP<=nROI
                out.min_TP = nROI+1;
                min_TP     = out.min_TP;
                disp('Minimum number of usable timepoints needed for Partial Correlation is at least #ROIs + 1 and has been reset to that value');
            end
            if use_parfor
                subs         = out.subs;
                num_rep_levs = out.num_rep_levs;
                if strcmpi(out.data_type,'Compute Conmats From:') || strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both')
                    ts = out.ts;
                else
                    ts = [];
                end
                if strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both')
                    decon_ts = out.deconv_ts;
                else
                    decon_ts = [];
                end
                data_type = out.data_type;
                parfor currsub = 1:length(subs)
                    logfile_fid = fopen(logfile_outname,'a');
                    if iscell(subs)
                        sub = subs{currsub};
                    else
                        sub = num2str(subs(currsub));
                    end
                    for rep_lev = 1:num_rep_levs
                        if ~isempty(ts) && (strcmpi(data_type,'Compute Conmats From:') || strcmpi(data_type,'Original timeseries') || strcmpi(data_type,'Both'))
                            if size(ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_ts                   = squeeze(ts{currsub,rep_lev});
                                    usevols                      = sum(isnan(currsub_ts),1)==0;
                                    outmats(:,:,currsub,rep_lev) = partialcorr(currsub_ts(:,usevols)');
                                else
                                    outmats(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable timepoints for computing partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                        if ~isempty(decon_ts) && (strcmpi(data_type,'Deconvolved timeseries') || strcmpi(data_type,'Both'))
                            if size(ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(decon_ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_deconv_ts                   = squeeze(decon_ts{currsub,rep_lev});
                                    deconv_usevols                      = sum(isnan(currsub_deconv_ts),1)==0;
                                    outmats_deconv(:,:,currsub,rep_lev) = partialcorr(currsub_deconv_ts(:,deconv_usevols)');
                                else
                                    outmats_deconv(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable deconvolved timepoints for computing partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats_deconv(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                    end
                    fclose(logfile_fid);
                end
            else
                logfile_fid = fopen(logfile_outname,'a');
                for currsub = 1:length(out.subs)
                    if iscell(out.subs)
                        sub = out.subs{currsub};
                    else
                        sub = num2str(out.subs(currsub));
                    end
                    for rep_lev = 1:out.num_rep_levs
                        if isfield(out,'ts') && (strcmpi(out.data_type,'Compute Conmats From:') || strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both'))
                            if size(out.ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(out.ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_ts                   = squeeze(out.ts{currsub,rep_lev});
                                    usevols                      = sum(isnan(currsub_ts),1)==0;
                                    outmats(:,:,currsub,rep_lev) = partialcorr(currsub_ts(:,usevols)');
                                else
                                    outmats(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable timepoints for computing partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                        if isfield(out,'deconv_ts') && (strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both'))
                            if size(out.deconv_ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(out.deconv_ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_deconv_ts                   = squeeze(out.deconv_ts{currsub,rep_lev});
                                    deconv_usevols                      = sum(isnan(currsub_deconv_ts),1)==0;
                                    outmats_deconv(:,:,currsub,rep_lev) = partialcorr(currsub_deconv_ts(:,deconv_usevols)');
                                else
                                    outmats_deconv(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable deconvolved timepoints for computing partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats_deconv(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                    end
                    prog = currsub/length(out.subs);
                    progressbar(prog)
                end
                fclose(logfile_fid);
            end
        case 'Regularized Partial Correlation'
            if use_parfor
                subs         = out.subs;
                num_rep_levs = out.num_rep_levs;
                if strcmpi(out.data_type,'Compute Conmats From:') || strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both')
                    ts = out.ts;
                else
                    ts = [];
                end
                if strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both')
                    decon_ts = out.deconv_ts;
                else
                    decon_ts = [];
                end
                data_type = out.data_type;
                parfor currsub = 1:length(subs)
                    logfile_fid = fopen(logfile_outname,'a');
                    if iscell(subs)
                        sub = subs{currsub};
                    else
                        sub = num2str(subs(currsub));
                    end
                    for rep_lev = 1:num_rep_levs
                        if ~isempty(ts) && (strcmpi(data_type,'Compute Conmats From:') || strcmpi(data_type,'Original timeseries') || strcmpi(data_type,'Both'))
                            if size(ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_ts                   = squeeze(ts{currsub,rep_lev});
                                    usevols                      = sum(isnan(currsub_ts),1)==0;
                                    [~,shrinkage]                = cov1para(currsub_ts(:,usevols)');
                                    shrunkcov                    = ((1-shrinkage)*cov(currsub_ts(:,usevols)'))+(shrinkage*eye(nROI));
                                    outmats(:,:,currsub,rep_lev) = 2*eye(nROI)+(-inv(shrunkcov)./sqrt(diag(inv(shrunkcov))*diag(inv(shrunkcov))'));
                                else
                                    outmats(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable timepoints for computing Ledoit-Wolf regularized partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                        if ~isempty(decon_ts) && (strcmpi(data_type,'Deconvolved timeseries') || strcmpi(data_type,'Both'))
                            if size(ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(decon_ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_deconv_ts                   = squeeze(decon_ts{currsub,rep_lev});
                                    deconv_usevols                      = sum(isnan(currsub_deconv_ts),1)==0;
                                    [~,shrinkage]                       = cov1para(currsub_deconv_ts(:,deconv_usevols)');
                                    shrunkcov                           = ((1-shrinkage)*cov(currsub_deconv_ts'))+(shrinkage*eye(nROI));
                                    outmats_deconv(:,:,currsub,rep_lev) = 2*eye(nROI)+(-inv(shrunkcov)./sqrt(diag(inv(shrunkcov))*diag(inv(shrunkcov))'));
                                else
                                    outmats_deconv(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable deconvolved timepoints for computing Ledoit-Wolf regularized partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats_deconv(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                    end
                    fclose(logfile_fid);
                end
            else
                logfile_fid = fopen(logfile_outname,'a');
                for currsub = 1:length(out.subs)
                    if iscell(out.subs)
                        sub = out.subs{currsub};
                    else
                        sub = num2str(out.subs(currsub));
                    end
                    for rep_lev = 1:out.num_rep_levs
                        if isfield(out,'ts') && (strcmpi(out.data_type,'Compute Conmats From:') || strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both'))
                            if size(out.ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(out.ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_ts                   = squeeze(out.ts{currsub,rep_lev});
                                    usevols                      = sum(isnan(currsub_ts),1)==0;
                                    [~,shrinkage]                = cov1para(currsub_ts(:,usevols)');
                                    shrunkcov                    = ((1-shrinkage)*cov(currsub_ts(:,usevols)'))+(shrinkage*eye(nROI));
                                    outmats(:,:,currsub,rep_lev) = 2*eye(nROI)+(-inv(shrunkcov)./sqrt(diag(inv(shrunkcov))*diag(inv(shrunkcov))'));
                                else
                                    outmats(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable timepoints for computing Ledoit-Wolf regularized partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                        if isfield(out,'deconv_ts') && (strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both'))
                            if size(out.deconv_ts{currsub},2)>1
                                if sum(sum(isnan(squeeze(out.deconv_ts{currsub,rep_lev})),1)==0)>=min_TP
                                    currsub_deconv_ts                   = squeeze(out.deconv_ts{currsub,rep_lev});
                                    deconv_usevols                      = sum(isnan(currsub_deconv_ts),1)==0;
                                    [~,shrinkage]                       = cov1para(currsub_deconv_ts(:,deconv_usevols)');
                                    shrunkcov                           = ((1-shrinkage)*cov(currsub_deconv_ts'))+(shrinkage*eye(nROI));
                                    outmats_deconv(:,:,currsub,rep_lev) = 2*eye(nROI)+(-inv(shrunkcov)./sqrt(diag(inv(shrunkcov))*diag(inv(shrunkcov))'));
                                else
                                    outmats_deconv(:,:,currsub,rep_lev) = NaN;
                                    fprintf(logfile_fid,'Less than i% usable deconvolved timepoints for computing Ledoit-Wolf regularized partial correlations for %s\n',min_TP,sub);
                                end
                            else
                                outmats_deconv(:,:,currsub,rep_lev) = NaN;
                            end
                        end
                    end
                    prog = currsub/length(out.subs);
                    progressbar(prog)
                end
                fclose(logfile_fid);
            end
    end
else
    switch out.connect_type
        case {'Pearson Correlation','Mutual Information','Robust Correlation','Kendall''s Tau','Spearman''s Rho','Gaussian Copula','t Copula','Frank Copula'}
            calc_cells  = logical(triu(ones(nROI),1));
            for currsub = 1:length(out.subs)
                for rep_lev = 1:out.num_rep_levs
                    currsub_ts = squeeze(out.beta_ts{currsub,rep_lev});
                    if size(currsub_ts,2)>1
                        for rowvar = 1:nROI
                            rowts = currsub_ts(rowvar,:)';
                            if use_parfor
                                connect_type = out.connect_type;
                                parfor colvar = 1:nROI
                                    if calc_cells(rowvar,colvar)==1
                                        switch connect_type
                                            case 'Pearson Correlation'
                                                outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts,currsub_ts(colvar,:)');
                                            case 'Mutual Information'
                                                outmats(rowvar,colvar,currsub,rep_lev) = kernelmi(rowts',currsub_ts(colvar,:));
                                            case 'Robust Correlation'
                                                outmats(rowvar,colvar,currsub,rep_lev) = bendcorr(rowts,currsub_ts(colvar,:)',0);
                                            case 'Kendall''s Tau'
                                                outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts,currsub_ts(colvar,:)','type','Kendall');
                                            case 'Spearman''s Rho'
                                                outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts,currsub_ts(colvar,:)','type','Spearman');
                                            case 'Gaussian Copula'
                                                rowts_kdens                            = ksdensity(rowts,rowts,'function','cdf');
                                                colts_kdens                            = ksdensity(currsub_ts(colvar,:)',currsub_ts(colvar,:)','function','cdf');
                                                tempcop                                = copulafit('Gaussian',[rowts_kdens,colts_kdens]);
                                                outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                            case 't Copula'
                                                rowts_kdens                            = ksdensity(rowts,rowts,'function','cdf');
                                                colts_kdens                            = ksdensity(currsub_ts(colvar,:)',currsub_ts(colvar,:)','function','cdf');
                                                tempcop                                = copulafit('t',[rowts_kdens,colts_kdens]);
                                                outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                            case 'Frank Copula'
                                                rowts_kdens                            = ksdensity(rowts,rowts,'function','cdf');
                                                colts_kdens                            = ksdensity(currsub_ts(colvar,:)',currsub_ts(colvar,:)','function','cdf');
                                                outmats(rowvar,colvar,currsub,rep_lev) = copulafit('Frank',[rowts_kdens,colts_kdens]);
                                        end
                                    end
                                end
                            else
                                for colvar = 1:nROI
                                    if calc_cells(rowvar,colvar)==1
                                        switch out.connect_type
                                            case 'Pearson Correlation'
                                                outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts,currsub_ts(colvar,:)');
                                            case 'Mutual Information'
                                                outmats(rowvar,colvar,currsub,rep_lev) = kernelmi(rowts',currsub_ts(colvar,:));
                                            case 'Robust Correlation'
                                                outmats(rowvar,colvar,currsub,rep_lev) = bendcorr(rowts,currsub_ts(colvar,:)',0);
                                            case 'Kendall''s Tau'
                                                outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts,currsub_ts(colvar,:)','type','Kendall');
                                            case 'Spearman''s Rho'
                                                outmats(rowvar,colvar,currsub,rep_lev) = corr(rowts,currsub_ts(colvar,:)','type','Spearman');
                                            case 'Gaussian Copula'
                                                rowts_kdens                            = ksdensity(rowts,rowts,'function','cdf');
                                                colts_kdens                            = ksdensity(currsub_ts(colvar,:)',currsub_ts(colvar,:)','function','cdf');
                                                tempcop                                = copulafit('Gaussian',[rowts_kdens,colts_kdens]);
                                                outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                            case 't Copula'
                                                rowts_kdens                            = ksdensity(rowts,rowts,'function','cdf');
                                                colts_kdens                            = ksdensity(currsub_ts(colvar,:)',currsub_ts(colvar,:)','function','cdf');
                                                tempcop                                = copulafit('t',[rowts_kdens,colts_kdens]);
                                                outmats(rowvar,colvar,currsub,rep_lev) = tempcop(1,2);
                                            case 'Frank Copula'
                                                rowts_kdens                            = ksdensity(rowts,rowts,'function','cdf');
                                                colts_kdens                            = ksdensity(currsub_ts(colvar,:)',currsub_ts(colvar,:)','function','cdf');
                                                outmats(rowvar,colvar,currsub,rep_lev) = copulafit('Frank',[rowts_kdens,colts_kdens]);
                                        end
                                    end
                                end
                            end
                        end
                        outmats(:,:,currsub,rep_lev) = outmats(:,:,currsub,rep_lev)+outmats(:,:,currsub,rep_lev)';
                    else
                        outmats(:,:,currsub,rep_lev) = NaN;
                    end
                end
                prog = currsub/length(out.subs);
                progressbar(prog)
            end
        case 'Partial Correlation'
            if use_parfor
                subs         = out.subs;
                num_rep_levs = out.num_rep_levs;
                ts           = out.ts;
                beta_ts      = out.beta_ts;
                parfor currsub = 1:length(subs)
                    for rep_lev = 1:num_rep_levs
                        if size(ts{currsub},2)>1
                            outmats(:,:,currsub,rep_lev) = partialcorr(squeeze(beta_ts{currsub,rep_lev})');
                        else
                            outmats(:,:,currsub,rep_lev) = NaN;
                        end
                    end
                end
            else
                for currsub = 1:length(out.subs)
                    for rep_lev = 1:out.num_rep_levs
                        if size(out.ts{currsub},2)>1
                            outmats(:,:,currsub,rep_lev) = partialcorr(squeeze(out.beta_ts{currsub,rep_lev})');
                        else
                            outmats(:,:,currsub,rep_lev) = NaN;
                        end
                    end
                    prog = currsub/length(out.subs);
                    progressbar(prog)
                end
            end
        case 'Regularized Partial Correlation'
            if use_parfor
                subs         = out.subs;
                num_rep_levs = out.num_rep_levs;
                ts           = out.ts;
                beta_ts      = out.beta_ts;
                parfor currsub = 1:length(subs)
                    for rep_lev = 1:num_rep_levs
                        if size(ts{currsub},2)>1
                            [~,shrinkage]                = cov1para(squeeze(beta_ts{currsub,rep_lev})');
                            shrunkcov                    = ((1-shrinkage)*cov(squeeze(beta_ts{currsub,rep_lev})'))+(shrinkage*eye(nROI));
                            outmats(:,:,currsub,rep_lev) = 2*eye(nROI)+(-inv(shrunkcov)./sqrt(diag(inv(shrunkcov))*diag(inv(shrunkcov))'));
                        else
                            outmats(:,:,currsub,rep_lev) = NaN;
                        end
                    end
                end
            else
                for currsub = 1:length(out.subs)
                    for rep_lev = 1:out.num_rep_levs
                        if size(out.ts{currsub},2)>1
                            [~,shrinkage]                = cov1para(squeeze(out.beta_ts{currsub,rep_lev})');
                            shrunkcov                    = ((1-shrinkage)*cov(squeeze(out.beta_ts{currsub,rep_lev})'))+(shrinkage*eye(nROI));
                            outmats(:,:,currsub,rep_lev) = 2*eye(nROI)+(-inv(shrunkcov)./sqrt(diag(inv(shrunkcov))*diag(inv(shrunkcov))'));
                        else
                            outmats(:,:,currsub,rep_lev) = NaN;
                        end
                    end
                    prog = currsub/length(out.subs);
                    progressbar(prog)
                end
            end
    end
end

if use_parfor
    try
        parpool close
    catch
        matlabpool close %#ok<DPOOL>
    end
end

if isfield(out,'deconv_ts_orig')
    out.deconv_ts_block = out.ts;
    out.deconv_ts       = out.deconv_ts_orig;
    out.ts              = out.ts_orig;
    out                 = rmfield(out,{'deconv_ts_orig','ts_orig'});
elseif isfield(out,'ts_orig')
    out.deconv_ts_block = out.ts;
    out.ts              = out.ts_orig;
    out                 = rmfield(out,'ts_orig');
end

if strcmpi(out.data_type,'Compute Conmats From:') || strcmpi(out.data_type,'Original timeseries') || strcmpi(out.data_type,'Both') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use node-specific mean HRF') || strcmp(out.div_by_cond,'Yes - Event-Related Design - Use canonical HRF')
    out.nonan   = sum(squeeze(sum(isnan(sum(outmats,4)),1)),1)'==0;
    out.allnan  = sum(squeeze(sum(isnan(sum(outmats,4)),1)),1)'==nROI^2;
    out.conmats = outmats;
    out.conmats(logical(repmat(eye(size(out.conmats,1)),[1,1,size(out.conmats,3),size(out.conmats,4),size(out.conmats,5)]))) = 0;                         % Set diagonal of connectivity matrices to 0
    save(mat_outname,'out');
end

if (strcmpi(out.data_type,'Deconvolved timeseries') || strcmpi(out.data_type,'Both')) && (~strcmp(out.div_by_cond,'Yes - Event-Related Design - Use node-specific mean HRF') && ~strcmp(out.div_by_cond,'Yes - Event-Related Design - Use canonical HRF'))
    out.nonan   = sum(squeeze(sum(isnan(sum(outmats_deconv,4)),1)),1)'==0;
    out.allnan  = sum(squeeze(sum(isnan(sum(outmats_deconv,4)),1)),1)'==nROI^2;
    out.conmats = outmats_deconv;
    out.conmats(logical(repmat(eye(size(out.conmats,1)),[1,1,size(out.conmats,3),size(out.conmats,4),size(out.conmats,5)]))) = 0;                         % Set diagonal of connectivity matrices to 0
    save(deconv_mat_outname,'out');
end