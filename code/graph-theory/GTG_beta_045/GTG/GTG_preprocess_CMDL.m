function GTG_preprocess_CMDL(config_filename)

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
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

load(config_filename);

% Check whether inputs have been specified
if ~isfield(out,'subs') %#ok<NODEF>
    msgbox('Enter a cell array of participant IDs','Error','error')
    return
end
if any(~cell2mat(cellfun(@ischar,out.subs,'UniformOutput',false)))
    out.subs = cellfun(@num2str,out.subs,'UniformOutput',false);
end
if ~isfield(out,'useavail') || strcmp(out.useavail,'Use files if available?')
    out.useavail = 'No, run all processing steps';
end
if ~isfield(out,'ROI_labels')
    msgbox('Enter a cell array of ROI labels','Error','error')
    return
end
if ~isfield(out,'ROIextord') || strcmp(out.ROIextord,'Extract from ROIs:')
    out.ROIextord = 'After all preprocessing';
end
if ~isfield(out,'saveintermed') || strcmp(out.saveintermed,'Save Intermediate Files?')
    out.saveintermed = 'No';
end
if ~isfield(out,'nTP')
    out.nTP = 120;
end
if ~isfield(out,'TR')
    out.TR = 2;
end
if ~isfield(out,'MinROIvox')
    out.MinROIvox = 5;
end
if ~isfield(out,'MinpercentROIvox')
    out.MinpercentROIvox = 50;
end
if ~isfield(out,'func_first_filename')
    msgbox('Enter the 4d EPI filename for the first participant','Error','error')
    return
end
if ~isfield(out,'first_funcmask_filename')
    msgbox('Enter the 3d brainmask filename for the first participant','Error','error')
    return
end
if ~isfield(out,'first_parc_filename')
    msgbox('Enter the filename for the 3d image containing ROI masks for the first participant','Error','error')
    return
end
if ~isfield(out,'parc_space') || strcmp(out.parc_space,'Image Space')
    out.parc_space = 'Functional';
end
if ~isfield(out,'ST_correct')
    out.ST_correct = 0;
end
if out.ST_correct==1 && (~isfield(out,'ST_ord') || strcmp(out.ST_ord,'Slice Collection Order'))
    msgbox('Please select the slice aquisition order','Error','error')
    return
elseif out.ST_correct==0 && ~isfield(out,'ST_ord')
    out.ST_ord = [];
end
if ~isfield(out,'detr')
    out.detr = 0;
end
if out.detr==1 && (~isfield(out,'detr_ord') || strcmp(out.detr_ord,'Polynomial Order'))
    out.detr_ord = 0;
elseif ~isfield(out,'detr_ord')
    out.detr_ord = [];
end
if ~isfield(out,'BP')
    out.BP = 0;
end
if out.BP==1 && (~isfield(out,'BP_type') || strcmp(out.BP_type,'Filter Type'))
    msgbox('Please select the bandpass filter type','Error','error')
    return
elseif out.BP==0 && ~isfield(out,'BP_type')
    out.BP_type = [];
end
if ~isfield(out,'despike')
    out.despike = 0;
end
if ~isfield(out,'BP_LP_cutoff')
    out.BP_LP_cutoff = 0.10;
elseif isempty(out.BP_LP_cutoff)
    out.BP_LP_cutoff = -1;
elseif out.BP_LP_cutoff<=0
    out.BP_LP_cutoff = -1;
end
if ~isfield(out,'BP_HP_cutoff')
    out.BP_HP_cutoff = 0.01;
elseif isempty(out.BP_HP_cutoff)
    out.BP_HP_cutoff = -1;
elseif out.BP_HP_cutoff<=0
    out.BP_HP_cutoff = -1;
end
if out.BP_HP_cutoff==-1 && out.BP_LP_cutoff==-1
    msgbox('Temporal filtering is selected but both lowpass and highpass values are invalid','Error','error')
    return
end
if isfield(out,'BP_type')
    if strcmp(out.BP_type,'Butterworth') && (~isfield(out,'butter_ord') || strcmp(out.butter_ord,'Butterworth Order'))
        out.butter_ord = 1;
    elseif ~strcmp(out.BP_type,'Butterworth') && (~isfield(out,'butter_ord') || strcmp(out.butter_ord,'Butterworth Order'))
        out.butter_ord = [];
    end
elseif isempty(out.BP_type) || ~isfield(out,'butter_ord')
    out.butter_ord = [];
end
if ~isfield(out,'MC')
    out.MC = 0;
end
if ~isfield(out,'motcens')
    out.motcens       = 0;
    out.motcens_logic = [];
elseif out.motcens==0
    out.motcens_logic = [];
elseif out.motcens==1
    out.motcens_logic = get(get(handles.ANDORmotcens_uipanel,'SelectedObject'),'String');
end
if ~isfield(out,'motcens_FD_cutoff')
    out.motcens_FD_cutoff = 0.3;
end
if ~isfield(out,'motcens_DVARS_cutoff')
    out.motcens_DVARS_cutoff = 2.5;
end
if ~isfield(out,'motpar_part')
    out.motpar_part = 0;
end
if ~isfield(out,'motpart1_part')
    out.motpart1_part = 0;
end
if ~isfield(out,'motparsqr_part')
    out.motparsqr_part = 0;
end
if ~isfield(out,'motparderiv_part')
    out.motparderiv_part = 0;
end
if ~isfield(out,'globsig_part')
    out.globsig_part = 0;
end
if ~isfield(out,'globsigderiv_part')
    out.globsigderiv_part = 0;
end
if ~isfield(out,'globsig_calcGCOR')
    out.globsig_calcGCOR = 0;
end
if ~isfield(out,'globsig_calcGNI')
    out.globsig_calcGNI = 0;
end
if ~isfield(out,'WMsig_part')
    out.WMsig_part = 0;
end
if ~isfield(out,'WMsigderiv_part')
    out.WMsigderiv_part = 0;
end
if out.WMsig_part==1 && ~isfield(out,'first_WM_filename')
    msgbox('Enter the filename for the white matter mask for the first participant','Error','error')
    return
end
if ~isfield(out,'WM_space') || strcmp(out.WM_space,'Image Space')
    out.WM_space = 'Functional';
end
if ~isfield(out,'WMmask_scope') || strcmp(out.WMmask_scope,'Scope of Mask')
    out.WMmask_scope = 'Entire Mask';
end
if ~isfield(out,'ventsig_part')
    out.ventsig_part = 0;
end
if ~isfield(out,'ventsigderiv_part')
    out.ventsigderiv_part = 0;
end
if out.ventsig_part==1 && ~isfield(out,'first_vent_filename')
    msgbox('Enter the filename for the ventricular mask for the first participant','Error','error')
    return
end
if ~isfield(out,'vent_space') || strcmp(out.vent_space,'Image Space')
    out.vent_space = 'Functional';
end
if ~isfield(out,'PCA_WMvent')
    out.PCA_WMvent = 0;
end
if ~isfield(out,'erode_masks')
    out.erode_masks = 0;
end
if out.erode_masks==1 && ~isfield(out,'erode_radius')
    out.erode_radius = 2;
end
if isfield(out,'erode_radius')
    out.erode_radius = round(out.erode_radius);
    if out.erode_radius<2
        out.erode_radius = 2;
    end
else
    out.erode_radius = [];
end
if ~isfield(out,'rem_AR')
    out.rem_AR = 0;
end
if out.rem_AR==1 && ~isfield(out,'num_hist_tp')
    out.num_hist_tp = 1;
end
if isfield(out,'num_hist_tp')
    out.num_hist_tp = round(out.num_hist_tp);
    if out.num_hist_tp<1
        out.num_hist_tp = 1;
    elseif out.num_hist_tp>round(out.nTP/2)
        out.num_hist_tp = round(out.nTP/2);
    end
else
    out.num_hist_tp = [];
end
if ~isfield(out,'deconv')
    out.deconv        = 0;
    out.deconv_method = [];
elseif out.deconv==0 && ~isfield(out,'deconv_method')
    out.deconv_method = [];
elseif out.deconv==1 && (~isfield(out,'deconv_method') || isempty(out.deconv_method))
    out.deconv_method = 'SPM';
end
if out.deconv==1 && ~isfield(out,'whendeconv')
    out.whendeconv = 'after';
elseif out.deconv==0 && ~isfield(out,'whendeconv')
    out.whendeconv = [];
end
if ~isfield(out,'ROI_summary_measure') || strcmp(out.ROI_summary_measure,'Timeseries Extraction Measure')
    out.ROI_summary_measure = 'Mean';
end
if strcmp(out.ROIextord,'After only slice timing/motion correction/detrending')
    out.despike          = 0;
    out.motcens          = 0;
    out.globsig_calcGCOR = 0;
    out.WMmask_scope     = 'Entire Mask';
end
if ~isfield(out,'outname')
    out.outname = [pwd,'/','out.mat'];
elseif isempty(strfind(out.outname,'/'))
    out.outname = [pwd,'/',out.outname];
elseif out.outname(end)=='/'
    out.outname = [out.outname,'out'];
end
if ~strcmpi(out.outname(end-3:end),'.mat')
    out.outname = [out.outname,'.mat'];
end

% Set base filenames
func_base_filename     = strrep(out.func_first_filename,out.subs{1},'SUBNUM');
funcmask_base_filename = strrep(out.first_funcmask_filename,out.subs{1},'SUBNUM');
parc_base_filename     = strrep(out.first_parc_filename,out.subs{1},'SUBNUM');
if isfield(out,'first_WM_filename')
    WM_base_filename = strrep(out.first_WM_filename,out.subs{1},'SUBNUM');
else
    WM_base_filename = ' ';
end
if isfield(out,'first_vent_filename')
    vent_base_filename = strrep(out.first_vent_filename,out.subs{1},'SUBNUM');
else
    vent_base_filename = ' ';
end

par_base_filename = out.par_base_filename;

if strcmp(out.parc_space,'Structural') || strcmp(out.WM_space,'Structural') || strcmp(out.vent_space,'Structural')
    [file,path]                      = uigetfile('*.mat','Select the structural-to-functional transformation matrix for the first participant');
    filename                         = [path,file];
    struct2func_regmat_base_filename = strrep(filename,out.subs{1},'SUBNUM');
else
    struct2func_regmat_base_filename = ' ';
end
if strcmp(out.parc_space,'Standard') || strcmp(out.WM_space,'Standard') || strcmp(out.vent_space,'Standard')
    [path,file]                        = uigetfile({'*.nii.gz';'*.nii'},'Select the standard-to-functional warp file for the first participant');
    filename                           = [path,file];
    standard2func_regmat_base_filename = strrep(filename,out.subs{1},'SUBNUM');
else
    standard2func_regmat_base_filename = ' ';
end

use_parfor = out.use_parfor;
num_par_workers = out.num_par_workers;

if use_parfor
    if isempty(gcp('nocreate'))
        if num_par_workers>feature('numCores')
            num_par_workers = feature('numCores');
        end
        if num_par_workers>1
            try
                parpool('open',num_par_workers);
            catch
                matlabpool('open',num_par_workers); %#ok<DPOOL>
            end
        else
            use_parfor = false;
        end
    end
end

out.preproc_outname = strrep(out.outname,'.mat','_preproc.mat');

% Write non-participant specific info to logfile
preproc_logfile_outname = strrep(out.outname,'.mat','_preproc_logfile.txt');
logfile_fid = fopen(preproc_logfile_outname,'w');
fprintf(logfile_fid,'Output filename = %s\n',out.preproc_outname);
fprintf(logfile_fid,'\nBase filenames:\n');
fprintf(logfile_fid,'Base filename for functional data = %s\n',func_base_filename);
fprintf(logfile_fid,'Base filename for functional mask = %s\n',funcmask_base_filename);
fprintf(logfile_fid,'Base filename for ROI mask = %s, which was input in %s space\n',parc_base_filename,out.parc_space);
if out.WMsig_part==1
    fprintf(logfile_fid,'Base filename for white matter mask = %s, which was input in %s space\n',WM_base_filename,out.WM_space);
end
if out.ventsig_part==1
    fprintf(logfile_fid,'Base filename for ventricle mask = %s, which was input in %s space\n',vent_base_filename,out.vent_space);
end
if strcmp(out.parc_space,'Structural') || strcmp(out.WM_space,'Structural') || strcmp(out.vent_space,'Structural')
    fprintf(logfile_fid,'Base filename for structural to functional warp = %s\n',struct2func_regmat_base_filename);
end
if strcmp(out.parc_space,'Standard') || strcmp(out.WM_space,'Standard') || strcmp(out.vent_space,'Standard')
    fprintf(logfile_fid,'Base filename for standard space to functional warp = %s\n',standard2func_regmat_base_filename);
end
fprintf(logfile_fid,'# of ROIs = %u\n',out.nROI);
fprintf(logfile_fid,'# of timepoints = %u\n',out.nTP);
fprintf(logfile_fid,'(User specified) TR = %u\n',out.TR);
fprintf(logfile_fid,'Minimum # of voxels required in masked ROI = %u\n',out.MinROIvox);
fprintf(logfile_fid,'Minimum percentage of ROI required to be retained after masking = %u\n',out.MinpercentROIvox);
if strcmp(out.ROIextord,'After all preprocessing')
    fprintf(logfile_fid,'The ROI timeseries were extracted after all preprocessing had occurred\n');
elseif strcmp(out.ROIextord,'After only slice timing/motion correction/detrending')
    fprintf(logfile_fid,'The ROI timeseries were extracted after only slice timing\motion correction/detrending (if these steps were performed), after which all preprocessing was performed on the extracted timeseries\n');
end
if strcmp(out.useavail,'Yes, do not run step if output file exists')
    fprintf(logfile_fid,'Existing files were used if available for a participant\n');
end
if out.ST_correct==1
    fprintf(logfile_fid,'Slicetiming correction was performed with %s (user specified) slice order\n',out.ST_ord);
end
if out.MC==1
    fprintf(logfile_fid,'Motion correction was performed\n');
end
if out.detr==1
    fprintf(logfile_fid,'The data were detrended for up to polynomials of order %u\n',out.detr_ord);
end
if out.despike==1
    fprintf(logfile_fid,'The data were wavelet despiked\n');
end
if out.BP==1 && strcmp(out.BP_type,'Ideal')
    if out.BP_LP_cutoff==-1
        fprintf(logfile_fid,'The data were highpass filtered with the "Ideal" filter (cutoff = %6.4f)\n',out.BP_HP_cutoff);
    elseif out.BP_HP_cutoff==-1
        fprintf(logfile_fid,'The data were lowpass filtered with the "Ideal" filter (cutoff = %6.4f)\n',out.BP_LP_cutoff);
    else
        fprintf(logfile_fid,'The data were bandpass filtered with the "Ideal" filter (lowpass cutoff = %6.4f, highpass cutoff = %6.4f)\n',out.BP_LP_cutoff,out.BP_HP_cutoff);
    end
elseif out.BP==1 && strcmp(out.BP_type,'Butterworth')
    if out.BP_LP_cutoff==-1
        fprintf(logfile_fid,'The data were highpass filtered with a Butterworth filter of order %u (cutoff = %6.4f)\n',out.butter_ord,out.BP_HP_cutoff);
    elseif out.BP_HP_cutoff==-1
        fprintf(logfile_fid,'The data were lowpass filtered with a Butterworth filter of order %u (cutoff = %6.4f)\n',out.butter_ord,out.BP_LP_cutoff);
    else
        fprintf(logfile_fid,'The data were bandpass filtered with a Butterworth filter of order %u (lowpass cutoff = %6.4f, highpass cutoff = %6.4f)\n',out.butter_ord,out.BP_HP_cutoff,out.BP_LP_cutoff);
    end
elseif out.BP==1 && strcmp(out.BP_type,'FSL')
    if out.BP_LP_cutoff==-1
        fprintf(logfile_fid,'The data were filtered with FSL''s non-linear highpass filter (cutoff = %6.4f)\n',out.BP_HP_cutoff);
    elseif out.BP_HP_cutoff==-1
        fprintf(logfile_fid,'The data were filtered with FSL''s Gaussian linear lowpass filter (cutoff = %6.4f)\n',out.BP_LP_cutoff);
    else
        fprintf(logfile_fid,'The data were filtered with FSL''s non-linear highpass and Gaussian linear lowpass filter (lowpass cutoff = %6.4f, highpass cutoff = %6.4f)\n',out.BP_LP_cutoff,out.BP_HP_cutoff);
    end
end
if out.motcens==1
    fprintf(logfile_fid,'Motion-censoring was performed with FD cutoff = %5.3fmm and DVARS cutoff = %6.3f\n',out.motcens_FD_cutoff,out.motcens_DVARS_cutoff);
    if strcmp(out.motcens_logic,'AND')
        fprintf(logfile_fid,'The criterion for censoring a timepoint was being above threshold for FD AND DVARS\n');
    elseif strcmp(out.motcens_logic,'OR')
        fprintf(logfile_fid,'The criterion for censoring a timepoint was being above threshold for FD OR DVARS\n');
    end
end
fprintf(logfile_fid,'\nData partialing:\n');
if out.motpar_part==1
    fprintf(logfile_fid,'Motion paramaters were partialed from the timeseries data\n');
    if out.motpart1_part==1
        fprintf(logfile_fid,'The t - 1 motion paramaters were partialed from the timeseries data\n');
    end
    if out.motparderiv_part==1
        fprintf(logfile_fid,'The 1st derivatives of the motion paramaters were partialed from the timeseries data\n');
    end
    if out.motparsqr_part==1
        fprintf(logfile_fid,'The squares of the motion paramaters were partialed from the timeseries data\n');
    end
end
if out.globsig_part==1 && out.globsig_calcGNI==0
    fprintf(logfile_fid,'Global signal was partialed from the timeseries data for all participants\n');
    if out.globsigderiv_part==1
        fprintf(logfile_fid,'The 1st derivative of the global signal was partialed from the timeseries data\n');
    end
elseif out.globsig_part==1 && out.globsig_calcGNI==1
    fprintf(logfile_fid,'Global signal was partialed from the timeseries data for participants with GNI < 3\n');
    if out.globsigderiv_part==1
        fprintf(logfile_fid,'The 1st derivative of the global signal was partialed from the timeseries data for participants with GNI < 3\n');
    end
end
if out.erode_masks==1 && out.WMsig_part==1
    fprintf(logfile_fid,'The white matter mask was (3d) eroded with a sphere of radius %3.0f voxels\n',out.erode_radius);
end
if out.WMsig_part==1 && out.PCA_WMvent==1
    fprintf(logfile_fid,'First 5 principal components of white matter signal (from entire WM mask) were partialed from the timeseries data\n');
elseif out.WMsig_part==1
    if strcmp(out.WMmask_scope,'Entire Mask')
        fprintf(logfile_fid,'White matter signal (from entire WM mask) was partialed from the timeseries data\n');
    elseif out.WMsig_part==1 && strcmp(out.WMmask_scope,'Local mask')
        fprintf(logfile_fid,'White matter signal (from WM within a ~45mm radius sphere around each voxel) was partialed from the timeseries data\n');
    end
end
if out.WMsig_part==1 && out.WMsigderiv_part==1
    fprintf(logfile_fid,'The 1st derivative of the white matter signal(s) was partialed from the timeseries data\n');
end
if out.erode_masks==1 && out.ventsig_part==1
    fprintf(logfile_fid,'The ventricle mask was (3d) eroded with a sphere of radius %3.0f voxels\n',out.erode_radius);
end
if out.ventsig_part==1
    if out.PCA_WMvent==1
        fprintf(logfile_fid,'First 5 principal components of ventricular signal were partialed from the timeseries data\n');
    else
        fprintf(logfile_fid,'Ventricular signal was partialed from the timeseries data\n');
    end
    if out.ventsigderiv_part==1
        fprintf(logfile_fid,'The 1st derivative of the ventricular signal(s) was partialed from the timeseries data\n');
    end
end
if out.rem_AR==1
    fprintf(logfile_fid,'Autocorrelations were removed from the timeseries (modelled using %3.0f past timepoints)\n',out.num_hist_tp);
end
if out.deconv==1
    if strcmp(out.whendeconv,'before')
        fprintf(logfile_fid,'Timeseries were deconvolved before extraction from ROIs using the %s method\n',out.deconv_method);
    elseif strcmp(out.whendeconv,'after')
        fprintf(logfile_fid,'Timeseries were deconvolved after extraction from ROIs using the %s method\n',out.deconv_method);
    end
end
fprintf(logfile_fid,'The summary measure used to extract the timeseries from each ROI was the %s\n',out.ROI_summary_measure);
fprintf(logfile_fid,'\nParticipant specific information:\n');

if use_parfor
    subs                 = out.subs;
    BP                   = out.BP;
    BP_HP_cutoff         = out.BP_HP_cutoff;
    BP_LP_cutoff         = out.BP_LP_cutoff;
    BP_type              = out.BP_type;
    butter_ord           = out.butter_ord;
    deconv               = out.deconv;
    deconv_method        = out.deconv_method;
    despike              = out.despike;
    detr                 = out.detr;
    detr_ord             = out.detr_ord;
    erode_masks          = out.erode_masks;
    erode_radius         = out.erode_radius;
    globsig_calcGCOR     = out.globsig_calcGCOR;
    globsig_calcGNI      = out.globsig_calcGNI;
    globsig_part         = out.globsig_part;
    globsigderiv_part    = out.globsigderiv_part;
    MC                   = out.MC;
    MinpercentROIvox     = out.MinpercentROIvox;
    MinROIvox            = out.MinROIvox;
    motcens              = out.motcens;
    motcens_DVARS_cutoff = out.motcens_DVARS_cutoff;
    motcens_FD_cutoff    = out.motcens_FD_cutoff;
    motcens_logic        = out.motcens_logic;
    motpar_part          = out.motpar_part;
    motparderiv_part     = out.motparderiv_part;
    motparsqr_part       = out.motparsqr_part;
    motpart1_part        = out.motpart1_part;
    nROI                 = out.nROI;
    nTP                  = out.nTP;
    num_hist_tp          = out.num_hist_tp;
    parc_space           = out.parc_space;
    PCA_WMvent           = out.PCA_WMvent;
    rem_AR               = out.rem_AR;
    ROI_num_labels       = out.ROI_num_labels;
    ROI_summary_measure  = out.ROI_summary_measure;
    ROIextord            = out.ROIextord;
    saveintermed         = out.saveintermed;
    ST_correct           = out.ST_correct;
    ST_ord               = out.ST_ord;
    TR                   = out.TR;
    useavail             = out.useavail;
    vent_space           = out.vent_space;
    ventsig_part         = out.ventsig_part;
    ventsigderiv_part    = out.ventsigderiv_part;
    whendeconv           = out.whendeconv;
    WM_space             = out.WM_space;
    WMmask_scope         = out.WMmask_scope;
    WMsig_part           = out.WMsig_part;
    WMsigderiv_part      = out.WMsigderiv_part;
    fclose(logfile_fid);
    parforprog           = ProgressBar(length(subs));
    parfor currsub = 1:length(subs)
        logfile_fid       = fopen(preproc_logfile_outname,'a');
        sub               = subs{currsub};
        func_filename     = strrep(func_base_filename,'SUBNUM',sub);
        parc_filename     = strrep(parc_base_filename,'SUBNUM',sub);
        funcmask_filename = strrep(funcmask_base_filename,'SUBNUM',sub);
        if strcmp('WM_base_filename',' ')
            WM_filename = ' ';
        else
            WM_filename = strrep(WM_base_filename,'SUBNUM',sub);
        end
        if strcmp('vent_base_filename',' ')
            vent_filename = ' ';
        else
            vent_filename = strrep(vent_base_filename,'SUBNUM',sub);
        end
        if MC==0 && (motpar_part==1 || motcens==1) && ~strcmp('par_base_filename',' ')
            par_filename = strrep(par_base_filename,'SUBNUM',sub);
        else
            par_filename = ' ';
        end
        if strcmp(parc_space,'Structural')
            parc_warp_filename = strrep(struct2func_regmat_base_filename,'SUBNUM',sub);
        elseif strcmp(parc_space,'Standard')
            parc_warp_filename = strrep(standard2func_regmat_base_filename,'SUBNUM',sub);
        else
            parc_warp_filename = ' ';
        end
        if strcmp(WM_space,'Structural')
            WM_warp_filename = strrep(struct2func_regmat_base_filename,'SUBNUM',sub);
        elseif strcmp(WM_space,'Standard')
            WM_warp_filename = strrep(standard2func_regmat_base_filename,'SUBNUM',sub);
        else
            WM_warp_filename = ' ';
        end
        if strcmp(vent_space,'Structural')
            vent_warp_filename = strrep(struct2func_regmat_base_filename,'SUBNUM',sub);
        elseif strcmp(vent_space,'Standard')
            vent_warp_filename = strrep(standard2func_regmat_base_filename,'SUBNUM',sub);
        else
            vent_warp_filename = ' ';
        end
        
        if exist(func_filename,'file') && exist(parc_filename,'file') && exist(funcmask_filename,'file')
            [ts{currsub,1},num_censored_vols(currsub,1),num_abovethresh_vols_FD(currsub,1),num_abovethresh_vols_DVARS(currsub,1),GCOR(currsub,1),mean_FD(currsub,1),mean_DVARS(currsub,1),deconv_ts{currsub,1}] = preprocess_for_graph(BP,BP_HP_cutoff,BP_LP_cutoff,BP_type,butter_ord,deconv,deconv_method,despike,detr,detr_ord,erode_masks,erode_radius,funcmask_filename,globsig_calcGCOR,globsig_calcGNI,globsig_part,globsigderiv_part,logfile_fid,MC,MinpercentROIvox,MinROIvox,motcens,motcens_DVARS_cutoff,motcens_FD_cutoff,motcens_logic,motpar_part,motparderiv_part,motparsqr_part,motpart1_part,nROI,nTP,num_hist_tp,func_filename,par_filename,parc_filename,parc_space,parc_warp_filename,PCA_WMvent,rem_AR,ROI_num_labels,ROI_summary_measure,ROIextord,saveintermed,ST_correct,ST_ord,sub,TR,useavail,vent_filename,vent_space,vent_warp_filename,ventsig_part,ventsigderiv_part,whendeconv,WM_filename,WM_space,WM_warp_filename,WMmask_scope,WMsig_part,WMsigderiv_part); %#ok<*PFOUS>
        else
            fprintf(logfile_fid,'Missing files for %s\n',sub);
            ts{currsub,1}                         = NaN;
            deconv_ts{currsub,1}                  = NaN;
            num_censored_vols(currsub,1)          = NaN;
            num_abovethresh_vols_FD(currsub,1)    = NaN;
            num_abovethresh_vols_DVARS(currsub,1) = NaN;
            GCOR(currsub,1)                       = NaN;
            mean_FD(currsub,1)                    = NaN;
            mean_DVARS(currsub,1)                 = NaN;
        end
        fclose(logfile_fid);
        parforprog.progress %#ok<PFBNS>
    end
    out.ts                         = ts;
    out.deconv_ts                  = deconv_ts;
    out.num_censored_vols          = num_censored_vols;
    out.num_abovethresh_vols_FD    = num_abovethresh_vols_FD;
    out.num_abovethresh_vols_DVARS = num_abovethresh_vols_DVARS;
    out.GCOR                       = GCOR;
    out.mean_FD                    = mean_FD;
    out.mean_DVARS                 = mean_DVARS;
    logfile_fid                    = fopen(preproc_logfile_outname,'a');
    parforprog.stop;
else
    progressbar('Total Progress')
    for currsub = 1:length(out.subs)
        sub               = out.subs{currsub};
        func_filename     = strrep(func_base_filename,'SUBNUM',sub);
        parc_filename     = strrep(parc_base_filename,'SUBNUM',sub);
        funcmask_filename = strrep(funcmask_base_filename,'SUBNUM',sub);
        if exist('WM_base_filename','var')
            WM_filename = strrep(WM_base_filename,'SUBNUM',sub);
        else
            WM_filename = ' ';
        end
        if exist('vent_base_filename','var')
            vent_filename = strrep(vent_base_filename,'SUBNUM',sub);
        else
            vent_filename = ' ';
        end
        if out.MC==0 && (out.motpar_part==1 || out.motcens==1)
            par_filename = strrep(par_base_filename,'SUBNUM',sub);
        else
            par_filename = ' ';
        end
        if strcmp(out.parc_space,'Structural')
            parc_warp_filename = strrep(struct2func_regmat_base_filename,'SUBNUM',sub);
        elseif strcmp(out.parc_space,'Standard')
            parc_warp_filename = strrep(standard2func_regmat_base_filename,'SUBNUM',sub);
        else
            parc_warp_filename = ' ';
        end
        if strcmp(out.WM_space,'Structural')
            WM_warp_filename = strrep(struct2func_regmat_base_filename,'SUBNUM',sub);
        elseif strcmp(out.WM_space,'Standard')
            WM_warp_filename = strrep(standard2func_regmat_base_filename,'SUBNUM',sub);
        else
            WM_warp_filename = ' ';
        end
        if strcmp(out.vent_space,'Structural')
            vent_warp_filename = strrep(struct2func_regmat_base_filename,'SUBNUM',sub);
        elseif strcmp(out.vent_space,'Standard')
            vent_warp_filename = strrep(standard2func_regmat_base_filename,'SUBNUM',sub);
        else
            vent_warp_filename = ' ';
        end
        if exist(func_filename,'file') && exist(parc_filename,'file') && exist(funcmask_filename,'file')
            [out.ts{currsub,1},out.num_censored_vols(currsub,1),out.num_abovethresh_vols_FD(currsub,1),out.num_abovethresh_vols_DVARS(currsub,1),out.GCOR(currsub,1),out.mean_FD(currsub,1),out.mean_DVARS(currsub,1),out.deconv_ts{currsub,1}] = preprocess_for_graph(out.BP,out.BP_HP_cutoff,out.BP_LP_cutoff,out.BP_type,out.butter_ord,out.deconv,out.deconv_method,out.despike,out.detr,out.detr_ord,out.erode_masks,out.erode_radius,funcmask_filename,out.globsig_calcGCOR,out.globsig_calcGNI,out.globsig_part,out.globsigderiv_part,logfile_fid,out.MC,out.MinpercentROIvox,out.MinROIvox,out.motcens,out.motcens_DVARS_cutoff,out.motcens_FD_cutoff,out.motcens_logic,out.motpar_part,out.motparderiv_part,out.motparsqr_part,out.motpart1_part,out.nROI,out.nTP,out.num_hist_tp,func_filename,par_filename,parc_filename,out.parc_space,parc_warp_filename,out.PCA_WMvent,out.rem_AR,out.ROI_num_labels,out.ROI_summary_measure,out.ROIextord,out.saveintermed,out.ST_correct,out.ST_ord,sub,out.TR,out.useavail,vent_filename,out.vent_space,vent_warp_filename,out.ventsig_part,out.ventsigderiv_part,out.whendeconv,WM_filename,out.WM_space,WM_warp_filename,out.WMmask_scope,out.WMsig_part,out.WMsigderiv_part);
        else
            fprintf(logfile_fid,'Missing files for %s\n',sub);
            out.ts{currsub,1}                         = NaN;
            out.deconv_ts{currsub,1}                  = NaN;
            out.num_censored_vols(currsub,1)          = NaN;
            out.num_abovethresh_vols_FD(currsub,1)    = NaN;
            out.num_abovethresh_vols_DVARS(currsub,1) = NaN;
            out.GCOR(currsub,1)                       = NaN;
            out.mean_FD(currsub,1)                    = NaN;
            out.mean_DVARS(currsub,1)                 = NaN;
        end
        
        save(out.preproc_outname,'out');
        prog = currsub/length(out.subs);
        progressbar(prog)
    end
end

if use_parfor
    try
        parpool close
    catch %#ok<CTCH>
        matlabpool close %#ok<DPOOL>
    end
end

fclose(logfile_fid);
save(out.preproc_outname,'out');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Conduct preprocessing for a participant
function [ts,num_censored_vols,num_abovethresh_vols_FD,num_abovethresh_vols_DVARS,GCOR,mean_FD,mean_DVARS,deconv_ts] = preprocess_for_graph(BP,BP_HP_cutoff,BP_LP_cutoff,BP_type,butter_ord,deconv,deconv_method,despike,detr,detr_ord,erode_masks,erode_radius,funcmask_filename,globsig_calcGCOR,globsig_calcGNI,globsig_part,globsigderiv_part,logfile_fid,MC,MinpercentROIvox,MinROIvox,motcens,motcens_DVARS_cutoff,motcens_FD_cutoff,motcens_logic,motpar_part,motparderiv_part,motparsqr_part,motpart1_part,nROI,nTP,num_hist_tp,orig_func_filename,par_filename,parc_filename,parc_space,parc_warp_filename,PCA_WMvent,rem_AR,ROI_num_labels,ROI_summary_measure,ROIextord,saveintermed,ST_correct,ST_ord,sub,TR,useavail,vent_filename,vent_space,vent_warp_filename,ventsig_part,ventsigderiv_part,whendeconv,WM_filename,WM_space,WM_warp_filename,WMmask_scope,WMsig_part,WMsigderiv_part)

% Set up parameters for temporal filter
if BP==1
    if strcmp(BP_type,'Butterworth')
        if BP_HP_cutoff==-1
            [butter_zeds,butter_poles,butter_gain] = butter(butter_ord,BP_LP_cutoff/((1/TR)/2),'low');
            [butter_sos,butter_gain]               = zp2sos(butter_zeds,butter_poles,butter_gain);
            butter_filt                            = dfilt.df2sos(butter_sos,butter_gain);
        elseif BP_LP_cutoff==-1
            [butter_zeds,butter_poles,butter_gain] = butter(butter_ord,BP_HP_cutoff/((1/TR)/2),'high');
            [butter_sos,butter_gain]               = zp2sos(butter_zeds,butter_poles,butter_gain);
            butter_filt                            = dfilt.df2sos(butter_sos,butter_gain);
        else
            [butter_zeds,butter_poles,butter_gain] = butter(butter_ord,[BP_HP_cutoff,BP_LP_cutoff]/((1/TR)/2),'bandpass');
            [butter_sos,butter_gain]               = zp2sos(butter_zeds,butter_poles,butter_gain);
            butter_filt                            = dfilt.df2sos(butter_sos,butter_gain);
        end
        if deconv==1
            [butter_zeds,butter_poles,butter_gain] = butter(butter_ord,BP_HP_cutoff/((1/TR)/2),'high');
            [butter_sos,butter_gain]               = zp2sos(butter_zeds,butter_poles,butter_gain);
            butter_filt_HP                         = dfilt.df2sos(butter_sos,butter_gain);
%             [butter_zeds,butter_poles,butter_gain] = butter(butter_ord,BP_LP_cutoff/((1/TR)/2),'low');
%             [butter_sos,butter_gain] = zp2sos(butter_zeds,butter_poles,butter_gain);
%             butter_filt_LP = dfilt.df2sos(butter_sos,butter_gain);
        end
    elseif strcmp(BP_type,'Ideal')
        filt_ind = calc_IdealFilter(nTP,TR,BP_HP_cutoff,BP_LP_cutoff);
        if deconv==1
            filt_ind_HP = calc_IdealFilter(nTP,TR,BP_HP_cutoff,0);
%             filt_ind_LP = calc_IdealFilter(nTP,TR,0,BP_LP_cutoff);
        end
    elseif strcmp(BP_type,'FSL')
        if BP_HP_cutoff==-1
            HP_sigma = -1;
        else
            % Method suggested on FSL listserve
            HP_sigma = ((1./BP_HP_cutoff)./sqrt(8*log(2)))./TR;
%             % Alternative method of calculating sigma
%             HP_sigma = (1/TR)./(2*pi*BP_HP_cutoff);
        end
        if BP_LP_cutoff==-1
            LP_sigma = -1;
        else
            % Method suggested on FSL listserve
            LP_sigma = ((1./BP_LP_cutoff)./sqrt(8*log(2)))./TR;
%             % Alternative method of calculating sigma
%             HP_sigma = (1/TR)./(2*pi*BP_LP_cutoff);
        end
    end
end

% Correct for slice-timing
if ST_correct==1
    switch ST_ord
        case 'Sequential Up'
            slice_order = '';
        case 'Sequential Down'
            slice_order = '--down';
        case 'Interleaved'
            slice_order = '--odd';
    end
    if strcmp(useavail,'Yes, do not run step if output file exists')
        if ~exist(strrep(orig_func_filename,'.nii','_st.nii'),'file')
            call_fsl_edit(['!slicetimer -i ',strrep(strrep(orig_func_filename,'.gz',''),'.nii',''),' -r ',num2str(TR),' ',slice_order]);
        end
    else
        call_fsl_edit(['!slicetimer -i ',strrep(strrep(orig_func_filename,'.gz',''),'.nii',''),' -r ',num2str(TR),' ',slice_order]);
    end
    orig_func_filename = strrep(orig_func_filename,'.nii','_st.nii');
    if ~exist(orig_func_filename,'file')
        if strcmp(orig_func_filename(end-2:end),'.gz')
            orig_func_filename(end-2:end) = [];
            if ~exist(orig_func_filename,'file')
                error('Output file from slicetiming correction not found');
            end
        elseif strcmp(orig_func_filename(end-3:end),'.nii')
            orig_func_filename = [orig_func_filename,'.gz'];
            if ~exist(orig_func_filename,'file')
                error('Output file from slicetiming correction not found');
            end
        else
            error('Output file from slicetiming correction not found');
        end
    end
end

% Motion correct
if MC==1
    if strcmp(useavail,'Yes, do not run step if output file exists')
        if ~exist(strrep(orig_func_filename,'.nii','_mcf.nii'),'file')
            call_fsl_edit(['!mcflirt -in ',strrep(strrep(orig_func_filename,'.gz',''),'.nii',''),' -plots']);
        end
    else
        call_fsl_edit(['!mcflirt -in ',strrep(strrep(orig_func_filename,'.gz',''),'.nii',''),' -plots']);
    end
    orig_func_filename = strrep(orig_func_filename,'.nii','_mcf.nii');
    if ~exist(orig_func_filename,'file')
        if strcmp(orig_func_filename(end-2:end),'.gz')
            orig_func_filename(end-2:end) = [];
            if ~exist(orig_func_filename,'file')
                error('Output file from motion correction not found');
            end
        elseif strcmp(orig_func_filename(end-3:end),'.nii')
            orig_func_filename = [orig_func_filename,'.gz'];
            if ~exist(orig_func_filename,'file')
                error('Output file from motion correction not found');
            end
        else
            error('Output file from motion correction not found');
        end;
    end
    par_filename = [strrep(strrep(orig_func_filename,'.gz',''),'.nii',''),'.par'];
end
func_filename = orig_func_filename;

% Detrend data
if detr==1
    if strcmp(useavail,'Yes, do not run step if output file exists')
        if ~exist(strrep(func_filename,'.nii',['_detr',num2str(detr_ord),'.nii']),'file')
            detrend_image(func_filename,funcmask_filename,detr_ord);
        end
    else
        detrend_image(func_filename,funcmask_filename,detr_ord);
    end
    func_filename = strrep(func_filename,'.nii',['_detr',num2str(detr_ord),'.nii']);
end

% Warp ROI atlas into functional space
if ~strcmp(parc_space,'Functional')
    out_parc_filename = strrep(parc_filename,'.nii','_warp2func.nii');
    if strcmp(useavail,'Yes, do not run step if output file exists')
        if ~exist(out_parc_filename,'file')
            if any(strfind(parc_warp_filename,'.mat'))
                call_fsl_edit(['!flirt -interp nearestneighbour -in ',parc_filename,' -ref ',funcmask_filename,' -out ',out_parc_filename,' -init ',parc_warp_filename,' -applyxfm']);
            elseif any(strfind(parc_warp_filename,'.nii'))
                call_fsl_edit(['!applywarp --ref=',funcmask_filename,' --in=',parc_filename,' --warp=',parc_warp_filename,' --out=',out_parc_filename,' --interp=nn']);
            end
        end
    else
        if any(strfind(parc_warp_filename,'.mat'))
            call_fsl_edit(['!flirt -interp nearestneighbour -in ',parc_filename,' -ref ',funcmask_filename,' -out ',out_parc_filename,' -init ',parc_warp_filename,' -applyxfm']);
        elseif any(strfind(parc_warp_filename,'.nii'))
            call_fsl_edit(['!applywarp --ref=',funcmask_filename,' --in=',parc_filename,' --warp=',parc_warp_filename,' --out=',out_parc_filename,' --interp=nn']);
        end
    end
    parc_filename = out_parc_filename;
    if ~exist(parc_filename,'file')
        if strcmp(parc_filename(end-2:end),'.gz')
            parc_filename(end-2:end) = [];
            if ~exist(parc_filename,'file')
                error('Warped ROI mask not found');
            end
        elseif strcmp(parc_filename(end-3:end),'.nii')
            parc_filename = [parc_filename,'.gz'];
            if ~exist(parc_filename,'file')
                error('Warped ROI mask not found');
            end
        else
            error('Warped ROI mask not found');
        end
    end
end

% Warp white matter mask into functional space
if ~strcmp(WM_space,'Functional')
    out_WM_filename = strrep(WM_filename,'.nii','_warp2func.nii');
    if strcmp(useavail,'Yes, do not run step if output file exists')
        if ~exist(out_WM_filename,'file')
            if any(strfind(WM_warp_filename,'.mat'))
                call_fsl_edit(['!flirt -interp nearestneighbour -in ',WM_filename,' -ref ',funcmask_filename,' -out ',out_WM_filename,' -init ',WM_warp_filename,' -applyxfm']);
            elseif any(strfind(vent_warp_filename,'.nii'))
                call_fsl_edit(['!applywarp --ref=',funcmask_filename,' --in=',WM_filename,' --warp=',WM_warp_filename,' --out=',out_WM_filename,' --interp=nn']);
            end
        end
    else
        if any(strfind(WM_warp_filename,'.mat'))
            call_fsl_edit(['!flirt -interp nearestneighbour -in ',WM_filename,' -ref ',funcmask_filename,' -out ',out_WM_filename,' -init ',WM_warp_filename,' -applyxfm']);
        elseif any(strfind(vent_warp_filename,'.nii'))
            call_fsl_edit(['!applywarp --ref=',funcmask_filename,' --in=',WM_filename,' --warp=',WM_warp_filename,' --out=',out_WM_filename,' --interp=nn']);
        end
    end
    WM_filename = out_WM_filename;
    if ~exist(WM_filename,'file')
        if strcmp(WM_filename(end-2:end),'.gz')
            WM_filename(end-2:end) = [];
            if ~exist(WM_filename,'file')
                error('Warped white matter mask not found');
            end
        elseif strcmp(WM_filename(end-3:end),'.nii')
            WM_filename = [WM_filename,'.gz'];
            if ~exist(WM_filename,'file')
                error('Warped white matter mask not found');
            end
        else
            error('Warped white matter mask not found');
        end
    end
end

% Warp ventricular mask into functional space
if ~strcmp(vent_space,'Functional')
    out_vent_filename = strrep(vent_filename,'.nii','_warp2func.nii');
    if strcmp(useavail,'Yes, do not run step if output file exists')
        if ~exist(out_vent_filename,'file')
            if any(strfind(vent_warp_filename,'.mat'))
                call_fsl_edit(['!flirt -interp nearestneighbour -in ',vent_filename,' -ref ',funcmask_filename,' -out ',out_vent_filename,' -init ',vent_warp_filename,' -applyxfm']);
            elseif any(strfind(vent_warp_filename,'.nii'))
                call_fsl_edit(['!applywarp --ref=',funcmask_filename,' --in=',vent_filename,' --warp=',vent_warp_filename,' --out=',out_vent_filename,' --interp=nn']);
            end
        end
    else
        if any(strfind(vent_warp_filename,'.mat'))
            call_fsl_edit(['!flirt -interp nearestneighbour -in ',vent_filename,' -ref ',funcmask_filename,' -out ',out_vent_filename,' -init ',vent_warp_filename,' -applyxfm']);
        elseif any(strfind(vent_warp_filename,'.nii'))
            call_fsl_edit(['!applywarp --ref=',funcmask_filename,' --in=',vent_filename,' --warp=',vent_warp_filename,' --out=',out_vent_filename,' --interp=nn']);
        end
    end
    vent_filename = out_vent_filename;
    if ~exist(vent_filename,'file')
        if strcmp(vent_filename(end-2:end),'.gz')
            vent_filename(end-2:end) = [];
            if ~exist(vent_filename,'file')
                error('Warped ventricular mask not found');
            end
        elseif strcmp(vent_filename(end-3:end),'.nii')
            vent_filename = [vent_filename,'.gz'];
            if ~exist(vent_filename,'file')
                error('Warped ventricular mask not found');
            end
        else
            error('Warped ventricular mask not found');
        end
    end
end

% Conduct preprocessing
if strcmp(ROIextord,'After all preprocessing')
    % Wavelet despiking
    if despike==1
        if strcmp(useavail,'Yes, do not run step if output file exists')
            if ~exist(strrep(func_filename,'.nii','_wds.nii'),'file')
                WaveletDespike(func_filename,strrep(strrep(func_filename,'.gz',''),'.nii',''),'verbose',0,'sp',0);
                func_filename = strrep(func_filename,'.nii','_wds.nii');
                orient        = evalc(['!fslorient -getorient ',func_filename]);
                if strcmp(orient(1:12),'NEUROLOGICAL')
                    call_fsl_edit(['!fslswapdim ',func_filename,' -x y z ',func_filename]);
                    call_fsl_edit(['!fslorient -swaporient ',func_filename]);
                end
            else
                func_filename = strrep(func_filename,'.nii','_wds.nii');
                orient = evalc(['!fslorient -getorient ',func_filename]);
                if strcmp(orient(1:12),'NEUROLOGICAL')
                    call_fsl_edit(['!fslswapdim ',func_filename,' -x y z ',func_filename]);
                    call_fsl_edit(['!fslorient -swaporient ',func_filename]);
                end
            end
        else
            WaveletDespike(func_filename,strrep(strrep(func_filename,'.gz',''),'.nii',''),'verbose',0,'sp',0);
            func_filename = strrep(func_filename,'.nii','_wds.nii');
            orient        = evalc(['!fslorient -getorient ',func_filename]);
            if strcmp(orient(1:12),'NEUROLOGICAL')
                call_fsl_edit(['!fslswapdim ',func_filename,' -x y z ',func_filename]);
                call_fsl_edit(['!fslorient -swaporient ',func_filename]);
            end
        end
        if ~exist(func_filename,'file')
            if strcmp(func_filename(end-2:end),'.gz')
                func_filename(end-2:end) = [];
                if ~exist(func_filename,'file')
                    error('Output file from wavelet despiking not found');
                end
            elseif strcmp(func_filename(end-3:end),'.nii')
                func_filename = [func_filename,'.gz'];
                if ~exist(func_filename,'file')
                    error('Output file from wavelet despiking not found');
                end
            else
                error('Output file from wavelet despiking not found');
            end
        end
    end
    
    % Create design matrix of nuisance signals
    desmat = [];
    
    % Extract signal from white matter
    if WMsig_part==1
        if erode_masks==1
            WM_eroded_filename = strrep(WM_filename,'.nii',['_eroded_',num2str(erode_radius),'.nii']);
            if strcmp(useavail,'Yes, do not run step if output file exists')
                if ~exist(WM_eroded_filename,'file')
                    tempmask       = load_nii_gz(WM_filename);
                    [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                    save_nii_gz(tempmask,WM_eroded_filename);
                end
            else
                tempmask       = load_nii_gz(WM_filename);
                [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                save_nii_gz(tempmask,WM_eroded_filename);
            end
            WM_filename = WM_eroded_filename;
        end
        if strcmp(WMmask_scope,'Entire Mask')
            if PCA_WMvent==1
                WM_extract_filename = [strrep(strrep(WM_filename,'.gz',''),'.nii',''),'_EIG.txt'];
                if ~exist(WM_extract_filename,'file')
                    call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename,' --eig --order=5']);
                else
                    vernum = 0;
                    while 1
                        if exist(strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                            vernum = vernum+1;
                        else
                            WM_extract_filename = strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                            break
                        end
                    end
                    call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename,' --eig --order=5']);
                end
                WM_preds = dlmread(WM_extract_filename);
                WM_preds = spm_detrend(WM_preds,detr_ord);
                desmat   = [desmat,WM_preds];
                
                if WMsigderiv_part==1
                    WM1stderiv_preds = WM_preds(3:end,:)-WM_preds(1:(end-2),:);
                    WM1stderiv_preds = [WM1stderiv_preds(1,:);WM1stderiv_preds;WM1stderiv_preds(end,:)];
                    WM1stderiv_preds = spm_detrend(WM1stderiv_preds,detr_ord);
                    desmat           = [desmat,WM1stderiv_preds];
                end
            else
                WM_extract_filename = [strrep(strrep(WM_filename,'.gz',''),'.nii',''),'.txt'];
                if ~exist(WM_extract_filename,'file')
                    call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename]);
                else
                    vernum = 0;
                    while 1
                        if exist(strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                            vernum = vernum+1;
                        else
                            WM_extract_filename = strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                            break
                        end
                    end
                    call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename]);
                end
                WM_pred = dlmread(WM_extract_filename);
                WM_pred = spm_detrend(WM_pred,detr_ord);
                desmat  = [desmat,WM_pred];
                
                if WMsigderiv_part==1
                    WM1stderiv_pred = WM_pred(3:end)-WM_pred(1:(end-2));
                    WM1stderiv_pred = [WM1stderiv_pred(1);WM1stderiv_pred;WM1stderiv_pred(end)];
                    WM1stderiv_pred = spm_detrend(WM1stderiv_pred,detr_ord);
                    desmat          = [desmat,WM1stderiv_pred];
                end
            end
        else
            full_WM_mask              = load_nii_gz(WM_filename);
            WM_fmask_inds             = find(full_WM_mask.img==1);
            WM_rad                    = round(45/mean(full_WM_mask.hdr.dime.pixdim(2:4)));
            WM_fmask_dims             = size(full_WM_mask.img);
            clear full_WM_mask
            [vx,vy,vz]                = meshgrid(-WM_rad:WM_rad);
            V                         = sqrt((vx.^2)+(vy.^2)+(vz.^2));
            V(V<=WM_rad)              = 1;
            V(V>WM_rad)               = 0;
            [WM_sp_x,WM_sp_y,WM_sp_z] = ind2sub(size(V),find(V==1));
            WM_sp_x                   = WM_sp_x-WM_rad-1;
            WM_sp_y                   = WM_sp_y-WM_rad-1;
            WM_sp_z                   = WM_sp_z-WM_rad-1;
        end
    end
    
    % Extract signal from ventricles
    if ventsig_part==1
        if erode_masks==1
            vent_eroded_filename = strrep(vent_filename,'.nii',['_eroded_',num2str(erode_radius),'.nii']);
            if strcmp(useavail,'Yes, do not run step if output file exists')
                if ~exist(vent_eroded_filename,'file')
                    tempmask       = load_nii_gz(vent_filename);
                    [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                    save_nii_gz(tempmask,vent_eroded_filename);
                end
            else
                tempmask       = load_nii_gz(vent_filename);
                [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                save_nii_gz(tempmask,vent_eroded_filename);
            end
            vent_filename = vent_eroded_filename;
        end
        if PCA_WMvent==1
            vent_extract_filename = [strrep(strrep(vent_filename,'.gz',''),'.nii',''),'_EIG.txt'];
            if ~exist(vent_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename,' --eig --order=5']);
            else
                vernum = 0;
                while 1
                    if exist(strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        vent_extract_filename = strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename,' --eig --order=5']);
            end
            vent_preds = dlmread(vent_extract_filename);
            vent_preds = spm_detrend(vent_preds,detr_ord);
            desmat     = [desmat,vent_preds];
            
            if ventsigderiv_part==1
                vent1stderiv_preds = vent_preds(3:end,:)-vent_preds(1:(end-2),:);
                vent1stderiv_preds = [vent1stderiv_preds(1,:);vent1stderiv_preds;vent1stderiv_preds(end,:)];
                vent1stderiv_preds = spm_detrend(vent1stderiv_preds,detr_ord);
                desmat             = [desmat,vent1stderiv_preds];
            end
        else
            vent_extract_filename = [strrep(strrep(vent_filename,'.gz',''),'.nii',''),'.txt'];
            if ~exist(vent_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename]);
            else
                vernum = 0;
                while 1
                    if exist(strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        vent_extract_filename = strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename]);
            end
            vent_pred = dlmread(vent_extract_filename);
            vent_pred = spm_detrend(vent_pred,detr_ord);
            desmat    = [desmat,vent_pred];
            
            if ventsigderiv_part==1
                vent1stderiv_pred = vent_pred(3:end)-vent_pred(1:(end-2));
                vent1stderiv_pred = [vent1stderiv_pred(1);vent1stderiv_pred;vent1stderiv_pred(end)];
                vent1stderiv_pred = spm_detrend(vent1stderiv_pred,detr_ord);
                desmat            = [desmat,vent1stderiv_pred];
            end
        end
    end
    
    % Extract global signal & calculate GNI
    if globsig_part==1
        if globsig_calcGNI==1
            GNI = Rest2GNI_edit(func_filename,funcmask_filename,100);
            fprintf(logfile_fid,'Original GNI for %s = %6.3f\n',sub,GNI);
        else
            GNI = 0;
        end
        if GNI<=3
            funcmask_extract_filename = [strrep(strrep(funcmask_filename,'.gz',''),'.nii',''),'.txt'];
            if ~exist(funcmask_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',funcmask_extract_filename,' -m ',funcmask_filename]);
            else
                vernum = 0;
                while 1
                    if exist(strrep(funcmask_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        funcmask_extract_filename = strrep(funcmask_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',funcmask_extract_filename,' -m ',funcmask_filename]);
            end
            global_pred = dlmread(funcmask_extract_filename);
            global_pred = spm_detrend(global_pred,detr_ord);
            desmat      = [desmat,global_pred];
            
            if globsigderiv_part==1
                global1stderiv_pred = global_pred(3:end)-global_pred(1:(end-2));
                global1stderiv_pred = [global1stderiv_pred(1);global1stderiv_pred;global1stderiv_pred(end)];
                global1stderiv_pred = spm_detrend(global1stderiv_pred,detr_ord);
                desmat              = [desmat,global1stderiv_pred];
            end
        end
    end
    
    % Extract motion parameters
    if motpar_part==1
        mot_pred = dlmread(par_filename);
        for par = 1:size(mot_pred,2)
            mot_pred(:,par) = spm_detrend(mot_pred(:,par),detr_ord);
        end
        desmat = [desmat,mot_pred];
        if motpart1_part==1
            mot1back_pred = [mot_pred(2:end,:);mot_pred(end,:)];
            for par = 1:size(mot_pred,2)
                mot1back_pred(:,par) = spm_detrend(mot1back_pred(:,par),detr_ord);
            end
            desmat = [desmat,mot1back_pred];
        end
        if motparderiv_part==1
            mot1stderiv_pred = mot_pred(3:end,:)-mot_pred(1:(end-2),:);
            mot1stderiv_pred = [mot1stderiv_pred(1,:);mot1stderiv_pred;mot1stderiv_pred(end,:)];
            for par = 1:size(mot_pred,2)
                mot1stderiv_pred(:,par) = spm_detrend(mot1stderiv_pred(:,par),detr_ord);
            end
            desmat = [desmat,mot1stderiv_pred];
        end
        if motparsqr_part==1
            motsquare_pred = mot_pred.^2;
            for par = 1:size(mot_pred,2)
                motsquare_pred(:,par) = spm_detrend(motsquare_pred(:,par),detr_ord);
            end
            desmat = [desmat,motsquare_pred];
            if motpart1_part==1
                mot1backsquare_pred = mot1back_pred.^2;
                for par = 1:size(mot_pred,2)
                    mot1backsquare_pred(:,par) = spm_detrend(mot1backsquare_pred(:,par),detr_ord);
                end
                desmat = [desmat,mot1backsquare_pred];
            end
        end
    end
    
    % Load functional data & mask
    func_data     = load_nii_gz(func_filename);
    funcmask      = load_nii_gz(funcmask_filename);
    func_data.img = double(func_data.img);
    funcmask.img  = double(funcmask.img);
    
    % Preallocate output data
    out_func_data.img = zeros(size(func_data.img));
    if deconv==1
        out_func_data_deconv.img = zeros(size(func_data.img));
    end
    
    % Partial nuisance variance then apply temporal filter
    if BP==1 && (WMsig_part==1 || ventsig_part==1 || (globsig_part==1 && GNI<=3) || motpar_part==1)
        desmat      = [ones(size(desmat,1),1),desmat];
        orig_desmat = desmat;
        if rank(desmat)~=min(size(desmat))
            error(['Initial desmat to partial out unwanted variance for ',sub,' is not full rank']);
        end
        func_filename = strrep(func_filename,'.nii','_remnuisance_bpfilt.nii');
        if deconv==1
            func_filename_deconv = strrep(func_filename,'bpfilt','hpfilt');
        end
        Q = desmat*pinv(full(desmat));
        for x = 1:size(func_data.img,1)
            for y = 1:size(func_data.img,2)
                for z = 1:size(func_data.img,3)
                    if funcmask.img(x,y,z)==1
                        if strcmp(WMmask_scope,'Local Mask')
                            curr_WM_sp_x = WM_sp_x+x;
                            curr_WM_sp_y = WM_sp_y+y;
                            curr_WM_sp_z = WM_sp_z+z;
                            inrangeinds  = curr_WM_sp_x>0 & curr_WM_sp_x<=size(func_data.img,1) & curr_WM_sp_y>0 & curr_WM_sp_y<=size(func_data.img,2) & curr_WM_sp_z>0 & curr_WM_sp_z<=size(func_data.img,3);
                            WM_sp_inds   = sub2ind(WM_fmask_dims,curr_WM_sp_x(inrangeinds),curr_WM_sp_y(inrangeinds),curr_WM_sp_z(inrangeinds));
                            loc_WM_inds  = intersect(WM_fmask_inds,WM_sp_inds);
                            if ~isempty(loc_WM_inds)
                                for tp = 1:size(func_data.img,4)
                                    tp_func           = func_data.img(:,:,:,tp);
                                    loc_WM_pred(tp,1) = sum(tp_func(loc_WM_inds))/length(loc_WM_inds); %#ok<*AGROW>
                                end
                                loc_WM_pred = spm_detrend(loc_WM_pred,detr_ord);
                                desmat      = [orig_desmat,loc_WM_pred];
                                if WMsigderiv_part==1
                                    loc_WM1stderiv_pred = loc_WM_pred(3:end)-loc_WM_pred(1:(end-2));
                                    loc_WM1stderiv_pred = [loc_WM1stderiv_pred(1);loc_WM1stderiv_pred;loc_WM1stderiv_pred(end)];
                                    loc_WM1stderiv_pred = spm_detrend(loc_WM1stderiv_pred,detr_ord);
                                    desmat              = [desmat,loc_WM1stderiv_pred];
                                end
                            end
                            Q = desmat*pinv(full(desmat));
                        end
                        tempdata = squeeze(func_data.img(x,y,z,:))-(Q*squeeze(func_data.img(x,y,z,:)));
                        if any(isnan(tempdata))
                            out_func_data.img(x,y,z,:) = NaN;
                            if deconv==1
                                out_func_data_deconv.img(x,y,z,:) = NaN;
                            end
                        else
                            tempdata = tempdata-mean(tempdata);
                            if deconv==1
                                tempdata_deconv = tempdata;
                            end
                            if exist('HP_sigma','var')
                                out_func_data.img(x,y,z,:) = fsl_temporal_filt(tempdata,HP_sigma,LP_sigma);
                                if deconv==1
                                    out_func_data_deconv.img(x,y,z,:) = fsl_temporal_filt(tempdata_deconv,HP_sigma,-1);
                                end
                            elseif exist('butter_filt','var')
                                tempdata_padded            = [zeros(15,size(tempdata,2));tempdata;zeros(15,size(tempdata,2))];
                                temp_butter                = filter(butter_filt,tempdata_padded);
                                temp_butter(end:-1:1)      = filter(butter_filt,temp_butter(end:-1:1));
                                out_func_data.img(x,y,z,:) = temp_butter(16:(end-15));
                                if deconv==1
                                    tempdata_padded                   = [zeros(15,size(tempdata_deconv,2));tempdata_deconv;zeros(15,size(tempdata_deconv,2))];
                                    temp_butter                       = filter(butter_filt_HP,tempdata_padded);
                                    temp_butter(end:-1:1)             = filter(butter_filt_HP,temp_butter(end:-1:1));
                                    out_func_data_deconv.img(x,y,z,:) = temp_butter(16:(end-15));
                                end
                            elseif exist('filt_ind','var')
                                padlength                  = rest_nextpow2_one35(nTP);
                                tempdata                   = tempdata-repmat(mean(tempdata),[size(tempdata,1),1]);
                                tempdata                   = [tempdata;zeros(padlength-nTP,size(tempdata,2))];
                                freq                       = fft(tempdata);
                                freq(filt_ind,:)           = 0;
                                tempdata                   = ifft(freq);
                                out_func_data.img(x,y,z,:) = tempdata(1:nTP,:);
                                if deconv==1
                                    tempdata_deconv                   = tempdata_deconv-repmat(mean(tempdata_deconv),[size(tempdata_deconv,1),1]);
                                    tempdata_deconv                   = [tempdata_deconv;zeros(padlength-nTP,size(tempdata_deconv,2))];
                                    freq                              = fft(tempdata_deconv);
                                    freq(filt_ind_HP,:)               = 0;
                                    tempdata_deconv                   = ifft(freq);
                                    out_func_data_deconv.img(x,y,z,:) = tempdata_deconv(1:nTP,:);
                                end
                            else
                                out_func_data.img(x,y,z,:) = tempdata;
                                if deconv==1
                                    out_func_data_deconv.img(x,y,z,:) = tempdata_deconv;
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif BP==1
        func_filename = strrep(func_filename,'.nii','_bpfilt.nii');
        if deconv==1
            func_filename_deconv = strrep(func_filename,'.nii','_hpfilt.nii');
        end
        for x = 1:size(func_data.img,1)
            for y = 1:size(func_data.img,2)
                for z = 1:size(func_data.img,3)
                    if funcmask.img(x,y,z)==1
                        tempdata = squeeze(func_data.img(x,y,z,:));
                        if any(isnan(tempdata))
                            out_func_data.img(x,y,z,:) = NaN;
                            if deconv==1
                                out_func_data_deconv.img(x,y,z,:) = NaN;
                            end
                        else
                            tempdata = tempdata-mean(tempdata);
                            if deconv==1
                                tempdata_deconv = tempdata;
                            end
                            if exist('HP_sigma','var')
                                out_func_data.img(x,y,z,:) = fsl_temporal_filt(tempdata,HP_sigma,LP_sigma);
                                if deconv==1
                                    out_func_data_deconv.img(x,y,z,:) = fsl_temporal_filt(tempdata_deconv,HP_sigma,-1);
                                end
                            elseif exist('butter_filt','var')
                                tempdata_padded            = [zeros(15,size(tempdata,2));tempdata;zeros(15,size(tempdata,2))];
                                temp_butter                = filter(butter_filt,tempdata_padded);
                                temp_butter(end:-1:1)      = filter(butter_filt,temp_butter(end:-1:1));
                                out_func_data.img(x,y,z,:) = temp_butter(16:(end-15));
                                if deconv==1
                                    tempdata_padded                   = [zeros(15,size(tempdata_deconv,2));tempdata_deconv;zeros(15,size(tempdata_deconv,2))];
                                    temp_butter                       = filter(butter_filt_HP,tempdata_padded);
                                    temp_butter(end:-1:1)             = filter(butter_filt_HP,temp_butter(end:-1:1));
                                    out_func_data_deconv.img(x,y,z,:) = temp_butter(16:(end-15));
                                end
                            elseif exist('filt_ind','var')
                                padlength                  = rest_nextpow2_one35(nTP);
                                tempdata                   = [tempdata;zeros(padlength-nTP,size(tempdata,2))];
                                freq                       = fft(tempdata);
                                freq(filt_ind,:)           = 0;
                                tempdata                   = ifft(freq);
                                out_func_data.img(x,y,z,:) = tempdata(1:nTP,:);
                                if deconv==1
                                    tempdata_deconv                   = [tempdata_deconv;zeros(padlength-nTP,size(tempdata_deconv,2))];
                                    freq                              = fft(tempdata_deconv);
                                    freq(filt_ind_HP,:)               = 0;
                                    tempdata_deconv                   = ifft(freq);
                                    out_func_data_deconv.img(x,y,z,:) = tempdata_deconv(1:nTP,:);
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif WMsig_part==1 || ventsig_part==1 || (globsig_part==1 && GNI<=3) || motpar_part==1
        desmat      = [ones(size(desmat,1),1),desmat];
        orig_desmat = desmat;
        if rank(desmat)~=min(size(desmat))
            error(['Initial desmat to partial out unwanted variance for ' sub ' is not full rank']);
        end
        func_filename = strrep(func_filename,'.nii','_remnuisance.nii');
        Q             = desmat*pinv(full(desmat));
        for x = 1:size(func_data.img,1)
            for y = 1:size(func_data.img,2)
                for z = 1:size(func_data.img,3)
                    if funcmask.img(x,y,z)==1
                        if strcmp(WMmask_scope,'Local Mask')
                            curr_WM_sp_x = WM_sp_x+x;
                            curr_WM_sp_y = WM_sp_y+y;
                            curr_WM_sp_z = WM_sp_z+z;
                            inrangeinds  = curr_WM_sp_x>0 & curr_WM_sp_x<=size(func_data.img,1) & curr_WM_sp_y>0 & curr_WM_sp_y<=size(func_data.img,2) & curr_WM_sp_z>0 & curr_WM_sp_z<=size(func_data.img,3);
                            WM_sp_inds   = sub2ind(WM_fmask_dims,curr_WM_sp_x(inrangeinds),curr_WM_sp_y(inrangeinds),curr_WM_sp_z(inrangeinds));
                            loc_WM_inds  = intersect(WM_fmask_inds,WM_sp_inds);
                            if ~isempty(loc_WM_inds)
                                for tp = 1:size(func_data.img,4)
                                    tp_func           = func_data.img(:,:,:,tp);
                                    loc_WM_pred(tp,1) = sum(tp_func(loc_WM_inds))/length(loc_WM_inds);
                                end
                                loc_WM_pred = spm_detrend(loc_WM_pred,detr_ord);
                                desmat      = [orig_desmat,loc_WM_pred];
                                if WMsigderiv_part==1
                                    loc_WM1stderiv_pred = loc_WM_pred(3:end)-loc_WM_pred(1:(end-2));
                                    loc_WM1stderiv_pred = [loc_WM1stderiv_pred(1);loc_WM1stderiv_pred;loc_WM1stderiv_pred(end)];
                                    loc_WM1stderiv_pred = spm_detrend(loc_WM1stderiv_pred,detr_ord);
                                    desmat              = [desmat,loc_WM1stderiv_pred];
                                end
                            end
                            Q = desmat*pinv(full(desmat));
                        end
                        out_func_data.img(x,y,z,:) = squeeze(func_data.img(x,y,z,:))-(Q*squeeze(func_data.img(x,y,z,:)));
                    end
                end
            end
        end
        if deconv==1
            func_filename_deconv     = func_filename;
            out_func_data_deconv.img = out_func_data.img;
        end
    else
        out_func_data.img = func_data.img;
        if deconv==1
            func_filename_deconv     = func_filename;
            out_func_data_deconv.img = out_func_data.img;
        end
    end
    
    % Save data at current step
    if strcmp(saveintermed,'Yes')
        out_func_data.hdr               = func_data.hdr;
        out_func_data.hdr.dime.dim(2:5) = size(out_func_data.img);
        out_func_data.hdr.dime.cal_max  = max(out_func_data.img(:));
        out_func_data.hdr.dime.cal_min  = min(out_func_data.img(out_func_data.img~=0));
        out_func_data.hdr.dime.datatype = 64;
        out_func_data.hdr.dime.bitpix   = 64;
        save_nii_gz(out_func_data,func_filename);
    end
    
    % Apply motion censoring, redoing all preprocessing afterward
    if motcens==1
        mot_pred = dlmread(par_filename);
        FD       = cat(1,0,sum(abs(diff([(5000-5000*cos(mot_pred(:,1:3))),mot_pred(:,4:6)])),2));
        mean_FD  = mean(FD);
        
        [~,DVARS]  = clct_DVARS(out_func_data.img,funcmask.img);
        DVARS      = cat(1,0,DVARS);
        mean_DVARS = mean(DVARS);
        
        num_abovethresh_vols_FD    = sum(FD>=motcens_FD_cutoff);
        num_abovethresh_vols_DVARS = sum(DVARS>=motcens_DVARS_cutoff);
        
        if strcmp(motcens_logic,'OR')
            temporal_mask = logical(FD<motcens_FD_cutoff & DVARS<motcens_DVARS_cutoff);
        elseif strcmp(motcens_logic,'AND')
            temporal_mask = logical(FD<motcens_FD_cutoff | DVARS<motcens_DVARS_cutoff);
        end
            
        censored_vols = [0;find(temporal_mask==0)];
        for c_vol = 1:(length(censored_vols)-1)
            seg_length = censored_vols(c_vol+1)-censored_vols(c_vol);
            if seg_length<=5
                temporal_mask(censored_vols(c_vol)+1:censored_vols(c_vol+1)) = 0;
            end
        end
        num_censored_vols = sum(temporal_mask==0);
        if num_censored_vols>0
            if (nTP-num_censored_vols>=50) && (num_censored_vols/nTP<=.75)
                func_filename = orig_func_filename;
                
                % Detrend data
                if detr==1
                    detrend_image(func_filename,funcmask_filename,detr_ord,temporal_mask);
                    func_filename = strrep(orig_func_filename,'.nii',['_motcensor_detr',num2str(detr_ord),'.nii']);
                else
                    orig                   = load_nii_gz(func_filename);
                    func_filename          = strrep(orig_func_filename,'.nii','_motcensor.nii');
                    orig.img               = orig.img(:,:,:,temporal_mask);
                    orig.hdr.dime.dim(2:5) = size(orig.img);
                    save_nii_gz(orig,func_filename);
                    clear orig
                end
                
                % Create design matrix of nuisance signals
                desmat = [];
                
                % Extract signal from white matter
                if WMsig_part==1 && strcmp(WMmask_scope,'Entire Mask')
                    if PCA_WMvent==1
                        WM_extract_filename = [strrep(strrep(WM_filename,'.gz',''),'.nii',''),'_EIG_motcensor.txt'];
                        if ~exist(WM_extract_filename,'file')
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename,' --eig --order=5']);
                        else
                            vernum = 0;
                            while 1
                                if exist(strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                                    vernum = vernum+1;
                                else
                                    WM_extract_filename = strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                                    break
                                end
                            end
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename,' --eig --order=5']);
                        end
                        WM_preds = dlmread(WM_extract_filename);
                        WM_preds = spm_detrend(WM_preds,detr_ord);
                        desmat   = [desmat,WM_preds];
                        
                        if WMsigderiv_part==1
                            WM1stderiv_preds = WM_preds(3:end,:)-WM_preds(1:(end-2),:);
                            WM1stderiv_preds = [WM1stderiv_preds(1,:);WM1stderiv_preds;WM1stderiv_preds(end,:)];
                            WM1stderiv_preds = spm_detrend(WM1stderiv_preds,detr_ord);
                            desmat           = [desmat,WM1stderiv_preds];
                        end
                    else
                        WM_extract_filename = [strrep(strrep(WM_filename,'.gz',''),'.nii',''),'_motcensor.txt'];
                        if ~exist(WM_extract_filename,'file')
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename]);
                        else
                            vernum = 0;
                            while 1
                                if exist(strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                                    vernum = vernum+1;
                                else
                                    WM_extract_filename = strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                                    break
                                end
                            end
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename]);
                        end
                        WM_pred = dlmread(WM_extract_filename);
                        WM_pred = spm_detrend(WM_pred,detr_ord);
                        desmat  = [desmat,WM_pred];
                        if WMsigderiv_part==1
                            WM1stderiv_pred = WM_pred(3:end)-WM_pred(1:(end-2));
                            WM1stderiv_pred = [WM1stderiv_pred(1);WM1stderiv_pred;WM1stderiv_pred(end)];
                            WM1stderiv_pred = spm_detrend(WM1stderiv_pred,detr_ord);
                            desmat          = [desmat,WM1stderiv_pred];
                        end
                    end
                end
                
                % Extract signal from ventricles
                if ventsig_part==1
                    if PCA_WMvent==1
                        vent_extract_filename = [strrep(strrep(vent_filename,'.gz',''),'.nii',''),'_EIG_motcensor.txt'];
                        if ~exist(vent_extract_filename,'file')
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename,' --eig --order=5']);
                        else
                            vernum = 0;
                            while 1
                                if exist(strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                                    vernum = vernum+1;
                                else
                                    vent_extract_filename = strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                                    break
                                end
                            end
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename,' --eig --order=5']);
                        end
                        vent_preds = dlmread(vent_extract_filename);
                        vent_preds = spm_detrend(vent_preds,detr_ord);
                        desmat     = [desmat,vent_preds];
                        
                        if ventsigderiv_part==1
                            vent1stderiv_preds = vent_preds(3:end,:)-vent_preds(1:(end-2),:);
                            vent1stderiv_preds = [vent1stderiv_preds(1,:);vent1stderiv_preds;vent1stderiv_preds(end,:)];
                            vent1stderiv_preds = spm_detrend(vent1stderiv_preds,detr_ord);
                            desmat             = [desmat,vent1stderiv_preds];
                        end
                    else
                        vent_extract_filename = [strrep(strrep(vent_filename,'.gz',''),'.nii',''),'_motcensor.txt'];
                        if ~exist(vent_extract_filename,'file')
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename]);
                        else
                            vernum = 0;
                            while 1
                                if exist(strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                                    vernum = vernum+1;
                                else
                                    vent_extract_filename = strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                                    break
                                end
                            end
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename]);
                        end
                        vent_pred = dlmread(vent_extract_filename);
                        vent_pred = spm_detrend(vent_pred,detr_ord);
                        desmat    = [desmat,vent_pred];
                        
                        if ventsigderiv_part==1
                            vent1stderiv_pred = vent_pred(3:end)-vent_pred(1:(end-2));
                            vent1stderiv_pred = [vent1stderiv_pred(1);vent1stderiv_pred;vent1stderiv_pred(end)];
                            vent1stderiv_pred = spm_detrend(vent1stderiv_pred,detr_ord);
                            desmat            = [desmat,vent1stderiv_pred];
                        end
                    end
                end
                
                % Extract global signal & calculate GNI
                if globsig_part==1
                    if globsig_calcGNI==1
                        GNI = Rest2GNI_edit(func_filename,funcmask_filename,1000); % If your machine freezes/slows down tremendously at this step b/c you have low RAM, you can change the 1000 value to 100 (or less) and it will use less RAM
                        fprintf(logfile_fid,'Motion-censored GNI for %s = %6.3f\n',sub,GNI);
                    else
                        GNI = 0;
                    end
                    if GNI<=3
                        funcmask_extract_filename = [strrep(strrep(funcmask_filename,'.gz',''),'.nii',''),'_motcensor.txt'];
                        if ~exist(funcmask_extract_filename,'file')
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',funcmask_extract_filename,' -m ',funcmask_filename]);
                        else
                            vernum = 0;
                            while 1
                                if exist(strrep(funcmask_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                                    vernum = vernum+1;
                                else
                                    funcmask_extract_filename = strrep(funcmask_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                                    break
                                end
                            end
                            call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',funcmask_extract_filename,' -m ',funcmask_filename]);
                        end
                        global_pred = dlmread(funcmask_extract_filename);
                        global_pred = spm_detrend(global_pred,detr_ord);
                        desmat      = [desmat,global_pred];
                        if globsigderiv_part==1
                            global1stderiv_pred = global_pred(3:end)-global_pred(1:(end-2));
                            global1stderiv_pred = [global1stderiv_pred(1);global1stderiv_pred;global1stderiv_pred(end)];
                            global1stderiv_pred = spm_detrend(global1stderiv_pred,detr_ord);
                            desmat              = [desmat,global1stderiv_pred];
                        end
                    end
                end
                
                % Extract motion parameters
                if motpar_part==1
                    mot_pred = dlmread(par_filename);
                    if motpart1_part==1
                        mot1back_pred = [mot_pred(2:end,:);mot_pred(end,:)];
                        mot1back_pred = mot1back_pred(temporal_mask,:);
                        for par = 1:size(mot_pred,2)
                            mot1back_pred(:,par) = spm_detrend(mot1back_pred(:,par),detr_ord);
                        end
                        desmat = [desmat,mot1back_pred];
                    end
                    if motparderiv_part==1
                        mot1stderiv_pred = mot_pred(3:end,:)-mot_pred(1:(end-2),:);
                        mot1stderiv_pred = [mot1stderiv_pred(1,:);mot1stderiv_pred;mot1stderiv_pred(end,:)];
                        mot1stderiv_pred = mot1stderiv_pred(temporal_mask,:);
                        for par = 1:size(mot_pred,2)
                            mot1stderiv_pred(:,par) = spm_detrend(mot1stderiv_pred(:,par),detr_ord);
                        end
                        desmat = [desmat,mot1stderiv_pred];
                    end
                    if motparsqr_part==1
                        motsquare_pred = mot_pred.^2;
                        motsquare_pred = motsquare_pred(temporal_mask,:);
                        for par = 1:size(mot_pred,2)
                            motsquare_pred(:,par) = spm_detrend(motsquare_pred(:,par),detr_ord);
                        end
                        desmat = [desmat,motsquare_pred];
                    end
                    mot_pred = mot_pred(temporal_mask,:);
                    for par = 1:size(mot_pred,2)
                        mot_pred(:,par) = spm_detrend(mot_pred(:,par),detr_ord);
                    end
                    desmat = [desmat,mot_pred];
                end
                
                % Load functional data
                func_data     = load_nii_gz(func_filename);
                func_data.img = double(func_data.img);
                
                % Preallocate output data
                out_func_data.img = zeros(size(func_data.img));
                if deconv==1
                    out_func_data_deconv.img = zeros(size(func_data.img,1),size(func_data.img,2),size(func_data.img,3),nTP);
                end
                
                % Partial nuisance variance
                if WMsig_part==1 || ventsig_part==1 || (globsig_part==1 && GNI<=3) || motpar_part==1
                    desmat      = [ones(size(desmat,1),1),desmat];
                    orig_desmat = desmat;
                    if rank(desmat)~=min(size(desmat))
                        error(['Motion censoring desmat to partial out unwanted variance for ' sub ' is not full rank']);
                    end
                    func_filename = strrep(func_filename,'.nii','_remnuisance.nii');
                    Q             = desmat*pinv(full(desmat));
                    for x = 1:size(func_data.img,1)
                        for y = 1:size(func_data.img,2)
                            for z = 1:size(func_data.img,3)
                                if funcmask.img(x,y,z)==1
                                    if strcmp(WMmask_scope,'Local Mask')
                                        curr_WM_sp_x = WM_sp_x+x;
                                        curr_WM_sp_y = WM_sp_y+y;
                                        curr_WM_sp_z = WM_sp_z+z;
                                        WM_sp_inds   = sub2ind(WM_fmask_dims,[curr_WM_sp_x,curr_WM_sp_y,curr_WM_sp_z]);
                                        loc_WM_inds  = intersect(WM_fmask_inds,WM_sp_inds);
                                        if ~isempty(loc_WM_inds)
                                            for tp = 1:size(func_data.img,4)
                                                tp_func           = func_data.img(:,:,:,tp);
                                                loc_WM_pred(tp,1) = mean(tp_func(loc_WM_inds));
                                            end
                                            loc_WM_pred = spm_detrend(loc_WM_pred,detr_ord);
                                            desmat      = [orig_desmat,loc_WM_pred];
                                            if WMsigderiv_part==1
                                                loc_WM1stderiv_pred = loc_WM_pred(3:end)-loc_WM_pred(1:(end-2));
                                                loc_WM1stderiv_pred = [loc_WM1stderiv_pred(1);loc_WM1stderiv_pred;loc_WM1stderiv_pred(end)];
                                                loc_WM1stderiv_pred = spm_detrend(loc_WM1stderiv_pred,detr_ord);
                                                desmat              = [desmat,loc_WM1stderiv_pred];
                                            end
                                        end
                                        Q = desmat*pinv(full(desmat));
                                    end
                                    out_func_data.img(x,y,z,:) = squeeze(func_data.img(x,y,z,:))-(Q*squeeze(func_data.img(x,y,z,:)));
                                end
                            end
                        end
                    end
                    out_func_data.hdr               = func_data.hdr;
                    out_func_data.hdr.dime.dim(2:5) = size(out_func_data.img);
                    out_func_data.hdr.dime.cal_max  = max(out_func_data.img(:));
                    out_func_data.hdr.dime.cal_min  = min(out_func_data.img(out_func_data.img~=0));
                    out_func_data.hdr.dime.datatype = 64;
                    out_func_data.hdr.dime.bitpix   = 64;
                    if strcmp(saveintermed,'Yes')
                        save_nii_gz(out_func_data,func_filename);
                    end
                else
                    out_func_data = func_data;
                end
                
                if BP==1
                    func_filename = strrep(func_filename,'.nii','_bpfilt.nii');
                    if deconv==1
                        func_filename_deconv = strrep(func_filename,'_bpfilt.nii','_hpfilt.nii');
                    end
                    orig_times       = (1:nTP);
                    uncensored_times = orig_times(temporal_mask);
                    freq_bins        = 1:1:nTP;
                    for bin = length(freq_bins):-1:1
                        theta               = (1/(2*freq_bins(bin)))*atan(sum(sin(2*freq_bins(bin)*uncensored_times))/sum(cos(2*freq_bins(bin)*uncensored_times)));
                        cosmult(:,bin)      = cos(freq_bins(bin).*(uncensored_times-theta))';
                        cosdiv(bin)         = sum(cos(freq_bins(bin).*(uncensored_times-theta)).^2);
                        sinmult(:,bin)      = sin(freq_bins(bin).*(uncensored_times-theta))';
                        sindiv(bin)         = sum(sin(freq_bins(bin).*(uncensored_times-theta)).^2);
                        interpmult(:,:,bin) = [cos(freq_bins(bin).*(orig_times-theta));sin(freq_bins(bin).*(orig_times-theta))];
                    end
                    out_func_data.img = out_func_data.img-repmat(mean(out_func_data.img,4),[1,1,1,size(out_func_data.img,4)]);
                    for x = 1:size(func_data.img,1)
                        for y = 1:size(func_data.img,2)
                            for z = 1:size(func_data.img,3)
                                if funcmask.img(x,y,z)==1
                                    for bin = length(freq_bins):-1:1
                                        cosval             = sum(squeeze(out_func_data.img(x,y,z,:)).*cosmult(:,bin))./cosdiv(bin);
                                        sinval             = sum(squeeze(out_func_data.img(x,y,z,:)).*sinmult(:,bin))./sindiv(bin);
                                        interp_data(:,bin) = [cosval,sinval]*interpmult(:,:,bin);
                                    end
                                    interp_data = sum(interp_data,2);
                                    if deconv==1
                                        interp_data_deconv = interp_data;
                                    end
                                    if any(isnan(interp_data))
                                        out_func_data.img(x,y,z,:) = NaN;
                                        if deconv==1
                                            out_func_data_deconv.img(x,y,z,:) = NaN;
                                        end
                                    else
                                        if exist('HP_sigma','var')
                                            temp_filt                  = fsl_temporal_filt(interp_data,HP_sigma,LP_sigma);
                                            out_func_data.img(x,y,z,:) = temp_filt(temporal_mask);
                                            if deconv==1
                                                out_func_data_deconv.img(x,y,z,:) = fsl_temporal_filt(interp_data_deconv,HP_sigma,-1);
                                            end
                                        elseif exist('butter_filt','var')
                                            interp_data_padded         = [zeros(15,size(interp_data,2));interp_data;zeros(15,size(interp_data,2))];
                                            temp_butter                = filter(butter_filt,interp_data_padded);
                                            temp_butter(end:-1:1)      = filter(butter_filt,temp_butter(end:-1:1));
                                            temp_butter                = temp_butter(16:(end-15));
                                            out_func_data.img(x,y,z,:) = temp_butter(temporal_mask,:);
                                            if deconv==1
                                                interp_data_padded                = [zeros(15,size(interp_data_deconv,2));interp_data_deconv;zeros(15,size(interp_data_deconv,2))];
                                                temp_butter                       = filter(butter_filt_HP,interp_data_padded);
                                                temp_butter(end:-1:1)             = filter(butter_filt_HP,temp_butter(end:-1:1));
                                                out_func_data_deconv.img(x,y,z,:) = temp_butter(16:(end-15));
                                            end
                                        elseif exist('filt_ind','var')
                                            padlength                  = rest_nextpow2_one35(nTP);
                                            interp_data                = [interp_data;zeros(padlength-nTP,size(interp_data,2))];
                                            freq                       = fft(interp_data);
                                            freq(filt_ind,:)           = 0;
                                            interp_data                = ifft(freq);
                                            interp_data                = interp_data(1:nTP,:);
                                            out_func_data.img(x,y,z,:) = interp_data(temporal_mask,:);
                                            if deconv==1
                                                interp_data_deconv                = [interp_data_deconv;zeros(padlength-nTP,size(interp_data_deconv,2))];
                                                freq                              = fft(interp_data_deconv);
                                                freq(filt_ind_HP,:)               = 0;
                                                interp_data_deconv                = ifft(freq);
                                                out_func_data_deconv.img(x,y,z,:) = interp_data_deconv(1:nTP,:);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    out_func_data.hdr               = func_data.hdr;
                    out_func_data.hdr.dime.dim(2:5) = size(out_func_data.img);
                    out_func_data.hdr.dime.cal_max  = max(out_func_data.img(:));
                    out_func_data.hdr.dime.cal_min  = min(out_func_data.img(out_func_data.img~=0));
                    out_func_data.hdr.dime.datatype = 64;
                    out_func_data.hdr.dime.bitpix   = 64;
                    if strcmp(saveintermed,'Yes')
                        save_nii_gz(out_func_data,func_filename);
                    end
                elseif deconv==1
                    func_filename_deconv = func_filename;
                    orig_times           = (1:nTP);
                    uncensored_times     = orig_times(temporal_mask);
                    freq_bins            = 1:1:nTP;
                    for bin = length(freq_bins):-1:1
                        theta               = (1/(2*freq_bins(bin)))*atan(sum(sin(2*freq_bins(bin)*uncensored_times))/sum(cos(2*freq_bins(bin)*uncensored_times)));
                        cosmult(:,bin)      = cos(freq_bins(bin).*(uncensored_times-theta))';
                        cosdiv(bin)         = sum(cos(freq_bins(bin).*(uncensored_times-theta)).^2);
                        sinmult(:,bin)      = sin(freq_bins(bin).*(uncensored_times-theta))';
                        sindiv(bin)         = sum(sin(freq_bins(bin).*(uncensored_times-theta)).^2);
                        interpmult(:,:,bin) = [cos(freq_bins(bin).*(orig_times-theta));sin(freq_bins(bin).*(orig_times-theta))];
                    end
                    out_func_data_for_deconv.img = out_func_data.img;
                    out_func_data_for_deconv.img = out_func_data_for_deconv.img-repmat(mean(out_func_data_for_deconv.img,4),[1,1,1,size(out_func_data_for_deconv.img,4)]);
                    for x = 1:size(func_data.img,1)
                        for y = 1:size(func_data.img,2)
                            for z = 1:size(func_data.img,3)
                                if funcmask.img(x,y,z)==1
                                    for bin = length(freq_bins):-1:1
                                        cosval             = sum(squeeze(out_func_data_for_deconv.img(x,y,z,:)).*cosmult(:,bin))./cosdiv(bin);
                                        sinval             = sum(squeeze(out_func_data_for_deconv.img(x,y,z,:)).*sinmult(:,bin))./sindiv(bin);
                                        interp_data(:,bin) = [cosval,sinval]*interpmult(:,:,bin);
                                    end
                                    out_func_data_deconv.img(x,y,z,:) = sum(interp_data,2);
                                end
                            end
                        end
                    end
                    clear out_func_data_for_deconv.img
                end
                data_accept = 1;
            else
                data_accept = 0;
            end
        else
            data_accept = 1;
        end
    else
        mean_FD                    = NaN;
        mean_DVARS                 = NaN;
        num_censored_vols          = 0;
        num_abovethresh_vols_FD    = 0;
        num_abovethresh_vols_DVARS = 0;
        data_accept                = 1;
    end
    
    % Save data at current step
    if strcmp(saveintermed,'Yes')
        out_func_data.hdr               = func_data.hdr;
        out_func_data.hdr.dime.dim(2:5) = size(out_func_data.img);
        out_func_data.hdr.dime.cal_max  = max(out_func_data.img(:));
        out_func_data.hdr.dime.cal_min  = min(out_func_data.img(out_func_data.img~=0));
        out_func_data.hdr.dime.datatype = 64;
        out_func_data.hdr.dime.bitpix   = 64;
        save_nii_gz(out_func_data,func_filename);
    end
    
    if deconv==1 && strcmp(whendeconv,'before')
        out_func_data_deconv_orig.img = out_func_data_deconv.img;
        out_func_data_deconv.img      = zeros(size(out_func_data_deconv.img));
        if strcmp(deconv_method,'SPM')
            curr_nTP = length(squeeze(out_func_data_deconv_orig.img(1,1,1,:)));
            dt  = TR/16;
            NT  = TR/dt;
            hrf = spm_hrf(dt);
            for x = 1:size(func_data.img,1)
                for y = 1:size(func_data.img,2)
                    for z = 1:size(func_data.img,3)
                        if funcmask.img(x,y,z)==1
                            curr_ts = squeeze(out_func_data_deconv_orig.img(x,y,z,:))-mean(squeeze(out_func_data_deconv_orig.img(x,y,z,:)));
                            if any(isnan(curr_ts))
                                out_func_data_deconv.img(x,y,z,:) = NaN;
                            else
                                P       = {};
                                xb      = spm_dctmtx(curr_nTP*NT+128,curr_nTP);
                                P{1}.X  = zeros(curr_nTP,curr_nTP);
                                for i = 1:curr_nTP
                                    Hx          = conv(xb(:,i),hrf);
                                    P{1}.X(:,i) = Hx((1:NT:curr_nTP*NT)+128);
                                end
                                xb       = xb(129:end,:);
                                P{1}.C   = speye(curr_nTP,curr_nTP)/4;
                                P{2}.X   = sparse(curr_nTP,1);
                                P{2}.C   = speye(curr_nTP,curr_nTP)*curr_nTP/trace(P{1}.X'*P{1}.X);
                                C        = spm_PEB(curr_ts,P);
                                ts_decon = xb*C{2}.E(1:curr_nTP);
                                
                                % Resample back to original scale
                                out_func_data_deconv.img(x,y,z,:) = resample(ts_decon,length(curr_ts),length(ts_decon));
                            end
                        end
                    end
                end
            end
        elseif strcmp(deconv_method,'spontaneous pseudo events')
            event_lag_max = round(10/TR);
            thr           = 1;
            for x = 1:size(func_data.img,1)
                for y = 1:size(func_data.img,2)
                    for z = 1:size(func_data.img,3)
                        if funcmask.img(x,y,z)==1
                            curr_ts = squeeze(out_func_data_deconv_orig.img(x,y,z,:));
                            if any(isnan(curr_ts))
                                out_func_data_deconv.img(x,y,z,:) = NaN;
                            else
                                out_func_data_deconv.img(x,y,z,:) = wgr_deconv_canonhrf_par_edit(zscore(curr_ts(:)),thr,event_lag_max,TR);
                            end
                        end
                    end
                end
            end
        elseif strcmp(deconv_method,'non-linear regression')
            nev_lr  = 0.01;
            epsilon = 0.005;
            hrf     = spm_hrf(TR);
            for x = 1:size(func_data.img,1)
                for y = 1:size(func_data.img,2)
                    for z = 1:size(func_data.img,3)
                        if funcmask.img(x,y,z)==1
                            curr_ts = squeeze(out_func_data_deconv_orig.img(x,y,z,:));
                            if any(isnan(curr_ts))
                                out_func_data_deconv.img(x,y,z,:) = NaN;
                            else
                                % Scale data:
                                curr_ts = rescale_ts(curr_ts,1,0);
                                
                                % Deconvolve:
                                temp = deconvolve_Bush_2011(curr_ts',hrf,nev_lr,epsilon);
                                
                                % Remove initial dummy timepoints:
                                out_func_data_deconv.img(x,y,z,:) = temp(length(hrf):end);
                            end
                        end
                    end
                end
            end
        end
        
        if BP==1
            func_filename_deconv = strrep(func_filename_deconv,'.nii','_deconv_bpfilt.nii');
            for x = 1:size(func_data.img,1)
                for y = 1:size(func_data.img,2)
                    for z = 1:size(func_data.img,3)
                        if funcmask.img(x,y,z)==1
                            tempdata = squeeze(out_func_data_deconv.img(x,y,z,:));
                            if any(isnan(tempdata))
                                out_func_data_deconv.img(x,y,z,:) = NaN;
                            else
                                tempdata = tempdata-mean(tempdata);
                                if exist('HP_sigma','var')
                                    out_func_data_deconv.img(x,y,z,:) = fsl_temporal_filt(tempdata,HP_sigma,LP_sigma);
                                elseif exist('butter_filt','var')
                                    tempdata_padded                   = [zeros(15,size(tempdata,2));tempdata;zeros(15,size(tempdata,2))];
                                    temp_butter                       = filter(butter_filt,tempdata_padded);
                                    temp_butter(end:-1:1)             = filter(butter_filt,temp_butter(end:-1:1));
                                    out_func_data_deconv.img(x,y,z,:) = temp_butter(16:(end-15));
                                elseif exist('filt_ind','var')
                                    padlength                         = rest_nextpow2_one35(nTP);
                                    tempdata                          = [tempdata;zeros(padlength-nTP,size(tempdata,2))];
                                    freq                              = fft(tempdata);
                                    freq(filt_ind,:)                  = 0;
                                    tempdata                          = ifft(freq);
                                    out_func_data_deconv.img(x,y,z,:) = tempdata(1:nTP,:);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Remove autocorrelation
    if rem_AR==1
        func_filename = strrep(func_filename,'.nii',['_remAR',num2str(num_hist_tp),'.nii']);
        
        if ~strcmp(useavail,'Yes, do not run step if output file exists') || ~exist(func_filename,'file')
            out_func_data_orig.img = out_func_data.img;
            [xdim,ydim,zdim,tdim]  = size(out_func_data.img);
            out_func_data.img      = zeros(xdim,ydim,zdim,(tdim-num_hist_tp));
            for x = 1:xdim
                for y = 1:ydim
                    for z = 1:zdim
                        if funcmask.img(x,y,z)==1
                            ar_desmat = [];
                            for htp = 1:num_hist_tp
                                ar_desmat = [ar_desmat,squeeze(out_func_data_orig.img(x,y,z,htp:(end-num_hist_tp+htp-1)))];
                            end
                            stats                      = regstats(squeeze(out_func_data_orig.img(x,y,z,(num_hist_tp+1):end)),ar_desmat,'linear','r');
                            out_func_data.img(x,y,z,:) = stats.r;
                        end
                    end
                end
            end
            clear out_func_data_orig.img
            if strcmp(saveintermed,'Yes')
                out_func_data.hdr               = func_data.hdr;
                out_func_data.hdr.dime.dim(2:5) = size(out_func_data.img);
                out_func_data.hdr.dime.cal_max  = max(out_func_data.img(:));
                out_func_data.hdr.dime.cal_min  = min(out_func_data.img(out_func_data.img~=0));
                out_func_data.hdr.dime.datatype = 64;
                out_func_data.hdr.dime.bitpix   = 64;
                save_nii_gz(out_func_data,func_filename);
            end
        end
        if deconv==1 && strcmp(whendeconv,'before')
            func_filename_deconv          = strrep(func_filename_deconv,'.nii',['_remAR',num2str(num_hist_tp),'.nii']);
            out_func_data_deconv_orig.img = out_func_data_deconv.img;
            [xdim,ydim,zdim,tdim]         = size(out_func_data_deconv.img);
            out_func_data_deconv.img      = zeros(xdim,ydim,zdim,(tdim-num_hist_tp));
            for x = 1:xdim
                for y = 1:ydim
                    for z = 1:zdim
                        if funcmask.img(x,y,z)==1
                            ar_desmat = [];
                            for htp = 1:num_hist_tp
                                ar_desmat = [ar_desmat,squeeze(out_func_data_deconv_orig.img(x,y,z,htp:(end-num_hist_tp+htp-1)))];
                            end
                            stats                             = regstats(squeeze(out_func_data_deconv_orig.img(x,y,z,(num_hist_tp+1):end)),ar_desmat,'linear','r');
                            out_func_data_deconv.img(x,y,z,:) = stats.r;
                        end
                    end
                end
            end
            clear out_func_data_deconv_orig.img
            if strcmp(saveintermed,'Yes')
                out_func_data_deconv.hdr               = func_data.hdr;
                out_func_data_deconv.hdr.dime.dim(2:5) = size(out_func_data_deconv.img);
                out_func_data_deconv.hdr.dime.cal_max  = max(out_func_data_deconv.img(:));
                out_func_data_deconv.hdr.dime.cal_min  = min(out_func_data_deconv.img(out_func_data_deconv.img~=0));
                out_func_data_deconv.hdr.dime.datatype = 64;
                out_func_data_deconv.hdr.dime.bitpix   = 64;
                save_nii_gz(out_func_data_deconv,func_filename_deconv);
            end
        end
    end
    
    % Extract timeseries for each ROI
    GCOR = NaN;
    if data_accept==1
        if globsig_calcGCOR==1
            GCOR_data = reshape(out_func_data.img,[numel(funcmask.img),size(out_func_data.img,4)]);
            GCOR_data = GCOR_data(logical(reshape(funcmask.img,[numel(funcmask.img),1])),:);
            GCOR_data = GCOR_data-(sum(GCOR_data,2)./size(GCOR_data,2))*ones(1,size(GCOR_data,2));
            GCOR      = norm(sum(GCOR_data./(sqrt(sum(GCOR_data.^2,2))*ones(1,size(GCOR_data,2))))./size(GCOR_data,1)).^2;
        end
        
        parc_data = load_nii_gz(parc_filename);
        nredu_TP  = size(out_func_data.img,4);
        ts        = zeros(nROI,nredu_TP);
        if deconv==1
            nredu_deconv_TP = size(out_func_data_deconv.img,4);
            deconv_ts       = zeros(nROI,nredu_deconv_TP);
        end
        
        for ROI = 1:nROI
            curr_label = ROI_num_labels(ROI);
            if iscell(curr_label)
                curr_label = curr_label{:};
            end
            ROImask = zeros(size(parc_data.img));
            if ~isempty(find(parc_data.img==curr_label,1))
                ROImask(logical(parc_data.img==curr_label)) = 1;
                nvox_orig                                   = sum(ROImask(:));
                ROImask                                     = ROImask.*funcmask.img;
                nvox_masked                                 = sum(ROImask(:));
                
                if nvox_masked>=MinROIvox
                    if (nvox_masked/nvox_orig)>=(MinpercentROIvox/100)
                        if strcmp(ROI_summary_measure,'Mean')
                            for t = 1:nredu_TP
                                ts(ROI,t) = sum(sum(sum(out_func_data.img(:,:,:,t).*ROImask)))./nvox_masked;
                            end
                            if deconv==1
                                for t = 1:nredu_deconv_TP
                                    deconv_ts(ROI,t) = sum(sum(sum(out_func_data_deconv.img(:,:,:,t).*ROImask)))./nvox_masked;
                                end
                            end
                        elseif strcmp(ROI_summary_measure,'Median')
                            for t = 1:nredu_TP
                                temp      = out_func_data.img(:,:,:,t);
                                ts(ROI,t) = median(temp(logical(ROImask)));
                            end
                            if deconv==1
                                for t = 1:nredu_deconv_TP
                                    temp             = out_func_data_deconv.img(:,:,:,t);
                                    deconv_ts(ROI,t) = median(temp(logical(ROImask)));
                                end
                            end
                        elseif strcmp(ROI_summary_measure,'Largest Principal Component')
                            ROImask_orig = ROImask;
                            for t = 1:nredu_TP
                                tempmeansig(t,1) = sum(sum(sum(out_func_data.img(:,:,:,t).*ROImask)))./nvox_masked;
                            end
                            ROImask          = logical(reshape(ROImask,[(size(ROImask,1)*size(ROImask,2)*size(ROImask,3)),1]));
                            temp             = reshape(out_func_data.img,[(size(ROImask,1)*size(ROImask,2)*size(ROImask,3)),nredu_TP]);
                            temp(~ROImask,:) = [];
                            try
                                [~,PCA_scores] = pca(temp','Algorithm','svd');
                            catch
                                [~,PCA_scores] = princomp(temp');
                            end
                            if sign(corr(tempmeansig,PCA_scores(:,1)))==(-1)
                                PCA_scores(:,1) = PCA_scores(:,1)*(-1);
                            end
                            ts(ROI,:) = PCA_scores(:,1);
                            clear tempmeansig
                            if deconv==1
                                ROImask = ROImask_orig;
                                for t = 1:nredu_deconv_TP
                                    tempmeansig(t,1) = sum(sum(sum(out_func_data_deconv.img(:,:,:,t).*ROImask)))./nvox_masked;
                                end
                                ROImask          = logical(reshape(ROImask,[(size(ROImask,1)*size(ROImask,2)*size(ROImask,3)),1]));
                                temp             = reshape(out_func_data_deconv.img,[(size(ROImask,1)*size(ROImask,2)*size(ROImask,3)),nredu_deconv_TP]);
                                temp(~ROImask,:) = [];
                                try
                                    [~,PCA_scores] = pca(temp','Algorithm','svd');
                                catch
                                    [~,PCA_scores] = princomp(temp');
                                end
                                if sign(corr(tempmeansig,PCA_scores(:,1)))==(-1)
                                    PCA_scores(:,1) = PCA_scores(:,1)*(-1);
                                end
                                deconv_ts(ROI,:) = PCA_scores(:,1);
                            end
                            clear tempmeansig
                        end
                        if sum(isnan(ts(ROI,:)))>0
                            fprintf(logfile_fid,'There is at least one NaN in the timeseries of ROI # %u for %s\n',curr_label,sub);
                        end
                        if deconv==1
                            if sum(isnan(deconv_ts(ROI,:)))>0
                                fprintf(logfile_fid,'There is at least one NaN in the deconvolved timeseries of ROI # %u for %s\n',curr_label,sub);
                            end
                        end
                    else
                        fprintf(logfile_fid,'Less than %u%% of the original voxels are in the masked version of ROI # %u for %s\n',MinpercentROIvox,curr_label,sub);
                        ts(ROI,:) = NaN;
                        if deconv==1
                            deconv_ts(ROI,:) = NaN;
                        end
                    end
                else
                    fprintf(logfile_fid,'Less than 5 voxels in masked ROI # %u for %s\n',curr_label,sub);
                    ts(ROI,:) = NaN;
                    if deconv==1
                        deconv_ts(ROI,:) = NaN;
                    end
                end
            else
                fprintf(logfile_fid,'No mask for ROI # %u for %s\n',curr_label,sub);
                ts(ROI,:) = NaN;
                if deconv==1
                    deconv_ts(ROI,:) = NaN;
                end
            end
        end
        
        if deconv==1 && strcmp(whendeconv,'after')
            deconv_ts_orig = deconv_ts;
            deconv_ts      = zeros(size(deconv_ts));
            if strcmp(deconv_method,'SPM')
                curr_nTP = size(deconv_ts,2);
                dt       = TR/16;
                NT       = TR/dt;
                hrf      = spm_hrf(dt);
                for ROI = 1:size(deconv_ts,1)
                    curr_ts = squeeze(deconv_ts_orig(ROI,:))-mean(squeeze(deconv_ts_orig(ROI,:)));
                    if any(isnan(curr_ts))
                        deconv_ts(ROI,:) = NaN;
                    else
                        P      = {};
                        xb     = spm_dctmtx(curr_nTP*NT+128,curr_nTP);
                        P{1}.X = zeros(curr_nTP,curr_nTP);
                        for i = 1:curr_nTP
                            Hx          = conv(xb(:,i),hrf);
                            P{1}.X(:,i) = Hx((1:NT:curr_nTP*NT)+128);
                        end
                        xb       = xb(129:end,:);
                        P{1}.C   = speye(curr_nTP,curr_nTP)/4;
                        P{2}.X   = sparse(curr_nTP,1);
                        P{2}.C   = speye(curr_nTP,curr_nTP)*curr_nTP/trace(P{1}.X'*P{1}.X);
                        C        = spm_PEB(curr_ts(:),P);
                        ts_decon = xb*C{2}.E(1:curr_nTP);
                        
                        % Resample back to original scale
                        deconv_ts(ROI,:) = resample(ts_decon,length(curr_ts),length(ts_decon));
                    end
                end
            elseif strcmp(deconv_method,'spontaneous pseudo events')
                event_lag_max = round(10/TR);
                thr           = 1;
                for ROI = 1:size(deconv_ts,1)
                    curr_ts = squeeze(deconv_ts_orig(ROI,:));
                    if any(isnan(curr_ts))
                        deconv_ts(ROI,:) = NaN;
                    else
                        deconv_ts(ROI,:) = wgr_deconv_canonhrf_par_edit(zscore(curr_ts(:)),thr,event_lag_max,TR);
                    end
                end
            elseif strcmp(deconv_method,'non-linear regression')
                nev_lr  = 0.01;
                epsilon = 0.005;
                hrf     = spm_hrf(TR);
                for ROI = 1:size(deconv_ts,1)
                    curr_ts = squeeze(deconv_ts_orig(ROI,:));
                    if any(isnan(curr_ts))
                        deconv_ts(ROI,:) = NaN;
                    else
                        % Scale data:
                        curr_ts = rescale_ts(curr_ts,1,0);
                        
                        % Deconvolve:
                        temp = deconvolve_Bush_2011(curr_ts,hrf,nev_lr,epsilon);
                        
                        % Remove initial dummy timepoints:
                        deconv_ts(ROI,:) = temp(length(hrf):end);
                    end
                end
            end
            
            if BP==1
                curr_nTP = size(deconv_ts,2);
                for ROI = 1:size(deconv_ts,1)
                    tempdata = deconv_ts(ROI,:);
                    if any(isnan(tempdata))
                        deconv_ts(ROI,:) = NaN;
                    else
                        tempdata = tempdata-mean(tempdata);
                        if exist('HP_sigma','var')
                            deconv_ts(ROI,:) = fsl_temporal_filt(tempdata,HP_sigma,LP_sigma);
                        elseif exist('butter_filt','var')
                            tempdata_padded       = [zeros(15,size(tempdata(:),2));tempdata(:);zeros(15,size(tempdata(:),2))];
                            temp_butter           = filter(butter_filt,tempdata_padded);
                            temp_butter(end:-1:1) = filter(butter_filt,temp_butter(end:-1:1));
                            deconv_ts(ROI,:)      = temp_butter(16:(end-15));
                        elseif exist('filt_ind','var')
                            padlength        = rest_nextpow2_one35(curr_nTP);
                            tempdata         = [tempdata(:);zeros(padlength-curr_nTP,size(tempdata(:),2))];
                            freq             = fft(tempdata);
                            freq(filt_ind,:) = 0;
                            tempdata         = ifft(freq);
                            deconv_ts(ROI,:) = tempdata(1:curr_nTP,:);
                        end
                    end
                end
            end
        end
    else
        ts(:,:) = NaN;
        if BP==1
            deconv_ts(:,:) = NaN;
        end
    end
    if deconv==0
        deconv_ts = NaN(size(ts));
    end
elseif strcmp(ROIextord,'After only slice timing/motion correction/detrending')
    % Load functional data/mask & ROI atlas
    parc_data     = load_nii_gz(parc_filename);
    func_data     = load_nii_gz(func_filename);
    funcmask      = load_nii_gz(funcmask_filename);
    func_data.img = double(func_data.img);
    funcmask.img  = double(funcmask.img);
    
    % Preallocate output data
    ts = zeros(nROI,nTP);
    
    % Extract timeseries for each ROI
    for ROI = 1:nROI
        curr_label = ROI_num_labels(ROI);
        if iscell(curr_label)
            curr_label = curr_label{:};
        end
        ROImask = zeros(size(parc_data.img));
        if ~isempty(find(parc_data.img==curr_label,1))
            ROImask(logical(parc_data.img==curr_label)) = 1;
            nvox_orig                                   = sum(ROImask(:));
            ROImask                                     = ROImask.*funcmask.img;
            nvox_masked                                 = sum(ROImask(:));
            
            if nvox_masked>=MinROIvox
                if (nvox_masked/nvox_orig)>=(MinpercentROIvox/100)
                    if strcmp(ROI_summary_measure,'Mean')
                        for t = 1:nTP
                            ts(ROI,t) = sum(sum(sum(func_data.img(:,:,:,t).*ROImask)))./nvox_masked;
                        end
                    elseif strcmp(ROI_summary_measure,'Median')
                        for t = 1:nTP
                            temp      = func_data.img(:,:,:,t);
                            ts(ROI,t) = median(temp(logical(ROImask)));
                        end
                    elseif strcmp(ROI_summary_measure,'Largest Principal Component')
                        for t = 1:nTP
                            tempmeansig(t,1) = sum(sum(sum(func_data.img(:,:,:,t).*ROImask)))./nvox_masked;
                        end
                        ROImask          = logical(reshape(ROImask,[(size(ROImask,1)*size(ROImask,2)*size(ROImask,3)),1]));
                        temp             = reshape(func_data.img,[(size(ROImask,1)*size(ROImask,2)*size(ROImask,3)),nTP]);
                        temp(~ROImask,:) = [];
                        try
                            [~,PCA_scores] = pca(temp','Algorithm','svd');
                        catch
                            [~,PCA_scores] = princomp(temp');
                        end
                        if sign(corr(tempmeansig,PCA_scores(:,1)))==(-1)
                            PCA_scores(:,1) = PCA_scores(:,1)*(-1);
                        end
                        ts(ROI,:) = PCA_scores(:,1);
                    end
                    if sum(isnan(ts(ROI,:)))>0
                        fprintf(logfile_fid,'There is at least one NaN in the timeseries of ROI # %u for %s\n',curr_label,sub);
                    end
                else
                    fprintf(logfile_fid,'Less than %u%% of the original voxels are in the masked version of ROI # %u for %s\n',MinpercentROIvox,curr_label,sub);
                    ts(ROI,:) = NaN;
                end
            else
                fprintf(logfile_fid,'Less than 5 voxels in masked ROI # %u for %s\n',curr_label,sub);
                ts(ROI,:) = NaN;
            end
        else
            fprintf(logfile_fid,'No mask for ROI # %u for %s\n',curr_label,sub);
            ts(ROI,:) = NaN;
        end
    end
    
    % Create design matrix of nuisance signals
    desmat = [];
    
    % Extract signal from white matter
    if WMsig_part==1
        if erode_masks==1
            WM_eroded_filename = strrep(WM_filename,'.nii',['_eroded_',num2str(erode_radius),'.nii']);
            if strcmp(useavail,'Yes, do not run step if output file exists')
                if ~exist(WM_eroded_filename,'file')
                    tempmask       = load_nii_gz(WM_filename);
                    [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                    save_nii_gz(tempmask,WM_eroded_filename);
                end
            else
                tempmask       = load_nii_gz(WM_filename);
                [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                save_nii_gz(tempmask,WM_eroded_filename);
            end
            WM_filename = WM_eroded_filename;
        end
        if PCA_WMvent==1
            WM_extract_filename = [strrep(strrep(WM_filename,'.gz',''),'.nii',''),'_EIG.txt'];
            if ~exist(WM_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename,' --eig --order=5']);
            else
                vernum = 0;
                while 1
                    if exist(strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        WM_extract_filename = strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename,' --eig --order=5']);
            end
            WM_preds = dlmread(WM_extract_filename);
            WM_preds = spm_detrend(WM_preds,detr_ord);
            desmat   = [desmat,WM_preds];
            
            if WMsigderiv_part==1
                WM1stderiv_preds = WM_preds(3:end,:)-WM_preds(1:(end-2),:);
                WM1stderiv_preds = [WM1stderiv_preds(1,:);WM1stderiv_preds;WM1stderiv_preds(end,:)];
                WM1stderiv_preds = spm_detrend(WM1stderiv_preds,detr_ord);
                desmat           = [desmat,WM1stderiv_preds];
            end
        else
            WM_extract_filename = [strrep(strrep(WM_filename,'.gz',''),'.nii',''),'.txt'];
            if ~exist(WM_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename]);
            else
                vernum = 0;
                while 1
                    if exist(strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        WM_extract_filename = strrep(WM_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',WM_extract_filename,' -m ',WM_filename]);
            end
            WM_pred = dlmread(WM_extract_filename);
            WM_pred = spm_detrend(WM_pred,detr_ord);
            desmat  = [desmat,WM_pred];
            
            if WMsigderiv_part==1
                WM1stderiv_pred = WM_pred(3:end)-WM_pred(1:(end-2));
                WM1stderiv_pred = [WM1stderiv_pred(1);WM1stderiv_pred;WM1stderiv_pred(end)];
                WM1stderiv_pred = spm_detrend(WM1stderiv_pred,detr_ord);
                desmat          = [desmat,WM1stderiv_pred];
            end
        end
    end
    
    % Extract signal from ventricles
    if ventsig_part==1
        if erode_masks==1
            vent_eroded_filename = strrep(vent_filename,'.nii',['_eroded_',num2str(erode_radius),'.nii']);
            if strcmp(useavail,'Yes, do not run step if output file exists')
                if ~exist(vent_eroded_filename,'file')
                    tempmask       = load_nii_gz(vent_filename);
                    [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                    save_nii_gz(tempmask,vent_eroded_filename);
                end
            else
                tempmask       = load_nii_gz(vent_filename);
                [tempmask.img] = erode_mask(tempmask.img,erode_radius);
                save_nii_gz(tempmask,vent_eroded_filename);
            end
            vent_filename = vent_eroded_filename;
        end
        if PCA_WMvent==1
            vent_extract_filename = [strrep(strrep(vent_filename,'.gz',''),'.nii',''),'_EIG.txt'];
            if ~exist(vent_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename,' --eig --order=5']);
            else
                vernum = 0;
                while 1
                    if exist(strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        vent_extract_filename = strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename,' --eig --order=5']);
            end
            vent_preds = dlmread(vent_extract_filename);
            vent_preds = spm_detrend(vent_preds,detr_ord);
            desmat     = [desmat,vent_preds];
            
            if ventsigderiv_part==1
                vent1stderiv_preds = vent_preds(3:end,:)-vent_preds(1:(end-2),:);
                vent1stderiv_preds = [vent1stderiv_preds(1,:);vent1stderiv_preds;vent1stderiv_preds(end,:)];
                vent1stderiv_preds = spm_detrend(vent1stderiv_preds,detr_ord);
                desmat             = [desmat,vent1stderiv_preds];
            end
        else
            vent_extract_filename = [strrep(strrep(vent_filename,'.gz',''),'.nii',''),'.txt'];
            if ~exist(vent_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename]);
            else
                vernum = 0;
                while 1
                    if exist(strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        vent_extract_filename = strrep(vent_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',vent_extract_filename,' -m ',vent_filename]);
            end
            vent_pred = dlmread(vent_extract_filename);
            vent_pred = spm_detrend(vent_pred,detr_ord);
            desmat    = [desmat,vent_pred];
            
            if ventsigderiv_part==1
                vent1stderiv_pred = vent_pred(3:end)-vent_pred(1:(end-2));
                vent1stderiv_pred = [vent1stderiv_pred(1);vent1stderiv_pred;vent1stderiv_pred(end)];
                vent1stderiv_pred = spm_detrend(vent1stderiv_pred,detr_ord);
                desmat            = [desmat,vent1stderiv_pred];
            end
        end
    end
    
    % Extract global signal & calculate GNI
    if globsig_part==1
        if globsig_calcGNI==1
            GNI = Rest2GNI_edit(func_filename,funcmask_filename,100);
            fprintf(logfile_fid,'Original GNI for %s = %6.3f\n',sub,GNI);
        else
            GNI = 0;
        end
        if GNI<=3
            funcmask_extract_filename = [strrep(strrep(funcmask_filename,'.gz',''),'.nii',''),'.txt'];
            if ~exist(funcmask_extract_filename,'file')
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',funcmask_extract_filename,' -m ',funcmask_filename]);
            else
                vernum = 0;
                while 1
                    if exist(strrep(funcmask_extract_filename,'.txt',['_',num2str(vernum),'.txt']),'file')
                        vernum = vernum+1;
                    else
                        funcmask_extract_filename = strrep(funcmask_extract_filename,'.txt',['_',num2str(vernum),'.txt']);
                        break
                    end
                end
                call_fsl_edit(['!fslmeants -i ',func_filename,' -o ',funcmask_extract_filename,' -m ',funcmask_filename]);
            end
            global_pred = dlmread(funcmask_extract_filename);
            global_pred = spm_detrend(global_pred,detr_ord);
            desmat      = [desmat,global_pred];
            if globsigderiv_part==1
                global1stderiv_pred = global_pred(3:end)-global_pred(1:(end-2));
                global1stderiv_pred = [global1stderiv_pred(1);global1stderiv_pred;global1stderiv_pred(end)];
                global1stderiv_pred = spm_detrend(global1stderiv_pred,detr_ord);
                desmat              = [desmat,global1stderiv_pred];
            end
        end
    end
    
    % Extract motion parameters
    if motpar_part==1
        mot_pred = dlmread(par_filename);
        for par = 1:size(mot_pred,2)
            mot_pred(:,par) = spm_detrend(mot_pred(:,par),detr_ord);
        end
        desmat = [desmat,mot_pred];
        if motpart1_part==1
            mot1back_pred = [mot_pred(2:end,:);mot_pred(end,:)];
            for par = 1:size(mot_pred,2)
                mot1back_pred(:,par) = spm_detrend(mot1back_pred(:,par),detr_ord);
            end
            desmat = [desmat,mot1back_pred];
        end
        if motparderiv_part==1
            mot1stderiv_pred = mot_pred(3:end,:)-mot_pred(1:(end-2),:);
            mot1stderiv_pred = [mot1stderiv_pred(1,:);mot1stderiv_pred;mot1stderiv_pred(end,:)];
            for par = 1:size(mot_pred,2)
                mot1stderiv_pred(:,par) = spm_detrend(mot1stderiv_pred(:,par),detr_ord);
            end
            desmat = [desmat,mot1stderiv_pred];
        end
        if motparsqr_part==1
            motsquare_pred = mot_pred.^2;
            for par = 1:size(mot_pred,2)
                motsquare_pred(:,par) = spm_detrend(motsquare_pred(:,par),detr_ord);
            end
            desmat = [desmat,motsquare_pred];
            if motpart1_part==1
                mot1backsquare_pred = mot1back_pred.^2;
                for par = 1:size(mot_pred,2)
                    mot1backsquare_pred(:,par) = spm_detrend(mot1backsquare_pred(:,par),detr_ord);
                end
                desmat = [desmat,mot1backsquare_pred];
            end
        end
    end
    
    % Partial nuisance variance then apply temporal filter
    if BP==1 && (WMsig_part==1 || ventsig_part==1 || (globsig_part==1 && GNI<=3) || motpar_part==1)
        desmat = [ones(size(desmat,1),1),desmat];
        if rank(desmat)~=min(size(desmat))
            error(['Initial desmat to partial out unwanted variance for ',sub,' is not full rank']);
        end
        Q = desmat*pinv(full(desmat));
        for ROI = 1:nROI
            if any(isnan(ts(ROI,:)))
                ts(ROI,:) = NaN;
            else
                tempdata = squeeze(ts(ROI,:))-(Q*squeeze(ts(ROI,:)));
                tempdata = tempdata-mean(tempdata);
                if exist('HP_sigma','var')
                    ts(ROI,:) = fsl_temporal_filt(tempdata,HP_sigma,LP_sigma);
                elseif exist('butter_filt','var')
                    tempdata_padded       = [zeros(15,size(tempdata,2));tempdata;zeros(15,size(tempdata,2))];
                    temp_butter           = filter(butter_filt,tempdata_padded);
                    temp_butter(end:-1:1) = filter(butter_filt,temp_butter(end:-1:1));
                    ts(ROI,:)             = temp_butter(16:(end-15));
                elseif exist('filt_ind','var')
                    padlength        = rest_nextpow2_one35(nTP);
                    tempdata         = [tempdata;zeros(padlength-nTP,size(tempdata,2))];
                    freq             = fft(tempdata);
                    freq(filt_ind,:) = 0;
                    tempdata         = ifft(freq);
                    ts(ROI,:)        = tempdata(1:nTP,:);
                end
            end
        end
    elseif BP==1
        for ROI = 1:nROI
            if any(isnan(ts(ROI,:)))
                ts(ROI,:) = NaN;
            else
                tempdata = squeeze(ts(ROI,:));
                tempdata = tempdata-mean(tempdata);
                if exist('HP_sigma','var')
                    ts(ROI,:) = fsl_temporal_filt(tempdata,HP_sigma,LP_sigma);
                elseif exist('butter_filt','var')
                    tempdata_padded       = [zeros(15,size(tempdata,2));tempdata;zeros(15,size(tempdata,2))];
                    temp_butter           = filter(butter_filt,tempdata_padded);
                    temp_butter(end:-1:1) = filter(butter_filt,temp_butter(end:-1:1));
                    ts(ROI,:)             = temp_butter(16:(end-15));
                else
                    padlength        = rest_nextpow2_one35(nTP);
                    tempdata         = [tempdata;zeros(padlength-nTP,size(tempdata,2))];
                    freq             = fft(tempdata);
                    freq(filt_ind,:) = 0;
                    tempdata         = ifft(freq);
                    ts(ROI,:)        = tempdata(1:nTP,:);
                end
            end
        end
    elseif WMsig_part==1 || ventsig_part==1 || (globsig_part==1 && GNI<=3) || motpar_part==1
        desmat = [ones(size(desmat,1),1),desmat];
        if rank(desmat)~=min(size(desmat))
            error(['Initial desmat to partial out unwanted variance for ' sub ' is not full rank']);
        end
        Q = desmat*pinv(full(desmat));
        for ROI = 1:nROI
            if sum(isnan(ts(ROI,:)))==0
                ts(ROI,:) = squeeze(ts(ROI,:))-(Q*squeeze(ts(ROI,:)));
            end
        end
    end
    
    % Remove autocorrelation
    if rem_AR==1
        ts_orig = ts;
        tdim    = size(ts,2);
        ts      = zeros(nROI,(tdim-num_hist_tp));
        for ROI = 1:nROI
            if sum(isnan(ts_orig(ROI,:)))==0
                ar_desmat = [];
                for htp = 1:num_hist_tp
                    ar_desmat = [ar_desmat,squeeze(ts_orig(ROI,htp:(end-num_hist_tp+htp-1)))];
                end
                stats     = regstats(squeeze(ts_orig(ROI,(num_hist_tp+1):end)),ar_desmat,'linear','r');
                ts(ROI,:) = stats.r;
            else
                ts(ROI,:) = NaN;
            end
        end
    end
    
    mean_FD                    = NaN;
    mean_DVARS                 = NaN;
    num_censored_vols          = 0;
    num_abovethresh_vols_FD    = 0;
    num_abovethresh_vols_DVARS = 0;
end



%%%% Calculate the "ideal" bandpass filter
function [filt_ind] = calc_IdealFilter(nTP,TR,LowCutoff_HighPass,HighCutoff_LowPass)
% Created using the REST toolbox's rest_IdealFilter.m function

if LowCutoff_HighPass<=0
    LowCutoff_HighPass = 0;
elseif HighCutoff_LowPass<=0
    HighCutoff_LowPass = 0;
end
    
sampleFreq 	 = 1/TR;
paddedLength = rest_nextpow2_one35(nTP);

% Get the frequency index
if (LowCutoff_HighPass>=sampleFreq/2) % All high stop
    idxLowCutoff_HighPass = paddedLength/2+1;
else % high pass, such as freq > 0.01 Hz
    idxLowCutoff_HighPass = ceil(LowCutoff_HighPass * paddedLength * TR+1);
end

if (HighCutoff_LowPass>=sampleFreq/2)||(HighCutoff_LowPass==0) % All low pass
    idxHighCutoff_LowPass = paddedLength/2+1;
else % Low pass, such as freq < 0.1 Hz
    idxHighCutoff_LowPass = fix(HighCutoff_LowPass * paddedLength * TR+1);
end

FrequencyMask                                                                                 = zeros(paddedLength,1);
FrequencyMask(idxLowCutoff_HighPass:idxHighCutoff_LowPass,1)                                  = 1;
FrequencyMask(paddedLength-idxLowCutoff_HighPass+2:-1:paddedLength-idxHighCutoff_LowPass+2,1) = 1;

filt_ind = find(FrequencyMask==0);



%%%% Detrend a timeseries for a defined order
function detrend_image(func_filename,funcmask_filename,detr_poly_ord,varargin)

orig = load_nii_gz(func_filename);
mask = load_nii_gz(funcmask_filename);
if ~isa(orig.img,'double')
    orig.img = double(orig.img);
end

if ~isempty(varargin)
    temporal_mask = varargin{1};
    out_filename  = strrep(func_filename,'.nii',['_motcensor_detr',num2str(detr_poly_ord),'.nii']);
else
    temporal_mask = true(size(orig.img,4),1);
    out_filename  = strrep(func_filename,'.nii',['_detr',num2str(detr_poly_ord),'.nii']);
end

orig.img    = orig.img(:,:,:,temporal_mask);
outdata.img = zeros(size(orig.img));

for x = 1:size(orig.img,1)
    for y = 1:size(orig.img,2)
        for z = 1:size(orig.img,3)
            if mask.img(x,y,z)==1
                outdata.img(x,y,z,:) = spm_detrend(squeeze(orig.img(x,y,z,:)),detr_poly_ord);
            end
        end
    end
end

outdata.hdr               = orig.hdr;
outdata.hdr.dime.dim(2:5) = size(outdata.img);
outdata.hdr.dime.datatype = 64;
outdata.hdr.dime.bitpix   = 64;
outdata.hdr.dime.cal_max  = max(outdata.img(:));
outdata.hdr.dime.cal_min  = min(outdata.img(outdata.img~=0));
save_nii_gz(outdata,out_filename);
