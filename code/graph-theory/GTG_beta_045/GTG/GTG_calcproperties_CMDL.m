function GTG_calcproperties_CMDL(config_filename)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: Beta 0.45 (03.30.16)
%
% Usage:    GTG_calcproperties_CMDL('config_filename')
% 
% Input:    config_filename: Name (with path if not in current directory) 
%                            of the .mat file saved by the user. Note that
%                            this file must have all the desired options
%                            selected
% 
% History:
% 03.30.16 - Beta 0.45 - initial release of command line version
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
if ~isfield(out,'conmats') %#ok<NODEF>
    msgbox('Enter assortativity matrices','Error','error')
    return
end
if ~isfield(out,'ROI_labels')
    msgbox('Enter a cell array of ROI labels','Error','error')
    return
end
if ~isfield(out,'outname')
    msgbox('Enter an output name','Error','error')
    return
end
if ~isfield(out,'type_weight_norm') || strcmp(out.type_weight_norm,'Weight Normalization')
    out.type_weight_norm = 'Divide by Mean';
end
if ~isfield(out,'use_zeros') || strcmp(out.use_zeros,'Include Zeros?')
    out.use_zeros = 'Yes';
end
if ~isfield(out,'repmeas_norm_type') || strcmp(out.repmeas_norm_type,'Normalization for Repeated Measures')
    out.repmeas_norm_type = 'Normalize Each Level Individually';
end
if ~isfield(out,'type_dens_thresh') || strcmp(out.type_dens_thresh,'Type of Density Thresholding')
    out.type_dens_thresh = 'Use different thresh (same density)';
end
if ~isfield(out,'denscalc_varmat')
    out.denscalc_varmat    = [];
    out.denscalc_var_names = {'intercept'};
end
if ~isfield(out,'denscalc_covarmat')
    out.denscalc_covarmat    = [];
    out.denscalc_covar_names = {};
end
if ~isfield(out,'partial_for_min_dens')
    out.partial_for_min_dens = 'No';
end
if ~isfield(out,'calcAUC_nodiscon')
    out.calcAUC_nodiscon = 0;
end
if ~isfield(out,'calcbinthresh')
    out.calcbinthresh = 0;
end
if ~isfield(out,'max_dens_pos')
    out.max_dens_pos = 0.6;
end
if ~isfield(out,'dens_step_pos')
    out.dens_step_pos = 0.01;
end
if ~isfield(out,'max_dens_neg')
    out.max_dens_neg = 0.6;
end
if ~isfield(out,'dens_step_neg')
    out.dens_step_neg = 0.01;
end
if ~isfield(out,'num_mod_runs')
    out.num_mod_runs = 1000;
end
if ~isfield(out,'weight_type')
    out.weight_type = 'Positive and Negative';
end
if ~isfield(out,'calc_max_club_size')
    out.calc_max_club_size = 0;
end
if ~isfield(out,'max_club_size')
    out.max_club_size = [];
end
if out.calc_max_club_size==1 && isfield(out,'max_rich_club_size')
    out = rmfield(out,'max_rich_club_size');
elseif out.calc_max_club_size==0 && ~isfield(out,'max_rich_club_size')
    out.max_rich_club_size = 30;
end
if ~isfield(out,'properties_calcd_fullmat') && ~isfield(out,'properties_calcd_thrmat')
    msgbox('Select at least one property to test','Error','error')
    return
elseif ~isfield(out,'properties_calcd_thrmat') && all(strcmp(out.properties_calcd_fullmat,'None')==1)
    msgbox('Select at least one property to test','Error','error')
    return
elseif ~isfield(out,'properties_calcd_fullmat') && all(strcmp(out.properties_calcd_thrmat,'None')==1)
    msgbox('Select at least one property to test','Error','error')
    return
elseif (isfield(out,'properties_calcd_fullmat') && all(strcmp(out.properties_calcd_fullmat,'None')==1)) && (isfield(out,'properties_calcd_thrmat') && all(strcmp(out.properties_calcd_thrmat,'None')==1))
    msgbox('Select at least one property to test','Error','error')
    return
end
if ~isfield(out,'properties_calcd_fullmat')
    out.properties_calcd_fullmat = {};
end
if ~isfield(out,'properties_calcd_thrmat')
    out.properties_calcd_thrmat  = {};
end

if isempty(strfind(out.outname,'/'))
    out.outname = [pwd,'/',out.outname];
elseif out.outname(end)=='/'
    out.outname = [out.outname,'out'];
end
if ~strcmpi(out.outname(end-3:end),'.mat')
    out.outname = [out.outname,'.mat'];
end

out.num_subs = size(out.conmats,3);                                                                                                    % Determine the number of participants
out.nROI     = length(out.ROI_labels);                                                                                         % Determine the number of nodes

% Set intial values to 0 (no testing)
out.calc_props_fullmat.assort            = 0;
out.calc_props_fullmat.bkg               = 0;
out.calc_props_fullmat.cpl               = 0;
out.calc_props_fullmat.close_cent        = 0;
out.calc_props_fullmat.clust_coef        = 0;
out.calc_props_fullmat.clust_coef_ZH     = 0;
out.calc_props_fullmat.clust_coef_signed = 0;
out.calc_props_fullmat.commn_cent        = 0;
out.calc_props_fullmat.div_coef          = 0;
out.calc_props_fullmat.edge_bet_cent     = 0;
out.calc_props_fullmat.eigvec_cent       = 0;
out.calc_props_fullmat.gate_coef         = 0;
out.calc_props_fullmat.glob_eff          = 0;
out.calc_props_fullmat.loc_assort        = 0;
out.calc_props_fullmat.loc_eff           = 0;
out.calc_props_fullmat.match             = 0;
out.calc_props_fullmat.node_bet_cent     = 0;
out.calc_props_fullmat.strength          = 0;
out.calc_props_fullmat.pagerank_cent     = 0;
out.calc_props_fullmat.part_coef         = 0;
out.calc_props_fullmat.rich_club         = 0;
out.calc_props_fullmat.swp               = 0;
out.calc_props_fullmat.trans             = 0;
out.calc_props_fullmat.mod_deg_z         = 0;

out.calc_props_thrmat.assort        = 0;
out.calc_props_thrmat.bkg           = 0;
out.calc_props_thrmat.cpl           = 0;
out.calc_props_thrmat.close_cent    = 0;
out.calc_props_thrmat.clust_coef    = 0;
out.calc_props_thrmat.clust_coef_ZH = 0;
out.calc_props_thrmat.commn_cent    = 0;
out.calc_props_thrmat.deg           = 0;
out.calc_props_thrmat.dens          = 0;
out.calc_props_thrmat.edge_bet_cent = 0;
out.calc_props_thrmat.eigvec_cent   = 0;
out.calc_props_thrmat.gate_coef     = 0;
out.calc_props_thrmat.glob_eff      = 0;
out.calc_props_thrmat.kcore_cent    = 0;
out.calc_props_thrmat.loc_assort    = 0;
out.calc_props_thrmat.loc_eff       = 0;
out.calc_props_thrmat.match         = 0;
out.calc_props_thrmat.node_bet_cent = 0;
out.calc_props_thrmat.pagerank_cent = 0;
out.calc_props_thrmat.part_coef     = 0;
out.calc_props_thrmat.rich_club     = 0;
out.calc_props_thrmat.small_world   = 0;
out.calc_props_thrmat.sub_cent      = 0;
out.calc_props_thrmat.trans         = 0;
out.calc_props_thrmat.mod_deg_z     = 0;

% Determine which properties the user selected for fully connected matrices
if ismember('Assortativity',out.properties_calcd_fullmat) && ~strcmp(out.weight_type,'Absolute Value')
    out.calc_props_fullmat.assort = 1;
elseif ismember('Assortativity',out.properties_calcd_fullmat)
    out.properties_calcd_fullmat(ismember(out.properties_calcd_fullmat,'Assortativity')) = [];
end
if ismember('Brokerage',out.properties_calcd_fullmat)
    out.calc_props_fullmat.bkg   = 1;
    out.properties_calcd_fullmat = [out.properties_calcd_fullmat;'Aggregate Brokerage'];
end
if ismember('Characteristic Path Length',out.properties_calcd_fullmat)
    out.calc_props_fullmat.cpl = 1;
end
if ismember('Closeness Centrality',out.properties_calcd_fullmat)
    out.calc_props_fullmat.close_cent = 1;
end
if ismember('Clustering Coefficient',out.properties_calcd_fullmat)
    out.calc_props_fullmat.clust_coef     = 1;
    out.calc_props_fullmat.clust_coef_tot = 1;
    out.properties_calcd_fullmat          = [out.properties_calcd_fullmat;'Mean Clustering Coefficient'];
end
if ismember('Clustering Coefficient (Zhang & Horvath)',out.properties_calcd_fullmat)
    out.calc_props_fullmat.clust_coef_ZH     = 1;
    out.calc_props_fullmat.clust_coef_ZH_tot = 1;
    out.properties_calcd_fullmat             = [out.properties_calcd_fullmat;'Mean Clustering Coefficient (Zhang & Horvath)'];
end
if ismember('Clustering Coefficient (sign incorporating)',out.properties_calcd_fullmat)
    out.calc_props_fullmat.clust_coef_signed     = 1;
    out.calc_props_fullmat.clust_coef_signed_tot = 1;
    out.properties_calcd_fullmat                 = [out.properties_calcd_fullmat;'Mean Clustering Coefficient (sign incorporating)'];
end
if ismember('Commn Centrality',out.properties_calcd_fullmat)
    out.calc_props_fullmat.commn_cent = 1;
end
if ismember('Diversity Coefficient',out.properties_calcd_fullmat)
    out.calc_props_fullmat.div_coef = 1;
end
if ismember('Edge Betweenness Centrality',out.properties_calcd_fullmat)
    out.calc_props_fullmat.edge_bet_cent = 1;
end
if ismember('Eigenvector Centrality',out.properties_calcd_fullmat)
    out.calc_props_fullmat.eigvec_cent = 1;
end
if ismember('Gateway Coefficient',out.properties_calcd_fullmat)
    out.calc_props_fullmat.gate_coef = 1;
end
if ismember('Global Efficiency',out.properties_calcd_fullmat)
    out.calc_props_fullmat.glob_eff = 1;
end
if ismember('Local Assortativity',out.properties_calcd_fullmat)
    out.calc_props_fullmat.loc_assort = 1;
end
if ismember('Local Efficiency',out.properties_calcd_fullmat)
    out.calc_props_fullmat.loc_eff     = 1;
    out.calc_props_fullmat.loc_eff_tot = 1;
    out.properties_calcd_fullmat       = [out.properties_calcd_fullmat;'Mean Local Efficiency'];
end
if ismember('Matching Index',out.properties_calcd_fullmat)
    out.calc_props_fullmat.match = 1;
end
if ismember('Node Betweenness Centrality',out.properties_calcd_fullmat)
    out.calc_props_fullmat.node_bet_cent = 1;
end
if ismember('Node Strength',out.properties_calcd_fullmat)
    out.calc_props_fullmat.strength     = 1;
    out.calc_props_fullmat.strength_tot = 1;
    out.properties_calcd_fullmat        = [out.properties_calcd_fullmat;'Total Node Strength'];
end
if ismember('PageRank Centrality',out.properties_calcd_fullmat)
    out.calc_props_fullmat.pagerank_cent = 1;
end
if ismember('Participation Coefficient',out.properties_calcd_fullmat)
    out.calc_props_fullmat.part_coef = 1;
end
if ismember('Rich Club Networks',out.properties_calcd_fullmat) && ~strcmp(out.weight_type,'Absolute Value')
    out.calc_props_fullmat.rich_club = 1;
elseif ismember('Rich Club Networks',out.properties_calcd_fullmat)
    out.properties_calcd_fullmat(ismember(out.properties_calcd_fullmat,'Rich Club Networks')) = [];
end
if ismember('Small World Propensity',out.properties_calcd_fullmat)
    out.calc_props_fullmat.swp = 1;
end
if ismember('Transitivity',out.properties_calcd_fullmat)
    out.calc_props_fullmat.trans = 1;
end
if ismember('Within-Module Degree Z-Score',out.properties_calcd_fullmat)
    out.calc_props_fullmat.mod_deg_z = 1;
end

% Determine which properties the user selected for thresholded matrices
if ismember('Assortativity',out.properties_calcd_thrmat)
    out.calc_props_thrmat.assort = 1;
end
if ismember('Brokerage',out.properties_calcd_thrmat)
    out.calc_props_thrmat.bkg   = 1;
    out.properties_calcd_thrmat = [out.properties_calcd_thrmat;'Aggregate Brokerage'];
end
if ismember('Characteristic Path Length',out.properties_calcd_thrmat)
    out.calc_props_thrmat.cpl = 1;
end
if ismember('Closeness Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.close_cent = 1;
end
if ismember('Clustering Coefficient',out.properties_calcd_thrmat)
    out.calc_props_thrmat.clust_coef     = 1;
    out.calc_props_thrmat.clust_coef_tot = 1;
    out.properties_calcd_thrmat          = [out.properties_calcd_thrmat;'Mean Clustering Coefficient'];
end
if ismember('Clustering Coefficient (Zhang & Horvath)',out.properties_calcd_thrmat)
    out.calc_props_thrmat.clust_coef_ZH     = 1;
    out.calc_props_thrmat.clust_coef_ZH_tot = 1;
    out.properties_calcd_thrmat             = [out.properties_calcd_thrmat;'Mean Clustering Coefficient (Zhang & Horvath)'];
end
if ismember('Commn Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.commn_cent = 1;
end
if ismember('Degree',out.properties_calcd_thrmat)
    out.calc_props_thrmat.deg = 1;
end
if ismember('Density',out.properties_calcd_thrmat)
    out.calc_props_thrmat.dens = 1;
end
if ismember('Edge Betweenness Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.edge_bet_cent = 1;
end
if ismember('Eigenvector Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.eigvec_cent = 1;
end
if ismember('Gateway Coefficient',out.properties_calcd_thrmat)
    out.calc_props_thrmat.gate_coef = 1;
end
if ismember('Global Efficiency',out.properties_calcd_thrmat)
    out.calc_props_thrmat.glob_eff = 1;
end
if ismember('K-Coreness Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.kcore_cent = 1;
end
if ismember('Local Assortativity',out.properties_calcd_thrmat)
    out.calc_props_thrmat.loc_assort = 1;
end
if ismember('Local Efficiency',out.properties_calcd_thrmat)
    out.calc_props_thrmat.loc_eff     = 1;
    out.calc_props_thrmat.loc_eff_tot = 1;
    out.properties_calcd_thrmat       = [out.properties_calcd_thrmat;'Mean Local Efficiency'];
end
if ismember('Matching Index',out.properties_calcd_thrmat)
    out.calc_props_thrmat.match = 1;
end
if ismember('Node Betweenness Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.node_bet_cent = 1;
end
if ismember('PageRank Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.pagerank_cent = 1;
end
if ismember('Participation Coefficient',out.properties_calcd_thrmat)
    out.calc_props_thrmat.part_coef = 1;
end
if ismember('Rich Club Networks',out.properties_calcd_thrmat)
    out.calc_props_thrmat.rich_club = 1;
end
if ismember('Small Worldness',out.properties_calcd_thrmat)
    out.calc_props_thrmat.small_world = 1;
end
if ismember('Subgraph Centrality',out.properties_calcd_thrmat)
    out.calc_props_thrmat.sub_cent = 1;
end
if ismember('Transitivity',out.properties_calcd_thrmat)
    out.calc_props_thrmat.trans = 1;
end
if ismember('Within-Module Degree Z-Score',out.properties_calcd_thrmat)
    out.calc_props_thrmat.mod_deg_z = 1;
end

out.conmats(logical(repmat(eye(size(out.conmats,1)),[1,1,size(out.conmats,3),size(out.conmats,4),size(out.conmats,5)]))) = 0;                         % Set diagonal of connectivity matrices to 0
out.conmats      = double(out.conmats);
out.outname      = strrep(out.outname,'.mat','_propcalc.mat');                                                                                                                                        % Create name for output file
out.num_rep_levs = size(out.conmats,4);                                                                                                                                                   % Determine # of repeated levels
use_parfor       = out.use_parfor;
num_par_workers  = out.num_par_workers;

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

if strcmp(out.weight_type,'Absolute Value')
    out.conmats = abs(out.conmats);
end

if out.calc_props_fullmat.assort                 ==1 || ...
        out.calc_props_fullmat.bkg               ==1 || ...
        out.calc_props_fullmat.cpl               ==1 || ...
        out.calc_props_fullmat.close_cent        ==1 || ...
        out.calc_props_fullmat.clust_coef        ==1 || ...
        out.calc_props_fullmat.clust_coef_ZH     ==1 || ...
        out.calc_props_fullmat.clust_coef_signed ==1 || ...
        out.calc_props_fullmat.commn_cent        ==1 || ...
        out.calc_props_fullmat.div_coef          ==1 || ...
        out.calc_props_fullmat.edge_bet_cent     ==1 || ...
        out.calc_props_fullmat.eigvec_cent       ==1 || ...
        out.calc_props_fullmat.gate_coef         ==1 || ...
        out.calc_props_fullmat.glob_eff          ==1 || ...
        out.calc_props_fullmat.loc_assort        ==1 || ...
        out.calc_props_fullmat.loc_eff           ==1 || ...
        out.calc_props_fullmat.match             ==1 || ...
        out.calc_props_fullmat.node_bet_cent     ==1 || ...
        out.calc_props_fullmat.rich_club         ==1 || ...
        out.calc_props_fullmat.strength          ==1 || ...
        out.calc_props_fullmat.pagerank_cent     ==1 || ...    
        out.calc_props_fullmat.part_coef         ==1 || ...
        out.calc_props_fullmat.swp               ==1 || ...
        out.calc_props_fullmat.trans             ==1 || ...
        out.calc_props_fullmat.mod_deg_z         ==1
    out.calcfullmat = 1;
else
    out.calcfullmat = 0;
end

if out.calc_props_thrmat.assort             ==1 || ...
        out.calc_props_thrmat.bkg           ==1 || ...
        out.calc_props_thrmat.cpl           ==1 || ...
        out.calc_props_thrmat.close_cent    ==1 || ...
        out.calc_props_thrmat.clust_coef    ==1 || ...
        out.calc_props_thrmat.clust_coef_ZH ==1 || ...
        out.calc_props_thrmat.commn_cent    ==1 || ...
        out.calc_props_thrmat.deg           ==1 || ...
        out.calc_props_thrmat.dens          ==1 || ...
        out.calc_props_thrmat.edge_bet_cent ==1 || ...
        out.calc_props_thrmat.eigvec_cent   ==1 || ...
        out.calc_props_thrmat.gate_coef     ==1 || ...
        out.calc_props_thrmat.glob_eff      ==1 || ...
        out.calc_props_thrmat.kcore_cent    ==1 || ...
        out.calc_props_thrmat.loc_assort    ==1 || ...
        out.calc_props_thrmat.loc_eff       ==1 || ...
        out.calc_props_thrmat.match         ==1 || ...
        out.calc_props_thrmat.node_bet_cent ==1 || ...
        out.calc_props_thrmat.pagerank_cent ==1 || ...
        out.calc_props_thrmat.part_coef     ==1 || ...
        out.calc_props_thrmat.rich_club     ==1 || ...
        out.calc_props_thrmat.small_world   ==1 || ...
        out.calc_props_thrmat.sub_cent      ==1 || ...
        out.calc_props_thrmat.trans         ==1 || ...
        out.calc_props_thrmat.mod_deg_z     ==1
    out.calcthrmat = 1;
else
    out.calcthrmat = 0;
end

if ~use_parfor
    if out.calcfullmat==1 && out.calcthrmat==1
        progressbar('Progress For Fully Connected Matrices','Progress For Thresholded Matrices','Progress For Calculating AUC for Thresholded Matrices')                                         % Initialize progress bars at zero
    elseif out.calcfullmat==1
        progressbar('Progress For Fully Connected Matrices')                                         % Initialize progress bars at zero
    else
        progressbar('Progress For Thresholded Matrices','Progress For Calculating AUC for Thresholded Matrices')                                         % Initialize progress bars at zero
    end
end

if out.calc_props_fullmat.commn_cent     ==1 || ...
        out.calc_props_fullmat.div_coef  ==1 || ...
        out.calc_props_fullmat.gate_coef ==1 || ...
        out.calc_props_fullmat.part_coef ==1 || ...
        out.calc_props_fullmat.mod_deg_z ==1 || ...
        out.calc_props_thrmat.commn_cent ==1 || ...
        out.calc_props_thrmat.gate_coef  ==1 || ...
        out.calc_props_thrmat.part_coef  ==1 || ...
        out.calc_props_thrmat.mod_deg_z  ==1            % If properties requiring a modular organization should be calculated
    if isfield(out,'mod_grps')
        out.mod_grps = out.mod_grps(:);
    else
        %%%% Create modularity based on mean network
        if size(out.conmats,4)>1
            [full_mean_conmat] = create_mean_conmats(reshape(out.conmats,[size(out.conmats,1),size(out.conmats,2),size(out.conmats,3)*size(out.conmats,4)]));                              % Calculate the mean connectivity matrix across participants; if there are repeated levels, concatenate along the participant dimension beforhand
        else
            [full_mean_conmat] = create_mean_conmats(out.conmats);                              % Calculate the mean connectivity matrix across participants; if there are repeated levels, concatenate along the participant dimension beforhand
        end
        full_mean_conmat(1:size(full_mean_conmat,1)+1:end) = 0;
        
        for run = out.num_mod_runs:-1:1                                                         % Loop on run
            M = community_louvain(full_mean_conmat,[],[],'negative_asym');                      % Calculate the initial organization
            
            % Initialize modularity values
            Q0 = -1;
            Q1 = 0;
            % Iteratively refine modularity
            while Q1-Q0>1e-5
                Q0     = Q1;
                [M,Q1] = community_louvain(full_mean_conmat,[],M,'negative_asym');
            end
            temp_grps(:,run) = M;
            Q(run)           = Q1;
        end
        out.mod_grps = temp_grps(:,find(Q==max(Q),1,'first'));
    end
end

out.conmats_orig = out.conmats;

%%%% Fully connected networks
if out.calcfullmat==1
    if out.calc_props_fullmat.clust_coef_signed==1
        out.conmats_CCsigned = out.conmats;
        switch out.type_weight_norm
            case 'Divide by Mean'
                switch out.use_zeros
                    case 'Yes'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(mean(mean(out.conmats_CCsigned(:,:,curr_sub,rep_lev))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(mean(mean(mean(out.conmats_CCsigned(:,:,curr_sub,:)))));
                                    end
                                end
                        end
                    case 'No'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        tempmat                                    = out.conmats_CCsigned(:,:,curr_sub,rep_lev);
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(mean(mean(tempmat(tempmat~=0))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        tempmat                                    = out.conmats_CCsigned(:,:,curr_sub,:);
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(mean(mean(mean(tempmat(tempmat~=0)))));
                                    end
                                end
                        end
                end
            case 'Divide by Median'
                switch out.use_zeros
                    case 'Yes'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(median(median(out.conmats_CCsigned(:,:,curr_sub,rep_lev))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(median(median(median(out.conmats_CCsigned(:,:,curr_sub,:)))));
                                    end
                                end
                        end
                    case 'No'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        tempmat                                    = out.conmats_CCsigned(:,:,curr_sub,rep_lev);
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(median(median(tempmat(tempmat~=0))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        tempmat                                    = out.conmats_CCsigned(:,:,curr_sub,:);
                                        out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/abs(median(median(median(tempmat(tempmat~=0)))));
                                    end
                                end
                        end
                end
                if any(isnan(out.conmats_CCsigned(:)))
                    disp('WARNING: NaN''s were input in the connectivity matrices and/or normalizing by the median resulted in at least 1 matrix being replaced with NaNs; check out.conmats_CCsigned to determine the extent of the problem')
                end
            case 'Divide by Max'
                switch out.repmeas_norm_type
                    case 'Normalize Each Level Individually'
                        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                            for curr_sub = 1:out.num_subs
                                out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/max(max(abs(out.conmats_CCsigned(:,:,curr_sub,rep_lev))));
                            end
                        end
                    case 'Normalize Across Levels'
                        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                            for curr_sub = 1:out.num_subs
                                out.conmats_CCsigned(:,:,curr_sub,rep_lev) = out.conmats_CCsigned(:,:,curr_sub,rep_lev)/max(max(max(abs(out.conmats_CCsigned(:,:,curr_sub,:)))));
                            end
                        end
                end
        end
        if max(abs(out.conmats_CCsigned(:)))>1
            out.conmats_CCsigned = out.conmats_CCsigned./max(abs(out.conmats_CCsigned(:)));
        end
    end
    
    if strcmp(out.weight_type,'Positive and Negative')
        posmats = out.conmats.*(out.conmats>0);
        negmats = out.conmats.*(out.conmats<0);
        switch out.type_weight_norm
            case 'Divide by Mean'
                switch out.use_zeros
                    case 'Yes'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(posmats(:,:,curr_sub,rep_lev)));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(mean(mean(negmats(:,:,curr_sub,rep_lev))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(mean(posmats(:,:,curr_sub,:))));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(mean(mean(mean(negmats(:,:,curr_sub,:)))));
                                    end
                                end
                        end
                    case 'No'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,rep_lev);
                                        tempneg                       = negmats(:,:,curr_sub,rep_lev);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(temppos(temppos>0)));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(mean(mean(tempneg(tempneg<0))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,:);
                                        tempneg                       = negmats(:,:,curr_sub,:);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(mean(temppos(temppos>0))));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(mean(mean(mean(tempneg(tempneg<0)))));
                                    end
                                end
                        end
                end
            case 'Divide by Median'
                switch out.use_zeros
                    case 'Yes'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(posmats(:,:,curr_sub,rep_lev)));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(median(median(negmats(:,:,curr_sub,rep_lev))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(median(posmats(:,:,curr_sub,:))));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(median(median(median(negmats(:,:,curr_sub,:)))));
                                    end
                                end
                        end
                    case 'No'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,rep_lev);
                                        tempneg                       = negmats(:,:,curr_sub,rep_lev);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(temppos(temppos>0)));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(median(median(tempneg(tempneg<0))));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,:);
                                        tempneg                       = negmats(:,:,curr_sub,:);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(median(temppos(temppos>0))));
                                        negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/abs(median(median(median(tempneg(tempneg<0)))));
                                    end
                                end
                        end
                end
                if any(isnan(posmats(:))) || any(isnan(negmats(:)))
                    disp('WARNING: NaN''s were input in the connectivity matrices and/or normalizing by the median resulted in at least 1 matrix being replaced with NaNs in the full matrices; check out.conmats_full_normed to determine the extent of the problem')
                end
            case 'Divide by Max'
                switch out.repmeas_norm_type
                    case 'Normalize Each Level Individually'
                        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                            for curr_sub = 1:out.num_subs
                                posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/max(max(posmats(:,:,curr_sub,rep_lev)));
                                negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/max(max(abs(negmats(:,:,curr_sub,rep_lev))));
                            end
                        end
                    case 'Normalize Across Levels'
                        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                            for curr_sub = 1:out.num_subs
                                posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/max(max(max(posmats(:,:,curr_sub,:))));
                                negmats(:,:,curr_sub,rep_lev) = negmats(:,:,curr_sub,rep_lev)/max(max(max(abs(negmats(:,:,curr_sub,:)))));
                            end
                        end
                end
        end
        out.conmats = posmats+negmats;
        clear posmats negmats
    else
        posmats = out.conmats.*(out.conmats>0);
        switch out.type_weight_norm
            case 'Divide by Mean'
                switch out.use_zeros
                    case 'Yes'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(posmats(:,:,curr_sub,rep_lev)));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(mean(posmats(:,:,curr_sub,:))));
                                    end
                                end
                        end
                    case 'No'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,rep_lev);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(temppos(temppos>0)));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,:);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/mean(mean(mean(temppos(temppos>0))));
                                    end
                                end
                        end
                end
            case 'Divide by Median'
                switch out.use_zeros
                    case 'Yes'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(posmats(:,:,curr_sub,rep_lev)));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(median(posmats(:,:,curr_sub,:))));
                                    end
                                end
                        end
                    case 'No'
                        switch out.repmeas_norm_type
                            case 'Normalize Each Level Individually'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,rep_lev);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(temppos(temppos>0)));
                                    end
                                end
                            case 'Normalize Across Levels'
                                for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                    for curr_sub = 1:out.num_subs
                                        temppos                       = posmats(:,:,curr_sub,:);
                                        posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/median(median(median(temppos(temppos>0))));
                                    end
                                end
                        end
                end
                if any(isnan(posmats(:)))
                    disp('WARNING: NaN''s were input in the connectivity matrices and/or normalizing by the median resulted in at least 1 matrix being replaced with NaNs in the full matrices; check out.conmats_full_normed to determine the extent of the problem')
                end
            case 'Divide by Max'
                switch out.repmeas_norm_type
                    case 'Normalize Each Level Individually'
                        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                            for curr_sub = 1:out.num_subs
                                posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/max(max(posmats(:,:,curr_sub,rep_lev)));
                            end
                        end
                    case 'Normalize Across Levels'
                        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                            for curr_sub = 1:out.num_subs
                                posmats(:,:,curr_sub,rep_lev) = posmats(:,:,curr_sub,rep_lev)/max(max(max(posmats(:,:,curr_sub,:))));
                            end
                        end
                end
        end
        out.conmats = posmats;
        clear posmats
    end
    
    if max(abs(out.conmats(:)))>1
        out.conmats = out.conmats./max(abs(out.conmats(:)));
    end
    
    fprintf('Calculating properties for fully connected matrices ...\n')
    
    %%%% Calculate network measures for each participant, for each repeated
    %%%% level
    if use_parfor
        calc_props_fullmat = out.calc_props_fullmat;
        num_subs           = out.num_subs;
        weight_type        = out.weight_type;
        if isfield(out,'mod_grps')
            mod_grps = out.mod_grps;
        end
        max_club_size = out.max_club_size;
        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
            conmats = squeeze(out.conmats(:,:,:,rep_lev));
            if calc_props_fullmat.clust_coef_signed==1 && ~strcmp(weight_type,'Absolute Value')                                                                                                        % If the clutering coefficient should be calculated
                conmats_CCsigned = squeeze(out.conmats_CCsigned(:,:,:,rep_lev));
            else
                conmats_CCsigned = nan(size(conmats));
            end
            
            parfor curr_sub = 1:num_subs                                                                                                                                                 % Loop through each participant
                curr_conmat = conmats(:,:,curr_sub);                                                                                                                    % Extract connectivity matrix for current participant
                
                if calc_props_fullmat.clust_coef_signed==1 && ~strcmp(weight_type,'Absolute Value')                                                                                                        % If the clutering coefficient should be calculated
                    [clust_coef_signed(curr_sub,:,rep_lev),~,clust_coef_signed_tot(curr_sub,rep_lev),~] = clustering_coef_wu_sign(conmats_CCsigned(:,:,curr_sub),3);                                        % Each output should be vector of size #ROIs
                end
                
                if strcmp(weight_type,'Positive and Negative')
                    if calc_props_fullmat.assort==1                                                                                                                                        % If assortativity should be calculated
                        [assort_pos(curr_sub,rep_lev),assort_neg(curr_sub,rep_lev)] = assortativity_wei_sign(curr_conmat);                                                                                      % Each output should be 1 number
                    end
                    
                    if calc_props_fullmat.bkg==1                                                                                                                                        % If brokerage should be calculated
                        [bkg_pos(curr_sub,:,rep_lev),bkg_neg(curr_sub,:,rep_lev),bkg_tot_pos(curr_sub,rep_lev),bkg_tot_neg(curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);
                    end
                    
                    if calc_props_fullmat.cpl==1 || calc_props_fullmat.glob_eff==1 || calc_props_fullmat.edge_bet_cent==1 || calc_props_fullmat.node_bet_cent==1                                                                                    % If either node or edge betweenness centrality should be calculated
                        length_mat = weight_conversion(curr_conmat,'lengths');                                                                                                                   % Convert weights to lengths
                        
                        if calc_props_fullmat.cpl==1 || calc_props_fullmat.glob_eff==1                                                                                                   % If either the characteristic path length or global efficiency should be calculated
                            if calc_props_fullmat.cpl==1 && calc_props_fullmat.glob_eff==1                                                                                               % If both should be calculated
                                [cpl_pos(curr_sub,rep_lev),cpl_neg(curr_sub,rep_lev),glob_eff_pos(curr_sub,rep_lev),glob_eff_neg(curr_sub,rep_lev)] = charpath_sign(length_mat);                                               % Each output should be 1 number
                            elseif calc_props_fullmat.cpl==1                                                                                                                                   % If only the characteristic path length should be calculated
                                [cpl_pos(curr_sub,rep_lev),cpl_neg(curr_sub,rep_lev)] = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            else                                                                                                                                                                     % If only the global efficiency should be calculated
                                [~,~,glob_eff_pos(curr_sub,rep_lev),glob_eff_neg(curr_sub,rep_lev)] = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            end
                        end
                        
                        if calc_props_fullmat.node_bet_cent==1                                                                                                                             % If node betweenness centrality  should be calculated
                            [node_bet_cent_pos(curr_sub,:,rep_lev),node_bet_cent_neg(curr_sub,:,rep_lev)] = betweenness_wei_sign(length_mat);                                                                            % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_fullmat.edge_bet_cent==1                                                                                                                             % If edge betweenness centrality  should be calculated
                            [edge_bet_cent_pos(curr_sub,:,:,rep_lev),edge_bet_cent_neg(curr_sub,:,:,rep_lev)] = edge_betweenness_wei_sign(length_mat);                                                                       % Each output should be square matrix of size #ROIs
                        end
                    end
                    
                    if calc_props_fullmat.close_cent==1                                                                                                        % If the clustering coefficient should be calculated
                        [close_cent_pos(curr_sub,:,rep_lev),close_cent_neg(curr_sub,:,rep_lev)] = closeness_cent_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.clust_coef==1                                                                                                        % If the clustering coefficient should be calculated
                        [clust_coef_pos(curr_sub,:,rep_lev),clust_coef_neg(curr_sub,:,rep_lev),clust_coef_tot_pos(curr_sub,rep_lev),clust_coef_tot_neg(curr_sub,rep_lev)] = clustering_coef_wu_sign(curr_conmat,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.clust_coef_ZH==1                                                                                                        % If the clustering coefficient should be calculated
                        [clust_coef_ZH_pos(curr_sub,:,rep_lev),clust_coef_ZH_neg(curr_sub,:,rep_lev),clust_coef_ZH_tot_pos(curr_sub,rep_lev),clust_coef_ZH_tot_neg(curr_sub,rep_lev)] = clustering_coef_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.commn_cent==1                                                                                                        % If the commn centrality should be calculated
                        [commn_cent_pos(curr_sub,:,rep_lev),commn_cent_neg(curr_sub,:,rep_lev)] = commn_cent_wu(curr_conmat,mod_grps);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.div_coef==1                                                                                                                                      % If the diversity coefficient should be calculated
                        [div_coef_pos(curr_sub,:,rep_lev),div_coef_neg(curr_sub,:,rep_lev)] = diversity_coef_sign(curr_conmat,mod_grps);       % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.eigvec_cent==1                                                                                                                                
                        [eigvec_cent_pos(curr_sub,:,rep_lev),eigvec_cent_neg(curr_sub,:,rep_lev)] = eigenvector_centrality_und_sign(curr_conmat);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.gate_coef==1                                                                                                        % If the gateway coefficient should be calculated
                        [gate_coef_pos(curr_sub,:,rep_lev),gate_coef_neg(curr_sub,:,rep_lev)] = gateway_coef_sign(curr_conmat,mod_grps,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.loc_assort==1
                        [loc_assort_pos(curr_sub,:,rep_lev),loc_assort_neg(curr_sub,:,rep_lev)] = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                        [loc_eff_pos(curr_sub,:,rep_lev),loc_eff_neg(curr_sub,:,rep_lev),loc_eff_tot_pos(curr_sub,rep_lev),loc_eff_tot_neg(curr_sub,rep_lev)] = efficiency_wei_sign(curr_conmat,1);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.match==1                                                                                                             % If the matching index should be calculated
                        [match_pos(curr_sub,:,:,rep_lev),match_neg(curr_sub,:,:,rep_lev)] = matching_ind_und_sign(curr_conmat);                                             % Each output should be square matrix of of size #ROIs
                    end
                    
                    if calc_props_fullmat.strength==1                                                                                                                                      % If node strength should be calculated
                        [strength_pos(curr_sub,:,rep_lev),strength_neg(curr_sub,:,rep_lev), ...
                            strength_tot_pos(curr_sub,rep_lev),strength_tot_neg(curr_sub,rep_lev)] = strengths_und_sign(curr_conmat);                % Produces four outputs, two describe entire network, two describe each node
                    end
                    
                    if calc_props_fullmat.pagerank_cent==1                                                                                                                               
                        [pagerank_cent_pos(curr_sub,:,rep_lev),pagerank_cent_neg(curr_sub,:,rep_lev)] = pagerank_centrality_sign(curr_conmat,0.85);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                        [part_coef_pos(curr_sub,:,rep_lev),part_coef_neg(curr_sub,:,rep_lev)] = participation_coef_sign(curr_conmat,mod_grps); % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                        if ~isempty(max_club_size)                                                                                                  % If the user has specified a maximum density
                            [rich_club_pos{curr_sub,rep_lev},rich_club_neg{curr_sub,rep_lev}] = rich_club_wu_sign(curr_conmat,max_club_size);                           % Each output should be vector of size max density
                        else                                                                                                                                 % If not
                            [rich_club_pos{curr_sub,rep_lev},rich_club_neg{curr_sub,rep_lev}] = rich_club_wu_sign(curr_conmat);                                             % Each output should be vector of size of max density
                        end
                    end
                    
                    if calc_props_fullmat.swp==1                                                                                                             % If small world propensity should be calculated
                        curr_conmat_pos           = curr_conmat.*(curr_conmat>0);
                        curr_conmat_neg           = -curr_conmat.*(curr_conmat<0);
                        swp_pos(curr_sub,rep_lev) = small_world_propensity(curr_conmat_pos);                                                  % Each output should be 1 number
                        swp_neg(curr_sub,rep_lev) = small_world_propensity(curr_conmat_neg);                                                  % Each output should be 1 number
                    end
                    
                    if calc_props_fullmat.trans==1                                                                                                             % If transitivity should be calculated
                        [trans_pos(curr_sub,rep_lev),trans_neg(curr_sub,rep_lev)] = transitivity_wu_sign(curr_conmat);                                                  % Each output should be 1 number
                    end
                    
                    if calc_props_fullmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                        [mod_deg_z_pos(curr_sub,:,rep_lev),mod_deg_z_neg(curr_sub,:,rep_lev)] = module_degree_zscore_sign(curr_conmat,mod_grps);                                                                 % Each output should be vector of size #ROIs
                    end
                else
                    if calc_props_fullmat.assort==1                                                                                                                                        % If assortativity should be calculated
                        assort_pos(curr_sub,rep_lev) = assortativity_wei_sign(curr_conmat);                                                                                      % Each output should be 1 number
                    end
                    
                    if calc_props_fullmat.bkg==1                                                                                                                                        % If brokerage should be calculated
                        [bkg_pos(curr_sub,:,rep_lev),~,bkg_tot_pos(curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);
                    end
                    
                    if calc_props_fullmat.cpl==1 || calc_props_fullmat.glob_eff==1 || calc_props_fullmat.edge_bet_cent==1 || calc_props_fullmat.node_bet_cent==1                                                                                    % If either node or edge betweenness centrality should be calculated
                        length_mat = weight_conversion(curr_conmat,'lengths');                                                                                                                   % Convert weights to lengths
                        
                        if calc_props_fullmat.cpl==1 || calc_props_fullmat.glob_eff==1                                                                                                   % If either the characteristic path length or global efficiency should be calculated
                            if calc_props_fullmat.cpl==1 && calc_props_fullmat.glob_eff==1                                                                                               % If both should be calculated
                                [cpl_pos(curr_sub,rep_lev),~,glob_eff_pos(curr_sub,rep_lev)] = charpath_sign(length_mat);                                               % Each output should be 1 number
                            elseif calc_props_fullmat.cpl==1                                                                                                                                   % If only the characteristic path length should be calculated
                                cpl_pos(curr_sub,rep_lev) = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            else                                                                                                                                                                     % If only the global efficiency should be calculated
                                [~,~,glob_eff_pos(curr_sub,rep_lev)] = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            end
                        end
                        
                        if calc_props_fullmat.node_bet_cent==1                                                                                                                             % If node betweenness centrality  should be calculated
                            node_bet_cent_pos(curr_sub,:,rep_lev) = betweenness_wei_sign(length_mat);                                                                            % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_fullmat.edge_bet_cent==1                                                                                                                             % If edge betweenness centrality  should be calculated
                            edge_bet_cent_pos(curr_sub,:,:,rep_lev) = edge_betweenness_wei_sign(length_mat);                                                                       % Each output should be square matrix of size #ROIs
                        end
                    end
                    
                    if calc_props_fullmat.close_cent==1                                                                                                 
                        close_cent_pos(curr_sub,:,rep_lev) = closeness_cent_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                        [clust_coef_pos(curr_sub,:,rep_lev),~,clust_coef_tot_pos(curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.clust_coef_ZH==1                                                                                                        % If the clutering coefficient should be calculated
                        [clust_coef_ZH_pos(curr_sub,:,rep_lev),~,clust_coef_ZH_tot_pos(curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.commn_cent==1                                                                                                        % If the commn centrality should be calculated
                        commn_cent_pos(curr_sub,:,rep_lev) = commn_cent_wu(curr_conmat,mod_grps);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.div_coef==1                                                                                                                                      % If the diversity coefficient should be calculated
                        div_coef_pos(curr_sub,:,rep_lev) = diversity_coef_sign(curr_conmat,mod_grps);       % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.eigvec_cent==1                                                                                                                                   % If the diversity coefficient should be calculated
                        eigvec_cent_pos(curr_sub,:,rep_lev) = eigenvector_centrality_und_sign(curr_conmat);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.gate_coef==1                                                                                                        % If the gateway coefficient should be calculated
                        gate_coef_pos(curr_sub,:,rep_lev) = gateway_coef_sign(curr_conmat,mod_grps,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.loc_assort==1                                                                                                           % If local efficiency should be calculated
                        loc_assort_pos(curr_sub,:,rep_lev) = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                        [loc_eff_pos(curr_sub,:,rep_lev),~,loc_eff_tot_pos(curr_sub,rep_lev),~] = efficiency_wei_sign(curr_conmat,1);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.match==1                                                                                                             % If the matching index should be calculated
                        match_pos(curr_sub,:,:,rep_lev) = matching_ind_und_sign(curr_conmat);                                             % Each output should be square matrix of of size #ROIs
                    end
                    
                    if calc_props_fullmat.strength==1                                                                                                                                      % If node strength should be calculated
                        [strength_pos(curr_sub,:,rep_lev),~, ...
                            strength_tot_pos(curr_sub,rep_lev)] = strengths_und_sign(curr_conmat);                % Produces four outputs, two describe entire network, two describe each node
                    end
                    
                    if calc_props_fullmat.pagerank_cent==1                                                                                                                               
                        pagerank_cent_pos(curr_sub,:,rep_lev) = pagerank_centrality_sign(curr_conmat,0.85);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                        part_coef_pos(curr_sub,:,rep_lev) = participation_coef_sign(curr_conmat,mod_grps); % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if calc_props_fullmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                        if ~isempty(max_club_size)                                                                                                  % If the user has specified a maximum density
                            rich_club_pos{curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat,max_club_size);                           % Each output should be vector of size max density
                        else                                                                                                                                 % If not
                            rich_club_pos{curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat);                                             % Each output should be vector of size of max density
                        end
                    end
                    
                    if calc_props_fullmat.swp==1                                                                                                             % If small world propensity should be calculated
                        curr_conmat_pos = curr_conmat.*(curr_conmat>0);
                        swp_pos(curr_sub,rep_lev) = small_world_propensity(curr_conmat_pos);                                                  % Each output should be 1 number
                    end
                    
                    if calc_props_fullmat.trans==1                                                                                                             % If transitivity should be calculated
                        trans_pos(curr_sub,rep_lev) = transitivity_wu_sign(curr_conmat);                                                  % Each output should be 1 number
                    end
                    
                    if calc_props_fullmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                        mod_deg_z_pos(curr_sub,:,rep_lev) = module_degree_zscore_sign(curr_conmat,mod_grps);                                                                 % Each output should be vector of size #ROIs
                    end
                end
            end
        end
        
        if strcmp(out.weight_type,'Positive and Negative')
            if calc_props_fullmat.assort==1                                                                                                                                        % If assortativity was calculated
                out.fullmat_graph_meas.assort_pos = assort_pos;
                out.fullmat_graph_meas.assort_neg = assort_neg;
                clear assort_pos assort_neg
            end
            if calc_props_fullmat.bkg==1
                out.fullmat_graph_meas.bkg_pos     = bkg_pos;
                out.fullmat_graph_meas.bkg_neg     = bkg_neg;
                out.fullmat_graph_meas.bkg_tot_pos = bkg_tot_pos;
                out.fullmat_graph_meas.bkg_tot_neg = bkg_tot_neg;
                clear bkg_pos bkg_neg bkg_tot_pos bkg_tot_neg
            end        
            if calc_props_fullmat.cpl==1                                                                                                                                   % If only the characteristic path length was calculated
                out.fullmat_graph_meas.cpl_pos = cpl_pos;
                out.fullmat_graph_meas.cpl_neg = cpl_neg;
                clear cpl_pos cpl_neg
            end
            if calc_props_fullmat.close_cent==1                                                
                out.fullmat_graph_meas.close_cent_pos = close_cent_pos;
                out.fullmat_graph_meas.close_cent_neg = close_cent_neg;
                clear close_cent_pos close_cent_neg
            end
            if calc_props_fullmat.clust_coef==1
                out.fullmat_graph_meas.clust_coef_pos     = clust_coef_pos;
                out.fullmat_graph_meas.clust_coef_neg     = clust_coef_neg;
                out.fullmat_graph_meas.clust_coef_tot_pos = clust_coef_tot_pos;
                out.fullmat_graph_meas.clust_coef_tot_neg = clust_coef_tot_neg;
                clear clust_coef_pos clust_coef_neg clust_coef_tot_pos clust_coef_tot_neg
            end
            if calc_props_fullmat.clust_coef_ZH==1                                             
                out.fullmat_graph_meas.clust_coef_ZH_pos     = clust_coef_ZH_pos;
                out.fullmat_graph_meas.clust_coef_ZH_neg     = clust_coef_ZH_neg;
                out.fullmat_graph_meas.clust_coef_ZH_tot_pos = clust_coef_ZH_tot_pos;
                out.fullmat_graph_meas.clust_coef_ZH_tot_neg = clust_coef_ZH_tot_neg;
                clear clust_coef_ZH_pos clust_coef_ZH_neg clust_coef_ZH_tot_pos clust_coef_ZH_tot_neg
            end
            if calc_props_fullmat.clust_coef_signed==1
                out.fullmat_graph_meas.clust_coef_signed     = clust_coef_signed;
                out.fullmat_graph_meas.clust_coef_signed_tot = clust_coef_signed_tot;
                clear clust_coef_signed clust_coef_signed_tot
            end
            if calc_props_fullmat.commn_cent==1                                                
                out.fullmat_graph_meas.commn_cent_pos = commn_cent_pos;
                out.fullmat_graph_meas.commn_cent_neg = commn_cent_neg;
                clear commn_cent_pos commn_cent_neg
            end
            if calc_props_fullmat.div_coef==1                                                                                                                                      % If the diversity coefficient was calculated
                out.fullmat_graph_meas.div_coef_pos = div_coef_pos;
                out.fullmat_graph_meas.div_coef_neg = div_coef_neg;
                clear div_coef_pos div_coef_neg
            end
            if calc_props_fullmat.edge_bet_cent==1                                                                                                                             % If edge betweenness centrality  was calculated
                out.fullmat_graph_meas.edge_bet_cent_pos = edge_bet_cent_pos;
                out.fullmat_graph_meas.edge_bet_cent_neg = edge_bet_cent_neg;
                clear edge_bet_cent_pos edge_bet_cent_neg
            end
            if calc_props_fullmat.eigvec_cent==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.eigvec_cent_pos = eigvec_cent_pos;
                out.fullmat_graph_meas.eigvec_cent_neg = eigvec_cent_neg;
                clear eigvec_cent_pos eigvec_cent_neg
            end
            if calc_props_fullmat.gate_coef==1                                                                                                                                         % If only the gateway coefficient was calculated
                out.fullmat_graph_meas.gate_coef_pos = gate_coef_pos;
                out.fullmat_graph_meas.gate_coef_neg = gate_coef_neg;
                clear gate_coef_pos gate_coef_neg
            end
            if calc_props_fullmat.glob_eff==1                                                                                                                                         % If only the global efficiency was calculated
                out.fullmat_graph_meas.glob_eff_pos = glob_eff_pos;
                out.fullmat_graph_meas.glob_eff_neg = glob_eff_neg;
                clear glob_eff_pos glob_eff_neg
            end
            if calc_props_fullmat.loc_assort==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.loc_assort_pos = loc_assort_pos;
                out.fullmat_graph_meas.loc_assort_neg = loc_assort_neg;
                clear loc_assort_pos loc_assort_neg
            end
            if calc_props_fullmat.loc_eff==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.loc_eff_pos     = loc_eff_pos;
                out.fullmat_graph_meas.loc_eff_neg     = loc_eff_neg;
                out.fullmat_graph_meas.loc_eff_tot_pos = loc_eff_tot_pos;
                out.fullmat_graph_meas.loc_eff_tot_neg = loc_eff_tot_neg;
                clear loc_eff_pos loc_eff_neg loc_eff_tot_pos loc_eff_tot_neg
            end
            if calc_props_fullmat.match==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.match_pos = match_pos;
                out.fullmat_graph_meas.match_neg = match_neg;
                clear match_pos match_neg
            end
            if calc_props_fullmat.node_bet_cent==1                                                                                                                             % If node betweenness centrality  was calculated
                out.fullmat_graph_meas.node_bet_cent_pos = node_bet_cent_pos;
                out.fullmat_graph_meas.node_bet_cent_neg = node_bet_cent_neg;
                clear node_bet_cent_pos node_bet_cent_neg
            end
            if calc_props_fullmat.strength==1                                                                                                                                      % If node strength was calculated
                out.fullmat_graph_meas.strength_pos     = strength_pos;
                out.fullmat_graph_meas.strength_neg     = strength_neg;
                out.fullmat_graph_meas.strength_tot_pos = strength_tot_pos;
                out.fullmat_graph_meas.strength_tot_neg = strength_tot_neg;
                clear strength_pos strength_neg strength_tot_pos strength_tot_neg
            end
            if calc_props_fullmat.pagerank_cent==1                                                                                                                                 % If the diversity coefficient was calculated
                out.fullmat_graph_meas.pagerank_cent_pos = pagerank_cent_pos;
                out.fullmat_graph_meas.pagerank_cent_neg = pagerank_cent_neg;
                clear pagerank_cent_pos pagerank_cent_neg
            end
            if calc_props_fullmat.part_coef==1                                                                                                                                     % If the participation coefficient was calculated
                out.fullmat_graph_meas.part_coef_pos = part_coef_pos;
                out.fullmat_graph_meas.part_coef_neg = part_coef_neg;
                clear part_coef_pos part_coef_neg
            end
            if calc_props_fullmat.rich_club==1                                                                                                                                     % If the participation coefficient was calculated
                out.fullmat_graph_meas.rich_club_pos = rich_club_pos;
                out.fullmat_graph_meas.rich_club_neg = rich_club_neg;
                clear rich_club_pos rich_club_neg
            end
            if calc_props_fullmat.swp==1                                                                                                                                     % If the small world propensity was calculated
                out.fullmat_graph_meas.swp_pos = swp_pos;
                out.fullmat_graph_meas.swp_neg = swp_neg;
                clear swp_pos swp_neg
            end
            if calc_props_fullmat.trans==1                                                                                                                                     % If the participation coefficient was calculated
                out.fullmat_graph_meas.trans_pos = trans_pos;
                out.fullmat_graph_meas.trans_neg = trans_neg;
                clear trans_pos trans_neg
            end
            if calc_props_fullmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score was calculated
                out.fullmat_graph_meas.mod_deg_z_pos = mod_deg_z_pos;
                out.fullmat_graph_meas.mod_deg_z_neg = mod_deg_z_neg;
                clear mod_deg_z_pos mod_deg_z_neg
            end
        else
            if calc_props_fullmat.assort==1                                                                                                                                        % If assortativity was calculated
                out.fullmat_graph_meas.assort_pos = assort_pos;
                clear assort_pos
            end
            if calc_props_fullmat.bkg==1
                out.fullmat_graph_meas.bkg_pos     = bkg_pos;
                out.fullmat_graph_meas.bkg_tot_pos = bkg_tot_pos;
                clear bkg_pos bkg_tot_pos
            end
            if calc_props_fullmat.cpl==1                                                                                                                                   % If only the characteristic path length was calculated
                out.fullmat_graph_meas.cpl_pos = cpl_pos;
                clear cpl_pos
            end
            if calc_props_fullmat.close_cent==1                                                                                                                                   % If only the characteristic path length was calculated
                out.fullmat_graph_meas.close_cent_pos = close_cent_pos;
                clear close_cent_pos
            end
            if calc_props_fullmat.clust_coef==1                                                                                                                                   % If only the characteristic path length was calculated
                out.fullmat_graph_meas.clust_coef_pos     = clust_coef_pos;
                out.fullmat_graph_meas.clust_coef_tot_pos = clust_coef_tot_pos;
                clear clust_coef_pos clust_coef_tot_pos
            end
            if calc_props_fullmat.clust_coef_ZH==1                                             
                out.fullmat_graph_meas.clust_coef_ZH_pos     = clust_coef_ZH_pos;
                out.fullmat_graph_meas.clust_coef_ZH_tot_pos = clust_coef_ZH_tot_pos;
                clear clust_coef_ZH_pos clust_coef_ZH_tot_pos
            end
            if calc_props_fullmat.clust_coef_signed==1 && ~strcmp(out.weight_type,'Absolute Value')  
                out.fullmat_graph_meas.clust_coef_signed     = clust_coef_signed;
                out.fullmat_graph_meas.clust_coef_signed_tot = clust_coef_signed_tot;
                clear clust_coef_signed clust_coef_signed_tot
            end
            if calc_props_fullmat.commn_cent==1                                                                                                                                   % If only the characteristic path length was calculated
                out.fullmat_graph_meas.commn_cent_pos = commn_cent_pos;
                clear commn_cent_pos
            end
            if calc_props_fullmat.div_coef==1                                                                                                                                      % If the diversity coefficient was calculated
                out.fullmat_graph_meas.div_coef_pos = div_coef_pos;
                clear div_coef_pos
            end
            if calc_props_fullmat.edge_bet_cent==1                                                                                                                             % If edge betweenness centrality  was calculated
                out.fullmat_graph_meas.edge_bet_cent_pos = edge_bet_cent_pos;
                clear edge_bet_cent_pos
            end
            if calc_props_fullmat.eigvec_cent==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.eigvec_cent_pos = eigvec_cent_pos;
                clear eigvec_cent_pos
            end
            if calc_props_fullmat.gate_coef==1                                                                                                                                         % If only the gateway coefficient was calculated
                out.fullmat_graph_meas.gate_coef_pos = gate_coef_pos;
                clear gate_coef_pos
            end
            if calc_props_fullmat.glob_eff==1                                                                                                                                         % If only the global efficiency was calculated
                out.fullmat_graph_meas.glob_eff_pos = glob_eff_pos;
                clear glob_eff_pos
            end
            if calc_props_fullmat.loc_assort==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.loc_assort_pos = loc_assort_pos;
                clear loc_assort_pos
            end
            if calc_props_fullmat.loc_eff==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.loc_eff_pos     = loc_eff_pos;
                out.fullmat_graph_meas.loc_eff_tot_pos = loc_eff_tot_pos;
                clear loc_eff_pos loc_eff_tot_pos
            end
            if calc_props_fullmat.match==1                                                                                                                                   % If the diversity coefficient was calculated
                out.fullmat_graph_meas.match_pos = match_pos;
                clear match_pos
            end
            if calc_props_fullmat.node_bet_cent==1                                                                                                                             % If node betweenness centrality  was calculated
                out.fullmat_graph_meas.node_bet_cent_pos = node_bet_cent_pos;
                clear node_bet_cent_pos
            end
            if calc_props_fullmat.strength==1                                                                                                                                      % If node strength was calculated
                out.fullmat_graph_meas.strength_pos     = strength_pos;
                out.fullmat_graph_meas.strength_tot_pos = strength_tot_pos;
                clear strength_pos strength_tot_pos
            end
            if calc_props_fullmat.pagerank_cent==1                                                                                                                                 % If the diversity coefficient was calculated
                out.fullmat_graph_meas.pagerank_cent_pos = pagerank_cent_pos;
                clear pagerank_cent_pos
            end
            if calc_props_fullmat.part_coef==1                                                                                                                                     % If the participation coefficient was calculated
                out.fullmat_graph_meas.part_coef_pos = part_coef_pos;
                clear part_coef_pos
            end
            if calc_props_fullmat.rich_club==1                                                                                                                                     % If the participation coefficient was calculated
                out.fullmat_graph_meas.rich_club_pos = rich_club_pos;
                clear rich_club_pos
            end
            if calc_props_fullmat.swp==1                                                                                                                                     % If the small world propensity was calculated
                out.fullmat_graph_meas.swp_pos = swp_pos;
                clear swp_pos
            end
            if calc_props_fullmat.trans==1                                                                                                                                     % If the participation coefficient was calculated
                out.fullmat_graph_meas.trans_pos = trans_pos;
                clear trans_pos
            end
            if calc_props_fullmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score was calculated
                out.fullmat_graph_meas.mod_deg_z_pos = mod_deg_z_pos;
                clear mod_deg_z_pos
            end
        end
    else
        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
            for curr_sub = 1:out.num_subs                                                                                                                                                 % Loop through each participant
                curr_conmat = squeeze(out.conmats(:,:,curr_sub,rep_lev));                                                                                                                    % Extract connectivity matrix for current participant
                
                if out.calc_props_fullmat.clust_coef_signed==1 && ~strcmp(out.weight_type,'Absolute Value')                                                                                                        % If the clutering coefficient should be calculated
                    [out.fullmat_graph_meas.clust_coef_signed(curr_sub,:,rep_lev),~,out.fullmat_graph_meas.clust_coef_signed_tot(curr_sub,rep_lev)] = clustering_coef_wu_sign(out.conmats_CCsigned(:,:,curr_sub,rep_lev),3);                                        % Each output should be vector of size #ROIs
                end
                
                if strcmp(out.weight_type,'Positive and Negative')
                    if out.calc_props_fullmat.assort==1                                                                                                                                        % If assortativity should be calculated
                        [out.fullmat_graph_meas.assort_pos(curr_sub,rep_lev),out.fullmat_graph_meas.assort_neg(curr_sub,rep_lev)] = assortativity_wei_sign(curr_conmat);                                                                                      % Each output should be 1 number
                    end
                    
                    if out.calc_props_fullmat.bkg==1                                                                                                                                        % If brokerage should be calculated
                        [out.fullmat_graph_meas.bkg_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.bkg_neg(curr_sub,:,rep_lev),out.fullmat_graph_meas.bkg_tot_pos(curr_sub,rep_lev),out.fullmat_graph_meas.bkg_tot_neg(curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);
                    end
                    
                    if out.calc_props_fullmat.cpl==1 || out.calc_props_fullmat.glob_eff==1 || out.calc_props_fullmat.edge_bet_cent==1 || out.calc_props_fullmat.node_bet_cent==1                                                                                    % If either node or edge betweenness centrality should be calculated
                        length_mat = weight_conversion(curr_conmat,'lengths');                                                                                                                   % Convert weights to lengths
                        
                        if out.calc_props_fullmat.cpl==1 || out.calc_props_fullmat.glob_eff==1                                                                                                   % If either the characteristic path length or global efficiency should be calculated
                            if out.calc_props_fullmat.cpl==1 && out.calc_props_fullmat.glob_eff==1                                                                                               % If both should be calculated
                                [out.fullmat_graph_meas.cpl_pos(curr_sub,rep_lev),out.fullmat_graph_meas.cpl_neg(curr_sub,rep_lev), ...
                                    out.fullmat_graph_meas.glob_eff_pos(curr_sub,rep_lev),out.fullmat_graph_meas.glob_eff_neg(curr_sub,rep_lev)] = charpath_sign(length_mat);                                               % Each output should be 1 number
                            elseif out.calc_props_fullmat.cpl==1                                                                                                                                   % If only characteristic path length should be calculated
                                [out.fullmat_graph_meas.cpl_pos(curr_sub,rep_lev),out.fullmat_graph_meas.cpl_neg(curr_sub,rep_lev)]        = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            else                                                                                                                                                                     % If only global efficiency should be calculated
                                [~,~,out.fullmat_graph_meas.glob_eff_pos(curr_sub,rep_lev),out.fullmat_graph_meas.glob_eff_neg(curr_sub,rep_lev)] = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            end
                        end
                        
                        if out.calc_props_fullmat.node_bet_cent==1                                                                                                                             % If node betweenness centrality  should be calculated
                            [out.fullmat_graph_meas.node_bet_cent_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.node_bet_cent_neg(curr_sub,:,rep_lev)]   = betweenness_wei_sign(length_mat);                                                                            % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_fullmat.edge_bet_cent==1                                                                                                                             % If edge betweenness centrality  should be calculated
                            [out.fullmat_graph_meas.edge_bet_cent_pos(curr_sub,:,:,rep_lev),out.fullmat_graph_meas.edge_bet_cent_neg(curr_sub,:,:,rep_lev)] = edge_betweenness_wei_sign(length_mat);                                                                       % Each output should be square matrix of size #ROIs
                        end
                    end
                    
                    if out.calc_props_fullmat.close_cent==1
                        [out.fullmat_graph_meas.close_cent_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.close_cent_neg(curr_sub,:,rep_lev)] = closeness_cent_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                        [out.fullmat_graph_meas.clust_coef_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.clust_coef_neg(curr_sub,:,rep_lev),out.fullmat_graph_meas.clust_coef_tot_pos(curr_sub,rep_lev),out.fullmat_graph_meas.clust_coef_tot_neg(curr_sub,rep_lev)] = clustering_coef_wu_sign(curr_conmat,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.clust_coef_ZH==1                                                                                                        % If the clutering coefficient should be calculated
                        [out.fullmat_graph_meas.clust_coef_ZH_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.clust_coef_ZH_neg(curr_sub,:,rep_lev),out.fullmat_graph_meas.clust_coef_ZH_tot_pos(curr_sub,rep_lev),out.fullmat_graph_meas.clust_coef_ZH_tot_neg(curr_sub,rep_lev)] = clustering_coef_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.commn_cent==1                                                                                                        % If the commn centrality should be calculated
                        [out.fullmat_graph_meas.commn_cent_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.commn_cent_neg(curr_sub,:,rep_lev)] = commn_cent_wu(curr_conmat,out.mod_grps);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.div_coef==1                                                                                                                                      % If the diversity coefficient should be calculated
                        [out.fullmat_graph_meas.div_coef_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.div_coef_neg(curr_sub,:,rep_lev)] = diversity_coef_sign(curr_conmat,out.mod_grps);       % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.eigvec_cent==1                                                                                                                                   % If the diversity coefficient should be calculated
                        [out.fullmat_graph_meas.eigvec_cent_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.eigvec_cent_neg(curr_sub,:,rep_lev)] = eigenvector_centrality_und_sign(curr_conmat);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.gate_coef==1                                                                                                        % If the gateway coefficient should be calculated
                        [out.fullmat_graph_meas.gate_coef_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.gate_coef_neg(curr_sub,:,rep_lev)] = gateway_coef_sign(curr_conmat,out.mod_grps,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.loc_assort==1
                        [out.fullmat_graph_meas.loc_assort_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.loc_assort_neg(curr_sub,:,rep_lev)] = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                        [out.fullmat_graph_meas.loc_eff_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.loc_eff_neg(curr_sub,:,rep_lev),out.fullmat_graph_meas.loc_eff_tot_pos(curr_sub,rep_lev),out.fullmat_graph_meas.loc_eff_tot_neg(curr_sub,rep_lev)] = efficiency_wei_sign(curr_conmat,1);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.match==1                                                                                                             % If the matching index should be calculated
                        [out.fullmat_graph_meas.match_pos(curr_sub,:,:,rep_lev),out.fullmat_graph_meas.match_neg(curr_sub,:,:,rep_lev)] = matching_ind_und_sign(curr_conmat);                                             % Each output should be square matrix of of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.strength==1                                                                                                                                      % If node strength should be calculated
                        [out.fullmat_graph_meas.strength_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.strength_neg(curr_sub,:,rep_lev), ...
                            out.fullmat_graph_meas.strength_tot_pos(curr_sub,rep_lev),out.fullmat_graph_meas.strength_tot_neg(curr_sub,rep_lev)] = strengths_und_sign(curr_conmat);                % Produces four outputs, two describe entire network, two describe each node
                    end
                    
                    if out.calc_props_fullmat.pagerank_cent==1                                                                                                                                 % If the diversity coefficient should be calculated
                        [out.fullmat_graph_meas.pagerank_cent_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.pagerank_cent_neg(curr_sub,:,rep_lev)] = pagerank_centrality_sign(curr_conmat,0.85);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                        [out.fullmat_graph_meas.part_coef_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.part_coef_neg(curr_sub,:,rep_lev)] = participation_coef_sign(curr_conmat,out.mod_grps); % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                        if ~isempty(out.max_club_size)
                            [out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev},out.fullmat_graph_meas.rich_club_neg{curr_sub,rep_lev}] = rich_club_wu_sign(curr_conmat,out.max_club_size);                           % Each output should be vector of size max density
                        else                                                                                                                                 % If not
                            [out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev},out.fullmat_graph_meas.rich_club_neg{curr_sub,rep_lev}] = rich_club_wu_sign(curr_conmat);                                             % Each output should be vector of size of max density
                        end
                    end
                    
                    if out.calc_props_fullmat.swp==1                                                                                                             % If small world propensity should be calculated
                        curr_conmat_pos                                  = abs(curr_conmat.*(curr_conmat>0));
                        curr_conmat_neg                                  = abs(-curr_conmat.*(curr_conmat<0));
                        out.fullmat_graph_meas.swp_pos(curr_sub,rep_lev) = small_world_propensity(curr_conmat_pos);                                                  % Each output should be 1 number
                        out.fullmat_graph_meas.swp_neg(curr_sub,rep_lev) = small_world_propensity(curr_conmat_neg);                                                  % Each output should be 1 number
                    end
                    
                    if out.calc_props_fullmat.trans==1                                                                                                             % If transitivity should be calculated
                        [out.fullmat_graph_meas.trans_pos(curr_sub,rep_lev),out.fullmat_graph_meas.trans_neg(curr_sub,rep_lev)] = transitivity_wu_sign(curr_conmat);                                                  % Each output should be 1 number
                    end
                    
                    if out.calc_props_fullmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                        [out.fullmat_graph_meas.mod_deg_z_pos(curr_sub,:,rep_lev),out.fullmat_graph_meas.mod_deg_z_neg(curr_sub,:,rep_lev)] = module_degree_zscore_sign(curr_conmat,out.mod_grps);                                                                 % Each output should be vector of size #ROIs
                    end
                    
                    prog = (curr_sub/out.num_subs)*(1-((rep_lev-1)/out.num_rep_levs));                                                                                                         % Calculate progress
                    if out.calcthrmat==1
                        progressbar(prog,[],[])                                                                                                                                                      % Update progress bar
                    else
                        progressbar(prog)                                                                                                                                                      % Update progress bar
                    end
                else
                    if out.calc_props_fullmat.assort==1                                                                                                                                        % If assortativity should be calculated
                        out.fullmat_graph_meas.assort_pos(curr_sub,rep_lev) = assortativity_wei_sign(curr_conmat);                                                                                      % Each output should be 1 number
                    end
                    
                    if out.calc_props_fullmat.bkg==1                                                                                                                                        % If brokerage should be calculated
                        [out.fullmat_graph_meas.bkg_pos(curr_sub,:,rep_lev),~,out.fullmat_graph_meas.bkg_tot_pos(curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);
                    end
                    
                    if out.calc_props_fullmat.cpl==1 || out.calc_props_fullmat.glob_eff==1 || out.calc_props_fullmat.edge_bet_cent==1 || out.calc_props_fullmat.node_bet_cent==1                                                                                    % If either node or edge betweenness centrality should be calculated
                        length_mat = weight_conversion(curr_conmat,'lengths');                                                                                                                   % Convert weights to lengths
                        
                        if out.calc_props_fullmat.cpl==1 || out.calc_props_fullmat.glob_eff==1                                                                                                   % If either the characteristic path length or global efficiency should be calculated
                            if out.calc_props_fullmat.cpl==1 && out.calc_props_fullmat.glob_eff==1                                                                                               % If both should be calculated
                                [out.fullmat_graph_meas.cpl_pos(curr_sub,rep_lev),~, ...
                                    out.fullmat_graph_meas.glob_eff_pos(curr_sub,rep_lev)] = charpath_sign(length_mat);                                               % Each output should be 1 number
                            elseif out.calc_props_fullmat.cpl==1                                                                                                                                   % If only the characteristic path length should be calculated
                                out.fullmat_graph_meas.cpl_pos(curr_sub,rep_lev) = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            else                                                                                                                                                                     % If only the global efficiency should be calculated
                                [~,~,out.fullmat_graph_meas.glob_eff_pos(curr_sub,rep_lev)] = charpath_sign(length_mat);                                                                                          % Each output should be 1 number
                            end
                        end
                        
                        if out.calc_props_fullmat.node_bet_cent==1                                                                                                                             % If node betweenness centrality  should be calculated
                            out.fullmat_graph_meas.node_bet_cent_pos(curr_sub,:,rep_lev) = betweenness_wei_sign(length_mat);                                                                            % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_fullmat.edge_bet_cent==1                                                                                                                             % If edge betweenness centrality  should be calculated
                            out.fullmat_graph_meas.edge_bet_cent_pos(curr_sub,:,:,rep_lev) = edge_betweenness_wei_sign(length_mat);                                                                       % Each output should be square matrix of size #ROIs
                        end
                    end
                    
                    if out.calc_props_fullmat.close_cent==1
                        out.fullmat_graph_meas.close_cent_pos(curr_sub,:,rep_lev) = closeness_cent_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                        [out.fullmat_graph_meas.clust_coef_pos(curr_sub,:,rep_lev),~,out.fullmat_graph_meas.clust_coef_tot_pos(curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.clust_coef_ZH==1                                                                                                        % If the clutering coefficient should be calculated
                        [out.fullmat_graph_meas.clust_coef_ZH_pos(curr_sub,:,rep_lev),~,out.fullmat_graph_meas.clust_coef_ZH_tot_pos(curr_sub,rep_lev)] = clustering_coef_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.commn_cent==1                                                                                                        % If the commn centrality should be calculated
                        out.fullmat_graph_meas.commn_cent_pos(curr_sub,:,rep_lev) = commn_cent_wu(curr_conmat,out.mod_grps);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.div_coef==1                                                                                                                                      % If the diversity coefficient should be calculated
                        out.fullmat_graph_meas.div_coef_pos(curr_sub,:,rep_lev) = diversity_coef_sign(curr_conmat,out.mod_grps);       % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.eigvec_cent==1                                                                                                                                   % If the diversity coefficient should be calculated
                        out.fullmat_graph_meas.eigvec_cent_pos(curr_sub,:,rep_lev) = eigenvector_centrality_und_sign(curr_conmat);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.gate_coef==1                                                                                                        % If the gateway coefficient should be calculated
                        out.fullmat_graph_meas.gate_coef_pos(curr_sub,:,rep_lev) = gateway_coef_sign(curr_conmat,out.mod_grps,1);                                        % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.loc_assort==1
                        out.fullmat_graph_meas.loc_assort_pos(curr_sub,:,rep_lev) = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                        [out.fullmat_graph_meas.loc_eff_pos(curr_sub,:,rep_lev),~,out.fullmat_graph_meas.loc_eff_tot_pos(curr_sub,rep_lev),~] = efficiency_wei_sign(curr_conmat,1);                                             % Each output should be vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.match==1                                                                                                             % If the matching index should be calculated
                        out.fullmat_graph_meas.match_pos(curr_sub,:,:,rep_lev) = matching_ind_und_sign(curr_conmat);                                             % Each output should be square matrix of of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.strength==1                                                                                                                                      % If node strength should be calculated
                        [out.fullmat_graph_meas.strength_pos(curr_sub,:,rep_lev),~, ...
                            out.fullmat_graph_meas.strength_tot_pos(curr_sub,rep_lev)] = strengths_und_sign(curr_conmat);                % Produces four outputs, two describe entire network, two describe each node
                    end
                    
                    if out.calc_props_fullmat.pagerank_cent==1                                                                                                                                 % If the diversity coefficient should be calculated
                        out.fullmat_graph_meas.pagerank_cent_pos(curr_sub,:,rep_lev) = pagerank_centrality_sign(curr_conmat,0.85);                                                                        % Produces vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                        out.fullmat_graph_meas.part_coef_pos(curr_sub,:,rep_lev) = participation_coef_sign(curr_conmat,out.mod_grps); % Produces two outputs, each a vector of size #ROIs
                    end
                    
                    if out.calc_props_fullmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                        if ~isempty(out.max_club_size)                                                                                                  % If the user has specified a maximum density
                            out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat,out.max_club_size);                           % Each output should be vector of size max density
                        else                                                                                                                                 % If not
                            out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat);                                             % Each output should be vector of size of max density
                        end
                    end
                    
                    if out.calc_props_fullmat.swp==1                                                                                                             % If small world propensity should be calculated
                        curr_conmat_pos                                  = abs(curr_conmat.*(curr_conmat>0));
                        out.fullmat_graph_meas.swp_pos(curr_sub,rep_lev) = small_world_propensity(curr_conmat_pos);                                                  % Each output should be 1 number
                    end
                    
                    if out.calc_props_fullmat.trans==1                                                                                                             % If transitivity should be calculated
                        out.fullmat_graph_meas.trans_pos(curr_sub,rep_lev) = transitivity_wu_sign(curr_conmat);                                                  % Each output should be 1 number
                    end
                    
                    if out.calc_props_fullmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                        out.fullmat_graph_meas.mod_deg_z_pos(curr_sub,:,rep_lev) = module_degree_zscore_sign(curr_conmat,out.mod_grps);                                                                 % Each output should be vector of size #ROIs
                    end
                    
                    prog = (curr_sub/out.num_subs)*(1-((rep_lev-1)/out.num_rep_levs));                                                                                                         % Calculate progress
                    if out.calcthrmat==1
                        progressbar(prog,[],[])                                                                                                                                                      % Update progress bar
                    else
                        progressbar(prog)                                                                                                                                                      % Update progress bar
                    end
                end
            end
        end
    end
    
    if out.calc_props_fullmat.rich_club==1                                                                                                            % If rich club networks were calculated
        %%%% Calculate maximum club size based on data if none provided
        if strcmp(out.weight_type,'Positive and Negative')
            if ~isempty(out.max_club_size)                                                                                                             % If max was provided
                out.max_club_size_full_pos = out.max_club_size;
                out.max_club_size_full_neg = out.max_club_size;
            else
                out.max_club_size_full_pos = 10000;                                                                                                                   % Set starting value as way higher than it could be
                out.max_club_size_full_neg = 10000;                                                                                                                   % Set starting value as way higher than it could be
                for rep_lev = 1:out.num_rep_levs                                                                                                           % Loop on repeated levels
                    for curr_sub = 1:out.num_subs                                                                                                        % Loop on participants
                        if length(out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev})<out.max_club_size_full_pos                                           % If this max is less than the current threshold
                            out.max_club_size_full_pos = length(out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev});                                         % Set new threshold
                        end
                        
                        if length(out.fullmat_graph_meas.rich_club_neg{curr_sub,rep_lev})<out.max_club_size_full_neg                                           % If this max is less than the current threshold
                            out.max_club_size_full_neg = length(out.fullmat_graph_meas.rich_club_neg{curr_sub,rep_lev});                                         % Set new threshold
                        end
                    end
                end
            end
            
            for rep_lev = 1:out.num_rep_levs                                                                                                           % Loop on repeated levels
                for curr_sub = out.num_subs:-1:1                                                                                                        % Loop on participants
                    temp_rcp(curr_sub,:,rep_lev) = out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev}(1:out.max_club_size_full_pos);
                    temp_rcn(curr_sub,:,rep_lev) = out.fullmat_graph_meas.rich_club_neg{curr_sub,rep_lev}(1:out.max_club_size_full_neg);
                end
            end
            
            out.fullmat_graph_meas.rich_club_pos = temp_rcp;
            out.fullmat_graph_meas.rich_club_neg = temp_rcn;
        else
            if ~isempty(out.max_club_size)                                                                                                             % If max was provided
                out.max_club_size_full_pos = out.max_club_size;
            else
                out.max_club_size_full_pos = 10000;                                                                                                                   % Set starting value as way higher than it could be
                for rep_lev = 1:out.num_rep_levs                                                                                                           % Loop on repeated levels
                    for curr_sub = 1:out.num_subs                                                                                                        % Loop on participants
                        if length(out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev})<out.max_club_size_full_pos                                           % If this max is less than the current threshold
                            out.max_club_size_full_pos = length(out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev});                                         % Set new threshold
                        end
                    end
                end
            end
            
            for rep_lev = 1:out.num_rep_levs                                                                                                           % Loop on repeated levels
                for curr_sub = out.num_subs:-1:1                                                                                                        % Loop on participants
                    temp_rcp(curr_sub,:,rep_lev) = out.fullmat_graph_meas.rich_club_pos{curr_sub,rep_lev}(1:out.max_club_size_full_pos);
                end
            end
            
            out.fullmat_graph_meas.rich_club_pos = temp_rcp;
        end
    end
    
    out.conmats_full_normed = out.conmats;
    out.conmats             = out.conmats_orig;
    
    fprintf('Done calculating properties for fully connected matrices!\n\n')                                                                                                             % Alert user
end




%%%%%%%%%%%%%%%%%%%%% Thresholded networks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if out.calcthrmat==1             % If any properties for thresholded networks should be calculated
    fprintf('Calculating minimum density ...\n')
    if isempty(out.denscalc_varmat)
        for rep_lev = out.num_rep_levs:-1:1                                          % For each repeated level
            full_mean_conmat_pos(:,:,rep_lev)            = create_mean_conmats(out.conmats(:,:,:,rep_lev)); % Calculate mean connectivity matrix
            full_mean_conmat_pos(full_mean_conmat_pos<0) = 0;
            full_density_pos(rep_lev)                    = find_min_graph_density(full_mean_conmat_pos(:,:,rep_lev)); % Find the minumum density at which the mean network remains connected
        end
        if any(isnan(full_density_pos(:)))
            fprintf('Warning: Minimum density could not be found for positive weights\n')
            out.min_dens_pos    = NaN;
            out.pos_mindens_nan = 1;
        else
            out.min_dens_pos    = ceil(100*max(full_density_pos(:)))/100;              % Note minimum density
            out.pos_mindens_nan = 0;
            if out.min_dens_pos>=out.max_dens_pos
                out.pos_mingreatermaxdens = 1;
                fprintf('Warning: Minimum density larger than specified maximum for positive weights\n')
            else
                out.pos_mingreatermaxdens = 0;
            end
        end
        
        if strcmp(out.weight_type,'Positive and Negative')
            for rep_lev = out.num_rep_levs:-1:1                                    % For each repeated level
                full_mean_conmat_neg(:,:,rep_lev)            = create_mean_conmats(-out.conmats(:,:,:,rep_lev)); % Calculate mean connectivity matrix
                full_mean_conmat_neg(full_mean_conmat_neg<0) = 0;
                full_density_neg(rep_lev)                    = find_min_graph_density(full_mean_conmat_neg(:,:,rep_lev)); % Find the minumum density at which the mean network remains connected
            end
            if any(isnan(full_density_neg(:)))
                fprintf('Warning: Minimum density could not be reached for negative weights\n')
                out.min_dens_neg    = NaN;
                out.neg_mindens_nan = 1;
            else
                out.min_dens_neg    = ceil(100*max(full_density_neg(:)))/100;          % Note minimum density
                out.neg_mindens_nan = 0;
                if out.min_dens_neg>=out.max_dens_neg
                    out.neg_mingreatermaxdens = 1;
                    fprintf('Warning: Minimum density larger than specified maximum for negative weights\n')
                else
                    out.neg_mingreatermaxdens = 0;
                end
            end
        end
    else
        if strcmp(out.partial_for_min_dens,'Yes')                                    % If the IVs should be partialed
            if rank([out.denscalc_varmat,out.denscalc_covarmat,ones(size(out.denscalc_varmat,1),1)])==(size([out.denscalc_varmat,out.denscalc_covarmat],2)+1) % If an intercept can be added without making the design matrix rank deficient
                out.denscalc_covarmat    = [out.denscalc_covarmat,ones(size(out.denscalc_varmat,1),1)];
                out.denscalc_covar_names = [out.denscalc_covar_names, 'intercept'];
            end
            min_dens_calc_vars = [];                                                 % Set empty matrix to start
            for var = size(out.denscalc_varmat,2):-1:1                               % For each IV (reverse indexed so min_dens_calc_vars is preallocated upon creation)
                des_part_ivs                    = [out.denscalc_covarmat,out.denscalc_varmat(:,setdiff((1:size(out.denscalc_varmat,2)),var))]; % Determine covariates to partial
                des_part_dv                     = out.denscalc_varmat(:,var);                            % Extract current IV from which to partial shared variance
                [~,~,min_dens_calc_vars(:,var)] = regress(des_part_dv,des_part_ivs); % Partial shared variance
            end
        else                                                                         % If IVs should not be partialed
            min_dens_calc_vars = out.denscalc_varmat;                                % Set matrix full variables
        end
                
        %%%% Find the minimum density. This is done by creating several
        %%%% mean networks, finding the minumum density at which each
        %%%% remains connected, taking the maximum value (which should be
        %%%% valid for all networks examined). Besides the overall mean 
        %%%% network, networks are created by stratifying each of the IVs
        %%%% selected above into groups (# og groups depending on the total
        %%%% N), and creating mean networks for each of these groups. This,
        %%%% stratification is done to avoid the mean density being invalid
        %%%% for a subset of participants (e.g., those who are high on some
        %%%% IV), which bias tests for that IV
        if out.num_subs<=20                           % If the total N is 15 or less
            min_dens_num_grps = 1;                      % Create 1 group for each IV
        elseif out.num_subs<=89                       % If the total N is between 20 & 89
            min_dens_num_grps = floor(out.num_subs/20); % Create 1 group for every ~20 participants
        else
            min_dens_num_grps = floor(out.num_subs/30); % Create 1 group for every 30 participants
        end
        
        if min_dens_num_grps==1
            for rep_lev = out.num_rep_levs:-1:1                                        % For each repeated level
                full_mean_conmat_pos(:,:,rep_lev)            = create_mean_conmats(out.conmats(:,:,:,rep_lev)); % Calculate mean connectivity matrix
                full_mean_conmat_pos(full_mean_conmat_pos<0) = 0;
                full_density_pos(rep_lev)                    = find_min_graph_density(full_mean_conmat_pos(:,:,rep_lev)); % Find the minumum density at which the mean network remains connected
            end
            if any(isnan(full_density_pos(:)))
                fprintf('Warning: Minimum density could not be reached for positive weights\n')
                out.min_dens_pos    = NaN;
                out.pos_mindens_nan = 1;
            else
                out.min_dens_pos    = ceil(100*max(full_density_pos(:)))/100;          % Note minimum density
                out.pos_mindens_nan = 0;
                if out.min_dens_pos>=out.max_dens_pos
                    out.pos_mingreatermaxdens = 1;
                    fprintf('Warning: Minimum density larger than specified maximum for positive weights\n')
                else
                    out.pos_mingreatermaxdens = 0;
                end
            end
            
            if strcmp(out.weight_type,'Positive and Negative')
                for rep_lev = out.num_rep_levs:-1:1                                    % For each repeated level
                    full_mean_conmat_neg(:,:,rep_lev)            = create_mean_conmats(-out.conmats(:,:,:,rep_lev)); % Calculate mean connectivity matrix
                    full_mean_conmat_neg(full_mean_conmat_neg<0) = 0;
                    full_density_neg(rep_lev)                    = find_min_graph_density(full_mean_conmat_neg(:,:,rep_lev)); % Find the minumum density at which the mean network remains connected
                end
                if any(isnan(full_density_neg(:)))
                    fprintf('Warning: Minimum density could not be reached for negative weights\n')
                    out.min_dens_neg    = NaN;
                    out.neg_mindens_nan = 1;
                else
                    out.min_dens_neg    = ceil(100*max(full_density_neg(:)))/100;          % Note minimum density
                    out.neg_mindens_nan = 0;
                    if out.min_dens_neg>=out.max_dens_neg
                        out.neg_mingreatermaxdens = 1;
                        fprintf('Warning: Minimum density larger than specified maximum for negative weights\n')
                    else
                        out.neg_mingreatermaxdens = 0;
                    end
                end
            end
        else
            for rep_lev = out.num_rep_levs:-1:1                                        % Loop on repeated levels
                for var = size(min_dens_calc_vars,2):-1:1                                % For each IV to use
                    [full_mean_conmat_pos(:,:,rep_lev),grp_mean_conmats_pos(:,:,:,rep_lev)] = create_mean_conmats(out.conmats(:,:,:,rep_lev),min_dens_calc_vars(:,var),1,'-num_grps',min_dens_num_grps); % Create an overall mean connectivity matrix and mean matrices for each subgrouping
                    full_mean_conmat_pos(full_mean_conmat_pos<0)                            = 0;
                    grp_mean_conmats_pos(grp_mean_conmats_pos<0)                            = 0;
                    full_density_pos(var,rep_lev)                                           = find_min_graph_density(full_mean_conmat_pos(:,:,rep_lev)); % Find minimum density for overall mean network
                    for grp = min_dens_num_grps:-1:1                                     % Loop on # of subgroups
                        grp_density_pos(var,grp,rep_lev) = find_min_graph_density(squeeze(grp_mean_conmats_pos(:,:,grp,rep_lev))); % Find minimum density for each subgroups network
                    end
                end
            end
            if any(isnan([full_density_pos(:),grp_density_pos(:)]))
                fprintf('Warning: Minimum density could not be reached for positive weights\n')
                out.min_dens_pos    = NaN;
                out.pos_mindens_nan = 1;
            else
                out.min_dens_pos    = ceil(100*max(max(full_density_pos(:)),max(grp_density_pos(:))))/100;                % Find maximum of all the minima
                out.pos_mindens_nan = 0;
                if out.min_dens_pos>=out.max_dens_pos
                    out.pos_mingreatermaxdens = 1;
                    fprintf('Warning: Minimum density larger than specified maximum for positive weights\n')
                else
                    out.pos_mingreatermaxdens = 0;
                end
            end
            
            if strcmp(out.weight_type,'Positive and Negative')
                [full_mean_conmat_neg(:,:,rep_lev),grp_mean_conmats_neg(:,:,:,rep_lev)] = create_mean_conmats(-out.conmats(:,:,:,rep_lev),min_dens_calc_vars(:,var),1,'-num_grps',min_dens_num_grps); % Create an overall mean connectivity matrix and mean matrices for each subgrouping
                full_mean_conmat_neg(full_mean_conmat_neg<0)                            = 0;
                grp_mean_conmats_neg(grp_mean_conmats_neg<0)                            = 0;
                for rep_lev = out.num_rep_levs:-1:1                                                              % Loop on repeated levels
                    for var = size(min_dens_calc_vars,2):-1:1                                                      % For each IV to use
                        full_density_neg(var,rep_lev) = find_min_graph_density(full_mean_conmat_neg(:,:,rep_lev)); % Find minimum density for overall mean network
                        for grp = min_dens_num_grps:-1:1                                                           % Loop on # of subgroups
                            grp_density_neg(var,grp,rep_lev) = find_min_graph_density(squeeze(grp_mean_conmats_neg(:,:,grp,rep_lev))); % Find minimum density for each subgroups network
                        end
                    end
                end
                if any(isnan([full_density_neg(:),grp_density_neg(:)]))
                    fprintf('Warning: Minimum density could not be reached for negative weights\n')
                    out.min_dens_neg    = NaN;
                    out.neg_mindens_nan = 1;
                else
                    out.min_dens_neg    = ceil(100*max(max(full_density_neg(:)),max(grp_density_neg(:))))/100;                % Find maximum of all the minima
                    out.neg_mindens_nan = 0;
                    if out.min_dens_neg>=out.max_dens_neg
                        out.neg_mingreatermaxdens = 1;
                        fprintf('Warning: Minimum density larger than specified maximum for negative weights\n')
                    else
                        out.neg_mingreatermaxdens = 0;
                    end
                end
            end
        end
    end
    
    for rep_lev = out.num_rep_levs:-1:1                                                                                                   % Loop on repeated levels
        max_dens_pos(rep_lev) = density_und(weight_conversion(threshold_absolute(full_mean_conmat_pos(:,:,rep_lev),0),'binarize'));
    end
    if out.max_dens_pos>min(max_dens_pos(:))
        out.max_dens_pos = floor(100*min(max_dens_pos(:)))/100;
        fprintf('Warning: Requested maximum density could not be reached for positive weights\n')
    end
    if mod((out.max_dens_pos-out.min_dens_pos),out.dens_step_pos)~=0
        out.dens_step_pos = (out.max_dens_pos-out.min_dens_pos)/floor((out.max_dens_pos-out.min_dens_pos)/out.dens_step_pos);                                                               % Find an appropriate step size close to that which they entered
    end
    out.dens_pos = out.min_dens_pos:out.dens_step_pos:out.max_dens_pos; % Calculate all the densities to use
    
    if strcmp(out.weight_type,'Positive and Negative') && out.neg_mindens_nan==0
        for rep_lev = out.num_rep_levs:-1:1                                                                                                   % Loop on repeated levels
            max_dens_neg(rep_lev) = density_und(weight_conversion(threshold_absolute(full_mean_conmat_neg(:,:,rep_lev),0),'binarize'));
        end
        if out.max_dens_neg>min(max_dens_neg(:))
            out.max_dens_neg = floor(100*min(max_dens_neg(:)))/100;
            fprintf('Warning: Requested maximum density could not be reached for negative weights\n')
        end
        if mod((out.max_dens_neg-out.min_dens_neg),out.dens_step_neg)~=0
            out.dens_step_neg = (out.max_dens_neg-out.min_dens_neg)/floor((out.max_dens_neg-out.min_dens_neg)/out.dens_step_neg);                                                               % Find an appropriate step size close to that which they entered
        end
        out.dens_neg = out.min_dens_neg:out.dens_step_neg:out.max_dens_neg; % Calculate all the densities to use
    end
    
    %%%% Threshold at different densities
    fprintf('Thresholding matrices ...\n')
    
    if out.pos_mindens_nan==0 && (out.max_dens_pos-out.min_dens_pos)>=out.dens_step_pos
        conmats_pos          = out.conmats.*(out.conmats>0);
        threshed_conmats_pos = zeros([size(out.conmats,1),size(out.conmats,2),size(out.conmats,3),length(out.dens_pos),size(out.conmats,4)]);                         % Preallocate
        for rep_lev = out.num_rep_levs:-1:1                                                                                                   % Loop on repeated levels
            for curr_dens = 1:length(out.dens_pos)                                                                                            % For each density level
                switch out.type_dens_thresh
                    case 'Use same thresh (densities will differ)'
                        temp_threshed                                 = threshold_proportional(full_mean_conmat_pos(:,:,rep_lev),out.dens_pos(curr_dens));
                        threshed_conmats_pos(:,:,:,curr_dens,rep_lev) = conmats_pos(:,:,:,rep_lev).*(conmats_pos(:,:,:,rep_lev)>=min(temp_threshed(temp_threshed~=0)));
                        for curr_sub = 1:out.num_subs
                            out.connected_nets_pos(curr_sub,curr_dens,rep_lev) = isempty(find(reachdist(weight_conversion(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev),'binarize'))==0,1));
                        end
                        clear temp_threshed
                    case 'Use different thresh (same density)'
                        for curr_sub = 1:out.num_subs
                            threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshold_proportional(conmats_pos(:,:,curr_sub,rep_lev),out.dens_pos(curr_dens));
                            out.connected_nets_pos(curr_sub,curr_dens,rep_lev)   = isempty(find(reachdist(weight_conversion(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev),'binarize'))==0,1));
                        end
                end
            end
        end
    end
    
    if strcmp(out.weight_type,'Positive and Negative') && out.neg_mindens_nan==0 && (out.max_dens_neg-out.min_dens_neg)>=out.dens_step_neg
        conmats_neg          = -out.conmats.*(out.conmats<0);
        threshed_conmats_neg = zeros([size(out.conmats,1),size(out.conmats,2),size(out.conmats,3),length(out.dens_neg),size(out.conmats,4)]);                         % Preallocate
        for rep_lev = out.num_rep_levs:-1:1                                                                                                   % Loop on repeated levels
            for curr_dens = 1:length(out.dens_neg)                                                                                            % For each density level
                switch out.type_dens_thresh
                    case 'Use same thresh (densities will differ)'
                        temp_threshed                                 = threshold_proportional(full_mean_conmat_neg(:,:,rep_lev),out.dens_neg(curr_dens));
                        threshed_conmats_neg(:,:,:,curr_dens,rep_lev) = conmats_neg(:,:,:,rep_lev).*(conmats_neg(:,:,:,rep_lev)>=min(temp_threshed(temp_threshed~=0)));
                        for curr_sub = 1:out.num_subs
                            out.connected_nets_neg(curr_sub,curr_dens,rep_lev) = isempty(find(reachdist(weight_conversion(threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev),'binarize'))==0,1));
                        end
                        clear temp_threshed
                    case 'Use different thresh (same density)'
                        for curr_sub = 1:out.num_subs
                            threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshold_proportional(conmats_neg(:,:,curr_sub,rep_lev),out.dens_neg(curr_dens));
                            out.connected_nets_neg(curr_sub,curr_dens,rep_lev)   = isempty(find(reachdist(weight_conversion(threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev),'binarize'))==0,1));
                        end
                end
            end
        end
    end
    
    if out.pos_mindens_nan==0 && (out.max_dens_pos-out.min_dens_pos)>=out.dens_step_pos
        if strcmp(out.weight_type,'Positive and Negative') && (out.neg_mindens_nan==1 || (out.max_dens_neg-out.min_dens_neg)<out.dens_step_neg)
            out.weight_type = 'Positive Only'; % If negative weight properties were requested, but not computable (i.e., no minimum density was found, not enough separation between min and max densities), only compute positive weights
        end
        
        if strcmp(out.weight_type,'Positive and Negative')
            switch out.type_weight_norm
                case 'Divide by Mean'
                    switch out.use_zeros
                        case 'Yes'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(mean(mean(threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev))));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(mean(threshed_conmats_pos(:,:,curr_sub,curr_dens,:))));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(mean(mean(mean(threshed_conmats_neg(:,:,curr_sub,curr_dens,:)))));
                                            end
                                        end
                                    end
                            end
                        case 'No'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(temppos(temppos>0)));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                tempneg                                              = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev);
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(mean(mean(tempneg(tempneg<0))));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,:);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(mean(temppos(temppos>0))));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                tempneg                                              = threshed_conmats_neg(:,:,curr_sub,curr_dens,:);
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(mean(mean(mean(tempneg(tempneg<0)))));
                                            end
                                        end
                                    end
                            end
                    end
                case 'Divide by Median'
                    switch out.use_zeros
                        case 'Yes'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(median(median(threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev))));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(median(threshed_conmats_pos(:,:,curr_sub,curr_dens,:))));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(median(median(median(threshed_conmats_neg(:,:,curr_sub,curr_dens,:)))));
                                            end
                                        end
                                    end
                            end
                        case 'No'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(temppos(temppos>0)));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                tempneg                                              = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev);
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(median(median(tempneg(tempneg<0))));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,:);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(median(temppos(temppos>0))));
                                            end
                                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                                tempneg                                              = threshed_conmats_neg(:,:,curr_sub,curr_dens,:);
                                                threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(median(median(median(tempneg(tempneg<0)))));
                                            end
                                        end
                                    end
                            end
                    end
                    if any(isnan(threshed_conmats_pos(:))) || any(isnan(threshed_conmats_neg(:)))
                        disp('WARNING: NaN''s were input in the connectivity matrices and/or normalizing by the median resulted in at least 1 matrix being replaced with NaNs in the thresholded matrices; investigate further')
                    end
                case 'Divide by Max'
                    switch out.repmeas_norm_type
                        case 'Normalize Each Level Individually'
                            for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                for curr_sub = 1:out.num_subs
                                    for curr_dens = 1:size(threshed_conmats_pos,4)
                                        threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/max(max(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)));
                                    end
                                    for curr_dens = 1:size(threshed_conmats_neg,4)
                                        threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(max(max(threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev))));
                                    end
                                end
                            end
                        case 'Normalize Across Levels'
                            for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                for curr_sub = 1:out.num_subs
                                    for curr_dens = 1:size(threshed_conmats_pos,4)
                                        threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/max(max(max(threshed_conmats_pos(:,:,curr_sub,curr_dens,:))));
                                    end
                                    for curr_dens = 1:size(threshed_conmats_neg,4)
                                        threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev)/abs(max(max(max(threshed_conmats_neg(:,:,curr_sub,curr_dens,:)))));
                                    end
                                end
                            end
                    end
            end
        else
            switch out.type_weight_norm
                case 'Divide by Mean'
                    switch out.use_zeros
                        case 'Yes'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(mean(threshed_conmats_pos(:,:,curr_sub,curr_dens,:))));
                                            end
                                        end
                                    end
                            end
                        case 'No'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(temppos(temppos>0)));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,:);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/mean(mean(mean(temppos(temppos>0))));
                                            end
                                        end
                                    end
                            end
                    end
                case 'Divide by Median'
                    switch out.use_zeros
                        case 'Yes'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(median(threshed_conmats_pos(:,:,curr_sub,curr_dens,:))));
                                            end
                                        end
                                    end
                            end
                        case 'No'
                            switch out.repmeas_norm_type
                                case 'Normalize Each Level Individually'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(temppos(temppos>0)));
                                            end
                                        end
                                    end
                                case 'Normalize Across Levels'
                                    for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                        for curr_sub = 1:out.num_subs
                                            for curr_dens = 1:size(threshed_conmats_pos,4)
                                                temppos                                              = threshed_conmats_pos(:,:,curr_sub,curr_dens,:);
                                                threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/median(median(median(temppos(temppos>0))));
                                            end
                                        end
                                    end
                            end
                    end
                    if any(isnan(threshed_conmats_pos(:)))
                        disp('WARNING: NaN''s were input in the connectivity matrices and/or normalizing by the median resulted in at least 1 matrix being replaced with NaNs in the thresholded matrices; investigate further')
                    end
                case 'Divide by Max'
                    switch out.repmeas_norm_type
                        case 'Normalize Each Level Individually'
                            for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                for curr_sub = 1:out.num_subs
                                    for curr_dens = 1:size(threshed_conmats_pos,4)
                                        threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/max(max(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)));
                                    end
                                end
                            end
                        case 'Normalize Across Levels'
                            for rep_lev = out.num_rep_levs:-1:1                                                                                                                                        % Loop through each repeated level
                                for curr_sub = 1:out.num_subs
                                    for curr_dens = 1:size(threshed_conmats_pos,4)
                                        threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev) = threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev)/max(max(max(threshed_conmats_pos(:,:,curr_sub,curr_dens,:))));
                                    end
                                end
                            end
                    end
            end
        end
        
        fprintf('Calculating properties for thresholded matrices ...\n')                                                                          % Let user know about progress
        
        %%%% Calculate network measures for each participant, for each
        %%%% threshold
        if use_parfor
            calc_props_thrmat = out.calc_props_thrmat;
            calcbinthresh     = out.calcbinthresh;
            num_subs          = out.num_subs;
            if isfield(out,'mod_grps')
                mod_grps = out.mod_grps;
            end
            max_club_size = out.max_club_size;
            
            if exist('threshed_conmats_pos','var')
                assort_pos             = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                assort_pos_bin         = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                bkg_tot_pos            = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                bkg_tot_pos_bin        = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                clust_coef_tot_pos     = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                clust_coef_tot_pos_bin = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                clust_coef_ZH_tot_pos  = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                cpl_pos                = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                cpl_pos_bin            = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                dens_pos               = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                glob_eff_pos           = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                glob_eff_pos_bin       = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                loc_eff_tot_pos        = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                loc_eff_tot_pos_bin    = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                rich_club_pos          = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                rich_club_pos_bin      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                small_world_pos        = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                trans_pos              = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                trans_pos_bin          = zeros(size(threshed_conmats_pos,4),out.num_subs,out.num_rep_levs);
                
                bkg_pos                = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                bkg_pos_bin            = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                close_cent_pos         = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                close_cent_pos_bin     = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                clust_coef_pos         = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                clust_coef_pos_bin     = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                clust_coef_ZH_pos      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                commn_cent_pos         = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                commn_cent_pos_bin     = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                deg_pos                = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                eigvec_cent_pos        = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                eigvec_cent_pos_bin    = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                gate_coef_pos          = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                gate_coef_pos_bin      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                kcore_cent_pos         = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                loc_assort_pos         = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                loc_eff_pos            = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                loc_eff_pos_bin        = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                mod_deg_z_pos          = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                mod_deg_z_pos_bin      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                node_bet_cent_pos      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                node_bet_cent_pos_bin  = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                pagerank_cent_pos      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                pagerank_cent_pos_bin  = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                part_coef_pos          = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                part_coef_pos_bin      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                subgraph_cent_pos      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.num_rep_levs);
                
                edge_bet_cent_pos      = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
                edge_bet_cent_pos_bin  = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
                match_pos              = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
                match_pos_bin          = zeros(size(threshed_conmats_pos,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
            end
            
            if exist('threshed_conmats_neg','var')
                assort_neg             = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                assort_neg_bin         = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                bkg_tot_neg            = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                bkg_tot_neg_bin        = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                clust_coef_tot_neg     = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                clust_coef_tot_neg_bin = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                clust_coef_ZH_tot_neg  = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                cpl_neg                = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                cpl_neg_bin            = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                dens_neg               = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                glob_eff_neg           = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                glob_eff_neg_bin       = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                loc_eff_tot_neg        = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                loc_eff_tot_neg_bin    = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                rich_club_neg          = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                rich_club_neg_bin      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                small_world_neg        = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                trans_neg              = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                trans_neg_bin          = zeros(size(threshed_conmats_neg,4),out.num_subs,out.num_rep_levs);
                
                bkg_neg                = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                bkg_neg_bin            = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                close_cent_neg         = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                close_cent_neg_bin     = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                clust_coef_neg         = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                clust_coef_neg_bin     = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                clust_coef_ZH_neg      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                commn_cent_neg         = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                commn_cent_neg_bin     = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                deg_neg                = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                eigvec_cent_neg        = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                eigvec_cent_neg_bin    = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                gate_coef_neg          = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                gate_coef_neg_bin      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                kcore_cent_neg         = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                loc_assort_neg         = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                loc_eff_neg            = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                loc_eff_neg_bin        = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                mod_deg_z_neg          = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                mod_deg_z_neg_bin      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                node_bet_cent_neg      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                node_bet_cent_neg_bin  = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                pagerank_cent_neg      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                pagerank_cent_neg_bin  = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                part_coef_neg          = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                part_coef_neg_bin      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                subgraph_cent_neg      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.num_rep_levs);
                
                edge_bet_cent_neg      = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
                edge_bet_cent_neg_bin  = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
                match_neg              = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
                match_neg_bin          = zeros(size(threshed_conmats_neg,4),out.num_subs,out.nROI,out.nROI,out.num_rep_levs);
            end
            
            for rep_lev = out.num_rep_levs:-1:1                                                                                                                % Loop on repeated levels
                for curr_dens = 1:size(threshed_conmats_pos,4)                                                                                                       % Loop on densities
                    parfor curr_sub = 1:num_subs                                                                                                                % Loop through each participant
                        curr_conmat = squeeze(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev));                                                                 % Extract connectivity matrix for current participant
                        
                        if calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                            assort_pos(curr_dens,curr_sub,rep_lev) = assortativity_wei(curr_conmat,0);                                             % Each output should be 1 number
                        end
                        
                        if calc_props_thrmat.bkg==1                                                                                                        % If the brokerage should be calculated
                            [bkg_pos(curr_dens,curr_sub,:,rep_lev),~,bkg_tot_pos(curr_dens,curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.cpl==1 || calc_props_thrmat.glob_eff==1                                                                               % If any properties requiring distance matrices should be calculated
                            dist_mat = distance_wei(weight_conversion(curr_conmat,'lengths'));                                                                   % Calculate distance matrix
                            if calc_props_thrmat.cpl==1 && calc_props_thrmat.glob_eff==1                                                                           % If both characteristic path length and global efficiency should be calculated
                                [cpl_pos(curr_dens,curr_sub,rep_lev),glob_eff_pos(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                            elseif calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                                [cpl_pos(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                        % Each output should be 1 number
                            else                                                                                                                                 % If only global efficiency should be calculated
                                [~,glob_eff_pos(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                 % Each output should be 1 number
                            end
                        end
                        
                        if calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                            close_cent_pos(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                            [clust_coef_pos(curr_dens,curr_sub,:,rep_lev),~,clust_coef_tot_pos(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,1);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.clust_coef_ZH==1                                                                                                        % If the clutering coefficient should be calculated
                            [clust_coef_ZH_pos(curr_dens,curr_sub,:,rep_lev),~,clust_coef_ZH_tot_pos(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                            commn_cent_pos(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(curr_conmat,mod_grps);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.deg==1 || calc_props_thrmat.dens==1 || calc_props_thrmat.kcore_cent==1 || calc_props_thrmat.sub_cent==1 || calc_props_thrmat.small_world==1 || calcbinthresh==1              % If any properties requiring binary matrices should be calculated
                            bin_conmat = weight_conversion(curr_conmat,'binarize');                                                                              % Calculate binary matrices
                            
                            if calc_props_thrmat.deg==1                                                                                                           % If degree should be calculated
                                deg_pos(curr_dens,curr_sub,:,rep_lev) = degrees_und(bin_conmat);                                                   % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.dens==1                                                                                                          % If density should be calculated
                                dens_pos(curr_dens,curr_sub,rep_lev) = density_und(bin_conmat);                                                    % Each output should be 1 number
                            end
                            
                            if calc_props_thrmat.kcore_cent==1                                                                                                          % If k-coreness centrality should be calculated
                                kcore_cent_pos(curr_dens,curr_sub,:,rep_lev) = kcoreness_centrality_bu(bin_conmat);                                                    % Each output should be 1 number
                            end
                            
                            if calc_props_thrmat.small_world==1                                                                                                   % If small-worldness should be calculated
                                small_world_pos(curr_dens,curr_sub,rep_lev) = HumphriesGurney_smallworldness_bu(bin_conmat);                       % Each output should be 1 number
                            end
                            
                            if calc_props_thrmat.sub_cent==1                                                                                                      % If subgraph centrality should be calculated
                                subgraph_cent_pos(curr_dens,curr_sub,:,rep_lev) = subgraph_centrality(bin_conmat);                                 % Each output should be vector of size #ROIs
                            end
                            
                            if calcbinthresh==1
                                if calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                                    assort_pos_bin(curr_dens,curr_sub,rep_lev) = assortativity_bin(bin_conmat,0);                                             % Each output should be 1 number
                                end
                                
                                if calc_props_thrmat.bkg==1                                                                                                        % If the brokerage should be calculated
                                    [bkg_pos_bin(curr_dens,curr_sub,:,rep_lev),bkg_tot_pos_bin(curr_dens,curr_sub,rep_lev)] = brokerage_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.cpl==1 || calc_props_thrmat.glob_eff==1                                                                               % If any properties requiring distance matrices should be calculated
                                    dist_mat = distance_bin(weight_conversion(bin_conmat,'lengths'));                                                                   % Calculate distance matrix
                                    if calc_props_thrmat.cpl==1 && calc_props_thrmat.glob_eff==1                                                                           % If both characteristic path length and global efficiency should be calculated
                                        [cpl_pos_bin(curr_dens,curr_sub,rep_lev),glob_eff_pos_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                                    elseif calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                                        [cpl_pos_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                        % Each output should be 1 number
                                    else                                                                                                                                 % If only global efficiency should be calculated
                                        [~,glob_eff_pos_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                 % Each output should be 1 number
                                    end
                                end
                                
                                if calc_props_thrmat.close_cent==1
                                    close_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(bin_conmat,2);                                        % Each output should be vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                                    clust_coef_pos_bin(curr_dens,curr_sub,:,rep_lev)   = clustering_coef_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                    clust_coef_tot_pos_bin(curr_dens,curr_sub,rep_lev) = mean(clust_coef_pos_bin(curr_dens,curr_sub,:,rep_lev));
                                end
                                
                                if calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                                    commn_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(bin_conmat,mod_grps);                                        % Each output should be vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.edge_bet_cent==1 || calc_props_thrmat.node_bet_cent==1                                                                % If any properties requiring length matrices should be calculated
                                    length_mat = weight_conversion(bin_conmat,'lengths');                                                                               % Calculate length matrix
                                    if calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                                        edge_bet_cent_pos_bin(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_bin(length_mat);                              % Each output should be square matrix of size #ROIs
                                    end
                                    if calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                                        node_bet_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = betweenness_bin(length_mat);                                     % Each output should be vector of size #ROIs
                                    end
                                end
                                
                                if calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                                    eigvec_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(bin_conmat);                               % Produces vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.gate_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                    gate_coef_pos_bin(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(bin_conmat,mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                                    loc_eff_pos_bin(curr_dens,curr_sub,:,rep_lev)   = efficiency_bin(bin_conmat,1);                                             % Each output should be vector of size #ROIs
                                    loc_eff_tot_pos_bin(curr_dens,curr_sub,rep_lev) = mean(loc_eff_pos_bin(curr_dens,curr_sub,:,rep_lev));
                                end
                                
                                if calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                                    match_pos_bin(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(bin_conmat);                                             % Each output should be square matrix of of size #ROIs
                                end
                                
                                if calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                                    pagerank_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(bin_conmat,0.85);                               % Produces vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                    part_coef_pos_bin(curr_dens,curr_sub,:,rep_lev) = participation_coef(bin_conmat,mod_grps); % Produces two outputs, each a vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                                    if ~isempty(max_club_size)                                                                                                  % If the user has specified a maximum density
                                        rich_club_pos_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat,max_club_size);                           % Each output should be vector of size max density
                                    else                                                                                                                                 % If not
                                        rich_club_pos_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat);                                             % Each output should be vector of size of max density
                                    end
                                end
                                
                                if calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                                    trans_pos_bin(curr_dens,curr_sub,rep_lev) = transitivity_bu(bin_conmat);                                                  % Each output should be 1 number
                                end
                                
                                if calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                                    mod_deg_z_pos_bin(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(bin_conmat,mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                                end
                            end
                        end
                        
                        if calc_props_thrmat.edge_bet_cent==1 || calc_props_thrmat.node_bet_cent==1                                                                % If any properties requiring length matrices should be calculated
                            length_mat = weight_conversion(curr_conmat,'lengths');                                                                               % Calculate length matrix
                            if calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                                edge_bet_cent_pos(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_wei(length_mat);                              % Each output should be square matrix of size #ROIs
                            end
                            if calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                                node_bet_cent_pos(curr_dens,curr_sub,:,rep_lev) = betweenness_wei(length_mat);                                     % Each output should be vector of size #ROIs
                            end
                        end
                        
                        if calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                            eigvec_cent_pos(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(curr_conmat);                               % Produces vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.gate_coef==1                                                                                                                                    % If the gateway coefficient should be calculated
                            gate_coef_pos(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(curr_conmat,mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.loc_assort==1                                                                                                           % If local efficiency should be calculated
                            loc_assort_pos(curr_dens,curr_sub,:,rep_lev) = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                            [loc_eff_pos(curr_dens,curr_sub,:,rep_lev),~,loc_eff_tot_pos(curr_dens,curr_sub,rep_lev),~] = efficiency_wei_sign(curr_conmat,1);                                             % Each output should be vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                            match_pos(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(curr_conmat);                                             % Each output should be square matrix of of size #ROIs
                        end
                        
                        if calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                            pagerank_cent_pos(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(curr_conmat,0.85);                               % Produces vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.part_coef==1                                                                                                                                    % If the participation coefficient should be calculated
                            part_coef_pos(curr_dens,curr_sub,:,rep_lev) = participation_coef(curr_conmat,mod_grps); % Produces two outputs, each a vector of size #ROIs
                        end
                        
                        if calc_props_thrmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                            if ~isempty(max_club_size)                                                                                                  % If the user has specified a maximum density
                                rich_club_pos{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat,max_club_size);                           % Each output should be vector of size max density
                            else                                                                                                                                 % If not
                                rich_club_pos{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat);                                             % Each output should be vector of size of max density
                            end
                        end
                        
                        if calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                            trans_pos(curr_dens,curr_sub,rep_lev) = transitivity_wu_sign(curr_conmat);                                                  % Each output should be 1 number
                        end
                        
                        if calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                            mod_deg_z_pos(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(curr_conmat,mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                        end
                    end
                end
                
                if calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                    thrmat_graph_meas.assort_pos = assort_pos;
                end
                if calc_props_thrmat.bkg==1                                                                                                        % If the clustering coefficient should be calculated
                    thrmat_graph_meas.bkg_pos     = bkg_pos;
                    thrmat_graph_meas.bkg_tot_pos = bkg_tot_pos;
                end
                if calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                    thrmat_graph_meas.cpl_pos = cpl_pos;
                end
                if calc_props_thrmat.glob_eff==1                                                                                                                        % If only global efficiency should be calculated
                    thrmat_graph_meas.glob_eff_pos = glob_eff_pos;
                end
                if calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                    thrmat_graph_meas.close_cent_pos = close_cent_pos;
                end
                if calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                    thrmat_graph_meas.clust_coef_pos     = clust_coef_pos;
                    thrmat_graph_meas.clust_coef_tot_pos = clust_coef_tot_pos;
                end
                if calc_props_thrmat.clust_coef_ZH==1                                                                                                        % If the clutering coefficient should be calculated
                    thrmat_graph_meas.clust_coef_ZH_pos     = clust_coef_ZH_pos;
                    thrmat_graph_meas.clust_coef_ZH_tot_pos = clust_coef_ZH_tot_pos;
                end
                if calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                    thrmat_graph_meas.commn_cent_pos = commn_cent_pos;
                end
                if calc_props_thrmat.deg==1                                                                                                           % If degree should be calculated
                    thrmat_graph_meas.deg_pos = deg_pos;
                end
                if calc_props_thrmat.dens==1                                                                                                          % If density should be calculated
                    thrmat_graph_meas.dens_pos = dens_pos;
                end
                if calc_props_thrmat.kcore_cent==1                                                                                                          % If k-coreness centrality should be calculated
                    thrmat_graph_meas.kcore_cent_pos = kcore_cent_pos;
                end
                if calc_props_thrmat.small_world==1                                                                                                   % If small-worldness should be calculated
                    thrmat_graph_meas.small_world_pos = small_world_pos;
                end
                if calc_props_thrmat.sub_cent==1                                                                                                      % If subgraph centrality should be calculated
                    thrmat_graph_meas.subgraph_cent_pos = subgraph_cent_pos;
                end
                if calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                    thrmat_graph_meas.edge_bet_cent_pos = edge_bet_cent_pos;
                end
                if calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                    thrmat_graph_meas.node_bet_cent_pos = node_bet_cent_pos;
                end
                if calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                    thrmat_graph_meas.eigvec_cent_pos = eigvec_cent_pos;
                end
                if calc_props_thrmat.gate_coef==1                                                                                                     % If gateway coefficient was calculated
                    thrmat_graph_meas.gate_coef_pos = gate_coef_pos;
                end
                if calc_props_thrmat.loc_assort==1                                                                                                           % If local efficiency should be calculated
                    thrmat_graph_meas.loc_assort_pos = loc_assort_pos;
                end
                if calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                    thrmat_graph_meas.loc_eff_pos     = loc_eff_pos;
                    thrmat_graph_meas.loc_eff_tot_pos = loc_eff_tot_pos;
                end
                if calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                    thrmat_graph_meas.match_pos = match_pos;
                end
                if calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                    thrmat_graph_meas.pagerank_cent_pos = pagerank_cent_pos;
                end
                if calc_props_thrmat.part_coef==1                                                                                                     % If pagerank centrality should be calculated
                    thrmat_graph_meas.part_coef_pos = part_coef_pos;
                end
                if calc_props_thrmat.rich_club==1
                    thrmat_graph_meas.rich_club_pos = rich_club_pos;
                end
                if calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                    thrmat_graph_meas.trans_pos = trans_pos;
                end
                if calc_props_thrmat.mod_deg_z==1                                                                                                             % If transitivity should be calculated
                    thrmat_graph_meas.mod_deg_z_pos = mod_deg_z_pos;
                end
                
                if calcbinthresh==1
                    if calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                        thrmat_graph_meas.assort_pos_bin = assort_pos_bin;
                    end
                    if calc_props_thrmat.bkg==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.bkg_pos_bin     = bkg_pos_bin;
                        thrmat_graph_meas.bkg_tot_pos_bin = bkg_tot_pos_bin;
                    end
                    if calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                        thrmat_graph_meas.cpl_pos_bin = cpl_pos_bin;
                    end
                    if calc_props_thrmat.glob_eff==1                                                                                                                        % If only global efficiency should be calculated
                        thrmat_graph_meas.glob_eff_pos_bin = glob_eff_pos_bin;
                    end
                    if calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.close_cent_pos_bin = close_cent_pos_bin;
                    end
                    if calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.clust_coef_pos_bin     = clust_coef_pos_bin;
                        thrmat_graph_meas.clust_coef_tot_pos_bin = clust_coef_tot_pos_bin;
                    end
                    if calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.commn_cent_pos_bin = commn_cent_pos_bin;
                    end
                    if calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                        thrmat_graph_meas.edge_bet_cent_pos_bin = edge_bet_cent_pos_bin;
                    end
                    if calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                        thrmat_graph_meas.node_bet_cent_pos_bin = node_bet_cent_pos_bin;
                    end
                    if calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                        thrmat_graph_meas.eigvec_cent_pos_bin = eigvec_cent_pos_bin;
                    end
                    if calc_props_thrmat.gate_coef==1                                                                                                     % If pagerank centrality should be calculated
                        thrmat_graph_meas.gate_coef_pos_bin = gate_coef_pos_bin;
                    end
                    if calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                        thrmat_graph_meas.loc_eff_pos_bin     = loc_eff_pos_bin;
                        thrmat_graph_meas.loc_eff_tot_pos_bin = loc_eff_tot_pos_bin;
                    end
                    if calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                        thrmat_graph_meas.match_pos_bin = match_pos_bin;
                    end
                    if calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                        thrmat_graph_meas.pagerank_cent_pos_bin = pagerank_cent_pos_bin;
                    end
                    if calc_props_thrmat.part_coef==1                                                                                                     % If pagerank centrality should be calculated
                        thrmat_graph_meas.part_coef_pos_bin = part_coef_pos_bin;
                    end
                    if calc_props_thrmat.rich_club==1
                        thrmat_graph_meas.rich_club_pos_bin = rich_club_pos_bin;
                    end
                    if calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                        thrmat_graph_meas.trans_pos_bin = trans_pos_bin;
                    end
                    if calc_props_thrmat.mod_deg_z==1                                                                                                             % If transitivity should be calculated
                        thrmat_graph_meas.mod_deg_z_pos_bin = mod_deg_z_pos_bin;
                    end
                end
                
                if exist('threshed_conmats_pos','var')
                    clear assort_pos assort_pos_bin bkg_tot_pos bkg_tot_pos_bin clust_coef_tot_pos clust_coef_tot_pos_bin clust_coef_ZH_tot_pos cpl_pos cpl_pos_bin dens_pos glob_eff_pos glob_eff_pos_bin kcore_cent_pos loc_eff_tot_pos loc_eff_tot_pos_bin rich_club_pos rich_club_pos_bin small_world_pos trans_pos trans_pos_bin bkg_pos bkg_pos_bin close_cent_pos close_cent_pos_bin clust_coef_pos clust_coef_pos_bin clust_coef_ZH_pos commn_cent_pos commn_cent_pos_bin deg_pos eigvec_cent_pos eigvec_cent_pos_bin gate_coef_pos gate_coef_pos_bin loc_assort_pos loc_eff_pos loc_eff_pos_bin mod_deg_z_pos mod_deg_z_pos_bin node_bet_cent_pos node_bet_cent_pos_bin pagerank_cent_pos pagerank_cent_pos_bin part_coef_pos part_coef_pos_bin subgraph_cent_pos edge_bet_cent_pos edge_bet_cent_pos_bin match_pos match_pos_bin
                end
                
                if strcmp(out.weight_type,'Positive and Negative')
                    for curr_dens = 1:size(threshed_conmats_neg,4)                                                                                        % Loop on densities
                        parfor curr_sub = 1:num_subs                                                                                                            % Loop through each participant
                            curr_conmat = squeeze(threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev));                                                  % Extract connectivity matrix for current participant
                            
                            if calc_props_thrmat.assort==1                                                                                                        % If assortativity should be calculated
                                assort_neg(curr_dens,curr_sub,rep_lev) = assortativity_wei(curr_conmat,0);                              % Each output should be 1 number
                            end
                            
                            if calc_props_thrmat.bkg==1                                                                                                        % If the brokerage should be calculated
                                [bkg_neg(curr_dens,curr_sub,:,rep_lev),~,bkg_tot_neg(curr_dens,curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);                                        % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.cpl==1 || calc_props_thrmat.glob_eff==1                                                                           % If any properties requiring distance matrices should be calculated
                                dist_mat = distance_wei(weight_conversion(curr_conmat,'lengths'));                                                               % Calculate distance matrix
                                if calc_props_thrmat.cpl==1 && calc_props_thrmat.glob_eff==1                                                                       % If both characteristic path length and global efficiency should be calculated
                                    [cpl_neg(curr_dens,curr_sub,rep_lev),glob_eff_neg(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                                elseif calc_props_thrmat.cpl==1                                                                                                   % If only characteristic path length should be calculated
                                    [cpl_neg(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                         % Each output should be 1 number
                                else                                                                                                                             % If only global efficiency should be calculated
                                    [~,glob_eff_neg(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                  % Each output should be 1 number
                                end
                            end
                            
                            if calc_props_thrmat.close_cent==1                                                                                                    % If the clutering coefficient should be calculated
                                close_cent_neg(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(curr_conmat,2);                         % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.clust_coef==1                                                                                                    % If the clutering coefficient should be calculated
                                [clust_coef_neg(curr_dens,curr_sub,:,rep_lev),~,clust_coef_tot_neg(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,1);                         % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.clust_coef_ZH==1                                                                                                    % If the clutering coefficient should be calculated
                                [clust_coef_ZH_neg(curr_dens,curr_sub,:,rep_lev),~,clust_coef_ZH_tot_neg(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,2);                         % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.commn_cent==1                                                                                                       
                                commn_cent_neg(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(curr_conmat,mod_grps);                                        % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.deg==1 || calc_props_thrmat.dens==1 || calc_props_thrmat.kcore_cent==1 || calc_props_thrmat.sub_cent==1 || calc_props_thrmat.small_world==1 || calcbinthresh==1            % If any properties requiring binary matrices should be calculated
                                bin_conmat = weight_conversion(curr_conmat,'binarize');                                                                          % Calculate binary matrices
                                
                                if calc_props_thrmat.deg==1                                                                                                       % If degree should be calculated
                                    deg_neg(curr_dens,curr_sub,:,rep_lev) = degrees_und(bin_conmat);                                    % Each output should be vector of size #ROIs
                                end
                                
                                if calc_props_thrmat.dens==1                                                                                                      % If density should be calculated
                                    dens_neg(curr_dens,curr_sub,rep_lev) = density_und(bin_conmat);                                     % Each output should be 1 number
                                end
                                
                                if calc_props_thrmat.kcore_cent==1                                                                                                          % If k-coreness centrality should be calculated
                                    kcore_cent_neg(curr_dens,curr_sub,:,rep_lev) = kcoreness_centrality_bu(bin_conmat);                                                    % Each output should be 1 number
                                end
                                
                                if calc_props_thrmat.small_world==1                                                                                               % If small-worldness should be calculated
                                    small_world_neg(curr_dens,curr_sub,rep_lev) = HumphriesGurney_smallworldness_bu(bin_conmat);        % Each output should be 1 number
                                end
                                
                                if calc_props_thrmat.sub_cent==1                                                                                                  % If subgraph centrality should be calculated
                                    subgraph_cent_neg(curr_dens,curr_sub,:,rep_lev) = subgraph_centrality(bin_conmat);                  % Each output should be vector of size #ROIs
                                end
                                
                                if calcbinthresh==1
                                    if calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                                        assort_neg_bin(curr_dens,curr_sub,rep_lev) = assortativity_bin(bin_conmat,0);                                             % Each output should be 1 number
                                    end
                                    
                                    if calc_props_thrmat.bkg==1                                                                                                        % If the brokerage should be calculated
                                        [bkg_neg_bin(curr_dens,curr_sub,:,rep_lev),bkg_tot_neg_bin(curr_dens,curr_sub,rep_lev)] = brokerage_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.cpl==1 || calc_props_thrmat.glob_eff==1                                                                               % If any properties requiring distance matrices should be calculated
                                        dist_mat = distance_bin(weight_conversion(bin_conmat,'lengths'));                                                                   % Calculate distance matrix
                                        if calc_props_thrmat.cpl==1 && calc_props_thrmat.glob_eff==1                                                                           % If both characteristic path length and global efficiency should be calculated
                                            [cpl_neg_bin(curr_dens,curr_sub,rep_lev),glob_eff_neg_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                                        elseif calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                                            [cpl_neg_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                        % Each output should be 1 number
                                        else                                                                                                                                 % If only global efficiency should be calculated
                                            [~,glob_eff_neg_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                 % Each output should be 1 number
                                        end
                                    end
                                    
                                    if calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                                        close_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(bin_conmat,2);                                        % Each output should be vector of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                                        clust_coef_neg_bin(curr_dens,curr_sub,:,rep_lev)   = clustering_coef_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                        clust_coef_tot_neg_bin(curr_dens,curr_sub,rep_lev) = mean(clust_coef_neg_bin(curr_dens,curr_sub,:,rep_lev));
                                    end
                                    
                                    if calc_props_thrmat.commn_cent==1
                                        commn_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(bin_conmat,mod_grps);                                        % Each output should be vector of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.edge_bet_cent==1 || calc_props_thrmat.node_bet_cent==1                                                                % If any properties requiring length matrices should be calculated
                                        length_mat = weight_conversion(bin_conmat,'lengths');                                                                               % Calculate length matrix
                                        if calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                                            edge_bet_cent_neg_bin(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_bin(length_mat);                              % Each output should be square matrix of size #ROIs
                                        end
                                        if calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                                            node_bet_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = betweenness_bin(length_mat);                                     % Each output should be vector of size #ROIs
                                        end
                                    end
                                    
                                    if calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                                        eigvec_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(bin_conmat);                               % Produces vector of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.gate_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                        gate_coef_neg_bin(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(bin_conmat,mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                                        loc_eff_neg_bin(curr_dens,curr_sub,:,rep_lev)   = efficiency_bin(bin_conmat,1);                                             % Each output should be vector of size #ROIs
                                        loc_eff_tot_neg_bin(curr_dens,curr_sub,rep_lev) = mean(loc_eff_neg_bin(curr_dens,curr_sub,:,rep_lev));
                                    end
                                    
                                    if calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                                        match_neg_bin(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(bin_conmat);                                             % Each output should be square matrix of of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                                        pagerank_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(bin_conmat,0.85);                               % Produces vector of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                        part_coef_neg_bin(curr_dens,curr_sub,:,rep_lev) = participation_coef(bin_conmat,mod_grps); % Produces two outputs, each a vector of size #ROIs
                                    end
                                    
                                    if calc_props_thrmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                                        if ~isempty(max_club_size)                                                                                                  % If the user has specified a maximum density
                                            rich_club_neg_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat,max_club_size);                           % Each output should be vector of size max density
                                        else                                                                                                                                 % If not
                                            rich_club_neg_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat);                                             % Each output should be vector of size of max density
                                        end
                                    end
                                    
                                    if calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                                        trans_neg_bin(curr_dens,curr_sub,rep_lev) = transitivity_bu(bin_conmat);                                                  % Each output should be 1 number
                                    end
                                    
                                    if calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                                        mod_deg_z_neg_bin(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(bin_conmat,mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                                    end
                                end
                            end
                            
                            if calc_props_thrmat.eigvec_cent==1                                                                                                   % If eigenvector centrality should be calculated
                                eigvec_cent_neg(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(curr_conmat);                % Produces vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.edge_bet_cent==1 || calc_props_thrmat.node_bet_cent==1                                                            % If any properties requiring length matrices should be calculated
                                length_mat = weight_conversion(curr_conmat,'lengths');                                                                           % Calculate length matrix
                                if calc_props_thrmat.edge_bet_cent==1                                                                                             % If edge betweenness centrality should be calculated
                                    edge_bet_cent_neg(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_wei(length_mat);               % Each output should be square matrix of size #ROIs
                                end
                                if calc_props_thrmat.node_bet_cent==1                                                                                             % If node betweenness centrality should be calculated
                                    node_bet_cent_neg(curr_dens,curr_sub,:,rep_lev) = betweenness_wei(length_mat);                      % Each output should be vector of size #ROIs
                                end
                            end
                            
                            if calc_props_thrmat.gate_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                gate_coef_neg(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(curr_conmat,mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.loc_assort==1                                                                                                           % If local efficiency should be calculated
                                loc_assort_neg(curr_dens,curr_sub,:,rep_lev) = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.loc_eff==1                                                                                                       % If local efficiency should be calculated
                                [loc_eff_neg(curr_dens,curr_sub,:,rep_lev),~,loc_eff_tot_neg(curr_dens,curr_sub,rep_lev),~] = efficiency_wei_sign(curr_conmat,1);                              % Each output should be vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.match==1                                                                                                         % If the matching index should be calculated
                                match_neg(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(curr_conmat);                              % Each output should be square matrix of of size #ROIs
                            end
                            
                            if calc_props_thrmat.pagerank_cent==1                                                                                                 % If pagerank centrality should be calculated
                                pagerank_cent_neg(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(curr_conmat,0.85);                % Produces vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                part_coef_neg(curr_dens,curr_sub,:,rep_lev) = participation_coef(curr_conmat,mod_grps); % Produces two outputs, each a vector of size #ROIs
                            end
                            
                            if calc_props_thrmat.rich_club==1                                                                                                     % If rich club networks should be calculated
                                if ~isempty(out.max_club_size)                                                                                              % If the user has specified a maximum density
                                    rich_club_neg{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat,max_club_size);            % Each output should be vector of size max density
                                else                                                                                                                             % If not
                                    rich_club_neg{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat);                              % Each output should be vector of size of max density
                                end
                            end
                            
                            if calc_props_thrmat.trans==1                                                                                                         % If transitivity should be calculated
                                trans_neg(curr_dens,curr_sub,rep_lev) = transitivity_wu_sign(curr_conmat);                                   % Each output should be 1 number
                            end
                            
                            if calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                                mod_deg_z_neg(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(curr_conmat,mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                            end
                        end
                    end
                    
                    if calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                        thrmat_graph_meas.assort_neg = assort_neg;
                    end
                    if calc_props_thrmat.bkg==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.bkg_neg     = bkg_neg;
                        thrmat_graph_meas.bkg_tot_neg = bkg_tot_neg;
                    end
                    if calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                        thrmat_graph_meas.cpl_neg = cpl_neg;
                    end
                    if calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.close_cent_neg = close_cent_neg;
                    end
                    if calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.clust_coef_neg     = clust_coef_neg;
                        thrmat_graph_meas.clust_coef_tot_neg = clust_coef_tot_neg;
                    end
                    if calc_props_thrmat.clust_coef_ZH==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.clust_coef_ZH_neg     = clust_coef_ZH_neg;
                        thrmat_graph_meas.clust_coef_ZH_tot_neg = clust_coef_ZH_tot_neg;
                    end
                    if calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                        thrmat_graph_meas.commn_cent_neg = commn_cent_neg;
                    end
                    if calc_props_thrmat.deg==1                                                                                                           % If degree should be calculated
                        thrmat_graph_meas.deg_neg = deg_neg;
                    end
                    if calc_props_thrmat.dens==1                                                                                                          % If density should be calculated
                        thrmat_graph_meas.dens_neg = dens_neg;
                    end
                    if calc_props_thrmat.glob_eff_neg==1                                                                                                                        % If only global efficiency should be calculated
                        thrmat_graph_meas.glob_eff_neg = glob_eff_neg;
                    end
                    if calc_props_thrmat.kcore_cent==1                                                                                                          % If k-coreness centrality should be calculated
                        thrmat_graph_meas.kcore_cent_neg = kcore_cent_neg;
                    end
                    if calc_props_thrmat.small_world==1                                                                                                   % If small-worldness should be calculated
                        thrmat_graph_meas.small_world_neg = small_world_neg;
                    end
                    if calc_props_thrmat.sub_cent==1                                                                                                      % If subgraph centrality should be calculated
                        thrmat_graph_meas.subgraph_cent_neg = subgraph_cent_neg;
                    end
                    if calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                        thrmat_graph_meas.edge_bet_cent_neg = edge_bet_cent_neg;
                    end
                    if calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                        thrmat_graph_meas.node_bet_cent_neg = node_bet_cent_neg;
                    end
                    if calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                        thrmat_graph_meas.eigvec_cent_neg = eigvec_cent_neg;
                    end
                    if calc_props_thrmat.gate_coef==1                                                                                                     % If pagerank centrality should be calculated
                        thrmat_graph_meas.gate_coef_neg = gate_coef_neg;
                    end
                    if calc_props_thrmat.loc_assort==1
                        thrmat_graph_meas.loc_assort_neg = loc_assort_neg;
                    end
                    if calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                        thrmat_graph_meas.loc_eff_neg     = loc_eff_neg;
                        thrmat_graph_meas.loc_eff_tot_neg = loc_eff_tot_neg;
                    end
                    if calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                        thrmat_graph_meas.match_neg = match_neg;
                    end
                    if calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                        thrmat_graph_meas.pagerank_cent_neg = pagerank_cent_neg;
                    end
                    if calc_props_thrmat.part_coef==1                                                                                                     % If pagerank centrality should be calculated
                        thrmat_graph_meas.part_coef_neg = part_coef_neg;
                    end
                    if calc_props_thrmat.rich_club==1
                        thrmat_graph_meas.rich_club_neg = rich_club_neg;
                    end
                    if calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                        thrmat_graph_meas.trans_neg = trans_neg;
                    end
                    if calc_props_thrmat.mod_deg_z==1                                                                                                             % If transitivity should be calculated
                        thrmat_graph_meas.mod_deg_z_neg = mod_deg_z_neg;
                    end
                    if calcbinthresh==1
                        if calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                            thrmat_graph_meas.assort_neg_bin = assort_neg_bin;
                        end
                        if calc_props_thrmat.bkg==1                                                                                                        % If the clutering coefficient should be calculated
                            thrmat_graph_meas.bkg_neg_bin     = bkg_neg_bin;
                            thrmat_graph_meas.bkg_tot_neg_bin = bkg_tot_neg_bin;
                        end
                        if calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                            thrmat_graph_meas.cpl_neg_bin = cpl_neg_bin;
                        end
                        if calc_props_thrmat.glob_eff==1                                                                                                                        % If only global efficiency should be calculated
                            thrmat_graph_meas.glob_eff_neg_bin = glob_eff_neg_bin;
                        end
                        if calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                            thrmat_graph_meas.close_cent_neg_bin = close_cent_neg_bin;
                        end
                        if calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                            thrmat_graph_meas.clust_coef_neg_bin     = clust_coef_neg_bin;
                            thrmat_graph_meas.clust_coef_tot_neg_bin = clust_coef_tot_neg_bin;
                        end
                        if calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                            thrmat_graph_meas.commn_cent_neg_bin = commn_cent_neg_bin;
                        end
                        if calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                            thrmat_graph_meas.edge_bet_cent_neg_bin = edge_bet_cent_neg_bin;
                        end
                        if calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                            thrmat_graph_meas.node_bet_cent_neg_bin = node_bet_cent_neg_bin;
                        end
                        if calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                            thrmat_graph_meas.eigvec_cent_neg_bin = eigvec_cent_neg_bin;
                        end
                        if calc_props_thrmat.gate_coef==1                                                                                                     % If pagerank centrality should be calculated
                            thrmat_graph_meas.gate_coef_neg_bin = gate_coef_neg_bin;
                        end
                        if calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                            thrmat_graph_meas.loc_eff_neg_bin     = loc_eff_neg_bin;
                            thrmat_graph_meas.loc_eff_tot_neg_bin = loc_eff_tot_neg_bin;
                        end
                        if calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                            thrmat_graph_meas.match_neg_bin = match_neg_bin;
                        end
                        if calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                            thrmat_graph_meas.pagerank_cent_neg_bin = pagerank_cent_neg_bin;
                        end
                        if calc_props_thrmat.part_coef==1                                                                                                     % If pagerank centrality should be calculated
                            thrmat_graph_meas.part_coef_neg_bin = part_coef_neg_bin;
                        end
                        if calc_props_thrmat.rich_club==1
                            thrmat_graph_meas.rich_club_neg_bin = rich_club_neg_bin;
                        end
                        if calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                            thrmat_graph_meas.trans_neg_bin = trans_neg_bin;
                        end
                        if calc_props_thrmat.mod_deg_z==1                                                                                                             % If transitivity should be calculated
                            thrmat_graph_meas.mod_deg_z_neg_bin = mod_deg_z_neg_bin;
                        end
                    end
                    
                    if exist('threshed_conmats_neg','var')
                        clear assort_neg assort_neg_bin bkg_tot_neg bkg_tot_neg_bin clust_coef_tot_neg clust_coef_tot_neg_bin clust_coef_ZH_tot_neg cpl_neg cpl_neg_bin dens_neg glob_eff_neg glob_eff_neg_bin kcore_cent_neg loc_eff_tot_neg loc_eff_tot_neg_bin rich_club_neg rich_club_neg_bin small_world_neg trans_neg trans_neg_bin bkg_neg bkg_neg_bin close_cent_neg close_cent_neg_bin clust_coef_neg clust_coef_neg_bin clust_coef_ZH_neg commn_cent_neg commn_cent_neg_bin deg_neg eigvec_cent_neg eigvec_cent_neg_bin gate_coef_neg gate_coef_neg_bin loc_assort_neg loc_eff_neg loc_eff_neg_bin mod_deg_z_neg mod_deg_z_neg_bin node_bet_cent_neg node_bet_cent_neg_bin pagerank_cent_neg pagerank_cent_neg_bin part_coef_neg part_coef_neg_bin subgraph_cent_neg edge_bet_cent_neg edge_bet_cent_neg_bin match_neg match_neg_bin
                    end
                end
            end
        else
            for rep_lev = out.num_rep_levs:-1:1                                                                                                                % Loop on repeated levels
                for curr_dens = 1:size(threshed_conmats_pos,4)                                                                                                       % Loop on densities
                    for curr_sub = 1:out.num_subs                                                                                                                % Loop through each participant
                        curr_conmat = squeeze(threshed_conmats_pos(:,:,curr_sub,curr_dens,rep_lev));                                                                 % Extract connectivity matrix for current participant
                        
                        if out.calc_props_thrmat.assort==1                                                                                                            % If assortativity should be calculated
                            thrmat_graph_meas.assort_pos(curr_dens,curr_sub,rep_lev) = assortativity_wei(curr_conmat,0);                                             % Each output should be 1 number
                        end
                        
                        if out.calc_props_thrmat.bkg == 1                                                                                                        % If the brokerage should be calculated
                            [thrmat_graph_meas.bkg_pos(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.bkg_tot_pos(curr_dens,curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.cpl==1 || out.calc_props_thrmat.glob_eff==1                                                                               % If any properties requiring distance matrices should be calculated
                            dist_mat = distance_wei(weight_conversion(curr_conmat,'lengths'));                                                                   % Calculate distance matrix
                            if out.calc_props_thrmat.cpl==1 && out.calc_props_thrmat.glob_eff==1                                                                           % If both characteristic path length and global efficiency should be calculated
                                [thrmat_graph_meas.cpl_pos(curr_dens,curr_sub,rep_lev),thrmat_graph_meas.glob_eff_pos(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                            elseif out.calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                                [thrmat_graph_meas.cpl_pos(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                        % Each output should be 1 number
                            else                                                                                                                                 % If only global efficiency should be calculated
                                [~,thrmat_graph_meas.glob_eff_pos(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                 % Each output should be 1 number
                            end
                        end
                        
                        if out.calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                            thrmat_graph_meas.close_cent_pos(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                            [thrmat_graph_meas.clust_coef_pos(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.clust_coef_tot_pos(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,1);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.clust_coef_ZH==1                                                                                                        % If the clutering coefficient should be calculated
                            [thrmat_graph_meas.clust_coef_ZH_pos(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.clust_coef_ZH_tot_pos(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,2);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                            thrmat_graph_meas.commn_cent_pos(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(curr_conmat,out.mod_grps);                                        % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.deg==1 || out.calc_props_thrmat.dens==1 || out.calc_props_thrmat.kcore_cent==1 || out.calc_props_thrmat.sub_cent==1 || out.calc_props_thrmat.small_world==1 || out.calcbinthresh==1               % If any properties requiring binary matrices should be calculated
                            bin_conmat = weight_conversion(curr_conmat,'binarize');                                                                              % Calcualte binary matrices
                            
                            if out.calc_props_thrmat.deg==1                                                                                                           % If degree should be calculated
                                thrmat_graph_meas.deg_pos(curr_dens,curr_sub,:,rep_lev) = degrees_und(bin_conmat);                                                   % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.dens==1                                                                                                          % If density should be calculated
                                thrmat_graph_meas.dens_pos(curr_dens,curr_sub,rep_lev) = density_und(bin_conmat);                                                    % Each output should be 1 number
                            end
                            
                            if out.calc_props_thrmat.kcore_cent==1                                                                                                          % If k-coreness centrality should be calculated
                                thrmat_graph_meas.kcore_cent_pos(curr_dens,curr_sub,:,rep_lev) = kcoreness_centrality_bu(bin_conmat);                                                    % Each output should be 1 number
                            end
                            
                            if out.calc_props_thrmat.small_world==1                                                                                                   % If small-worldness should be calculated
                                thrmat_graph_meas.small_world_pos(curr_dens,curr_sub,rep_lev) = HumphriesGurney_smallworldness_bu(bin_conmat);                       % Each output should be 1 number
                            end
                            
                            if out.calc_props_thrmat.sub_cent==1                                                                                                      % If subgraph centrality should be calculated
                                thrmat_graph_meas.subgraph_cent_pos(curr_dens,curr_sub,:,rep_lev) = subgraph_centrality(bin_conmat);                                 % Each output should be vector of size #ROIs
                            end
                            
                            if out.calcbinthresh==1
                                if out.calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                                    thrmat_graph_meas.assort_pos_bin(curr_dens,curr_sub,rep_lev) = assortativity_bin(bin_conmat,0);                                             % Each output should be 1 number
                                end
                                
                                if out.calc_props_thrmat.bkg == 1                                                                                                        % If the brokerage should be calculated
                                    [thrmat_graph_meas.bkg_pos_bin(curr_dens,curr_sub,:,rep_lev),thrmat_graph_meas.bkg_tot_pos_bin(curr_dens,curr_sub,rep_lev)] = brokerage_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.cpl==1 || out.calc_props_thrmat.glob_eff==1                                                                               % If any properties requiring distance matrices should be calculated
                                    dist_mat = distance_bin(weight_conversion(bin_conmat,'lengths'));                                                                   % Calculate distance matrix
                                    if out.calc_props_thrmat.cpl==1 && out.calc_props_thrmat.glob_eff==1                                                                           % If both characteristic path length and global efficiency should be calculated
                                        [thrmat_graph_meas.cpl_pos_bin(curr_dens,curr_sub,rep_lev),thrmat_graph_meas.glob_eff_pos_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                                    elseif out.calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                                        [thrmat_graph_meas.cpl_pos_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                        % Each output should be 1 number
                                    else                                                                                                                                 % If only global efficiency should be calculated
                                        [~,thrmat_graph_meas.glob_eff_pos_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                 % Each output should be 1 number
                                    end
                                end
                                
                                if out.calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                                    thrmat_graph_meas.close_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(bin_conmat,2);                                        % Each output should be vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                                    thrmat_graph_meas.clust_coef_pos_bin(curr_dens,curr_sub,:,rep_lev)   = clustering_coef_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                    thrmat_graph_meas.clust_coef_tot_pos_bin(curr_dens,curr_sub,rep_lev) = mean(thrmat_graph_meas.clust_coef_pos_bin(curr_dens,curr_sub,:,rep_lev));
                                end
                                
                                if out.calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                                    thrmat_graph_meas.commn_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(bin_conmat,out.mod_grps);                                        % Each output should be vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.edge_bet_cent==1 || out.calc_props_thrmat.node_bet_cent==1                                                                % If any properties requiring length matrices should be calculated
                                    length_mat = weight_conversion(bin_conmat,'lengths');                                                                               % Calculate length matrix
                                    if out.calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                                        thrmat_graph_meas.edge_bet_cent_pos_bin(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_bin(length_mat);                              % Each output should be square matrix of size #ROIs
                                    end
                                    if out.calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                                        thrmat_graph_meas.node_bet_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = betweenness_bin(length_mat);                                     % Each output should be vector of size #ROIs
                                    end
                                end
                                
                                if out.calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                                    thrmat_graph_meas.eigvec_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(bin_conmat);                               % Produces vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.gate_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                    thrmat_graph_meas.gate_coef_pos_bin(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(bin_conmat,out.mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                                    thrmat_graph_meas.loc_eff_pos_bin(curr_dens,curr_sub,:,rep_lev)   = efficiency_bin(bin_conmat,1);                                             % Each output should be vector of size #ROIs
                                    thrmat_graph_meas.loc_eff_tot_pos_bin(curr_dens,curr_sub,rep_lev) = mean(thrmat_graph_meas.loc_eff_pos_bin(curr_dens,curr_sub,:,rep_lev));
                                end
                                
                                if out.calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                                    thrmat_graph_meas.match_pos_bin(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(bin_conmat);                                             % Each output should be square matrix of of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                                    thrmat_graph_meas.pagerank_cent_pos_bin(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(bin_conmat,0.85);                               % Produces vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                    thrmat_graph_meas.part_coef_pos_bin(curr_dens,curr_sub,:,rep_lev) = participation_coef(bin_conmat,out.mod_grps); % Produces two outputs, each a vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                                    if ~isempty(out.max_club_size)                                                                                                  % If the user has specified a maximum density
                                        thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat,out.max_club_size);                           % Each output should be vector of size max density
                                    else                                                                                                                                 % If not
                                        thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat);                                             % Each output should be vector of size of max density
                                    end
                                end
                                
                                if out.calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                                    thrmat_graph_meas.trans_pos_bin(curr_dens,curr_sub,rep_lev) = transitivity_bu(bin_conmat);                                                  % Each output should be 1 number
                                end
                                
                                if out.calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                                    thrmat_graph_meas.mod_deg_z_pos_bin(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(bin_conmat,out.mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                                end
                            end
                        end
                        
                        if out.calc_props_thrmat.edge_bet_cent==1 || out.calc_props_thrmat.node_bet_cent==1                                                                % If any properties requiring length matrices should be calculated
                            length_mat = weight_conversion(curr_conmat,'lengths');                                                                               % Calculate length matrix
                            if out.calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                                thrmat_graph_meas.edge_bet_cent_pos(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_wei(length_mat);                              % Each output should be square matrix of size #ROIs
                            end
                            if out.calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                                thrmat_graph_meas.node_bet_cent_pos(curr_dens,curr_sub,:,rep_lev) = betweenness_wei(length_mat);                                     % Each output should be vector of size #ROIs
                            end
                        end
                        
                        if out.calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                            thrmat_graph_meas.eigvec_cent_pos(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(curr_conmat);                               % Produces vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.gate_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                            thrmat_graph_meas.gate_coef_pos(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(curr_conmat,out.mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.loc_assort==1                                                                                                           % If local efficiency should be calculated
                            thrmat_graph_meas.loc_assort_pos(curr_dens,curr_sub,:,rep_lev) = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                            [thrmat_graph_meas.loc_eff_pos(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.loc_eff_tot_pos(curr_dens,curr_sub,rep_lev),~] = efficiency_wei_sign(curr_conmat,1);                                             % Each output should be vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                            thrmat_graph_meas.match_pos(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(curr_conmat);                                             % Each output should be square matrix of of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                            thrmat_graph_meas.pagerank_cent_pos(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(curr_conmat,0.85);                               % Produces vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                            thrmat_graph_meas.part_coef_pos(curr_dens,curr_sub,:,rep_lev) = participation_coef(curr_conmat,out.mod_grps); % Produces two outputs, each a vector of size #ROIs
                        end
                        
                        if out.calc_props_thrmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                            if ~isempty(out.max_club_size)                                                                                                  % If the user has specified a maximum density
                                thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat,out.max_club_size);                           % Each output should be vector of size max density
                            else                                                                                                                                 % If not
                                thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat);                                             % Each output should be vector of size of max density
                            end
                        end
                        
                        if out.calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                            thrmat_graph_meas.trans_pos(curr_dens,curr_sub,rep_lev) = transitivity_wu_sign(curr_conmat);                                                  % Each output should be 1 number
                        end
                        
                        if out.calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                            thrmat_graph_meas.mod_deg_z_pos(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(curr_conmat,out.mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                        end
                        
                        prog = (((curr_sub/out.num_subs)/size(threshed_conmats_pos,4))/out.num_rep_levs)+(((curr_dens-1)/size(threshed_conmats_pos,4))/out.num_rep_levs)+(1-(rep_lev/out.num_rep_levs));                                % Calculate progress
                        if strcmp(out.weight_type,'Positive and Negative')
                            prog = prog/2;
                        end
                        if out.calcfullmat==1
                            progressbar([],prog,[])                                                                                                                  % Update user
                        else
                            progressbar(prog,[])                                                                                                                                                      % Update progress bar
                        end
                    end
                end
                
                if strcmp(out.weight_type,'Positive and Negative')
                    for curr_dens = 1:size(threshed_conmats_neg,4)                                                                                        % Loop on densities
                        for curr_sub = 1:out.num_subs                                                                                                            % Loop through each participant
                            curr_conmat = squeeze(threshed_conmats_neg(:,:,curr_sub,curr_dens,rep_lev));                                                  % Extract connectivity matrix for current participant
                            
                            if out.calc_props_thrmat.assort==1                                                                                                        % If assortativity should be calculated
                                thrmat_graph_meas.assort_neg(curr_dens,curr_sub,rep_lev) = assortativity_wei(curr_conmat,0);                              % Each output should be 1 number
                            end
                            
                            if out.calc_props_thrmat.bkg == 1                                                                                                        % If the brokerage should be calculated
                                [thrmat_graph_meas.bkg_neg(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.bkg_tot_neg(curr_dens,curr_sub,rep_lev)] = brokerage_wu_sign(curr_conmat);                                        % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.cpl==1 || out.calc_props_thrmat.glob_eff==1                                                                           % If any properties requiring distance matrices should be calculated
                                dist_mat = distance_wei(weight_conversion(curr_conmat,'lengths'));                                                               % Calculate distance matrix
                                if out.calc_props_thrmat.cpl==1 && out.calc_props_thrmat.glob_eff==1                                                                       % If both characteristic path length and global efficiency should be calculated
                                    [thrmat_graph_meas.cpl_neg(curr_dens,curr_sub,rep_lev),thrmat_graph_meas.glob_eff_neg(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                                elseif out.calc_props_thrmat.cpl==1                                                                                                   % If only characteristic path length should be calculated
                                    [thrmat_graph_meas.cpl_neg(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                         % Each output should be 1 number
                                else                                                                                                                             % If only global efficiency should be calculated
                                    [~,thrmat_graph_meas.glob_eff_neg(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                  % Each output should be 1 number
                                end
                            end
                            
                            if out.calc_props_thrmat.close_cent==1                                                                                                    % If the clutering coefficient should be calculated
                                thrmat_graph_meas.close_cent_neg(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(curr_conmat,2);                         % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.clust_coef==1                                                                                                    % If the clutering coefficient should be calculated
                                [thrmat_graph_meas.clust_coef_neg(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.clust_coef_tot_neg(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,1);                         % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.clust_coef_ZH==1                                                                                                    % If the clutering coefficient should be calculated
                                [thrmat_graph_meas.clust_coef_ZH_neg(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.clust_coef_ZH_tot_neg(curr_dens,curr_sub,rep_lev),~] = clustering_coef_wu_sign(curr_conmat,2);                         % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                                thrmat_graph_meas.commn_cent_neg(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(curr_conmat,out.mod_grps);                                        % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.deg==1 || out.calc_props_thrmat.dens==1 || out.calc_props_thrmat.kcore_cent==1 || out.calc_props_thrmat.sub_cent==1 || out.calc_props_thrmat.small_world==1 || out.calcbinthresh==1           % If any properties requiring binary matrices should be calculated
                                bin_conmat = weight_conversion(curr_conmat,'binarize');                                                                          % Calcualte binary matrices
                                
                                if out.calc_props_thrmat.deg==1                                                                                                       % If degree should be calculated
                                    thrmat_graph_meas.deg_neg(curr_dens,curr_sub,:,rep_lev) = degrees_und(bin_conmat);                                    % Each output should be vector of size #ROIs
                                end
                                
                                if out.calc_props_thrmat.dens==1                                                                                                      % If density should be calculated
                                    thrmat_graph_meas.dens_neg(curr_dens,curr_sub,rep_lev) = density_und(bin_conmat);                                     % Each output should be 1 number
                                end
                                
                                if out.calc_props_thrmat.kcore_cent==1                                                                                                          % If k-coreness centrality should be calculated
                                    thrmat_graph_meas.kcore_cent_neg(curr_dens,curr_sub,:,rep_lev) = kcoreness_centrality_bu(bin_conmat);                                                    % Each output should be 1 number
                                end
                                
                                if out.calc_props_thrmat.small_world==1                                                                                               % If small-worldness should be calculated
                                    thrmat_graph_meas.small_world_neg(curr_dens,curr_sub,rep_lev) = HumphriesGurney_smallworldness_bu(bin_conmat);        % Each output should be 1 number
                                end
                                
                                if out.calc_props_thrmat.sub_cent==1                                                                                                  % If subgraph centrality should be calculated
                                    thrmat_graph_meas.subgraph_cent_neg(curr_dens,curr_sub,:,rep_lev) = subgraph_centrality(bin_conmat);                  % Each output should be vector of size #ROIs
                                end
                                
                                if out.calcbinthresh==1
                                    if out.calc_props_thrmat.assort==1                                                                                                            %#ok<*PFBNS> % If assortativity should be calculated
                                        thrmat_graph_meas.assort_neg_bin(curr_dens,curr_sub,rep_lev) = assortativity_bin(bin_conmat,0);                                             % Each output should be 1 number
                                    end
                                    
                                    if out.calc_props_thrmat.bkg == 1                                                                                                        % If the brokerage should be calculated
                                        [thrmat_graph_meas.bkg_neg_bin(curr_dens,curr_sub,:,rep_lev),thrmat_graph_meas.bkg_tot_neg_bin(curr_dens,curr_sub,rep_lev)] = brokerage_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.cpl==1 || out.calc_props_thrmat.glob_eff==1                                                                               % If any properties requiring distance matrices should be calculated
                                        dist_mat = distance_bin(weight_conversion(bin_conmat,'lengths'));                                                                   % Calculate distance matrix
                                        if out.calc_props_thrmat.cpl==1 && out.calc_props_thrmat.glob_eff==1                                                                           % If both characteristic path length and global efficiency should be calculated
                                            [thrmat_graph_meas.cpl_neg_bin(curr_dens,curr_sub,rep_lev),thrmat_graph_meas.glob_eff_neg_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat); % Each output should be 1 number
                                        elseif out.calc_props_thrmat.cpl==1                                                                                                       % If only characteristic path length should be calculated
                                            [thrmat_graph_meas.cpl_neg_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                        % Each output should be 1 number
                                        else                                                                                                                                 % If only global efficiency should be calculated
                                            [~,thrmat_graph_meas.glob_eff_neg_bin(curr_dens,curr_sub,rep_lev)] = charpath(dist_mat);                                                 % Each output should be 1 number
                                        end
                                    end
                                    
                                    if out.calc_props_thrmat.close_cent==1                                                                                                        % If the clutering coefficient should be calculated
                                        thrmat_graph_meas.close_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = closeness_cent_wu_sign(bin_conmat,2);                                        % Each output should be vector of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.clust_coef==1                                                                                                        % If the clutering coefficient should be calculated
                                        thrmat_graph_meas.clust_coef_neg_bin(curr_dens,curr_sub,:,rep_lev)   = clustering_coef_bu(bin_conmat);                                        % Each output should be vector of size #ROIs
                                        thrmat_graph_meas.clust_coef_tot_neg_bin(curr_dens,curr_sub,rep_lev) = mean(thrmat_graph_meas.clust_coef_neg_bin(curr_dens,curr_sub,:,rep_lev));
                                    end
                                    
                                    if out.calc_props_thrmat.commn_cent==1                                                                                                        % If the clutering coefficient should be calculated
                                        thrmat_graph_meas.commn_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = commn_cent_wu(bin_conmat,out.mod_grps);                                        % Each output should be vector of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.edge_bet_cent==1 || out.calc_props_thrmat.node_bet_cent==1                                                                % If any properties requiring length matrices should be calculated
                                        length_mat = weight_conversion(bin_conmat,'lengths');                                                                               % Calculate length matrix
                                        if out.calc_props_thrmat.edge_bet_cent==1                                                                                                 % If edge betweenness centrality should be calculated
                                            thrmat_graph_meas.edge_bet_cent_neg_bin(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_bin(length_mat);                              % Each output should be square matrix of size #ROIs
                                        end
                                        if out.calc_props_thrmat.node_bet_cent==1                                                                                                 % If node betweenness centrality should be calculated
                                            thrmat_graph_meas.node_bet_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = betweenness_bin(length_mat);                                     % Each output should be vector of size #ROIs
                                        end
                                    end
                                    
                                    if out.calc_props_thrmat.eigvec_cent==1                                                                                                       % If eigenvector centrality should be calculated
                                        thrmat_graph_meas.eigvec_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(bin_conmat);                               % Produces vector of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.gate_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                        thrmat_graph_meas.gate_coef_neg_bin(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(bin_conmat,out.mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.loc_eff==1                                                                                                           % If local efficiency should be calculated
                                        thrmat_graph_meas.loc_eff_neg_bin(curr_dens,curr_sub,:,rep_lev)   = efficiency_bin(bin_conmat,1);                                             % Each output should be vector of size #ROIs
                                        thrmat_graph_meas.loc_eff_tot_neg_bin(curr_dens,curr_sub,rep_lev) = mean(thrmat_graph_meas.loc_eff_neg_bin(curr_dens,curr_sub,:,rep_lev));
                                    end
                                    
                                    if out.calc_props_thrmat.match==1                                                                                                             % If the matching index should be calculated
                                        thrmat_graph_meas.match_neg_bin(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(bin_conmat);                                             % Each output should be square matrix of of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.pagerank_cent==1                                                                                                     % If pagerank centrality should be calculated
                                        thrmat_graph_meas.pagerank_cent_neg_bin(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(bin_conmat,0.85);                               % Produces vector of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                        thrmat_graph_meas.part_coef_neg_bin(curr_dens,curr_sub,:,rep_lev) = participation_coef(bin_conmat,out.mod_grps); % Produces two outputs, each a vector of size #ROIs
                                    end
                                    
                                    if out.calc_props_thrmat.rich_club==1                                                                                                         % If rich club networks should be calculated
                                        if ~isempty(out.max_club_size)                                                                                                  % If the user has specified a maximum density
                                            thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat,out.max_club_size);                           % Each output should be vector of size max density
                                        else                                                                                                                                 % If not
                                            thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev} = rich_club_bu(bin_conmat);                                             % Each output should be vector of size of max density
                                        end
                                    end
                                    
                                    if out.calc_props_thrmat.trans==1                                                                                                             % If transitivity should be calculated
                                        thrmat_graph_meas.trans_neg_bin(curr_dens,curr_sub,rep_lev) = transitivity_bu(bin_conmat);                                                  % Each output should be 1 number
                                    end
                                    
                                    if out.calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                                        thrmat_graph_meas.mod_deg_z_neg_bin(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(bin_conmat,out.mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                                    end
                                end
                            end
                            
                            if out.calc_props_thrmat.eigvec_cent==1                                                                                                   % If eigenvector centrality should be calculated
                                thrmat_graph_meas.eigvec_cent_neg(curr_dens,curr_sub,:,rep_lev) = eigenvector_centrality_und(curr_conmat);                % Produces vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.edge_bet_cent==1 || out.calc_props_thrmat.node_bet_cent==1                                                            % If any properties requiring length matrices should be calculated
                                length_mat = weight_conversion(curr_conmat,'lengths');                                                                           % Calculate length matrix
                                if out.calc_props_thrmat.edge_bet_cent==1                                                                                             % If edge betweenness centrality should be calculated
                                    thrmat_graph_meas.edge_bet_cent_neg(curr_dens,curr_sub,:,:,rep_lev) = edge_betweenness_wei(length_mat);               % Each output should be square matrix of size #ROIs
                                end
                                if out.calc_props_thrmat.node_bet_cent==1                                                                                             % If node betweenness centrality should be calculated
                                    thrmat_graph_meas.node_bet_cent_neg(curr_dens,curr_sub,:,rep_lev) = betweenness_wei(length_mat);                      % Each output should be vector of size #ROIs
                                end
                            end
                            
                            if out.calc_props_thrmat.gate_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                thrmat_graph_meas.gate_coef_neg(curr_dens,curr_sub,:,rep_lev) = gateway_coef_sign(curr_conmat,out.mod_grps,1); % Produces two outputs, each a vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.loc_assort==1                                                                                                           % If local efficiency should be calculated
                                thrmat_graph_meas.loc_assort_neg(curr_dens,curr_sub,:,rep_lev) = local_assortativity_wu_sign(curr_conmat);                                             % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.loc_eff==1                                                                                                       % If local efficiency should be calculated
                                [thrmat_graph_meas.loc_eff_neg(curr_dens,curr_sub,:,rep_lev),~,thrmat_graph_meas.loc_eff_tot_neg(curr_dens,curr_sub,rep_lev),~] = efficiency_wei_sign(curr_conmat,1);                              % Each output should be vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.match==1                                                                                                         % If the matching index should be calculated
                                thrmat_graph_meas.match_neg(curr_dens,curr_sub,:,:,rep_lev) = matching_ind_und(curr_conmat);                              % Each output should be square matrix of of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.pagerank_cent==1                                                                                                 % If pagerank centrality should be calculated
                                thrmat_graph_meas.pagerank_cent_neg(curr_dens,curr_sub,:,rep_lev) = pagerank_centrality_sign(curr_conmat,0.85);                % Produces vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.part_coef==1                                                                                                                                     % If the participation coefficient should be calculated
                                thrmat_graph_meas.part_coef_neg(curr_dens,curr_sub,:,rep_lev) = participation_coef(curr_conmat,out.mod_grps); % Produces two outputs, each a vector of size #ROIs
                            end
                            
                            if out.calc_props_thrmat.rich_club==1                                                                                                     % If rich club networks should be calculated
                                if ~isempty(out.max_club_size)                                                                                              % If the user has specified a maximum density
                                    thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat,out.max_club_size);            % Each output should be vector of size max density
                                else                                                                                                                             % If not
                                    thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev} = rich_club_wu_sign(curr_conmat);                              % Each output should be vector of size of max density
                                end
                            end
                            
                            if out.calc_props_thrmat.trans==1                                                                                                         % If transitivity should be calculated
                                thrmat_graph_meas.trans_neg(curr_dens,curr_sub,rep_lev) = transitivity_wu_sign(curr_conmat);                                   % Each output should be 1 number
                            end
                            
                            if out.calc_props_thrmat.mod_deg_z==1                                                                                                                                     % If the within-module degree z-score should be calculated
                                thrmat_graph_meas.mod_deg_z_neg(curr_dens,curr_sub,:,rep_lev) = module_degree_zscore(curr_conmat,out.mod_grps,0);                                                                 % Each output should be vector of size #ROIs
                            end
                            
                            prog = (((((curr_sub/out.num_subs)/size(threshed_conmats_neg,4))/out.num_rep_levs)+(((curr_dens-1)/size(threshed_conmats_neg,4))/out.num_rep_levs)+(1-(rep_lev/out.num_rep_levs)))*0.5)+0.5;                                % Calculate progress
                            if out.calcfullmat==1
                                progressbar([],prog,[])                                                                                                                  % Update user
                            else
                                progressbar(prog,[])                                                                                                                                                      % Update progress bar
                            end
                        end
                    end
                end
            end
        end
        
        if out.calc_props_thrmat.rich_club==1                                                                                                              % If rich club networks were calculated
            %%%% Calculate maximum club size based on data if none provided
            if ~isempty(out.max_club_size)                                                                                                             % If no max was provided
                out.max_club_size_thr_pos = out.max_club_size;
                out.max_club_size_thr_neg = out.max_club_size;
                if out.calcbinthresh==1
                    out.max_club_size_thr_pos_bin = out.max_club_size;
                    out.max_club_size_thr_neg_bin = out.max_club_size;
                end
            else
                out.max_club_size_thr_pos = 10000;                                                                                                                   % Set starting value as way higher than it could be
                for rep_lev = 1:out.num_rep_levs                                                                                                           % Loop on repeated levels
                    for curr_dens = 1:size(threshed_conmats_pos,4)                                                                                               % Loop on densities
                        for curr_sub = 1:out.num_subs                                                                                                        % Loop on participants
                            if length(thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_pos                                           % If this max is less than the current threshold
                                out.max_club_size_thr_pos = length(thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev});                                         % Set new threshold
                            end
                        end
                    end
                end
                
                if out.calcbinthresh==1
                    out.max_club_size_thr_pos_bin = 10000;                                                                                                                   % Set starting value as way higher than it could be
                    for rep_lev = 1:out.num_rep_levs                                                                                                           % Loop on repeated levels
                        for curr_dens = 1:size(threshed_conmats_pos,4)                                                                                               % Loop on densities
                            for curr_sub = 1:out.num_subs                                                                                                        % Loop on participants
                                if length(thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_pos_bin                                   % If this max is less than the current threshold
                                    out.max_club_size_thr_pos_bin = length(thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev});                                         % Set new threshold
                                end
                            end
                        end
                    end
                end
                
                if strcmp(out.weight_type,'Positive and Negative')
                    out.max_club_size_thr_neg = 10000;                                                                                                    % Set starting value as way higher than it could be
                    for rep_lev = 1:out.num_rep_levs                                                                                                       % Loop on repeated levels
                        for curr_dens = 1:size(threshed_conmats_neg,4)                                                                                % Loop on densities
                            for curr_sub = 1:out.num_subs                                                                                                    % Loop on participants
                                if length(thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_neg                 % If this max is less than the current threshold
                                    out.max_club_size_thr_neg = length(thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev});               % Set new threshold
                                end
                            end
                        end
                    end
                    
                    if out.calcbinthresh==1
                        out.max_club_size_thr_neg_bin = 10000;                                                                                                    % Set starting value as way higher than it could be
                        for rep_lev = 1:out.num_rep_levs                                                                                                       % Loop on repeated levels
                            for curr_dens = 1:size(threshed_conmats_neg,4)                                                                                % Loop on densities
                                for curr_sub = 1:out.num_subs                                                                                                    % Loop on participants
                                    if length(thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_neg_bin                 % If this max is less than the current threshold
                                        out.max_club_size_thr_neg_bin = length(thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev});               % Set new threshold
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('Done calculating properties for thresholded matrices!\n\n')                                                                                 % Alert user
        
        
        %%%% Calculate area under the curve for each measure, for each
        %%%% participant (ignoring NaNs)
        fprintf('Calculating AUC for threshmat properties ...\n')
        if use_parfor
            progressbar('Progress For Calculating AUC for Thresholded Matrices')                                         % Initialize progress bar at zero
        end
        
        for rep_lev = out.num_rep_levs:-1:1                                                                                                                                                                 % Loop on repeated levels
            for curr_sub = 1:out.num_subs                                                                                                                                                                     % Loop on participant
                % Assortativity
                if out.calc_props_thrmat.assort==1                                                                                                                                                                 % If assortativity was calculated
                    accept_vals                                                       = isfinite(thrmat_graph_meas.assort_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.assort_pos(:,curr_sub,rep_lev))==0);                                                           % Determine which values are finite and real
                    out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);                                                                                                         % Determine how many acceptable values there are
                    if out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC(curr_sub,rep_lev)>1                                                                                                                      % If there is more than 1 acceptable value
                        out.AUC_thrmat_graph_meas.assort_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                    elseif out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC(curr_sub,rep_lev)==1                                                                                                                 % If there is only one acceptable value
                        out.AUC_thrmat_graph_meas.assort_pos(curr_sub,rep_lev) = thrmat_graph_meas.assort_pos(accept_vals,curr_sub,rep_lev);                                                                          % Set the AUC as that value
                    else                                                                                                                                                                                      % If no acceptable values are found
                        out.AUC_thrmat_graph_meas.assort_pos(curr_sub,rep_lev) = NaN;                                                                                                                             % Set the AUC as NaN
                    end
                    
                    if out.calcbinthresh==1
                        accept_vals_bin                                                       = isfinite(thrmat_graph_meas.assort_pos_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.assort_pos_bin(:,curr_sub,rep_lev))==0);                                                   % Determine which values are finite and real
                        out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);                                                                                                         % Determine how many acceptable values there are
                        if out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC(curr_sub,rep_lev)>1                                                                                                                      % If there is more than 1 acceptable value
                            out.AUC_thrmat_graph_meas.assort_pos_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                        elseif out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC(curr_sub,rep_lev)==1                                                                                                                 % If there is only one acceptable value
                            out.AUC_thrmat_graph_meas.assort_pos_bin(curr_sub,rep_lev) = thrmat_graph_meas.assort_pos_bin(accept_vals_bin,curr_sub,rep_lev);                                                                  % Set the AUC as that value
                        else                                                                                                                                                                                              % If no acceptable values are found
                            out.AUC_thrmat_graph_meas.assort_pos_bin(curr_sub,rep_lev) = NaN;                                                                                                                             % Set the AUC as NaN
                        end
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                                = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);                                                                                                         % Determine how many acceptable values there are
                        if out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1                                                                                                                      % If there is more than 1 acceptable value
                            out.AUC_thrmat_graph_meas.assort_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                        elseif out.AUC_thrmat_graph_meas.assort_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1                                                                                                                 % If there is only one acceptable value
                            out.AUC_thrmat_graph_meas.assort_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.assort_pos(accept_vals,curr_sub,rep_lev);                                                                      % Set the AUC as that value
                        else                                                                                                                                                                                                   % If no acceptable values are found
                            out.AUC_thrmat_graph_meas.assort_pos_nodiscon(curr_sub,rep_lev) = NaN;                                                                                                                             % Set the AUC as NaN
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                                = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);                                                                                                         % Determine how many acceptable values there are
                            if out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1                                                                                                                      % If there is more than 1 acceptable value
                                out.AUC_thrmat_graph_meas.assort_pos_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                            elseif out.AUC_thrmat_graph_meas.assort_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1                                                                                                                 % If there is only one acceptable value
                                out.AUC_thrmat_graph_meas.assort_pos_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.assort_pos_bin(accept_vals_bin,curr_sub,rep_lev);                                                                  % Set the AUC as that value
                            else                                                                                                                                                                                                       % If no acceptable values are found
                                out.AUC_thrmat_graph_meas.assort_pos_bin_nodiscon(curr_sub,rep_lev) = NaN;                                                                                                                             % Set the AUC as NaN
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                       = isfinite(thrmat_graph_meas.assort_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.assort_neg(:,curr_sub,rep_lev))==0);                                        % Determine which values are finite and real
                        out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);                                                                                          % Determine how many acceptable values there are
                        if out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC(curr_sub,rep_lev)>1                                                                                                       % If there is more than 1 acceptable value
                            out.AUC_thrmat_graph_meas.assort_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                        elseif out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC(curr_sub,rep_lev)==1                                                                                                  % If there is only one acceptable value
                            out.AUC_thrmat_graph_meas.assort_neg(curr_sub,rep_lev) = thrmat_graph_meas.assort_neg(accept_vals,curr_sub,rep_lev);                                                       % Set the AUC as that value
                        else                                                                                                                                                                           % If no acceptable values are found
                            out.AUC_thrmat_graph_meas.assort_neg(curr_sub,rep_lev) = NaN;                                                                                                              % Set the AUC as NaN
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                       = isfinite(thrmat_graph_meas.assort_neg_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.assort_neg_bin(:,curr_sub,rep_lev))==0);                                    % Determine which values are finite and real
                            out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);                                                                                          % Determine how many acceptable values there are
                            if out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC(curr_sub,rep_lev)>1                                                                                                       % If there is more than 1 acceptable value
                                out.AUC_thrmat_graph_meas.assort_neg_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                            elseif out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC(curr_sub,rep_lev)==1                                                                                                  % If there is only one acceptable value
                                out.AUC_thrmat_graph_meas.assort_neg_bin(curr_sub,rep_lev) = thrmat_graph_meas.assort_neg_bin(accept_vals_bin,curr_sub,rep_lev);                                                   % Set the AUC as that value
                            else                                                                                                                                                                               % If no acceptable values are found
                                out.AUC_thrmat_graph_meas.assort_neg_bin(curr_sub,rep_lev) = NaN;                                                                                                              % Set the AUC as NaN
                            end
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);                                                                                                         % Determine how many acceptable values there are
                            if out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1                                                                                                                      % If there is more than 1 acceptable value
                                out.AUC_thrmat_graph_meas.assort_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                            elseif out.AUC_thrmat_graph_meas.assort_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1                                                                                                                 % If there is only one acceptable value
                                out.AUC_thrmat_graph_meas.assort_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.assort_neg(accept_vals,curr_sub,rep_lev);                                                                      % Set the AUC as that value
                            else                                                                                                                                                                                                   % If no acceptable values are found
                                out.AUC_thrmat_graph_meas.assort_neg_nodiscon(curr_sub,rep_lev) = NaN;                                                                                                                             % Set the AUC as NaN
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);                                                                                                         % Determine how many acceptable values there are
                                if out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1                                                                                                                      % If there is more than 1 acceptable value
                                    out.AUC_thrmat_graph_meas.assort_neg_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.assort_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1); % Calculate the AUC, normalized for the number of acceptable values
                                elseif out.AUC_thrmat_graph_meas.assort_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1                                                                                                                 % If there is only one acceptable value
                                    out.AUC_thrmat_graph_meas.assort_neg_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.assort_neg_bin(accept_vals_bin,curr_sub,rep_lev);                                                                  % Set the AUC as that value
                                else                                                                                                                                                                                                       % If no acceptable values are found
                                    out.AUC_thrmat_graph_meas.assort_neg_bin_nodiscon(curr_sub,rep_lev) = NaN;                                                                                                                             % Set the AUC as NaN
                                end
                            end
                        end
                    end
                end
                
                % Brokerage (total)
                if out.calc_props_thrmat.bkg==1
                    accept_vals                                                        = isfinite(thrmat_graph_meas.bkg_tot_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.bkg_tot_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.bkg_tot_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.bkg_tot_pos(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.bkg_tot_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcbinthresh==1
                        accept_vals_bin                                                        = isfinite(thrmat_graph_meas.bkg_tot_pos_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.bkg_tot_pos_bin(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                        if out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.bkg_tot_pos_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.bkg_tot_pos_bin(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.bkg_tot_pos_bin(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                                 = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.bkg_tot_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.bkg_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.bkg_tot_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.bkg_tot_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                                 = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.bkg_tot_pos_bin_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                        = isfinite(thrmat_graph_meas.bkg_tot_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.bkg_tot_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.bkg_tot_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.bkg_tot_neg(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.bkg_tot_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                        = isfinite(thrmat_graph_meas.bkg_tot_neg_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.bkg_tot_neg_bin(:,curr_sub,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.bkg_tot_neg_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.bkg_tot_neg_bin(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.bkg_tot_neg_bin(curr_sub,rep_lev) = NaN;
                            end
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                 = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.bkg_tot_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.bkg_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.bkg_tot_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.bkg_tot_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                 = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.bkg_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.bkg_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.bkg_tot_neg_bin_nodiscon(curr_sub,rep_lev) = NaN;
                                end
                            end
                        end
                    end
                end
                
                % Characteristic Path Length
                if out.calc_props_thrmat.cpl==1
                    accept_vals                                                    = isfinite(thrmat_graph_meas.cpl_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.cpl_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.cpl_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.cpl_pos(curr_sub,rep_lev) = thrmat_graph_meas.cpl_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.cpl_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcbinthresh==1
                        accept_vals_bin                                                    = isfinite(thrmat_graph_meas.cpl_pos_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.cpl_pos_bin(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                        if out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.cpl_pos_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.cpl_pos_bin(curr_sub,rep_lev) = thrmat_graph_meas.cpl_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.cpl_pos_bin(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                             = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.cpl_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.cpl_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.cpl_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.cpl_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.cpl_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                             = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.cpl_pos_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.cpl_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.cpl_pos_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.cpl_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.cpl_pos_bin_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                    = isfinite(thrmat_graph_meas.cpl_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.cpl_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.cpl_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.cpl_neg(curr_sub,rep_lev) = thrmat_graph_meas.cpl_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.cpl_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                    = isfinite(thrmat_graph_meas.cpl_neg_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.cpl_neg_bin(:,curr_sub,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.cpl_neg_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.cpl_neg_bin(curr_sub,rep_lev) = thrmat_graph_meas.cpl_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.cpl_neg_bin(curr_sub,rep_lev) = NaN;
                            end
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                             = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.cpl_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.cpl_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.cpl_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.cpl_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.cpl_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                             = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.cpl_neg_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.cpl_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.cpl_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.cpl_neg_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.cpl_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.cpl_neg_bin_nodiscon(curr_sub,rep_lev) = NaN;
                                end
                            end
                        end
                    end
                end
                
                % Clustering Coefficient (total)
                if out.calc_props_thrmat.clust_coef==1
                    accept_vals                                                               = isfinite(thrmat_graph_meas.clust_coef_tot_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_tot_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.clust_coef_tot_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.clust_coef_tot_pos(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.clust_coef_tot_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcbinthresh==1
                        accept_vals_bin                                                               = isfinite(thrmat_graph_meas.clust_coef_tot_pos_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_tot_pos_bin(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                        if out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                                        = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.clust_coef_tot_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.clust_coef_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.clust_coef_tot_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.clust_coef_tot_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                                        = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.clust_coef_tot_pos_bin_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                               = isfinite(thrmat_graph_meas.clust_coef_tot_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_tot_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.clust_coef_tot_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.clust_coef_tot_neg(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.clust_coef_tot_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                               = isfinite(thrmat_graph_meas.clust_coef_tot_neg_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_tot_neg_bin(:,curr_sub,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin(curr_sub,rep_lev) = NaN;
                            end
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                        = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.clust_coef_tot_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.clust_coef_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.clust_coef_tot_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.clust_coef_tot_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                        = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.clust_coef_tot_neg_bin_nodiscon(curr_sub,rep_lev) = NaN;
                                end
                            end
                        end
                    end
                end
                
                % Zhang & Horvath Clustering Coefficient (total)
                if out.calc_props_thrmat.clust_coef_ZH==1
                    accept_vals                                                                  = isfinite(thrmat_graph_meas.clust_coef_ZH_tot_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_ZH_tot_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_ZH_tot_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                                           = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_ZH_tot_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                                  = isfinite(thrmat_graph_meas.clust_coef_ZH_tot_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_ZH_tot_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_ZH_tot_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                           = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.clust_coef_ZH_tot_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_tot_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                end
                
                % Local Efficiency (total)
                if out.calc_props_thrmat.loc_eff==1
                    accept_vals                                                            = isfinite(thrmat_graph_meas.loc_eff_tot_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_tot_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.loc_eff_tot_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.loc_eff_tot_pos(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.loc_eff_tot_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcbinthresh==1
                        accept_vals_bin                                                            = isfinite(thrmat_graph_meas.loc_eff_tot_pos_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_tot_pos_bin(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                        if out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                                     = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.loc_eff_tot_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.loc_eff_tot_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.loc_eff_tot_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.loc_eff_tot_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                                     = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.loc_eff_tot_pos_bin_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                            = isfinite(thrmat_graph_meas.loc_eff_tot_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_tot_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.loc_eff_tot_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.loc_eff_tot_neg(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.loc_eff_tot_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                            = isfinite(thrmat_graph_meas.loc_eff_tot_neg_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_tot_neg_bin(:,curr_sub,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin(curr_sub,rep_lev) = NaN;
                            end
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                     = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.loc_eff_tot_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.loc_eff_tot_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.loc_eff_tot_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.loc_eff_tot_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                     = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.loc_eff_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.loc_eff_tot_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.loc_eff_tot_neg_bin_nodiscon(curr_sub,rep_lev) = NaN;
                                end
                            end
                        end
                    end
                end
                
                if out.calc_props_thrmat.bkg==1 || out.calc_props_thrmat.clust_coef==1 || out.calc_props_thrmat.clust_coef_ZH==1 || out.calc_props_thrmat.deg==1 || out.calc_props_thrmat.eigvec_cent==1|| out.calc_props_thrmat.gate_coef==1 || out.calc_props_thrmat.kcore_cent==1 || out.calc_props_thrmat.loc_eff==1 || out.calc_props_thrmat.node_bet_cent==1 || out.calc_props_thrmat.pagerank_cent==1 || out.calc_props_thrmat.part_coef==1 || out.calc_props_thrmat.sub_cent==1 || out.calc_props_thrmat.mod_deg_z==1
                    for curr_ROI = 1:out.nROI
                        % Brokerage (per node)
                        if out.calc_props_thrmat.bkg==1
                            accept_vals                                                             = isfinite(thrmat_graph_meas.bkg_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.bkg_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.bkg_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.bkg_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.bkg_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                             = isfinite(thrmat_graph_meas.bkg_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.bkg_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.bkg_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.bkg_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.bkg_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                      = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.bkg_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.bkg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.bkg_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.bkg_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                      = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.bkg_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.bkg_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.bkg_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.bkg_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                             = isfinite(thrmat_graph_meas.bkg_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.bkg_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.bkg_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.bkg_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.bkg_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                             = isfinite(thrmat_graph_meas.bkg_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.bkg_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.bkg_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.bkg_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.bkg_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                      = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.bkg_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.bkg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.bkg_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.bkg_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                      = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.bkg_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.bkg_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.bkg_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.bkg_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.bkg_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.bkg_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Closeness Centrality (per node)
                        if out.calc_props_thrmat.close_cent==1
                            accept_vals                                                                    = isfinite(thrmat_graph_meas.close_cent_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.close_cent_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.close_cent_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.close_cent_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.close_cent_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                    = isfinite(thrmat_graph_meas.close_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.close_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.close_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.close_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.close_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                             = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.close_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.close_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.close_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.close_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                             = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.close_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.close_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.close_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.close_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                    = isfinite(thrmat_graph_meas.close_cent_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.close_cent_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.close_cent_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.close_cent_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.close_cent_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                    = isfinite(thrmat_graph_meas.close_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.close_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.close_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.close_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.close_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                             = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.close_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.close_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.close_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.close_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                             = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.close_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.close_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.close_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.close_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.close_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.close_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Clustering Coefficient (per node)
                        if out.calc_props_thrmat.clust_coef==1
                            accept_vals                                                                    = isfinite(thrmat_graph_meas.clust_coef_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.clust_coef_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.clust_coef_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.clust_coef_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                    = isfinite(thrmat_graph_meas.clust_coef_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.clust_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.clust_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.clust_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                             = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.clust_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.clust_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.clust_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.clust_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                             = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.clust_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.clust_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.clust_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.clust_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                    = isfinite(thrmat_graph_meas.clust_coef_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.clust_coef_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.clust_coef_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.clust_coef_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                    = isfinite(thrmat_graph_meas.clust_coef_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.clust_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.clust_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.clust_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                             = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.clust_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.clust_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.clust_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.clust_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                             = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.clust_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.clust_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.clust_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.clust_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Zhang & Horvath Clustering Coefficient (per node)
                        if out.calc_props_thrmat.clust_coef_ZH==1
                            accept_vals                                                                       = isfinite(thrmat_graph_meas.clust_coef_ZH_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_ZH_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_ZH_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                                = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_ZH_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.clust_coef_ZH_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                       = isfinite(thrmat_graph_meas.clust_coef_ZH_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.clust_coef_ZH_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.clust_coef_ZH_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.clust_coef_ZH_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_ZH_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.clust_coef_ZH_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                                = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.clust_coef_ZH_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.clust_coef_ZH_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.clust_coef_ZH_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                        end
                        
                        % Commn Centrality (per node)
                        if out.calc_props_thrmat.commn_cent==1
                            accept_vals                                                                    = isfinite(thrmat_graph_meas.commn_cent_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.commn_cent_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.commn_cent_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.commn_cent_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.commn_cent_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                    = isfinite(thrmat_graph_meas.commn_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.commn_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.commn_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.commn_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.commn_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                             = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.commn_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.commn_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.commn_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.commn_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                             = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.commn_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.commn_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.commn_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.commn_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                    = isfinite(thrmat_graph_meas.commn_cent_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.commn_cent_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.commn_cent_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.commn_cent_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.commn_cent_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                    = isfinite(thrmat_graph_meas.commn_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.commn_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.commn_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.commn_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.commn_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                             = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.commn_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.commn_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.commn_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.commn_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                             = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.commn_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.commn_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.commn_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.commn_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.commn_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.commn_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Degree
                        if out.calc_props_thrmat.deg==1
                            accept_vals                                                             = isfinite(thrmat_graph_meas.deg_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.deg_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.deg_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.deg_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.deg_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.deg_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.deg_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                      = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.deg_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.deg_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.deg_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.deg_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.deg_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.deg_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                             = isfinite(thrmat_graph_meas.deg_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.deg_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.deg_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.deg_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.deg_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.deg_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.deg_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                      = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.deg_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.deg_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.deg_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.deg_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.deg_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.deg_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                        end
                        
                        % Eigenvector Centrality
                        if out.calc_props_thrmat.eigvec_cent==1
                            accept_vals                                                                     = isfinite(thrmat_graph_meas.eigvec_cent_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.eigvec_cent_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.eigvec_cent_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.eigvec_cent_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.eigvec_cent_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                     = isfinite(thrmat_graph_meas.eigvec_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.eigvec_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                              = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.eigvec_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.eigvec_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.eigvec_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.eigvec_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                              = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.eigvec_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                     = isfinite(thrmat_graph_meas.eigvec_cent_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.eigvec_cent_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.eigvec_cent_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.eigvec_cent_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.eigvec_cent_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                     = isfinite(thrmat_graph_meas.eigvec_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.eigvec_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                              = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.eigvec_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.eigvec_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.eigvec_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.eigvec_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                              = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.eigvec_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.eigvec_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.eigvec_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % K-Coreness
                        if out.calc_props_thrmat.kcore_cent==1
                            accept_vals                                                                    = isfinite(thrmat_graph_meas.kcore_cent_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.kcore_cent_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.kcore_cent_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.kcore_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.kcore_cent_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.kcore_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.kcore_cent_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                             = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.kcore_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.kcore_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.kcore_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.kcore_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.kcore_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.kcore_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                    = isfinite(thrmat_graph_meas.kcore_cent_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.kcore_cent_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.kcore_cent_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.kcore_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.kcore_cent_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.kcore_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.kcore_cent_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                             = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.kcore_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.kcore_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.kcore_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.kcore_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.kcore_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.kcore_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                        end
                        
                        % Gateway Coefficient
                        if out.calc_props_thrmat.gate_coef==1
                            accept_vals                                                                   = isfinite(thrmat_graph_meas.gate_coef_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.gate_coef_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.gate_coef_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.gate_coef_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.gate_coef_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                   = isfinite(thrmat_graph_meas.gate_coef_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.gate_coef_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.gate_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.gate_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.gate_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                            = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.gate_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.gate_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.gate_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.gate_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                            = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.gate_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.gate_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.gate_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.gate_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                   = isfinite(thrmat_graph_meas.gate_coef_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.gate_coef_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.gate_coef_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.gate_coef_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.gate_coef_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                   = isfinite(thrmat_graph_meas.gate_coef_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.gate_coef_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.gate_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.gate_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.gate_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                            = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.gate_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.gate_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.gate_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.gate_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                            = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.gate_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.gate_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.gate_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.gate_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.gate_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.gate_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Local Assortativity
                        if out.calc_props_thrmat.loc_assort==1
                            accept_vals                                                                    = isfinite(thrmat_graph_meas.loc_assort_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.loc_assort_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.loc_assort_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_assort_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.loc_assort_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_assort_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.loc_assort_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                             = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.loc_assort_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_assort_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.loc_assort_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.loc_assort_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_assort_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.loc_assort_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                    = isfinite(thrmat_graph_meas.loc_assort_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.loc_assort_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.loc_assort_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_assort_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.loc_assort_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_assort_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.loc_assort_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                             = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.loc_assort_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_assort_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.loc_assort_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.loc_assort_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_assort_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.loc_assort_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                        end
                        
                        % Local Efficiency (per node)
                        if out.calc_props_thrmat.loc_eff==1
                            accept_vals                                                                 = isfinite(thrmat_graph_meas.loc_eff_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.loc_eff_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.loc_eff_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.loc_eff_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                 = isfinite(thrmat_graph_meas.loc_eff_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.loc_eff_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.loc_eff_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.loc_eff_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                          = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.loc_eff_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.loc_eff_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.loc_eff_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.loc_eff_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                          = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.loc_eff_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.loc_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.loc_eff_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.loc_eff_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                 = isfinite(thrmat_graph_meas.loc_eff_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.loc_eff_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.loc_eff_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.loc_eff_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                 = isfinite(thrmat_graph_meas.loc_eff_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.loc_eff_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.loc_eff_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.loc_eff_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.loc_eff_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                          = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.loc_eff_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.loc_eff_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.loc_eff_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.loc_eff_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                          = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.loc_eff_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.loc_eff_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.loc_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.loc_eff_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.loc_eff_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.loc_eff_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Node betweenness
                        if out.calc_props_thrmat.node_bet_cent==1
                            accept_vals                                                                       = isfinite(thrmat_graph_meas.node_bet_cent_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.node_bet_cent_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.node_bet_cent_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.node_bet_cent_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.node_bet_cent_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                       = isfinite(thrmat_graph_meas.node_bet_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.node_bet_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                                = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.node_bet_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.node_bet_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.node_bet_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.node_bet_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                                = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.node_bet_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                       = isfinite(thrmat_graph_meas.node_bet_cent_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.node_bet_cent_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.node_bet_cent_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.node_bet_cent_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.node_bet_cent_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                       = isfinite(thrmat_graph_meas.node_bet_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.node_bet_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                                = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.node_bet_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.node_bet_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.node_bet_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.node_bet_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                                = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.node_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.node_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.node_bet_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Pagerank Centrality
                        if out.calc_props_thrmat.pagerank_cent==1
                            accept_vals                                                                       = isfinite(thrmat_graph_meas.pagerank_cent_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.pagerank_cent_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.pagerank_cent_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.pagerank_cent_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.pagerank_cent_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                       = isfinite(thrmat_graph_meas.pagerank_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.pagerank_cent_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                                = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.pagerank_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.pagerank_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.pagerank_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.pagerank_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                                = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.pagerank_cent_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                       = isfinite(thrmat_graph_meas.pagerank_cent_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.pagerank_cent_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.pagerank_cent_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.pagerank_cent_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.pagerank_cent_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                       = isfinite(thrmat_graph_meas.pagerank_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.pagerank_cent_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                                = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.pagerank_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.pagerank_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.pagerank_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.pagerank_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                                = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.pagerank_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.pagerank_cent_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.pagerank_cent_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Participation Coefficient
                        if out.calc_props_thrmat.part_coef==1
                            accept_vals                                                                   = isfinite(thrmat_graph_meas.part_coef_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.part_coef_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.part_coef_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.part_coef_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.part_coef_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                   = isfinite(thrmat_graph_meas.part_coef_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.part_coef_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.part_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.part_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.part_coef_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                            = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.part_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.part_coef_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.part_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.part_coef_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                            = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.part_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.part_coef_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.part_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.part_coef_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                   = isfinite(thrmat_graph_meas.part_coef_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.part_coef_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.part_coef_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.part_coef_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.part_coef_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                   = isfinite(thrmat_graph_meas.part_coef_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.part_coef_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.part_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.part_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.part_coef_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                            = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.part_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.part_coef_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.part_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.part_coef_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                            = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.part_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.part_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.part_coef_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.part_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.part_coef_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.part_coef_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Subgraph Centrality
                        if out.calc_props_thrmat.sub_cent==1
                            accept_vals                                                                       = isfinite(thrmat_graph_meas.subgraph_cent_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.subgraph_cent_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.subgraph_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.subgraph_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.subgraph_cent_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.subgraph_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.subgraph_cent_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.subgraph_cent_numvalsAUC_pos(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.subgraph_cent_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.subgraph_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.subgraph_cent_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                                = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.subgraph_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.subgraph_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.subgraph_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.subgraph_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.subgraph_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.subgraph_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.subgraph_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.subgraph_cent_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.subgraph_cent_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                       = isfinite(thrmat_graph_meas.subgraph_cent_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.subgraph_cent_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.subgraph_cent_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.subgraph_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.subgraph_cent_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.subgraph_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.subgraph_cent_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                                = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.subgraph_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.subgraph_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.subgraph_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.subgraph_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.subgraph_cent_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.subgraph_cent_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                        end
                        
                        % Within Module Degree z-score
                        if out.calc_props_thrmat.mod_deg_z==1
                            accept_vals                                                                   = isfinite(thrmat_graph_meas.mod_deg_z_pos(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.mod_deg_z_pos(:,curr_sub,curr_ROI,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                out.AUC_thrmat_graph_meas.mod_deg_z_pos(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                out.AUC_thrmat_graph_meas.mod_deg_z_pos(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.mod_deg_z_pos(curr_sub,curr_ROI,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                   = isfinite(thrmat_graph_meas.mod_deg_z_pos_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.mod_deg_z_pos_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                            = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.mod_deg_z_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_pos(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.mod_deg_z_pos_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.mod_deg_z_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_pos(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.mod_deg_z_pos_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                            = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_pos_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.mod_deg_z_pos_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                            end
                            
                            if strcmp(out.weight_type,'Positive and Negative')
                                accept_vals                                                                   = isfinite(thrmat_graph_meas.mod_deg_z_neg(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.mod_deg_z_neg(:,curr_sub,curr_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.mod_deg_z_neg(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.mod_deg_z_neg(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.mod_deg_z_neg(curr_sub,curr_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                   = isfinite(thrmat_graph_meas.mod_deg_z_neg_bin(:,curr_sub,curr_ROI,rep_lev)) & (imag(thrmat_graph_meas.mod_deg_z_neg_bin(:,curr_sub,curr_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                            = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.mod_deg_z_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_neg(accept_vals,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.mod_deg_z_neg_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.mod_deg_z_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_neg(accept_vals,curr_sub,curr_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.mod_deg_z_neg_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                            = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = trapz(thrmat_graph_meas.mod_deg_z_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = thrmat_graph_meas.mod_deg_z_neg_bin(accept_vals_bin,curr_sub,curr_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.mod_deg_z_neg_bin_nodiscon(curr_sub,curr_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Density
                if out.calc_props_thrmat.dens==1
                    accept_vals                                                     = isfinite(thrmat_graph_meas.dens_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.dens_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.dens_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.dens_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.dens_pos(curr_sub,rep_lev) = thrmat_graph_meas.dens_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.dens_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                              = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.dens_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.dens_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.dens_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.dens_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.dens_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.dens_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                     = isfinite(thrmat_graph_meas.dens_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.dens_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.dens_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.dens_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.dens_neg(curr_sub,rep_lev) = thrmat_graph_meas.dens_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.dens_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                              = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.dens_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.dens_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.dens_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.dens_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.dens_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.dens_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                end
                
                if out.calc_props_thrmat.edge_bet_cent==1 || out.calc_props_thrmat.match==1
                    for curr_row_ROI = 1:out.nROI
                        for curr_col_ROI = 1:out.nROI
                            % Edge Betweenness
                            if out.calc_props_thrmat.edge_bet_cent==1
                                accept_vals                                                                                        = isfinite(thrmat_graph_meas.edge_bet_cent_pos(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.edge_bet_cent_pos(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.edge_bet_cent_pos(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.edge_bet_cent_pos(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.edge_bet_cent_pos(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                                        = isfinite(thrmat_graph_meas.edge_bet_cent_pos_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.edge_bet_cent_pos_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                                                 = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_pos_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.edge_bet_cent_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_pos_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_pos_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                                                 = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_pos_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                                
                                if strcmp(out.weight_type,'Positive and Negative')
                                    accept_vals                                                                                        = isfinite(thrmat_graph_meas.edge_bet_cent_neg(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.edge_bet_cent_neg(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_neg(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_neg(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_neg(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                                        = isfinite(thrmat_graph_meas.edge_bet_cent_neg_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.edge_bet_cent_neg_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                        end
                                    end
                                    
                                    if out.calcAUC_nodiscon==1
                                        accept_vals                                                                                                 = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                        if out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_neg_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.edge_bet_cent_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_neg_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_neg_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                        end
                                        
                                        if out.calcbinthresh==1
                                            accept_vals_bin                                                                                                 = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                            out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                            if out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                                out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.edge_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                            elseif out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                                out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.edge_bet_cent_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                            else
                                                out.AUC_thrmat_graph_meas.edge_bet_cent_neg_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                            end
                                        end
                                    end
                                end
                            end
                            
                            % Matching Index
                            if out.calc_props_thrmat.match==1
                                accept_vals                                                                                = isfinite(thrmat_graph_meas.match_pos(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.match_pos(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                out.AUC_thrmat_graph_meas.match_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.match_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.match_pos(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.match_pos_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.match_pos(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.match_pos(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                end
                                
                                if out.calcbinthresh==1
                                    accept_vals_bin                                                                                = isfinite(thrmat_graph_meas.match_pos_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.match_pos_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.match_pos_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.match_pos_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.match_pos_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                    end
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals                                                                                         = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.match_pos_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.match_pos_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.match_pos_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_pos_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.match_pos_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.match_pos_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_pos(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.match_pos_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                                         = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.match_pos_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.match_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.match_pos_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_pos_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.match_pos_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                        end
                                    end
                                end
                                
                                if strcmp(out.weight_type,'Positive and Negative')
                                    accept_vals                                                                                = isfinite(thrmat_graph_meas.match_neg(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.match_neg(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                    out.AUC_thrmat_graph_meas.match_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                    if out.AUC_thrmat_graph_meas.match_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.match_neg(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.match_neg_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.match_neg(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                    else
                                        out.AUC_thrmat_graph_meas.match_neg(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                    end
                                    
                                    if out.calcbinthresh==1
                                        accept_vals_bin                                                                                = isfinite(thrmat_graph_meas.match_neg_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)) & (imag(thrmat_graph_meas.match_neg_bin(:,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))==0);
                                        out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                        if out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.match_neg_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.match_neg_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.match_neg_bin(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                        end
                                    end
                                    
                                    if out.calcAUC_nodiscon==1
                                        accept_vals                                                                                         = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                        out.AUC_thrmat_graph_meas.match_neg_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals);
                                        if out.AUC_thrmat_graph_meas.match_neg_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                            out.AUC_thrmat_graph_meas.match_neg_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_neg_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                        elseif out.AUC_thrmat_graph_meas.match_neg_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                            out.AUC_thrmat_graph_meas.match_neg_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_neg(accept_vals,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                        else
                                            out.AUC_thrmat_graph_meas.match_neg_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                        end
                                        
                                        if out.calcbinthresh==1
                                            accept_vals_bin                                                                                         = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                            out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = sum(accept_vals_bin);
                                            if out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)>1
                                                out.AUC_thrmat_graph_meas.match_neg_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = trapz(thrmat_graph_meas.match_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev))/(out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)-1);
                                            elseif out.AUC_thrmat_graph_meas.match_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev)==1
                                                out.AUC_thrmat_graph_meas.match_neg_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = thrmat_graph_meas.match_neg_bin(accept_vals_bin,curr_sub,curr_row_ROI,curr_col_ROI,rep_lev);
                                            else
                                                out.AUC_thrmat_graph_meas.match_neg_bin_nodiscon(curr_sub,curr_row_ROI,curr_col_ROI,rep_lev) = NaN;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Global Efficiency
                if out.calc_props_thrmat.glob_eff==1
                    accept_vals                                                         = isfinite(thrmat_graph_meas.glob_eff_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.glob_eff_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.glob_eff_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.glob_eff_pos(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.glob_eff_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcbinthresh==1
                        accept_vals_bin                                                         = isfinite(thrmat_graph_meas.glob_eff_pos_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.glob_eff_pos_bin(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                        if out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.glob_eff_pos_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.glob_eff_pos_bin(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.glob_eff_pos_bin(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                                  = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.glob_eff_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.glob_eff_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.glob_eff_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.glob_eff_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                                  = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.glob_eff_pos_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.glob_eff_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.glob_eff_pos_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.glob_eff_pos_bin_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                         = isfinite(thrmat_graph_meas.glob_eff_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.glob_eff_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.glob_eff_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.glob_eff_neg(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.glob_eff_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                         = isfinite(thrmat_graph_meas.glob_eff_neg_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.glob_eff_neg_bin(:,curr_sub,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.glob_eff_neg_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.glob_eff_neg_bin(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.glob_eff_neg_bin(curr_sub,rep_lev) = NaN;
                            end
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                  = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.glob_eff_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.glob_eff_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.glob_eff_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.glob_eff_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                                  = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.glob_eff_neg_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.glob_eff_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.glob_eff_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.glob_eff_neg_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.glob_eff_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.glob_eff_neg_bin_nodiscon(curr_sub,rep_lev) = NaN;
                                end
                            end
                        end
                    end
                end
                
                % Rich Club
                if out.calc_props_thrmat.rich_club==1
                    curr_club_stats = zeros(size(threshed_conmats_pos,4),out.max_club_size_thr_pos);
                    for curr_dens = 1:size(threshed_conmats_pos,4)
                        if length(thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_pos
                            curr_club_stats(curr_dens,1:length(thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev}))       = thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev}(1:end);
                            curr_club_stats(curr_dens,(length(thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev})+1):end) = NaN;
                        else
                            curr_club_stats(curr_dens,:) = thrmat_graph_meas.rich_club_pos{curr_dens,curr_sub,rep_lev}(1:out.max_club_size_thr_pos);
                        end
                    end
                    for curr_club_size = 1:out.max_club_size_thr_pos
                        accept_vals                                                                         = isfinite(curr_club_stats(:,curr_club_size)) & (imag(curr_club_stats(:,curr_club_size))==0);
                        out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC(curr_sub,curr_club_size,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC(curr_sub,curr_club_size,rep_lev)>1
                            out.AUC_thrmat_graph_meas.rich_club_pos(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats(accept_vals,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC(curr_sub,curr_club_size,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC(curr_sub,curr_club_size,rep_lev)==1
                            out.AUC_thrmat_graph_meas.rich_club_pos(curr_sub,curr_club_size,rep_lev) = curr_club_stats(accept_vals,curr_club_size);
                        else
                            out.AUC_thrmat_graph_meas.rich_club_pos(curr_sub,curr_club_size,rep_lev) = NaN;
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                                  = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)>1
                                out.AUC_thrmat_graph_meas.rich_club_pos_nodiscon(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats(accept_vals,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.rich_club_pos_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)==1
                                out.AUC_thrmat_graph_meas.rich_club_pos_nodiscon(curr_sub,curr_club_size,rep_lev) = curr_club_stats(accept_vals,curr_club_size);
                            else
                                out.AUC_thrmat_graph_meas.rich_club_pos_nodiscon(curr_sub,curr_club_size,rep_lev) = NaN;
                            end
                        end
                    end
                    
                    if out.calcbinthresh==1
                        curr_club_stats_bin = zeros(size(threshed_conmats_pos,4),out.max_club_size_thr_pos_bin);
                        for curr_dens = 1:size(threshed_conmats_pos,4)
                            if length(thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_pos_bin
                                curr_club_stats_bin(curr_dens,1:length(thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev}))       = thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev}(1:end);
                                curr_club_stats_bin(curr_dens,(length(thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev})+1):end) = NaN;
                            else
                                curr_club_stats_bin(curr_dens,:) = thrmat_graph_meas.rich_club_pos_bin{curr_dens,curr_sub,rep_lev}(1:out.max_club_size_thr_pos_bin);
                            end
                        end
                        for curr_club_size = 1:out.max_club_size_thr_pos_bin
                            accept_vals_bin                                                                         = isfinite(curr_club_stats_bin(:,curr_club_size)) & (imag(curr_club_stats_bin(:,curr_club_size))==0);
                            out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev)>1
                                out.AUC_thrmat_graph_meas.rich_club_pos_bin(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats_bin(accept_vals_bin,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev)==1
                                out.AUC_thrmat_graph_meas.rich_club_pos_bin(curr_sub,curr_club_size,rep_lev) = curr_club_stats_bin(accept_vals_bin,curr_club_size);
                            else
                                out.AUC_thrmat_graph_meas.rich_club_pos_bin(curr_sub,curr_club_size,rep_lev) = NaN;
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals_bin                                                                                  = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.rich_club_pos_bin_nodiscon(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats_bin(accept_vals_bin,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.rich_club_pos_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.rich_club_pos_bin_nodiscon(curr_sub,curr_club_size,rep_lev) = curr_club_stats_bin(accept_vals_bin,curr_club_size);
                                else
                                    out.AUC_thrmat_graph_meas.rich_club_pos_bin_nodiscon(curr_sub,curr_club_size,rep_lev) = NaN;
                                end
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        curr_club_stats = zeros(size(threshed_conmats_neg,4),out.max_club_size_thr_neg);
                        for curr_dens = 1:size(threshed_conmats_neg,4)
                            if length(thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_neg
                                curr_club_stats(curr_dens,1:length(thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev}))       = thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev}(1:end);
                                curr_club_stats(curr_dens,(length(thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev})+1):end) = NaN;
                            else
                                curr_club_stats(curr_dens,:) = thrmat_graph_meas.rich_club_neg{curr_dens,curr_sub,rep_lev}(1:out.max_club_size_thr_neg);
                            end
                        end
                        for curr_club_size = 1:out.max_club_size_thr_neg
                            accept_vals                                                                         = isfinite(curr_club_stats(:,curr_club_size)) & (imag(curr_club_stats(:,curr_club_size))==0);
                            out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC(curr_sub,curr_club_size,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC(curr_sub,curr_club_size,rep_lev)>1
                                out.AUC_thrmat_graph_meas.rich_club_neg(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats(accept_vals,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC(curr_sub,curr_club_size,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC(curr_sub,curr_club_size,rep_lev)==1
                                out.AUC_thrmat_graph_meas.rich_club_neg(curr_sub,curr_club_size,rep_lev) = curr_club_stats(accept_vals,curr_club_size);
                            else
                                out.AUC_thrmat_graph_meas.rich_club_neg(curr_sub,curr_club_size,rep_lev) = NaN;
                            end
                            
                            if out.calcAUC_nodiscon==1
                                accept_vals                                                                                  = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev) = sum(accept_vals);
                                if out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.rich_club_neg_nodiscon(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats(accept_vals,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.rich_club_neg_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.rich_club_neg_nodiscon(curr_sub,curr_club_size,rep_lev) = curr_club_stats(accept_vals,curr_club_size);
                                else
                                    out.AUC_thrmat_graph_meas.rich_club_neg_nodiscon(curr_sub,curr_club_size,rep_lev) = NaN;
                                end
                            end
                        end
                        
                        if out.calcbinthresh==1
                            curr_club_stats_bin = zeros(size(threshed_conmats_neg,4),out.max_club_size_thr_neg_bin);
                            for curr_dens = 1:size(threshed_conmats_neg,4)
                                if length(thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev})<out.max_club_size_thr_neg_bin
                                    curr_club_stats_bin(curr_dens,1:length(thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev}))       = thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev}(1:end);
                                    curr_club_stats_bin(curr_dens,(length(thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev})+1):end) = NaN;
                                else
                                    curr_club_stats_bin(curr_dens,:) = thrmat_graph_meas.rich_club_neg_bin{curr_dens,curr_sub,rep_lev}(1:out.max_club_size_thr_neg_bin);
                                end
                            end
                            for curr_club_size = 1:out.max_club_size_thr_neg_bin
                                accept_vals_bin                                                                         = isfinite(curr_club_stats_bin(:,curr_club_size)) & (imag(curr_club_stats_bin(:,curr_club_size))==0);
                                out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.rich_club_neg_bin(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats_bin(accept_vals_bin,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC(curr_sub,curr_club_size,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.rich_club_neg_bin(curr_sub,curr_club_size,rep_lev) = curr_club_stats_bin(accept_vals_bin,curr_club_size);
                                else
                                    out.AUC_thrmat_graph_meas.rich_club_neg_bin(curr_sub,curr_club_size,rep_lev) = NaN;
                                end
                                
                                if out.calcAUC_nodiscon==1
                                    accept_vals_bin                                                                                  = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                    out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev) = sum(accept_vals_bin);
                                    if out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)>1
                                        out.AUC_thrmat_graph_meas.rich_club_neg_bin_nodiscon(curr_sub,curr_club_size,rep_lev) = trapz(curr_club_stats_bin(accept_vals_bin,curr_club_size))/(out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)-1);
                                    elseif out.AUC_thrmat_graph_meas.rich_club_neg_bin_numvalsAUC_nodiscon(curr_sub,curr_club_size,rep_lev)==1
                                        out.AUC_thrmat_graph_meas.rich_club_neg_bin_nodiscon(curr_sub,curr_club_size,rep_lev) = curr_club_stats_bin(accept_vals_bin,curr_club_size);
                                    else
                                        out.AUC_thrmat_graph_meas.rich_club_neg_bin_nodiscon(curr_sub,curr_club_size,rep_lev) = NaN;
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Small Worldness
                if out.calc_props_thrmat.small_world==1
                    accept_vals                                                            = isfinite(thrmat_graph_meas.small_world_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.small_world_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.small_world_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.small_world_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.small_world_pos(curr_sub,rep_lev) = thrmat_graph_meas.small_world_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.small_world_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                                     = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.small_world_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.small_world_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.small_world_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.small_world_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.small_world_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.small_world_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                            = isfinite(thrmat_graph_meas.small_world_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.small_world_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.small_world_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.small_world_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.small_world_neg(curr_sub,rep_lev) = thrmat_graph_meas.small_world_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.small_world_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                                     = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.small_world_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.small_world_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.small_world_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.small_world_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.small_world_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.small_world_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                end
                
                % Transitivity
                if out.calc_props_thrmat.trans==1
                    accept_vals                                                      = isfinite(thrmat_graph_meas.trans_pos(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.trans_pos(:,curr_sub,rep_lev))==0);
                    out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                    if out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC(curr_sub,rep_lev)>1
                        out.AUC_thrmat_graph_meas.trans_pos(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC(curr_sub,rep_lev)-1);
                    elseif out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC(curr_sub,rep_lev)==1
                        out.AUC_thrmat_graph_meas.trans_pos(curr_sub,rep_lev) = thrmat_graph_meas.trans_pos(accept_vals,curr_sub,rep_lev);
                    else
                        out.AUC_thrmat_graph_meas.trans_pos(curr_sub,rep_lev) = NaN;
                    end
                    
                    if out.calcbinthresh==1
                        accept_vals_bin                                                      = isfinite(thrmat_graph_meas.trans_pos_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.trans_pos_bin(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                        if out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.trans_pos_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.trans_pos_bin(curr_sub,rep_lev) = thrmat_graph_meas.trans_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.trans_pos_bin(curr_sub,rep_lev) = NaN;
                        end
                    end
                    
                    if out.calcAUC_nodiscon==1
                        accept_vals                                                               = logical(accept_vals.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                        out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.trans_pos_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_pos(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.trans_pos_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.trans_pos_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.trans_pos(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.trans_pos_nodiscon(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                               = logical(accept_vals_bin.*out.connected_nets_pos(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.trans_pos_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_pos_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.trans_pos_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.trans_pos_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.trans_pos_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.trans_pos_bin_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                        end
                    end
                    
                    if strcmp(out.weight_type,'Positive and Negative')
                        accept_vals                                                      = isfinite(thrmat_graph_meas.trans_neg(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.trans_neg(:,curr_sub,rep_lev))==0);
                        out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals);
                        if out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC(curr_sub,rep_lev)>1
                            out.AUC_thrmat_graph_meas.trans_neg(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC(curr_sub,rep_lev)-1);
                        elseif out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC(curr_sub,rep_lev)==1
                            out.AUC_thrmat_graph_meas.trans_neg(curr_sub,rep_lev) = thrmat_graph_meas.trans_neg(accept_vals,curr_sub,rep_lev);
                        else
                            out.AUC_thrmat_graph_meas.trans_neg(curr_sub,rep_lev) = NaN;
                        end
                        
                        if out.calcbinthresh==1
                            accept_vals_bin                                                      = isfinite(thrmat_graph_meas.trans_neg_bin(:,curr_sub,rep_lev)) & (imag(thrmat_graph_meas.trans_neg_bin(:,curr_sub,rep_lev))==0);
                            out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC(curr_sub,rep_lev) = sum(accept_vals_bin);
                            if out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.trans_neg_bin(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.trans_neg_bin(curr_sub,rep_lev) = thrmat_graph_meas.trans_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.trans_neg_bin(curr_sub,rep_lev) = NaN;
                            end
                        end
                        
                        if out.calcAUC_nodiscon==1
                            accept_vals                                                               = logical(accept_vals.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                            out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals);
                            if out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                out.AUC_thrmat_graph_meas.trans_neg_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_neg(accept_vals,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                            elseif out.AUC_thrmat_graph_meas.trans_neg_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                out.AUC_thrmat_graph_meas.trans_neg_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.trans_neg(accept_vals,curr_sub,rep_lev);
                            else
                                out.AUC_thrmat_graph_meas.trans_neg_nodiscon(curr_sub,rep_lev) = NaN;
                            end
                            
                            if out.calcbinthresh==1
                                accept_vals_bin                                                               = logical(accept_vals_bin.*out.connected_nets_neg(curr_sub,:,rep_lev)');
                                out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev) = sum(accept_vals_bin);
                                if out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)>1
                                    out.AUC_thrmat_graph_meas.trans_neg_bin_nodiscon(curr_sub,rep_lev) = trapz(thrmat_graph_meas.trans_neg_bin(accept_vals_bin,curr_sub,rep_lev))/(out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)-1);
                                elseif out.AUC_thrmat_graph_meas.trans_neg_bin_numvalsAUC_nodiscon(curr_sub,rep_lev)==1
                                    out.AUC_thrmat_graph_meas.trans_neg_bin_nodiscon(curr_sub,rep_lev) = thrmat_graph_meas.trans_neg_bin(accept_vals_bin,curr_sub,rep_lev);
                                else
                                    out.AUC_thrmat_graph_meas.trans_neg_bin_nodiscon(curr_sub,rep_lev) = NaN;
                                end
                            end
                        end
                    end
                end
                
                prog = (curr_sub/out.num_subs)*(1-((rep_lev-1)/out.num_rep_levs));
                if ~use_parfor
                    if out.calcfullmat==1
                        progressbar([],[],prog)                                                                                                                  % Update user
                    else
                        progressbar([],prog)                                                                                                                                                      % Update progress bar
                    end
                else
                    progressbar(prog)
                end
            end
        end
        out.conmats_thr_normed = out.conmats;
        out.conmats            = out.conmats_orig;
    else
        fprintf('Warning: Thresholded properties not computed (see earlier warnings)\n')
    end
end

clear out.conmats_orig

if use_parfor
    try
        parpool close
    catch %#ok<CTCH>
        try
            matlabpool close %#ok<DPOOL>
        catch
            delete(gcp('nocreate'))
        end
    end
end
%%%% Save data
save(out.outname,'out');
fprintf('Done calculating properties!!\n\n')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Create mean contrast matrices

function [full_mean_conmat,grp_mean_conmat] = create_mean_conmats(Z,varargin)

if nargin>1
    grouping_var = varargin{1};
    contin       = varargin{2};
    if nargin>3
        switch varargin{3};
            case '-num_grps'
                num_grps = varargin{4};
            case '-sub_per_grp'
                subs_per_grp = varargin{4};
        end
    else
        subs_per_grp = 50;
    end
end

fish_z                                             = atanh(Z);
full_mean_conmat                                   = atan(mean(fish_z,3));
full_mean_conmat(1:size(full_mean_conmat,1)+1:end) = 0;

if exist('grouping_var','var')
    if contin==0
        group_mem = grouping_var;
    else
        if exist('num_grps','var')
            subs_per_grp = round(size(Z,3)/num_grps);
        else
            num_grps     = round(size(Z,3)/subs_per_grp);
            subs_per_grp = round(size(Z,3)/num_grps);
        end
        
        [~,ind] = sort(grouping_var);
        fish_z  = fish_z(:,:,ind);
        
        for grp = 1:num_grps
            if grp~=num_grps
                group_mem((((grp-1)*subs_per_grp)+1):(grp*subs_per_grp)) = grp;
            else
                group_mem((((grp-1)*subs_per_grp)+1):size(Z,3)) = grp;
            end
        end
    end
    grps            = unique(group_mem);
    grp_mean_conmat = zeros([size(full_mean_conmat),length(grps)]);
    for grp = 1:length(grps)
        temp_conmat                              = atan(mean(fish_z(:,:,logical(group_mem==grps(grp))),3));
        temp_conmat(1:size(temp_conmat,1)+1:end) = 0;
        grp_mean_conmat(:,:,grp)                 = temp_conmat;
    end
end


%%%% Find the minimum density for which a particular set of graphs remain
%%%% connected
function [density,thr,step] = find_min_graph_density(mean_conmat,varargin)

min_step  = 0.00001;
thr       = 0.5;
step      = 0.5;
num_loops = 0;
max_loops = (1/min_step);
if nargin>1
    min_step = varargin{1};
    thr      = varargin{2};
end

while step>min_step && num_loops<=max_loops
    bin_mat = weight_conversion(threshold_absolute(mean_conmat,thr),'binarize');
    R       = reachdist(bin_mat);
    if isempty(find(R==0,1))
        num_in_loops = 0;
        while isempty(find(R==0,1)) && num_in_loops<=(max_loops/100)
            thr          = thr+step;
            bin_mat      = weight_conversion(threshold_absolute(mean_conmat,thr),'binarize');
            R            = reachdist(bin_mat);
            num_in_loops = num_in_loops+1;
        end
    else
        num_in_loops = 0;
        while ~isempty(find(R==0,1)) && num_in_loops<=(max_loops/100)
            thr          = thr-step;
            bin_mat      = weight_conversion(threshold_absolute(mean_conmat,thr),'binarize');
            R            = reachdist(bin_mat);
            num_in_loops = num_in_loops+1;
        end
    end
    step = step/2;
    if thr<0
        density = NaN;
        thr     = NaN;
        return;
    end
    num_loops = num_loops+1;
end

if num_loops==max_loops
    density = NaN;
else
    if ~isempty(find(R==0,1))
        thr = thr-2*step;
    end
    density = density_und(weight_conversion(threshold_absolute(mean_conmat,thr),'binarize'));
end
