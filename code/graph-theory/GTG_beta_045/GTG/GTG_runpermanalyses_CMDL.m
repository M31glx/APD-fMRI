function GTG_runpermanalyses_CMDL(config_filename)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: Beta 0.45 (03.30.16)
%
% Usage:    GTG_runpermanalyses_CMDL('config_filename')
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

% Remove potentially conflicting fields
if isfield(out,'full') %#ok<NODEF>
    out = rmfield(out,'full');
end
if isfield(out,'thrAUC')
    out = rmfield(out,'thrAUC');
end
if isfield(out,'sig_find')
    out = rmfield(out,'sig_find');
end

% Check whether inputs have been specified
if ~exist('out','var')
    msgbox('Enter output from previous stage (structure with graph properties)','Error','error')
    return
end
if isfield(handles.outtemp,'alpha')
    out.alpha = handles.outtemp.alpha;
else
    out.alpha = 0.05;
end
if isfield(handles.outtemp,'outname')
    out.outname = handles.outtemp.outname;
end
if isfield(handles.outtemp,'num_perms')
    out.num_perms = handles.outtemp.num_perms;
else
    out.num_perms = 5000;
end
if isfield(handles.outtemp,'HLA_reg_type')
    out.HLA_reg_type = handles.outtemp.HLA_reg_type;
else
    out.HLA_reg_type = 'OLS';
end
if isfield(handles.outtemp,'weightdirec')
    out.weightdirec = handles.outtemp.weightdirec;
else
    out.weightdirec = 'Positive Only';
end
if isfield(handles.outtemp,'MC_corr')
    out.MC_corr = handles.outtemp.MC_corr;
else
    out.MC_corr = 0;
end
if isfield(handles.outtemp,'same_desmat') || ~strcmp(handles.outtemp.same_desmat,'Use Stage 3 Design Matrix?') 
    out.same_desmat = handles.outtemp.same_desmat;
else
    msgbox('Indicate whether you want to use the design matrix from the previous stage','Error','error')
    return
end
if isfield(handles.outtemp,'desmat')
    out.desmat = handles.outtemp.desmat;
else
    msgbox('Enter a design matrix','Error','error')
    return
end
if isfield(handles.outtemp,'Contrast_or_F')
    out.Contrast_or_F = handles.outtemp.Contrast_or_F;
else
    out.Contrast_or_F = 'Contrasts';
end
if isfield(handles.outtemp,'num_contrasts')
    out.num_contrasts = handles.outtemp.num_contrasts;
else
    msgbox('Enter the number of contrasts desired','Error','error')
    return
end
if isfield(handles.outtemp,'contrasts')
    out.contrasts = handles.outtemp.contrasts;
else
    msgbox('Enter the contrast weights','Error','error')
    return
end
if any(sum(abs(out.contrasts),2)==0)
    msgbox('At least 1 contrast has no non-zero weights','Error','error')
    return
end
if isfield(handles.outtemp,'properties_tested_fullmat')
    out.properties_tested_fullmat = handles.outtemp.properties_tested_fullmat;
end
if isfield(handles.outtemp,'properties_tested_thrmat')
    out.properties_tested_thrmat = handles.outtemp.properties_tested_thrmat;
end
if ~isfield(out,'properties_tested_fullmat') && ~isfield(out,'properties_tested_thrmat')
    msgbox('Select at least one property to test','Error','error')
    return
elseif ~isfield(out,'properties_tested_thrmat') && all(strcmp(out.properties_tested_fullmat,'None')==1)
    msgbox('Select at least one property to test','Error','error')
    return
elseif ~isfield(out,'properties_tested_fullmat') && all(strcmp(out.properties_tested_thrmat,'None')==1)
    msgbox('Select at least one property to test','Error','error')
    return
elseif (isfield(out,'properties_tested_fullmat') && all(strcmp(out.properties_tested_fullmat,'None')==1)) && (isfield(out,'properties_tested_thrmat') && all(strcmp(out.properties_tested_thrmat,'None')==1))
    msgbox('Select at least one property to test','Error','error')
    return
elseif ~isfield(out,'properties_tested_fullmat')
    out.properties_tested_fullmat = {};
elseif ~isfield(out,'properties_tested_thrmat')
    out.properties_tested_thrmat  = {};
end
if isempty(out.denscalc_varmat)
    out.denscalc_varmat = ones(out.num_subs,1);
end
if isfield(handles.outtemp,'IV_names')
    out.IV_names = handles.outtemp.IV_names;
end
if isfield(handles.outtemp,'test_all_nodes')
    out.test_all_nodes = handles.outtemp.test_all_nodes;
end
if isfield(handles.outtemp,'node_mask')
    out.node_mask = handles.outtemp.node_mask;
end
if isfield(handles.outtemp,'test_all_edges')
    out.test_all_edges = handles.outtemp.test_all_edges;
end
if isfield(handles.outtemp,'edge_mask')
    out.edge_mask = handles.outtemp.edge_mask;
end

% Set intial values to 0 (no testing)
out.test_props_fullmat.assort                = 0;
out.test_props_fullmat.bkg                   = 0;
out.test_props_fullmat.bkg_tot               = 0;
out.test_props_fullmat.cpl                   = 0;
out.test_props_fullmat.close_cent            = 0;
out.test_props_fullmat.clust_coef            = 0;
out.test_props_fullmat.clust_coef_tot        = 0;
out.test_props_fullmat.clust_coef_ZH         = 0;
out.test_props_fullmat.clust_coef_ZH_tot     = 0;
out.test_props_fullmat.clust_coef_signed     = 0;
out.test_props_fullmat.clust_coef_signed_tot = 0;
out.test_props_fullmat.commn_cent            = 0;
out.test_props_fullmat.div_coef              = 0;
out.test_props_fullmat.edge_bet_cent         = 0;
out.test_props_fullmat.eigvec_cent           = 0;
out.test_props_fullmat.gate_coef             = 0;
out.test_props_fullmat.glob_eff              = 0;
out.test_props_fullmat.loc_assort            = 0;
out.test_props_fullmat.loc_eff               = 0;
out.test_props_fullmat.loc_eff_tot           = 0;
out.test_props_fullmat.match                 = 0;
out.test_props_fullmat.node_bet_cent         = 0;
out.test_props_fullmat.rich_club             = 0;
out.test_props_fullmat.strength              = 0;
out.test_props_fullmat.strength_tot          = 0;
out.test_props_fullmat.pagerank_cent         = 0;
out.test_props_fullmat.part_coef             = 0;
out.test_props_fullmat.swp                   = 0;
out.test_props_fullmat.trans                 = 0;
out.test_props_fullmat.mod_deg_z             = 0;

out.test_props_thrmat.assort            = 0;
out.test_props_thrmat.bkg               = 0;
out.test_props_thrmat.bkg_tot           = 0;
out.test_props_thrmat.cpl               = 0;
out.test_props_thrmat.close_cent        = 0;
out.test_props_thrmat.clust_coef        = 0;
out.test_props_thrmat.clust_coef_tot    = 0;
out.test_props_thrmat.clust_coef_ZH     = 0;
out.test_props_thrmat.clust_coef_ZH_tot = 0;
out.test_props_thrmat.commn_cent        = 0;
out.test_props_thrmat.deg               = 0;
out.test_props_thrmat.dens              = 0;
out.test_props_thrmat.edge_bet_cent     = 0;
out.test_props_thrmat.eigvec_cent       = 0;
out.test_props_thrmat.gate_coef         = 0;
out.test_props_thrmat.glob_eff          = 0;
out.test_props_thrmat.kcore_cent        = 0;
out.test_props_thrmat.loc_assort        = 0;
out.test_props_thrmat.loc_eff           = 0;
out.test_props_thrmat.loc_eff_tot       = 0;
out.test_props_thrmat.match             = 0;
out.test_props_thrmat.node_bet_cent     = 0;
out.test_props_thrmat.pagerank_cent     = 0;
out.test_props_thrmat.part_coef         = 0;
out.test_props_thrmat.rich_club         = 0;
out.test_props_thrmat.small_world       = 0;
out.test_props_thrmat.sub_cent          = 0;
out.test_props_thrmat.trans             = 0;
out.test_props_thrmat.mod_deg_z         = 0;

net_full_short_names  = {};
net_full_long_names   = {};
node_full_short_names = {};
node_full_long_names  = {};
edge_full_short_names = {};
edge_full_long_names  = {};

% Determine which properties the user selected for fully connected matrices
if ismember('Assortativity',out.properties_tested_fullmat)
    out.test_props_fullmat.assort = 1;
    net_full_short_names          = [net_full_short_names;'assort'];
    net_full_long_names           = [net_full_long_names;'Assortativity'];
end
if ismember('Brokerage',out.properties_tested_fullmat)
    out.test_props_fullmat.bkg = 1;
    node_full_short_names      = [node_full_short_names;'bkg'];
    node_full_long_names       = [node_full_long_names;'Brokerage'];
end
if ismember('Aggregate Brokerage',out.properties_tested_fullmat)
    out.test_props_fullmat.bkg_tot = 1;
    net_full_short_names           = [net_full_short_names;'bkg_tot'];
    net_full_long_names            = [net_full_long_names;'Aggregate Brokerage'];
end
if ismember('Characteristic Path Length',out.properties_tested_fullmat)
    out.test_props_fullmat.cpl = 1;
    net_full_short_names       = [net_full_short_names;'cpl'];
    net_full_long_names        = [net_full_long_names;'Characteristic Path Length'];
end
if ismember('Closeness Centrality',out.properties_tested_fullmat)
    out.test_props_fullmat.close_cent = 1;
    node_full_short_names             = [node_full_short_names;'close_cent'];
    node_full_long_names              = [node_full_long_names;'Closeness Centrality'];
end
if ismember('Clustering Coefficient',out.properties_tested_fullmat)
    out.test_props_fullmat.clust_coef = 1;
    node_full_short_names             = [node_full_short_names;'clust_coef'];
    node_full_long_names              = [node_full_long_names;'Clustering Coefficient'];
end
if ismember('Commn Centrality',out.properties_tested_fullmat)
    out.test_props_fullmat.commn_cent = 1;
    node_full_short_names             = [node_full_short_names;'commn_cent'];
    node_full_long_names              = [node_full_long_names;'Commn Centrality'];
end
if ismember('Mean Clustering Coefficient',out.properties_tested_fullmat)
    out.test_props_fullmat.clust_coef_tot = 1;
    net_full_short_names                  = [net_full_short_names;'clust_coef_tot'];
    net_full_long_names                   = [net_full_long_names;'Mean Clustering Coefficient'];
end
if ismember('Clustering Coefficient (Zhang & Horvath)',out.properties_tested_fullmat)
    out.test_props_fullmat.clust_coef_ZH = 1;
    node_full_short_names                = [node_full_short_names;'clust_coef_ZH'];
    node_full_long_names                 = [node_full_long_names;'Clustering Coefficient (Zhang & Horvath)'];
end
if ismember('Mean Clustering Coefficient (Zhang & Horvath)',out.properties_tested_fullmat)
    out.test_props_fullmat.clust_coef_ZH_tot = 1;
    net_full_short_names                     = [net_full_short_names;'clust_coef_ZH_tot'];
    net_full_long_names                      = [net_full_long_names;'Mean Clustering Coefficient (Zhang & Horvath)'];
end
if ismember('Clustering Coefficient (sign incorporating)',out.properties_tested_fullmat)
    out.test_props_fullmat.clust_coef_signed = 1;
    node_full_short_names                    = [node_full_short_names;'clust_coef_signed'];
    node_full_long_names                     = [node_full_long_names;'Clustering Coefficient (sign incorporating)'];
end
if ismember('Mean Clustering Coefficient (sign incorporating)',out.properties_tested_fullmat)
    out.test_props_fullmat.clust_coef_signed_tot = 1;
    net_full_short_names                         = [net_full_short_names;'clust_coef_signed_tot'];
    net_full_long_names                          = [net_full_long_names;'Mean Clustering Coefficient (sign incorporating)'];
end
if ismember('Diversity Coefficient',out.properties_tested_fullmat)
    out.test_props_fullmat.div_coef = 1;
    node_full_short_names           = [node_full_short_names;'div_coef'];
    node_full_long_names            = [node_full_long_names;'Diversity Coefficient'];
end
if ismember('Edge Betweenness Centrality',out.properties_tested_fullmat) || ismember('Edge Betweeness Centrality',out.properties_tested_fullmat)
    out.test_props_fullmat.edge_bet_cent = 1;
    edge_full_short_names                = [edge_full_short_names;'edge_bet_cent'];
    edge_full_long_names                 = [edge_full_long_names;'Edge Betweenness Centrality'];
end
if ismember('Eigenvector Centrality',out.properties_tested_fullmat)
    out.test_props_fullmat.eigvec_cent = 1;
    node_full_short_names              = [node_full_short_names;'eigvec_cent'];
    node_full_long_names               = [node_full_long_names;'Eigenvector Centrality'];
end
if ismember('Global Efficiency',out.properties_tested_fullmat)
    out.test_props_fullmat.glob_eff = 1;
    net_full_short_names            = [net_full_short_names;'glob_eff'];
    net_full_long_names             = [net_full_long_names;'Global Efficiency'];
end
if ismember('Gateway Coefficient',out.properties_tested_fullmat)
    out.test_props_fullmat.gate_coef = 1;
    node_full_short_names            = [node_full_short_names;'gate_coef'];
    node_full_long_names             = [node_full_long_names;'Gateway Coefficient'];
end
if ismember('Local Assortativity',out.properties_tested_fullmat)
    out.test_props_fullmat.loc_assort = 1;
    node_full_short_names             = [node_full_short_names;'loc_assort'];
    node_full_long_names              = [node_full_long_names;'Local Assortativity'];
end
if ismember('Local Efficiency',out.properties_tested_fullmat)
    out.test_props_fullmat.loc_eff = 1;
    node_full_short_names          = [node_full_short_names;'loc_eff'];
    node_full_long_names           = [node_full_long_names;'Local Efficiency'];
end
if ismember('Mean Local Efficiency',out.properties_tested_fullmat)
    out.test_props_fullmat.loc_eff_tot = 1;
    net_full_short_names               = [net_full_short_names;'loc_eff_tot'];
    net_full_long_names                = [net_full_long_names;'Mean Local Efficiency'];
end
if ismember('Matching Index',out.properties_tested_fullmat)
    out.test_props_fullmat.match = 1;
    edge_full_short_names        = [edge_full_short_names;'match'];
    edge_full_long_names         = [edge_full_long_names;'Matching Index'];
end
if ismember('Node Betweenness Centrality',out.properties_tested_fullmat) || ismember('Node Betweeness Centrality',out.properties_tested_fullmat)
    out.test_props_fullmat.node_bet_cent = 1;
    node_full_short_names                = [node_full_short_names;'node_bet_cent'];
    node_full_long_names                 = [node_full_long_names;'Node Betweenness Centrality'];
end
if ismember('Node Strength',out.properties_tested_fullmat)
    out.test_props_fullmat.strength = 1;
    node_full_short_names           = [node_full_short_names;'strength'];
    node_full_long_names            = [node_full_long_names;'Node Strength'];
end
if ismember('Total Node Strength',out.properties_tested_fullmat)
    out.test_props_fullmat.strength_tot = 1;
    net_full_short_names                = [net_full_short_names;'strength_tot'];
    net_full_long_names                 = [net_full_long_names;'Total Node Strength'];
end
if ismember('PageRank Centrality',out.properties_tested_fullmat)
    out.test_props_fullmat.pagerank_cent = 1;
    node_full_short_names                = [node_full_short_names;'pagerank_cent'];
    node_full_long_names                 = [node_full_long_names;'PageRank Centrality'];
end
if ismember('Participation Coefficient',out.properties_tested_fullmat)
    out.test_props_fullmat.part_coef = 1;
    node_full_short_names            = [node_full_short_names;'part_coef'];
    node_full_long_names             = [node_full_long_names;'Participation Coefficient'];
end
if ismember('Rich Club Networks',out.properties_tested_fullmat)
    out.test_props_fullmat.rich_club = 1;
end
if ismember('Small World Propensity',out.properties_tested_fullmat)
    out.test_props_fullmat.swp = 1;
    net_full_short_names       = [net_full_short_names;'swp'];
    net_full_long_names        = [net_full_long_names;'Small World Propensity'];
end
if ismember('Transitivity',out.properties_tested_fullmat)
    out.test_props_fullmat.trans = 1;
    net_full_short_names         = [net_full_short_names;'trans'];
    net_full_long_names          = [net_full_long_names;'Transitivity'];
end
if ismember('Within-Module Degree Z-Score',out.properties_tested_fullmat)
    out.test_props_fullmat.mod_deg_z = 1;
    node_full_short_names            = [node_full_short_names;'mod_deg_z'];
    node_full_long_names             = [node_full_long_names;'Within-Module Degree Z-Score'];
end

net_thr_short_names  = {};
net_thr_long_names   = {};
node_thr_short_names = {};
node_thr_long_names  = {};
edge_thr_short_names = {};
edge_thr_long_names  = {};

% Determine which properties the user selected for thresholded matrices
if ismember('Assortativity',out.properties_tested_thrmat)
    out.test_props_thrmat.assort = 1;
    net_thr_short_names          = [net_thr_short_names;'assort'];
    net_thr_long_names           = [net_thr_long_names;'Assortativity'];
end
if ismember('Brokerage',out.properties_tested_thrmat)
    out.test_props_thrmat.bkg = 1;
    node_thr_short_names      = [node_thr_short_names;'bkg'];
    node_thr_long_names       = [node_thr_long_names;'Brokerage'];
end
if ismember('Aggregate Brokerage',out.properties_tested_thrmat)
    out.test_props_thrmat.bkg_tot = 1;
    net_thr_short_names           = [net_thr_short_names;'bkg_tot'];
    net_thr_long_names            = [net_thr_long_names;'Aggregate Brokerage'];
end
if ismember('Characteristic Path Length',out.properties_tested_thrmat)
    out.test_props_thrmat.cpl = 1;
    net_thr_short_names       = [net_thr_short_names;'cpl'];
    net_thr_long_names        = [net_thr_long_names;'Characteristic Path Length'];
end
if ismember('Closeness Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.close_cent = 1;
    node_thr_short_names             = [node_thr_short_names;'close_cent'];
    node_thr_long_names              = [node_thr_long_names;'Closeness Centrality'];
end
if ismember('Clustering Coefficient',out.properties_tested_thrmat)
    out.test_props_thrmat.clust_coef = 1;
    node_thr_short_names             = [node_thr_short_names;'clust_coef'];
    node_thr_long_names              = [node_thr_long_names;'Clustering Coefficient'];
end
if ismember('Commn Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.commn_cent = 1;
    node_thr_short_names             = [node_thr_short_names;'commn_cent'];
    node_thr_long_names              = [node_thr_long_names;'Commn Centrality'];
end
if ismember('Mean Clustering Coefficient',out.properties_tested_thrmat)
    out.test_props_thrmat.clust_coef_tot = 1;
    net_thr_short_names                  = [net_thr_short_names;'clust_coef_tot'];
    net_thr_long_names                   = [net_thr_long_names;'Mean Clustering Coefficient'];
end
if ismember('Clustering Coefficient (Zhang & Horvath)',out.properties_tested_thrmat)
    out.test_props_thrmat.clust_coef_ZH = 1;
    node_thr_short_names                = [node_thr_short_names;'clust_coef_ZH'];
    node_thr_long_names                 = [node_thr_long_names;'Clustering Coefficient (Zhang & Horvath)'];
end
if ismember('Mean Clustering Coefficient (Zhang & Horvath)',out.properties_tested_thrmat)
    out.test_props_thrmat.clust_coef_ZH_tot = 1;
    net_thr_short_names                     = [net_thr_short_names;'clust_coef_ZH_tot'];
    net_thr_long_names                      = [net_thr_long_names;'Mean Clustering Coefficient (Zhang & Horvath)'];
end
if ismember('Degree',out.properties_tested_thrmat)
    out.test_props_thrmat.deg = 1;
    node_thr_short_names      = [node_thr_short_names;'deg'];
    node_thr_long_names       = [node_thr_long_names;'Degree'];
end
if ismember('Density',out.properties_tested_thrmat)
    out.test_props_thrmat.dens = 1;
    net_thr_short_names        = [net_thr_short_names;'dens'];
    net_thr_long_names         = [net_thr_long_names;'Density'];
end
if ismember('Edge Betweenness Centrality',out.properties_tested_thrmat) || ismember('Edge Betweeness Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.edge_bet_cent = 1;
    edge_thr_short_names                = [edge_thr_short_names;'edge_bet_cent'];
    edge_thr_long_names                 = [edge_thr_long_names;'Edge Betweenness Centrality'];
end
if ismember('Eigenvector Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.eigvec_cent = 1;
    node_thr_short_names              = [node_thr_short_names;'eigvec_cent'];
    node_thr_long_names               = [node_thr_long_names;'Eigenvector Centrality'];
end
if ismember('Global Efficiency',out.properties_tested_thrmat)
    out.test_props_thrmat.glob_eff = 1;
    net_thr_short_names            = [net_thr_short_names;'glob_eff'];
    net_thr_long_names             = [net_thr_long_names;'Global Efficiency'];
end
if ismember('K-Coreness Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.kcore_cent = 1;
    node_thr_short_names             = [node_thr_short_names;'kcore_cent'];
    node_thr_long_names              = [node_thr_long_names;'K-Coreness Centrality'];
end
if ismember('Gateway Coefficient',out.properties_tested_thrmat)
    out.test_props_thrmat.gate_coef = 1;
    node_thr_short_names            = [node_thr_short_names;'gate_coef'];
    node_thr_long_names             = [node_thr_long_names;'Gateway Coefficient'];
end
if ismember('Local Assortativity',out.properties_tested_thrmat)
    out.test_props_thrmat.loc_assort = 1;
    node_thr_short_names             = [node_thr_short_names;'loc_assort'];
    node_thr_long_names              = [node_thr_long_names;'Local Assortativity'];
end
if ismember('Local Efficiency',out.properties_tested_thrmat)
    out.test_props_thrmat.loc_eff = 1;
    node_thr_short_names          = [node_thr_short_names;'loc_eff'];
    node_thr_long_names           = [node_thr_long_names;'Local Efficiency'];
end
if ismember('Mean Local Efficiency',out.properties_tested_thrmat)
    out.test_props_thrmat.loc_eff_tot = 1;
    net_thr_short_names               = [net_thr_short_names;'loc_eff_tot'];
    net_thr_long_names                = [net_thr_long_names;'Mean Local Efficiency'];
end
if ismember('Matching Index',out.properties_tested_thrmat)
    out.test_props_thrmat.match = 1;
    edge_thr_short_names        = [edge_thr_short_names;'match'];
    edge_thr_long_names         = [edge_thr_long_names;'Matching Index'];
end
if ismember('Node Betweenness Centrality',out.properties_tested_thrmat) || ismember('Node Betweeness Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.node_bet_cent = 1;
    node_thr_short_names                = [node_thr_short_names;'node_bet_cent'];
    node_thr_long_names                 = [node_thr_long_names;'Node Betweenness Centrality'];
end
if ismember('PageRank Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.pagerank_cent = 1;
    node_thr_short_names                = [node_thr_short_names;'pagerank_cent'];
    node_thr_long_names                 = [node_thr_long_names;'PageRank Centrality'];
end
if ismember('Participation Coefficient',out.properties_tested_thrmat)
    out.test_props_thrmat.part_coef = 1;
    node_thr_short_names            = [node_thr_short_names;'part_coef'];
    node_thr_long_names             = [node_thr_long_names;'Participation Coefficient'];
end
if ismember('Rich Club Networks',out.properties_tested_thrmat)
    out.test_props_thrmat.rich_club = 1;
end
if ismember('Small Worldness',out.properties_tested_thrmat)
    out.test_props_thrmat.small_world = 1;
    net_thr_short_names               = [net_thr_short_names;'small_world'];
    net_thr_long_names                = [net_thr_long_names;'Small Worldness'];
end
if ismember('Subgraph Centrality',out.properties_tested_thrmat)
    out.test_props_thrmat.sub_cent = 1;
    node_thr_short_names           = [node_thr_short_names;'subgraph_cent'];
    node_thr_long_names            = [node_thr_long_names;'Subgraph Centrality'];
end
if ismember('Transitivity',out.properties_tested_thrmat)
    out.test_props_thrmat.trans = 1;
    net_thr_short_names         = [net_thr_short_names;'trans'];
    net_thr_long_names          = [net_thr_long_names;'Transitivity'];
end
if ismember('Within-Module Degree Z-Score',out.properties_tested_thrmat)
    out.test_props_thrmat.mod_deg_z = 1;
    node_thr_short_names            = [node_thr_short_names;'mod_deg_z'];
    node_thr_long_names             = [node_thr_long_names;'Within-Module Degree Z-Score'];
end

if out.test_props_fullmat.bkg                   ==1 || ...
        out.test_props_fullmat.close_cent       ==1 || ...
        out.test_props_fullmat.clust_coef       ==1 || ...
        out.test_props_fullmat.clust_coef_ZH    ==1 || ...
        out.test_props_fullmat.clust_coef_signed==1 || ...
        out.test_props_fullmat.commn_cent       ==1 || ...
        out.test_props_fullmat.div_coef         ==1 || ...
        out.test_props_fullmat.eigvec_cent      ==1 || ...
        out.test_props_fullmat.loc_eff          ==1 || ...
        out.test_props_fullmat.gate_coef        ==1 || ...
        out.test_props_fullmat.node_bet_cent    ==1 || ...
        out.test_props_fullmat.strength         ==1 || ...
        out.test_props_fullmat.pagerank_cent    ==1 || ...
        out.test_props_fullmat.part_coef        ==1 || ...
        out.test_props_fullmat.mod_deg_z        ==1 || ...
        out.test_props_thrmat.close_cent        ==1 || ...
        out.test_props_thrmat.clust_coef        ==1 || ...
        out.test_props_thrmat.clust_coef_ZH     ==1 || ...
        out.test_props_thrmat.commn_cent        ==1 || ...
        out.test_props_thrmat.deg               ==1 || ...
        out.test_props_thrmat.eigvec_cent       ==1 || ...
        out.test_props_thrmat.kcore_cent        ==1 || ...
        out.test_props_thrmat.loc_eff           ==1 || ...
        out.test_props_thrmat.node_bet_cent     ==1 || ...
        out.test_props_thrmat.pagerank_cent     ==1 || ...
        out.test_props_thrmat.part_coef         ==1 || ...
        out.test_props_thrmat.sub_cent          ==1 || ...
        out.test_props_thrmat.mod_deg_z         ==1
    
    if ~isfield(out,'test_all_nodes') && ~isfield(out,'node_mask')
        msgbox('Indicate whether to test all nodes','Error','error')
        return
    elseif out.test_all_nodes==0 && ~isfield(out,'node_mask')
        msgbox('Select which nodes to test','Error','error')
        return
    elseif out.test_all_nodes==0 && sum(out.node_mask(:))==0
        msgbox('Select at least one node to test','Error','error')
        return
    end
end
if out.test_props_fullmat.edge_bet_cent    ==1 || ...
        out.test_props_fullmat.match       ==1 || ...
        out.test_props_thrmat.edge_bet_cent==1 || ...
        out.test_props_thrmat.match        ==1
    
    if ~isfield(out,'test_all_edges') && ~isfield(out,'edge_mask')
        msgbox('Indicate whether to test all edges','Error','error')
        return
    elseif out.test_all_edges==0 && ~isfield(out,'edge_mask')
        msgbox('Select which edges to test','Error','error')
        return
    elseif out.test_all_edges==0 && sum(out.edge_mask(:))==0
        msgbox('Select at least one edge to test','Error','error')
        return
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

if size(out.desmat,2)==1 && unique(out.desmat)==1
    perms        = ones(out.num_subs,out.num_perms,2);
    perms(:,:,2) = -1;
    for curr_perm = out.num_perms:-1:1
        for curr_sub = out.num_subs:-1:1
            perms(curr_sub,curr_perm,:) = perms(curr_sub,curr_perm,randperm(2));
        end
    end
    perms = perms(:,:,1);
else
    perms = zeros(out.num_subs,out.num_perms);
    for curr_perm = out.num_perms:-1:1
        perms(:,curr_perm) = randperm(out.num_subs);
    end
end

if out.num_rep_levs>1
    within_perms = zeros(out.num_subs,out.num_rep_levs,out.num_perms);
    orig_ord     = zeros(out.num_subs,out.num_rep_levs);
    for lev = 1:out.num_rep_levs
        orig_ord(:,lev) = ((((lev*out.num_subs)-out.num_subs)+1):(lev*out.num_subs))';
    end
    for curr_perm = out.num_perms:-1:1
        for curr_sub = out.num_subs:-1:1
            within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
        end
    end
    out.HLA_reg_type = 'OLS';
else
    within_perms = [];
end

if out.MC_corr==1
    MC_permdis = NaN(out.num_perms,1);
end

if out.test_props_fullmat.rich_club==1 || out.test_props_thrmat.rich_club==1
    out.rich_club_min_n = str2double(cell2mat(inputdlg(sprintf('Rich clubs often cannot be calculated for a subset of participants.\nThis data is automatically excluded in analyses, leaving only a\nsubset of participants available for analysis.\n\nEnter the minimum number of participants you will allow:\n'),'Min n for rich club',2)));
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
                parpool('open',num_par_workers);
            catch
                matlabpool('open',num_par_workers); %#ok<DPOOL>
            end
        else
            use_parfor = false;
        end
    end
end

rng('default');
rng('shuffle');

out.outname    = strrep(out.outname,'.mat','_permoutput.mat');
sigeffects_fid = fopen(strrep(out.outname,'.mat','_sig_analyses.txt'),'w');
desmat_vars    = 'Design matrix = [';
for var = 1:length(out.IV_names)
    desmat_vars = [desmat_vars,' ',out.IV_names{1,var}];
end
desmat_vars = [desmat_vars,']'];
fprintf(sigeffects_fid,'%s\n\n',desmat_vars);

for con = 1:size(out.contrasts,1)
    fprintf(sigeffects_fid,'Contrast %s:  [%s]\n',num2str(con),num2str(out.contrasts(con,:)));
end
fprintf(sigeffects_fid,'\n\n');

if out.num_rep_levs>2
    NaNout = ones((out.num_rep_levs+1),1).*NaN;
else
    NaNout = ones(out.num_rep_levs,1).*NaN;
end

for con = 1:size(out.contrasts,1)
    fprintf('Running permutation analyses for contrast/F-test %s ...\n',num2str(con))
    curr_con = out.contrasts(con,:)';
    progressbar(['Progress For Contrast/F-Test ',num2str(con)])                                         % Initialize progress bars at zero
    
    tot_props_to_calc = out.test_props_fullmat.assort+ ...
        out.test_props_fullmat.bkg+ ...
        out.test_props_fullmat.bkg_tot+ ...
        out.test_props_fullmat.cpl+ ...
        out.test_props_fullmat.close_cent+ ...
        out.test_props_fullmat.clust_coef+ ...
        out.test_props_fullmat.clust_coef_tot+ ...
        out.test_props_fullmat.clust_coef_ZH+ ...
        out.test_props_fullmat.clust_coef_ZH_tot+ ...
        out.test_props_fullmat.clust_coef_signed+ ...
        out.test_props_fullmat.clust_coef_signed_tot+ ...
        out.test_props_fullmat.commn_cent+ ...
        out.test_props_fullmat.div_coef+ ...
        out.test_props_fullmat.edge_bet_cent+ ...
        out.test_props_fullmat.eigvec_cent+ ...
        out.test_props_fullmat.glob_eff+ ...
        out.test_props_fullmat.gate_coef+ ...
        out.test_props_fullmat.loc_assort+ ...
        out.test_props_fullmat.loc_eff+ ...
        out.test_props_fullmat.loc_eff_tot+ ...
        out.test_props_fullmat.match+ ...
        out.test_props_fullmat.node_bet_cent+ ...
        out.test_props_fullmat.rich_club+ ...
        out.test_props_fullmat.strength+ ...
        out.test_props_fullmat.strength_tot+ ...
        out.test_props_fullmat.pagerank_cent+ ...
        out.test_props_fullmat.part_coef+ ...
        out.test_props_fullmat.swp+ ...
        out.test_props_fullmat.trans+ ...
        out.test_props_fullmat.mod_deg_z+ ...
        out.test_props_thrmat.assort+ ...
        out.test_props_thrmat.bkg+ ...
        out.test_props_thrmat.bkg_tot+ ...
        out.test_props_thrmat.cpl+ ...
        out.test_props_thrmat.close_cent+ ...
        out.test_props_thrmat.clust_coef+ ...
        out.test_props_thrmat.clust_coef_tot+ ...
        out.test_props_thrmat.clust_coef_ZH+ ...
        out.test_props_thrmat.clust_coef_ZH_tot+ ...
        out.test_props_thrmat.commn_cent+ ...
        out.test_props_thrmat.deg+ ...
        out.test_props_thrmat.dens+ ...
        out.test_props_thrmat.edge_bet_cent+ ...
        out.test_props_thrmat.eigvec_cent+ ...
        out.test_props_thrmat.glob_eff+ ...
        out.test_props_thrmat.kcore_cent+ ...
        out.test_props_thrmat.gate_coef+ ...
        out.test_props_thrmat.loc_assort+ ...
        out.test_props_thrmat.loc_eff+ ...
        out.test_props_thrmat.loc_eff_tot+ ...
        out.test_props_thrmat.match+ ...
        out.test_props_thrmat.node_bet_cent+ ...
        out.test_props_thrmat.pagerank_cent+ ...
        out.test_props_thrmat.part_coef+ ...
        out.test_props_thrmat.rich_club+ ...
        out.test_props_thrmat.small_world+ ...
        out.test_props_thrmat.sub_cent+ ...
        out.test_props_thrmat.trans+ ...
        out.test_props_thrmat.mod_deg_z;
    curr_props_calcd = 0;
    
    net_full_props_to_calc = out.test_props_fullmat.assort+ ...
        out.test_props_fullmat.bkg_tot+ ...
        out.test_props_fullmat.cpl+ ...
        out.test_props_fullmat.clust_coef_tot+ ...
        out.test_props_fullmat.clust_coef_ZH_tot+ ...
        out.test_props_fullmat.clust_coef_signed_tot+ ...
        out.test_props_fullmat.glob_eff+ ...
        out.test_props_fullmat.loc_eff_tot+ ...
        out.test_props_fullmat.strength_tot+ ...
        out.test_props_fullmat.swp+ ...
        out.test_props_fullmat.trans;
    
    net_thr_props_to_calc = out.test_props_thrmat.assort+ ...
        out.test_props_thrmat.bkg_tot+ ...
        out.test_props_thrmat.cpl+ ...
        out.test_props_thrmat.clust_coef_tot+ ...
        out.test_props_thrmat.clust_coef_ZH_tot+ ...
        out.test_props_thrmat.dens+ ...
        out.test_props_thrmat.glob_eff+ ...
        out.test_props_thrmat.loc_eff_tot+ ...
        out.test_props_thrmat.small_world+ ...
        out.test_props_thrmat.trans;
    
    node_full_props_to_calc = out.test_props_fullmat.bkg+ ...
        out.test_props_fullmat.close_cent+ ...
        out.test_props_fullmat.clust_coef+ ...
        out.test_props_fullmat.clust_coef_ZH+ ...
        out.test_props_fullmat.clust_coef_signed+ ...
        out.test_props_fullmat.commn_cent+ ...
        out.test_props_fullmat.div_coef+ ...
        out.test_props_fullmat.eigvec_cent+ ...
        out.test_props_fullmat.gate_coef+ ...
        out.test_props_fullmat.loc_assort+ ...
        out.test_props_fullmat.loc_eff+ ...
        out.test_props_fullmat.node_bet_cent+ ...
        out.test_props_fullmat.strength+ ...
        out.test_props_fullmat.pagerank_cent+ ...
        out.test_props_fullmat.part_coef+ ...
        out.test_props_fullmat.mod_deg_z;
    
    node_thr_props_to_calc = out.test_props_thrmat.bkg+ ...
        out.test_props_thrmat.close_cent+ ...
        out.test_props_thrmat.clust_coef+ ...
        out.test_props_thrmat.clust_coef_ZH+ ...
        out.test_props_thrmat.commn_cent+ ...
        out.test_props_thrmat.deg+ ...
        out.test_props_thrmat.eigvec_cent+ ...
        out.test_props_thrmat.kcore_cent+ ...
        out.test_props_thrmat.gate_coef+ ...
        out.test_props_thrmat.loc_assort+ ...
        out.test_props_thrmat.loc_eff+ ...
        out.test_props_thrmat.node_bet_cent+ ...
        out.test_props_thrmat.pagerank_cent+ ...
        out.test_props_thrmat.part_coef+ ...
        out.test_props_thrmat.sub_cent+ ...
        out.test_props_thrmat.mod_deg_z;
    
    edge_full_props_to_calc = out.test_props_fullmat.edge_bet_cent+ ...
        out.test_props_fullmat.match;
    
    edge_thr_props_to_calc = out.test_props_thrmat.edge_bet_cent+ ...
        out.test_props_thrmat.match;
    
    %%%% Less efficient, but more parsable, code for creating contrast
    %%%% predictor
    % Q = inv(out.desmat'*out.desmat);
    % F1 = inv(curr_con'*Q*curr_con);
    % currpred_desmat = out.desmat*Q*curr_con*F1;
    % Pc = curr_con*inv(curr_con'*Q*curr_con)*curr_con'*Q;
    % orth_con = null(curr_con');
    % c3 = orth_con-Pc*orth_con;
    % F3 = inv(c3'*Q*c3);
    % currcovars_desmat = out.desmat*Q*c3*F3;
    
    %%%% Create contrast predictor (and orthogonal predictors for rest of
    %%%% design matrix)
    if size(out.desmat,2)==1 && unique(out.desmat)==1
        currcovars_desmat = [];
        curr_full_desmat  = out.desmat;
    else
        if strcmp(out.Contrast_or_F,'Contrasts')
            currpred_desmat   = out.desmat/(out.desmat'*out.desmat)*curr_con/(curr_con'/(out.desmat'*out.desmat)*curr_con);
            Pc                = curr_con/(curr_con'/(out.desmat'*out.desmat)*curr_con)*curr_con'/(out.desmat'*out.desmat);
            orth_con          = null(curr_con');
            c3                = orth_con-Pc*orth_con;
            currcovars_desmat = out.desmat/(out.desmat'*out.desmat)*c3/(c3'/(out.desmat'*out.desmat)*c3);
            curr_full_desmat  = [currpred_desmat,currcovars_desmat];
        else
            currpred_desmat   = out.desmat(:,logical(curr_con));
            currcovars_desmat = out.desmat(:,~logical(curr_con));
            curr_full_desmat  = [currpred_desmat,currcovars_desmat];
        end
    end
    
    %%%% Methods for calculating p-values
    % Method 1
    % 1. Calc p for neg t-value or 1-p for pos t-value
    % 2. Double that value
    % This method may be problematic if the distribution is extremely
    % skewed
    %
    % Method 2
    % 1. Calc p for both pos and neg values of t
    % 2. Subtract p val for neg t-value from p val for pos t-value
    % 3. Subtract that value from 1
    %
    % To get an accurate 1-tailed p-value, halve the p from method 1
    
    %%%% Run permutation analyses
    % Network-wide properties
    % Full matrices
    for net_prop = 1:net_full_props_to_calc
        curr_short_name = net_full_short_names{net_prop};
        if isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1']) && ~strcmp(curr_short_name,'clust_coef_signed_tot')
            eval(['[out.full.',curr_short_name,'_pos.beta(con,:),out.full.',curr_short_name,'_pos.test_stat(con,:),out.full.',curr_short_name,'_pos.crit_val(con,:),out.full.',curr_short_name,'_pos.p_1tail(con,:),out.full.',curr_short_name,'_pos.p_2tail(con,:),out.full.',curr_short_name,'_pos.hadNaN(con,:),out.full.',curr_short_name,'_pos.hadimag(con,:),out.full.',curr_short_name,'_pos.hadInf(con,:),permdis,out.full.',curr_short_name,'_pos.nonperm_p(con,:)] = GTG_GLM(out.fullmat_graph_meas.',curr_short_name,'_pos,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
            if out.MC_corr==1
                MC_permdis = min(MC_permdis,permdis);
            end
            if strcmp(out.weightdirec,'Positive and Negative')
                eval(['[out.full.',curr_short_name,'_neg.beta(con,:),out.full.',curr_short_name,'_neg.test_stat(con,:),out.full.',curr_short_name,'_neg.crit_val(con,:),out.full.',curr_short_name,'_neg.p_1tail(con,:),out.full.',curr_short_name,'_neg.p_2tail(con,:),out.full.',curr_short_name,'_neg.hadNaN(con,:),out.full.',curr_short_name,'_neg.hadimag(con,:),out.full.',curr_short_name,'_neg.hadInf(con,:),permdis,out.full.',curr_short_name,'_neg.nonperm_p(con,:)] = GTG_GLM(out.fullmat_graph_meas.',curr_short_name,'_neg,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                if out.MC_corr==1
                    MC_permdis = min(MC_permdis,permdis);
                end
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        elseif isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1'])
            eval(['[out.full.',curr_short_name,'.beta(con,:),out.full.',curr_short_name,'.test_stat(con,:),out.full.',curr_short_name,'.crit_val(con,:),out.full.',curr_short_name,'.p_1tail(con,:),out.full.',curr_short_name,'.p_2tail(con,:),out.full.',curr_short_name,'.hadNaN(con,:),out.full.',curr_short_name,'.hadimag(con,:),out.full.',curr_short_name,'.hadInf(con,:),permdis,out.full.',curr_short_name,'.nonperm_p(con,:)] = GTG_GLM(out.fullmat_graph_meas.',curr_short_name,',curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
            if out.MC_corr==1
                MC_permdis = min(MC_permdis,permdis);
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        end
    end
    
    % Thresholded matrices
    for net_prop = 1:net_thr_props_to_calc
        curr_short_name = net_thr_short_names{net_prop};
        if isfield(out.test_props_thrmat,curr_short_name) && eval(['out.test_props_thrmat.',curr_short_name,'==1'])
            eval(['[out.thrAUC.',curr_short_name,'_pos.beta(con,:),out.thrAUC.',curr_short_name,'_pos.test_stat(con,:),out.thrAUC.',curr_short_name,'_pos.crit_val(con,:),out.thrAUC.',curr_short_name,'_pos.p_1tail(con,:),out.thrAUC.',curr_short_name,'_pos.p_2tail(con,:),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,:),out.thrAUC.',curr_short_name,'_pos.hadimag(con,:),out.thrAUC.',curr_short_name,'_pos.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_pos.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
            if out.MC_corr==1
                MC_permdis = min(MC_permdis,permdis);
            end
            if out.calcbinthresh==1 && ~strcmp(curr_short_name,'clust_coef_ZH_tot')
                eval(['[out.thrAUC.',curr_short_name,'_pos_bin.beta(con,:),out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,:),out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,:),out.thrAUC.',curr_short_name,'_pos_bin.p_1tail(con,:),out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,:),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,:),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,:),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_pos_bin.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_bin,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                if out.MC_corr==1
                    MC_permdis = min(MC_permdis,permdis);
                end
            end
            if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                eval(['[out.thrAUC.',curr_short_name,'_neg.beta(con,:),out.thrAUC.',curr_short_name,'_neg.test_stat(con,:),out.thrAUC.',curr_short_name,'_neg.crit_val(con,:),out.thrAUC.',curr_short_name,'_neg.p_1tail(con,:),out.thrAUC.',curr_short_name,'_neg.p_2tail(con,:),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,:),out.thrAUC.',curr_short_name,'_neg.hadimag(con,:),out.thrAUC.',curr_short_name,'_neg.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_neg.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                if out.MC_corr==1
                    MC_permdis = min(MC_permdis,permdis);
                end
                if out.calcbinthresh==1 && ~strcmp(curr_short_name,'clust_coef_ZH_tot')
                    eval(['[out.thrAUC.',curr_short_name,'_neg_bin.beta(con,:),out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,:),out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,:),out.thrAUC.',curr_short_name,'_neg_bin.p_1tail(con,:),out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,:),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,:),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,:),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_neg_bin.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_bin,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                    if out.MC_corr==1
                        MC_permdis = min(MC_permdis,permdis);
                    end
                end
            end
            if out.calcAUC_nodiscon==1
                eval(['[out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_1tail(con,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_pos_nodiscon.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_nodiscon,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                if out.MC_corr==1
                    MC_permdis = min(MC_permdis,permdis);
                end
                if out.calcbinthresh==1 && ~strcmp(curr_short_name,'clust_coef_ZH_tot')
                    eval(['[out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_1tail(con,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_bin_nodiscon,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                    if out.MC_corr==1
                        MC_permdis = min(MC_permdis,permdis);
                    end
                end
                if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                    eval(['[out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_1tail(con,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_neg_nodiscon.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_nodiscon,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                    if out.MC_corr==1
                        MC_permdis = min(MC_permdis,permdis);
                    end
                    if out.calcbinthresh==1 && ~strcmp(curr_short_name,'clust_coef_ZH_tot')
                        eval(['[out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_1tail(con,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,:),permdis,out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.nonperm_p(con,:)] = GTG_GLM(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_bin_nodiscon,curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                        if out.MC_corr==1
                            MC_permdis = min(MC_permdis,permdis);
                        end
                    end
                end
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        end
    end
    
    % Node properties
    % Full matrices
    for node_prop = 1:node_full_props_to_calc
        curr_short_name = node_full_short_names{node_prop};
        if isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1']) && ~strcmp(curr_short_name,'clust_coef_signed')
            for curr_ROI = out.nROI:-1:1
                if out.node_mask(curr_ROI)==1
                    eval(['[out.full.',curr_short_name,'_pos.beta(con,curr_ROI,:),out.full.',curr_short_name,'_pos.test_stat(con,curr_ROI,:),out.full.',curr_short_name,'_pos.crit_val(con,curr_ROI,:),out.full.',curr_short_name,'_pos.p_1tail(con,curr_ROI,:),out.full.',curr_short_name,'_pos.p_2tail(con,curr_ROI,:),out.full.',curr_short_name,'_pos.hadNaN(con,curr_ROI,:),out.full.',curr_short_name,'_pos.hadimag(con,curr_ROI,:),out.full.',curr_short_name,'_pos.hadInf(con,curr_ROI,:),permdis,out.full.',curr_short_name,'_pos.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.fullmat_graph_meas.',curr_short_name,'_pos(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                    if out.MC_corr==1
                        MC_permdis = min(MC_permdis,permdis);
                    end
                    if strcmp(out.weightdirec,'Positive and Negative')
                        eval(['[out.full.',curr_short_name,'_neg.beta(con,curr_ROI,:),out.full.',curr_short_name,'_neg.test_stat(con,curr_ROI,:),out.full.',curr_short_name,'_neg.crit_val(con,curr_ROI,:),out.full.',curr_short_name,'_neg.p_1tail(con,curr_ROI,:),out.full.',curr_short_name,'_neg.p_2tail(con,curr_ROI,:),out.full.',curr_short_name,'_neg.hadNaN(con,curr_ROI,:),out.full.',curr_short_name,'_neg.hadimag(con,curr_ROI,:),out.full.',curr_short_name,'_neg.hadInf(con,curr_ROI,:),permdis,out.full.',curr_short_name,'_neg.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.fullmat_graph_meas.',curr_short_name,'_neg(:,curr_ROI,:)), curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                        if out.MC_corr==1
                            MC_permdis = min(MC_permdis,permdis);
                        end
                    end
                else
                    eval(['out.full.',curr_short_name,'_pos.hadNaN(con,curr_ROI,:)    = NaN;']);
                    eval(['out.full.',curr_short_name,'_pos.hadimag(con,curr_ROI,:)   = NaN;']);
                    eval(['out.full.',curr_short_name,'_pos.hadInf(con,curr_ROI,:)    = NaN;']);
                    eval(['out.full.',curr_short_name,'_pos.beta(con,curr_ROI,:)      = NaNout;']);
                    eval(['out.full.',curr_short_name,'_pos.test_stat(con,curr_ROI,:) = NaNout;']);
                    eval(['out.full.',curr_short_name,'_pos.crit_val(con,curr_ROI,:)  = NaNout;']);
                    eval(['out.full.',curr_short_name,'_pos.p_1tail(con,curr_ROI,:)   = NaNout;']);
                    eval(['out.full.',curr_short_name,'_pos.p_2tail(con,curr_ROI,:)   = NaNout;']);
                    if strcmp(out.weightdirec,'Positive and Negative')
                        eval(['out.full.',curr_short_name,'_neg.hadNaN(con,curr_ROI,:)    = NaN;']);
                        eval(['out.full.',curr_short_name,'_neg.hadimag(con,curr_ROI,:)   = NaN;']);
                        eval(['out.full.',curr_short_name,'_neg.hadInf(con,curr_ROI,:)    = NaN;']);
                        eval(['out.full.',curr_short_name,'_neg.beta(con,curr_ROI,:)      = NaNout;']);
                        eval(['out.full.',curr_short_name,'_neg.test_stat(con,curr_ROI,:) = NaNout;']);
                        eval(['out.full.',curr_short_name,'_neg.crit_val(con,curr_ROI,:)  = NaNout;']);
                        eval(['out.full.',curr_short_name,'_neg.p_1tail(con,curr_ROI,:)   = NaNout;']);
                        eval(['out.full.',curr_short_name,'_neg.p_2tail(con,curr_ROI,:)   = NaNout;']);
                    end
                end
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        elseif isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1'])
            for curr_ROI = out.nROI:-1:1
                if out.node_mask(curr_ROI)==1
                    eval(['[out.full.',curr_short_name,'.beta(con,curr_ROI,:),out.full.',curr_short_name,'.test_stat(con,curr_ROI,:),out.full.',curr_short_name,'.crit_val(con,curr_ROI,:),out.full.',curr_short_name,'.p_1tail(con,curr_ROI,:),out.full.',curr_short_name,'.p_2tail(con,curr_ROI,:),out.full.',curr_short_name,'.hadNaN(con,curr_ROI,:),out.full.',curr_short_name,'.hadimag(con,curr_ROI,:),out.full.',curr_short_name,'.hadInf(con,curr_ROI,:),permdis,out.full.',curr_short_name,'.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.fullmat_graph_meas.',curr_short_name,'(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                    if out.MC_corr==1
                        MC_permdis = min(MC_permdis,permdis);
                    end
                else
                    eval(['out.full.',curr_short_name,'.hadNaN(con,curr_ROI,:)    = NaN;']);
                    eval(['out.full.',curr_short_name,'.hadimag(con,curr_ROI,:)   = NaN;']);
                    eval(['out.full.',curr_short_name,'.hadInf(con,curr_ROI,:)    = NaN;']);
                    eval(['out.full.',curr_short_name,'.beta(con,curr_ROI,:)      = NaNout;']);
                    eval(['out.full.',curr_short_name,'.test_stat(con,curr_ROI,:) = NaNout;']);
                    eval(['out.full.',curr_short_name,'.crit_val(con,curr_ROI,:)  = NaNout;']);
                    eval(['out.full.',curr_short_name,'.p_1tail(con,curr_ROI,:)   = NaNout;']);
                    eval(['out.full.',curr_short_name,'.p_2tail(con,curr_ROI,:)   = NaNout;']);
                end
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        end
    end
    
    % Thresholded matrices
    for node_prop = 1:node_thr_props_to_calc
        curr_short_name = node_thr_short_names{node_prop};
        if isfield(out.test_props_thrmat,curr_short_name) && eval(['out.test_props_thrmat.',curr_short_name,'==1'])
            for curr_ROI = out.nROI:-1:1
                if out.node_mask(curr_ROI)==1
                    eval(['[out.thrAUC.',curr_short_name,'_pos.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos.test_stat(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                    if out.MC_corr==1
                        MC_permdis = min(MC_permdis,permdis);
                    end
                    
                    if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                        eval(['[out.thrAUC.',curr_short_name,'_pos_bin.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos_bin.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_bin(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                        if out.MC_corr==1
                            MC_permdis = min(MC_permdis,permdis);
                        end
                    end
                else
                    eval(['out.thrAUC.',curr_short_name,'_pos.hadNaN(con,curr_ROI,:)    = NaN;']);
                    eval(['out.thrAUC.',curr_short_name,'_pos.hadimag(con,curr_ROI,:)   = NaN;']);
                    eval(['out.thrAUC.',curr_short_name,'_pos.hadInf(con,curr_ROI,:)    = NaN;']);
                    eval(['out.thrAUC.',curr_short_name,'_pos.beta(con,curr_ROI,:)      = NaNout;']);
                    eval(['out.thrAUC.',curr_short_name,'_pos.test_stat(con,curr_ROI,:) = NaNout;']);
                    eval(['out.thrAUC.',curr_short_name,'_pos.crit_val(con,curr_ROI,:)  = NaNout;']);
                    eval(['out.thrAUC.',curr_short_name,'_pos.p_1tail(con,curr_ROI,:)   = NaNout;']);
                    eval(['out.thrAUC.',curr_short_name,'_pos.p_2tail(con,curr_ROI,:)   = NaNout;']);
                    
                    if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,curr_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,curr_ROI,:)   = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,curr_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.beta(con,curr_ROI,:)      = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,curr_ROI,:) = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,curr_ROI,:)  = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.p_1tail(con,curr_ROI,:)   = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,curr_ROI,:)   = NaNout;']);
                    end
                end
                if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                    if out.node_mask(curr_ROI)==1
                        eval(['[out.thrAUC.',curr_short_name,'_neg.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg.test_stat(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                        if out.MC_corr==1
                            MC_permdis = min(MC_permdis,permdis);
                        end
                        
                        if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                            eval(['[out.thrAUC.',curr_short_name,'_neg_bin.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg_bin.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_bin(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                            if out.MC_corr==1
                                MC_permdis = min(MC_permdis,permdis);
                            end
                        end
                    else
                        eval(['out.thrAUC.',curr_short_name,'_neg.hadNaN(con,curr_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_neg.hadimag(con,curr_ROI,:)   = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_neg.hadInf(con,curr_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_neg.beta(con,curr_ROI,:)      = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_neg.test_stat(con,curr_ROI,:) = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_neg.crit_val(con,curr_ROI,:)  = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_neg.p_1tail(con,curr_ROI,:)   = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_neg.p_2tail(con,curr_ROI,:)   = NaNout;']);
                        
                        if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,curr_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,curr_ROI,:)   = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,curr_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.beta(con,curr_ROI,:)      = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,curr_ROI,:) = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,curr_ROI,:)  = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.p_1tail(con,curr_ROI,:)   = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,curr_ROI,:)   = NaNout;']);
                        end
                    end
                end
                if out.calcAUC_nodiscon==1
                    if out.node_mask(curr_ROI)==1
                        eval(['[out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos_nodiscon.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_nodiscon(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                        if out.MC_corr==1
                            MC_permdis = min(MC_permdis,permdis);
                        end
                        if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                            eval(['[out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_bin_nodiscon(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                            if out.MC_corr==1
                                MC_permdis = min(MC_permdis,permdis);
                            end
                        end
                    else
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,curr_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,curr_ROI,:)   = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,curr_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,curr_ROI,:)      = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,curr_ROI,:) = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,curr_ROI,:)  = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.p_1tail(con,curr_ROI,:)   = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,curr_ROI,:)   = NaNout;']);
                        
                        if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,curr_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,curr_ROI,:)   = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,curr_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,curr_ROI,:)      = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,curr_ROI,:) = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,curr_ROI,:)  = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_1tail(con,curr_ROI,:)   = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,curr_ROI,:)   = NaNout;']);
                        end
                    end
                    if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                        if out.node_mask(curr_ROI)==1
                            eval(['[out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg_discon.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_nodiscon(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                            if out.MC_corr==1
                                MC_permdis = min(MC_permdis,permdis);
                            end
                            if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                                eval(['[out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_1tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,curr_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,curr_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.nonperm_p(con,curr_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_bin_nodiscon(:,curr_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                                if out.MC_corr==1
                                    MC_permdis = min(MC_permdis,permdis);
                                end
                            end
                        else
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,curr_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,curr_ROI,:)   = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,curr_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,curr_ROI,:)      = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,curr_ROI,:) = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,curr_ROI,:)  = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.p_1tail(con,curr_ROI,:)   = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,curr_ROI,:)   = NaNout;']);
                            
                            if out.calcbinthresh==1 && ~any(strcmp(curr_short_name,{'clust_coef_ZH','loc_assort'}))
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,curr_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,curr_ROI,:)   = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,curr_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,curr_ROI,:)      = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,curr_ROI,:) = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,curr_ROI,:)  = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_1tail(con,curr_ROI,:)   = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,curr_ROI,:)   = NaNout;']);
                            end
                        end
                    end
                end
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        end
    end
    
    % Edge properties
    % Full matrices
    for edge_prop = 1:edge_full_props_to_calc
        curr_short_name = edge_full_short_names{edge_prop};
        curr_long_name  = edge_full_long_names{edge_prop};
        if isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1'])
            for curr_row_ROI = out.nROI:-1:1
                for curr_col_ROI = out.nROI:-1:1
                    if eval(['out.edge_mask(curr_row_ROI,curr_col_ROI)==1 && any(any(squeeze(out.fullmat_graph_meas.',curr_short_name,'_pos(:,curr_row_ROI,curr_col_ROI,:))))']);
                        eval(['[out.full.',curr_short_name,'_pos.beta(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_pos.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_pos.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_pos.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_pos.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_pos.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_pos.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_pos.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.full.',curr_short_name,'_pos.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.fullmat_graph_meas.',curr_short_name,'_pos(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                        if out.MC_corr==1
                            MC_permdis = min(MC_permdis,permdis);
                        end
                    else
                        eval(['out.full.',curr_short_name,'_pos.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                        eval(['out.full.',curr_short_name,'_pos.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                        eval(['out.full.',curr_short_name,'_pos.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                        eval(['out.full.',curr_short_name,'_pos.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                        eval(['out.full.',curr_short_name,'_pos.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                        eval(['out.full.',curr_short_name,'_pos.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                        eval(['out.full.',curr_short_name,'_pos.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                        eval(['out.full.',curr_short_name,'_pos.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                    end
                    if strcmp(out.weightdirec,'Positive and Negative')
                        if eval(['out.edge_mask(curr_row_ROI,curr_col_ROI)==1 && any(any(squeeze(out.fullmat_graph_meas.',curr_short_name,'_neg(:,curr_row_ROI,curr_col_ROI,:))))']);
                            eval(['[out.full.',curr_short_name,'_neg.beta(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_neg.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_neg.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_neg.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_neg.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_neg.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_neg.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.full.',curr_short_name,'_neg.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.full.',curr_short_name,'_neg.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.fullmat_graph_meas.',curr_short_name,'_neg(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                            if out.MC_corr==1
                                MC_permdis = min(MC_permdis,permdis);
                            end
                        else
                            eval(['out.full.',curr_short_name,'_neg.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.full.',curr_short_name,'_neg.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                            eval(['out.full.',curr_short_name,'_neg.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.full.',curr_short_name,'_neg.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                            eval(['out.full.',curr_short_name,'_neg.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                            eval(['out.full.',curr_short_name,'_neg.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                            eval(['out.full.',curr_short_name,'_neg.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            eval(['out.full.',curr_short_name,'_neg.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                        end
                    end
                end
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        end
    end
    % Thresholded matrices
    for edge_prop = 1:edge_thr_props_to_calc
        curr_short_name = edge_thr_short_names{edge_prop};
        curr_long_name  = edge_thr_long_names{edge_prop};
        if isfield(out.test_props_thrmat,curr_short_name) && eval(['out.test_props_thrmat.',curr_short_name,'==1'])
            for curr_row_ROI = out.nROI:-1:1
                for curr_col_ROI = out.nROI:-1:1
                    if eval(['out.edge_mask(curr_row_ROI,curr_col_ROI)==1 && any(any(squeeze(out.fullmat_graph_meas.',curr_short_name,'(:,curr_row_ROI,curr_col_ROI,:))))']);
                        eval(['[out.thrAUC.',curr_short_name,'_pos.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                        if out.MC_corr==1
                            MC_permdis = min(MC_permdis,permdis);
                        end
                        if out.calcbinthresh==1
                            eval(['[out.thrAUC.',curr_short_name,'_pos_bin.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos_bin.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_bin(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                            if out.MC_corr==1
                                MC_permdis = min(MC_permdis,permdis);
                            end
                        end
                    else
                        eval(['out.thrAUC.',curr_short_name,'_pos.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                        eval(['out.thrAUC.',curr_short_name,'_pos.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                        if out.calcbinthresh==1
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                        end
                    end
                    if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                        if eval(['out.edge_mask(curr_row_ROI,curr_col_ROI)==1 && any(any(squeeze(out.fullmat_graph_meas.',curr_short_name,'(:,curr_row_ROI,curr_col_ROI,:))))']);
                            eval(['[out.thrAUC.',curr_short_name,'_neg.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                            if out.MC_corr==1
                                MC_permdis = min(MC_permdis,permdis);
                            end
                            if out.calcbinthresh==1
                                eval(['[out.thrAUC.',curr_short_name,'_neg_bin.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg_bin.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_bin(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                                if out.MC_corr==1
                                    MC_permdis = min(MC_permdis,permdis);
                                end
                            end
                        else
                            eval(['out.thrAUC.',curr_short_name,'_neg.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_neg.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            if out.calcbinthresh==1
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            end
                        end
                    end
                    if out.calcAUC_nodiscon==1
                        if eval(['out.edge_mask(curr_row_ROI,curr_col_ROI)==1 && any(any(squeeze(out.fullmat_graph_meas.',curr_short_name,'(:,curr_row_ROI,curr_col_ROI,:))))']);
                            eval(['[out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos_nodiscon.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_nodiscon(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                            if out.MC_corr==1
                                MC_permdis = min(MC_permdis,permdis);
                            end
                            
                            if out.calcbinthresh==1
                                eval(['[out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_pos_bin_nodiscon(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                                if out.MC_corr==1
                                    MC_permdis = min(MC_permdis,permdis);
                                end
                            end
                        else
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            if out.calcbinthresh==1
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                            end
                        end
                        if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                            if eval(['out.edge_mask(curr_row_ROI,curr_col_ROI)==1 && any(any(squeeze(out.fullmat_graph_meas.',curr_short_name,'(:,curr_row_ROI,curr_col_ROI,:))))']);
                                eval(['[out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg_discon.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_nodiscon(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                                if out.MC_corr==1
                                    MC_permdis = min(MC_permdis,permdis);
                                end
                                if out.calcbinthresh==1
                                    eval(['[out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:),permdis,out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.nonperm_p(con,curr_row_ROI,curr_col_ROI,:)] = GTG_GLM(squeeze(out.AUC_thrmat_graph_meas.',curr_short_name,'_neg_bin_nodiscon(:,curr_row_ROI,curr_col_ROI,:)),curr_full_desmat,currcovars_desmat,out.Contrast_or_F,perms,out.HLA_reg_type,use_parfor,within_perms,out.alpha,out.MC_corr);']);
                                    if out.MC_corr==1
                                        MC_permdis = min(MC_permdis,permdis);
                                    end
                                end
                            else
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                                eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                                if out.calcbinthresh==1
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,curr_row_ROI,curr_col_ROI,:)   = NaN;']);
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,curr_row_ROI,curr_col_ROI,:)    = NaN;']);
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,curr_row_ROI,curr_col_ROI,:)      = NaNout;']);
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,curr_row_ROI,curr_col_ROI,:) = NaNout;']);
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,curr_row_ROI,curr_col_ROI,:)  = NaNout;']);
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_1tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                                    eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,curr_row_ROI,curr_col_ROI,:)   = NaNout;']);
                                end
                            end
                        end
                    end
                end
            end
            curr_props_calcd = curr_props_calcd+1;
            prog             = curr_props_calcd/tot_props_to_calc;
            progressbar(prog)
        end
    end
    
    % Rich Club Networks
    if out.test_props_fullmat.rich_club==1
        for curr_size = out.max_club_size_full_pos:-1:1
            tempdata = squeeze(out.fullmat_graph_meas.rich_club_pos(:,curr_size,:));
            if any(isnan(tempdata(:)))
                out.full.rich_club_pos.hadNaN(con,curr_size) = 1;
            else
                out.full.rich_club_pos.hadNaN(con,curr_size) = 0;
            end
            notnan                                     = ~sum(isnan(tempdata),2);
            out.full.rich_club_pos.num_subs(curr_size) = sum(notnan);
            if out.full.rich_club_pos.num_subs(curr_size)<out.rich_club_min_n
                out.full.rich_club_pos.hadimag(con,curr_size,:)   = NaN;
                out.full.rich_club_pos.hadInf(con,curr_size,:)    = NaN;
                out.full.rich_club_pos.beta(con,curr_size,:)      = NaNout;
                out.full.rich_club_pos.test_stat(con,curr_size,:) = NaNout;
                out.full.rich_club_pos.crit_val(con,curr_size,:)  = NaNout;
                out.full.rich_club_pos.p_1tail(con,curr_size,:)   = NaNout;
                out.full.rich_club_pos.p_2tail(con,curr_size,:)   = NaNout;
            else
                temp_perms = perms(logical(perms<=out.full.rich_club_pos.num_subs(curr_size)));
                temp_perms = reshape(temp_perms,out.full.rich_club_pos.num_subs(curr_size),(length(temp_perms)/out.full.rich_club_pos.num_subs(curr_size)));
                if out.num_rep_levs>1
                    temp_within_perms = zeros(out.full.rich_club_pos.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                    orig_ord          = [1:out.full.rich_club_pos.num_subs(curr_size);(out.full.rich_club_pos.num_subs(curr_size)+1):(out.full.rich_club_pos.num_subs(curr_size)*2)]';
                    for curr_perm = out.num_perms:-1:1
                        for curr_sub = out.full.rich_club_pos.num_subs(curr_size):-1:1
                            temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                        end
                    end
                else
                    temp_within_perms = [];
                end
                [out.full.rich_club_pos.beta(con,curr_size,:), ...
                    out.full.rich_club_pos.test_stat(con,curr_size,:), ...
                    out.full.rich_club_pos.crit_val(con,curr_size,:), ...
                    out.full.rich_club_pos.p_1tail(con,curr_size,:), ...
                    out.full.rich_club_pos.p_2tail(con,curr_size,:), ...
                    ~, ...
                    out.full.rich_club_pos.hadimag(con,curr_size,:), ...
                    out.full.rich_club_pos.hadInf(con,curr_size,:), ...
                    permdis, ...
                    out.full.rich_club_pos.nonperm_p(con,curr_size,:)] = ...
                    GTG_GLM(tempdata(notnan,:), ...
                    curr_full_desmat(notnan,:), ...
                    currcovars_desmat(notnan,:), ...
                    out.Contrast_or_F, ...
                    temp_perms, ...
                    out.HLA_reg_type, ...
                    use_parfor, ...
                    temp_within_perms, ...
                    out.alpha, ...
                    out.MC_corr);
%                 if out.MC_corr==1
%                     MC_permdis = min(MC_permdis,permdis);
%                 end
            end
        end
        if strcmp(out.weightdirec,'Positive and Negative')
            for curr_size = out.max_club_size_full_neg:-1:1
                tempdata = squeeze(out.fullmat_graph_meas.rich_club_neg(:,curr_size,:));
                if any(isnan(tempdata(:)))
                    out.full.rich_club_neg.hadNaN(con,curr_size) = 1;
                else
                    out.full.rich_club_neg.hadNaN(con,curr_size) = 0;
                end
                notnan                                     = ~sum(isnan(tempdata),2);
                out.full.rich_club_neg.num_subs(curr_size) = sum(notnan);
                if out.full.rich_club_neg.num_subs(curr_size)<out.rich_club_min_n
                    out.full.rich_club_neg.hadimag(con,curr_size,:)   = NaN;
                    out.full.rich_club_neg.hadInf(con,curr_size,:)    = NaN;
                    out.full.rich_club_neg.beta(con,curr_size,:)      = NaNout;
                    out.full.rich_club_neg.test_stat(con,curr_size,:) = NaNout;
                    out.full.rich_club_neg.crit_val(con,curr_size,:)  = NaNout;
                    out.full.rich_club_neg.p_1tail(con,curr_size,:)   = NaNout;
                    out.full.rich_club_neg.p_2tail(con,curr_size,:)   = NaNout;
                else
                    temp_perms = perms(logical(perms<=out.full.rich_club_neg.num_subs(curr_size)));
                    temp_perms = reshape(temp_perms,out.full.rich_club_neg.num_subs(curr_size),(length(temp_perms)/out.full.rich_club_neg.num_subs(curr_size)));
                    if out.num_rep_levs>1
                        temp_within_perms = zeros(out.full.rich_club_neg.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                        orig_ord          = [1:out.full.rich_club_neg.num_subs(curr_size);(out.full.rich_club_neg.num_subs(curr_size)+1):(out.full.rich_club_neg.num_subs(curr_size)*2)]';
                        for curr_perm = out.num_perms:-1:1
                            for curr_sub = out.full.rich_club_neg.num_subs(curr_size):-1:1
                                temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                            end
                        end
                    else
                        temp_within_perms = [];
                    end
                    [out.full.rich_club_neg.beta(con,curr_size,:), ...
                        out.full.rich_club_neg.test_stat(con,curr_size,:), ...
                        out.full.rich_club_neg.crit_val(con,curr_size,:), ...
                        out.full.rich_club_neg.p_1tail(con,curr_size,:), ...
                        out.full.rich_club_neg.p_2tail(con,curr_size,:), ...
                        ~, ...
                        out.full.rich_club_neg.hadimag(con,curr_size,:), ...
                        out.full.rich_club_neg.hadInf(con,curr_size,:), ...
                        permdis, ...
                        out.full.rich_club_neg.nonperm_p(con,curr_size,:)] = ...
                        GTG_GLM(tempdata(notnan,:), ...
                        curr_full_desmat(notnan,:), ...
                        currcovars_desmat(notnan,:), ...
                        out.Contrast_or_F, ...
                        temp_perms, ...
                        out.HLA_reg_type, ...
                        use_parfor, ...
                        temp_within_perms, ...
                        out.alpha, ...
                        out.MC_corr);
%                     if out.MC_corr==1
%                         MC_permdis = min(MC_permdis,permdis);
%                     end
                end
            end
        end
        curr_props_calcd = curr_props_calcd+1;
        prog             = curr_props_calcd/tot_props_to_calc;
        progressbar(prog)
    end
    if out.test_props_thrmat.rich_club==1
        for curr_size = out.max_club_size_thr_pos:-1:1
            tempdata = squeeze(out.AUC_thrmat_graph_meas.rich_club_pos(:,curr_size,:));
            if any(isnan(tempdata(:)))
                out.thrAUC.rich_club_pos.hadNaN(con,curr_size) = 1;
            else
                out.thrAUC.rich_club_pos.hadNaN(con,curr_size) = 0;
            end
            notnan                                       = ~sum(isnan(tempdata),2);
            out.thrAUC.rich_club_pos.num_subs(curr_size) = sum(notnan);
            if out.thrAUC.rich_club_pos.num_subs(curr_size)<out.rich_club_min_n
                out.thrAUC.rich_club_pos.hadimag(con,curr_size,:)   = NaN;
                out.thrAUC.rich_club_pos.hadInf(con,curr_size,:)    = NaN;
                out.thrAUC.rich_club_pos.beta(con,curr_size,:)      = NaNout;
                out.thrAUC.rich_club_pos.test_stat(con,curr_size,:) = NaNout;
                out.thrAUC.rich_club_pos.crit_val(con,curr_size,:)  = NaNout;
                out.thrAUC.rich_club_pos.p_1tail(con,curr_size,:)   = NaNout;
                out.thrAUC.rich_club_pos.p_2tail(con,curr_size,:)   = NaNout;
            else
                temp_perms = perms(logical(perms<=out.thrAUC.rich_club_pos.num_subs(curr_size)));
                temp_perms = reshape(temp_perms,out.thrAUC.rich_club_pos.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_pos.num_subs(curr_size)));
                if out.num_rep_levs>1
                    temp_within_perms = zeros(out.thrAUC.rich_club_pos.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                    orig_ord          = [1:out.thrAUC.rich_club_pos.num_subs(curr_size);(out.thrAUC.rich_club_pos.num_subs(curr_size)+1):(out.thrAUC.rich_club_pos.num_subs(curr_size)*2)]';
                    for curr_perm = out.num_perms:-1:1
                        for curr_sub = out.thrAUC.rich_club_pos.num_subs(curr_size):-1:1
                            temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                        end
                    end
                else
                    temp_within_perms = [];
                end
                [out.thrAUC.rich_club_pos.beta(con,curr_size,:), ...
                    out.thrAUC.rich_club_pos.test_stat(con,curr_size,:), ...
                    out.thrAUC.rich_club_pos.crit_val(con,curr_size,:), ...
                    out.thrAUC.rich_club_pos.p_1tail(con,curr_size,:), ...
                    out.thrAUC.rich_club_pos.p_2tail(con,curr_size,:), ...
                    ~, ...
                    out.thrAUC.rich_club_pos.hadimag(con,curr_size,:), ...
                    out.thrAUC.rich_club_pos.hadInf(con,curr_size,:), ...
                    permdis, ...
                    out.thrAUC.rich_club_pos.nonperm_p(con,curr_size,:)] = ...
                    GTG_GLM(tempdata(notnan,:), ...
                    curr_full_desmat(notnan,:), ...
                    currcovars_desmat(notnan,:), ...
                    out.Contrast_or_F, ...
                    temp_perms, ...
                    out.HLA_reg_type, ...
                    use_parfor, ...
                    temp_within_perms, ...
                    out.alpha, ...
                    out.MC_corr);
%                 if out.MC_corr==1
%                     MC_permdis = min(MC_permdis,permdis);
%                 end
            end
        end
        if out.calcbinthresh==1
            for curr_size = out.max_club_size_thr_pos_bin:-1:1
                tempdata = squeeze(out.AUC_thrmat_graph_meas.rich_club_pos_bin(:,curr_size,:));
                if any(isnan(tempdata(:)))
                    out.thrAUC.rich_club_pos_bin.hadNaN(con,curr_size) = 1;
                else
                    out.thrAUC.rich_club_pos_bin.hadNaN(con,curr_size) = 0;
                end
                notnan                                           = ~sum(isnan(tempdata),2);
                out.thrAUC.rich_club_pos_bin.num_subs(curr_size) = sum(notnan);
                if out.thrAUC.rich_club_pos_bin.num_subs(curr_size)<out.rich_club_min_n
                    out.thrAUC.rich_club_pos_bin.hadimag(con,curr_size,:)   = NaN;
                    out.thrAUC.rich_club_pos_bin.hadInf(con,curr_size,:)    = NaN;
                    out.thrAUC.rich_club_pos_bin.beta(con,curr_size,:)      = NaNout;
                    out.thrAUC.rich_club_pos_bin.test_stat(con,curr_size,:) = NaNout;
                    out.thrAUC.rich_club_pos_bin.crit_val(con,curr_size,:)  = NaNout;
                    out.thrAUC.rich_club_pos_bin.p_1tail(con,curr_size,:)   = NaNout;
                    out.thrAUC.rich_club_pos_bin.p_2tail(con,curr_size,:)   = NaNout;
                else
                    temp_perms = perms(logical(perms<=out.thrAUC.rich_club_pos_bin.num_subs(curr_size)));
                    temp_perms = reshape(temp_perms,out.thrAUC.rich_club_pos_bin.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_pos_bin.num_subs(curr_size)));
                    if out.num_rep_levs>1
                        temp_within_perms = zeros(out.thrAUC.rich_club_pos_bin.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                        orig_ord          = [1:out.thrAUC.rich_club_pos_bin.num_subs(curr_size);(out.thrAUC.rich_club_pos_bin.num_subs(curr_size)+1):(out.thrAUC.rich_club_pos_bin.num_subs(curr_size)*2)]';
                        for curr_perm = out.num_perms:-1:1
                            for curr_sub = out.thrAUC.rich_club_pos_bin.num_subs(curr_size):-1:1
                                temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                            end
                        end
                    else
                        temp_within_perms = [];
                    end
                    [out.thrAUC.rich_club_pos_bin.beta(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_bin.test_stat(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_bin.crit_val(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_bin.p_1tail(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_bin.p_2tail(con,curr_size,:), ...
                        ~, ...
                        out.thrAUC.rich_club_pos_bin.hadimag(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_bin.hadInf(con,curr_size,:), ...
                        permdis, ...
                        out.thrAUC.rich_club_pos_bin.nonperm_p(con,curr_size,:)] = ...
                        GTG_GLM(tempdata(notnan,:), ...
                        curr_full_desmat(notnan,:), ...
                        currcovars_desmat(notnan,:), ...
                        out.Contrast_or_F, ...
                        temp_perms, ...
                        out.HLA_reg_type, ...
                        use_parfor, ...
                        temp_within_perms, ...
                        out.alpha, ...
                        out.MC_corr);
%                     if out.MC_corr==1
%                         MC_permdis = min(MC_permdis,permdis);
%                     end
                end
            end
        end
        if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
            for curr_size = out.max_club_size_thr_neg:-1:1
                tempdata = squeeze(out.AUC_thrmat_graph_meas.rich_club_neg(:,curr_size,:));
                if any(isnan(tempdata(:)))
                    out.thrAUC.rich_club_neg.hadNaN(con,curr_size) = 1;
                else
                    out.thrAUC.rich_club_neg.hadNaN(con,curr_size) = 0;
                end
                notnan                                       = ~sum(isnan(tempdata),2);
                out.thrAUC.rich_club_neg.num_subs(curr_size) = sum(notnan);
                if out.thrAUC.rich_club_neg.num_subs(curr_size)<out.rich_club_min_n
                    out.thrAUC.rich_club_neg.hadimag(con,curr_size,:)   = NaN;
                    out.thrAUC.rich_club_neg.hadInf(con,curr_size,:)    = NaN;
                    out.thrAUC.rich_club_neg.beta(con,curr_size,:)      = NaNout;
                    out.thrAUC.rich_club_neg.test_stat(con,curr_size,:) = NaNout;
                    out.thrAUC.rich_club_neg.crit_val(con,curr_size,:)  = NaNout;
                    out.thrAUC.rich_club_neg.p_1tail(con,curr_size,:)   = NaNout;
                    out.thrAUC.rich_club_neg.p_2tail(con,curr_size,:)   = NaNout;
                else
                    temp_perms = perms(logical(perms<=out.thrAUC.rich_club_neg.num_subs(curr_size)));
                    temp_perms = reshape(temp_perms,out.thrAUC.rich_club_neg.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_neg.num_subs(curr_size)));
                    if out.num_rep_levs>1
                        temp_within_perms = zeros(out.thrAUC.rich_club_neg.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                        orig_ord          = [1:out.thrAUC.rich_club_neg.num_subs(curr_size);(out.thrAUC.rich_club_neg.num_subs(curr_size)+1):(out.thrAUC.rich_club_neg.num_subs(curr_size)*2)]';
                        for curr_perm = out.num_perms:-1:1
                            for curr_sub = out.thrAUC.rich_club_neg.num_subs(curr_size):-1:1
                                temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                            end
                        end
                    else
                        temp_within_perms = [];
                    end
                    [out.thrAUC.rich_club_neg.beta(con,curr_size,:), ...
                        out.thrAUC.rich_club_neg.test_stat(con,curr_size,:), ...
                        out.thrAUC.rich_club_neg.crit_val(con,curr_size,:), ...
                        out.thrAUC.rich_club_neg.p_1tail(con,curr_size,:), ...
                        out.thrAUC.rich_club_neg.p_2tail(con,curr_size,:), ...
                        ~, ...
                        out.thrAUC.rich_club_neg.hadimag(con,curr_size,:), ...
                        out.thrAUC.rich_club_neg.hadInf(con,curr_size,:), ...
                        permdis, ...
                        out.thrAUC.rich_club_neg.nonperm_p(con,curr_size,:)] = ...
                        GTG_GLM(tempdata(notnan,:), ...
                        curr_full_desmat(notnan,:), ...
                        currcovars_desmat(notnan,:), ...
                        out.Contrast_or_F, ...
                        temp_perms, ...
                        out.HLA_reg_type, ...
                        use_parfor, ...
                        temp_within_perms, ...
                        out.alpha, ...
                        out.MC_corr);
%                     if out.MC_corr==1
%                         MC_permdis = min(MC_permdis,permdis);
%                     end
                end
            end
        end
        if out.calcbinthresh==1
            if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                for curr_size = out.max_club_size_thr_neg_bin:-1:1
                    tempdata = squeeze(out.AUC_thrmat_graph_meas.rich_club_neg_bin(:,curr_size,:));
                    if any(isnan(tempdata(:)))
                        out.thrAUC.rich_club_neg_bin.hadNaN(con,curr_size) = 1;
                    else
                        out.thrAUC.rich_club_neg_bin.hadNaN(con,curr_size) = 0;
                    end
                    notnan                                           = ~sum(isnan(tempdata),2);
                    out.thrAUC.rich_club_neg_bin.num_subs(curr_size) = sum(notnan);
                    if out.thrAUC.rich_club_neg_bin.num_subs(curr_size)<out.rich_club_min_n
                        out.thrAUC.rich_club_neg_bin.hadimag(con,curr_size,:)   = NaN;
                        out.thrAUC.rich_club_neg_bin.hadInf(con,curr_size,:)    = NaN;
                        out.thrAUC.rich_club_neg_bin.beta(con,curr_size,:)      = NaNout;
                        out.thrAUC.rich_club_neg_bin.test_stat(con,curr_size,:) = NaNout;
                        out.thrAUC.rich_club_neg_bin.crit_val(con,curr_size,:)  = NaNout;
                        out.thrAUC.rich_club_neg_bin.p_1tail(con,curr_size,:)   = NaNout;
                        out.thrAUC.rich_club_neg_bin.p_2tail(con,curr_size,:)   = NaNout;
                    else
                        temp_perms = perms(logical(perms<=out.thrAUC.rich_club_neg_bin.num_subs(curr_size)));
                        temp_perms = reshape(temp_perms,out.thrAUC.rich_club_neg_bin.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_neg_bin.num_subs(curr_size)));
                        if out.num_rep_levs>1
                            temp_within_perms = zeros(out.thrAUC.rich_club_neg_bin.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                            orig_ord          = [1:out.thrAUC.rich_club_neg_bin.num_subs(curr_size);(out.thrAUC.rich_club_neg_bin.num_subs(curr_size)+1):(out.thrAUC.rich_club_neg_bin.num_subs(curr_size)*2)]';
                            for curr_perm = out.num_perms:-1:1
                                for curr_sub = out.thrAUC.rich_club_neg_bin.num_subs(curr_size):-1:1
                                    temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                                end
                            end
                        else
                            temp_within_perms = [];
                        end
                        [out.thrAUC.rich_club_neg_bin.beta(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_bin.test_stat(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_bin.crit_val(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_bin.p_1tail(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_bin.p_2tail(con,curr_size,:), ...
                            ~, ...
                            out.thrAUC.rich_club_neg_bin.hadimag(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_bin.hadInf(con,curr_size,:), ...
                            permdis, ...
                            out.thrAUC.rich_club_neg_bin.nonperm_p(con,curr_size,:)] = ...
                            GTG_GLM(tempdata(notnan,:), ...
                            curr_full_desmat(notnan,:), ...
                            currcovars_desmat(notnan,:), ...
                            out.Contrast_or_F, ...
                            temp_perms, ...
                            out.HLA_reg_type, ...
                            use_parfor, ...
                            temp_within_perms, ...
                            out.alpha, ...
                            out.MC_corr);
%                         if out.MC_corr==1
%                             MC_permdis = min(MC_permdis,permdis);
%                         end
                    end
                end
            end
        end
        if out.calcAUC_nodiscon==1
            for curr_size = out.max_club_size_thr_pos:-1:1
                tempdata  = squeeze(out.AUC_thrmat_graph_meas.rich_club_pos_nodiscon(:,curr_size,:));
                if any(isnan(tempdata(:)))
                    out.thrAUC.rich_club_pos_nodiscon.hadNaN(con,curr_size) = 1;
                else
                    out.thrAUC.rich_club_pos_nodiscon.hadNaN(con,curr_size) = 0;
                end
                notnan                                                = ~sum(isnan(tempdata),2);
                out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size) = sum(notnan);
                if out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size)<out.rich_club_min_n
                    out.thrAUC.rich_club_pos_nodiscon.hadimag(con,curr_size,:)   = NaN;
                    out.thrAUC.rich_club_pos_nodiscon.hadInf(con,curr_size,:)    = NaN;
                    out.thrAUC.rich_club_pos_nodiscon.beta(con,curr_size,:)      = NaNout;
                    out.thrAUC.rich_club_pos_nodiscon.test_stat(con,curr_size,:) = NaNout;
                    out.thrAUC.rich_club_pos_nodiscon.crit_val(con,curr_size,:)  = NaNout;
                    out.thrAUC.rich_club_pos_nodiscon.p_1tail(con,curr_size,:)   = NaNout;
                    out.thrAUC.rich_club_pos_nodiscon.p_2tail(con,curr_size,:)   = NaNout;
                else
                    temp_perms = perms(logical(perms<=out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size)));
                    temp_perms = reshape(temp_perms,out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size)));
                    if out.num_rep_levs>1
                        temp_within_perms = zeros(out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                        orig_ord          = [1:out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size);(out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size)+1):(out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size)*2)]';
                        for curr_perm = out.num_perms:-1:1
                            for curr_sub = out.thrAUC.rich_club_pos_nodiscon.num_subs(curr_size):-1:1
                                temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                            end
                        end
                    else
                        temp_within_perms = [];
                    end
                    [out.thrAUC.rich_club_pos_nodiscon.beta(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_nodiscon.test_stat(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_nodiscon.crit_val(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_nodiscon.p_1tail(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_nodiscon.p_2tail(con,curr_size,:), ...
                        ~, ...
                        out.thrAUC.rich_club_pos_nodiscon.hadimag(con,curr_size,:), ...
                        out.thrAUC.rich_club_pos_nodiscon.hadInf(con,curr_size,:), ...
                        permdis, ...
                        out.thrAUC.rich_club_pos_nodiscon.nonperm_p(con,curr_size,:)] = ...
                        GTG_GLM(tempdata(notnan,:), ...
                        curr_full_desmat(notnan,:), ...
                        currcovars_desmat(notnan,:), ...
                        out.Contrast_or_F, ...
                        temp_perms, ...
                        out.HLA_reg_type, ...
                        use_parfor, ...
                        temp_within_perms, ...
                        out.alpha,out.MC_corr);
%                     if out.MC_corr==1
%                         MC_permdis = min(MC_permdis,permdis);
%                     end
                end
            end
            if out.calcbinthresh==1
                for curr_size = out.max_club_size_thr_pos_bin:-1:1
                    tempdata  = squeeze(out.AUC_thrmat_graph_meas.rich_club_pos_bin_nodiscon(:,curr_size,:));
                    if any(isnan(tempdata(:)))
                        out.thrAUC.rich_club_pos_bin_nodiscon.hadNaN(con,curr_size) = 1;
                    else
                        out.thrAUC.rich_club_pos_bin_nodiscon.hadNaN(con,curr_size) = 0;
                    end
                    notnan                                                    = ~sum(isnan(tempdata),2);
                    out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size) = sum(notnan);
                    if out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size)<out.rich_club_min_n
                        out.thrAUC.rich_club_pos_bin_nodiscon.hadimag(con,curr_size,:)   = NaN;
                        out.thrAUC.rich_club_pos_bin_nodiscon.hadInf(con,curr_size,:)    = NaN;
                        out.thrAUC.rich_club_pos_bin_nodiscon.beta(con,curr_size,:)      = NaNout;
                        out.thrAUC.rich_club_pos_bin_nodiscon.test_stat(con,curr_size,:) = NaNout;
                        out.thrAUC.rich_club_pos_bin_nodiscon.crit_val(con,curr_size,:)  = NaNout;
                        out.thrAUC.rich_club_pos_bin_nodiscon.p_1tail(con,curr_size,:)   = NaNout;
                        out.thrAUC.rich_club_pos_bin_nodiscon.p_2tail(con,curr_size,:)   = NaNout;
                    else
                        temp_perms = perms(logical(perms<=out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size)));
                        temp_perms = reshape(temp_perms,out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size)));
                        if out.num_rep_levs>1
                            temp_within_perms = zeros(out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                            orig_ord = [1:out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size);(out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size)+1):(out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size)*2)]';
                            for curr_perm = out.num_perms:-1:1
                                for curr_sub = out.thrAUC.rich_club_pos_bin_nodiscon.num_subs(curr_size):-1:1
                                    temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                                end
                            end
                        else
                            temp_within_perms = [];
                        end
                        [out.thrAUC.rich_club_pos_bin_nodiscon.beta(con,curr_size,:), ...
                            out.thrAUC.rich_club_pos_bin_nodiscon.test_stat(con,curr_size,:), ...
                            out.thrAUC.rich_club_pos_bin_nodiscon.crit_val(con,curr_size,:), ...
                            out.thrAUC.rich_club_pos_bin_nodiscon.p_1tail(con,curr_size,:), ...
                            out.thrAUC.rich_club_pos_bin_nodiscon.p_2tail(con,curr_size,:), ...
                            ~, ...
                            out.thrAUC.rich_club_pos_bin_nodiscon.hadimag(con,curr_size,:), ...
                            out.thrAUC.rich_club_pos_bin_nodiscon.hadInf(con,curr_size,:), ...
                            permdis, ...
                            out.thrAUC.rich_club_pos_bin_nodiscon.nonperm_p(con,curr_size,:)] = ...
                            GTG_GLM(tempdata(notnan,:), ...
                            curr_full_desmat(notnan,:), ...
                            currcovars_desmat(notnan,:), ...
                            out.Contrast_or_F, ...
                            temp_perms, ...
                            out.HLA_reg_type, ...
                            use_parfor, ...
                            temp_within_perms, ...
                            out.alpha, ...
                            out.MC_corr);
%                         if out.MC_corr==1
%                             MC_permdis = min(MC_permdis,permdis);
%                         end
                    end
                end
            end
            if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                for curr_size = out.max_club_size_thr_neg:-1:1
                    tempdata = squeeze(out.AUC_thrmat_graph_meas.rich_club_neg_nodiscon(:,curr_size,:));
                    if any(isnan(tempdata(:)))
                        out.thrAUC.rich_club_neg_nodiscon.hadNaN(con,curr_size) = 1;
                    else
                        out.thrAUC.rich_club_neg_nodiscon.hadNaN(con,curr_size) = 0;
                    end
                    notnan                                                = ~sum(isnan(tempdata),2);
                    out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size) = sum(notnan);
                    if out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size)<out.rich_club_min_n
                        out.thrAUC.rich_club_neg_nodiscon.hadimag(con,curr_size,:)   = NaN;
                        out.thrAUC.rich_club_neg_nodiscon.hadInf(con,curr_size,:)    = NaN;
                        out.thrAUC.rich_club_neg_nodiscon.beta(con,curr_size,:)      = NaNout;
                        out.thrAUC.rich_club_neg_nodiscon.test_stat(con,curr_size,:) = NaNout;
                        out.thrAUC.rich_club_neg_nodiscon.crit_val(con,curr_size,:)  = NaNout;
                        out.thrAUC.rich_club_neg_nodiscon.p_1tail(con,curr_size,:)   = NaNout;
                        out.thrAUC.rich_club_neg_nodiscon.p_2tail(con,curr_size,:)   = NaNout;
                    else
                        temp_perms = perms(logical(perms<=out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size)));
                        temp_perms = reshape(temp_perms,out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size)));
                        if out.num_rep_levs>1
                            temp_within_perms = zeros(out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                            orig_ord          = [1:out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size);(out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size)+1):(out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size)*2)]';
                            for curr_perm = out.num_perms:-1:1
                                for curr_sub = out.thrAUC.rich_club_neg_nodiscon.num_subs(curr_size):-1:1
                                    temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                                end
                            end
                        else
                            temp_within_perms = [];
                        end
                        [out.thrAUC.rich_club_neg_nodiscon.beta(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_nodiscon.test_stat(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_nodiscon.crit_val(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_nodiscon.p_1tail(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_nodiscon.p_2tail(con,curr_size,:), ...
                            ~, ...
                            out.thrAUC.rich_club_neg_nodiscon.hadimag(con,curr_size,:), ...
                            out.thrAUC.rich_club_neg_nodiscon.hadInf(con,curr_size,:), ...
                            permdis, ...
                            out.thrAUC.rich_club_neg_nodiscon.nonperm_p(con,curr_size,:)] = ...
                            GTG_GLM(tempdata(notnan,:), ...
                            curr_full_desmat(notnan,:), ...
                            currcovars_desmat(notnan,:), ...
                            out.Contrast_or_F, ...
                            temp_perms, ...
                            out.HLA_reg_type, ...
                            use_parfor, ...
                            temp_within_perms, ...
                            out.alpha, ...
                            out.MC_corr);
%                         if out.MC_corr==1
%                             MC_permdis = min(MC_permdis,permdis);
%                         end
                    end
                end
                if out.calcbinthresh==1
                    for curr_size = out.max_club_size_thr_neg_bin:-1:1
                        tempdata  = squeeze(out.AUC_thrmat_graph_meas.rich_club_neg_bin_nodiscon(:,curr_size,:));
                        if any(isnan(tempdata(:)))
                            out.thrAUC.rich_club_neg_bin_nodiscon.hadNaN(con,curr_size) = 1;
                        else
                            out.thrAUC.rich_club_neg_bin_nodiscon.hadNaN(con,curr_size) = 0;
                        end
                        notnan                                                    = ~sum(isnan(tempdata),2);
                        out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size) = sum(notnan);
                        if out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size)<out.rich_club_min_n
                            out.thrAUC.rich_club_neg_bin_nodiscon.hadimag(con,curr_size,:)   = NaN;
                            out.thrAUC.rich_club_neg_bin_nodiscon.hadInf(con,curr_size,:)    = NaN;
                            out.thrAUC.rich_club_neg_bin_nodiscon.beta(con,curr_size,:)      = NaNout;
                            out.thrAUC.rich_club_neg_bin_nodiscon.test_stat(con,curr_size,:) = NaNout;
                            out.thrAUC.rich_club_neg_bin_nodiscon.crit_val(con,curr_size,:)  = NaNout;
                            out.thrAUC.rich_club_neg_bin_nodiscon.p_1tail(con,curr_size,:)   = NaNout;
                            out.thrAUC.rich_club_neg_bin_nodiscon.p_2tail(con,curr_size,:)   = NaNout;
                        else
                            temp_perms = perms(logical(perms<=out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size)));
                            temp_perms = reshape(temp_perms,out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size),(length(temp_perms)/out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size)));
                            if out.num_rep_levs>1
                                temp_within_perms = zeros(out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size),out.num_rep_levs,out.num_perms);
                                orig_ord          = [1:out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size);(out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size)+1):(out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size)*2)]';
                                for curr_perm = out.num_perms:-1:1
                                    for curr_sub = out.thrAUC.rich_club_neg_bin_nodiscon.num_subs(curr_size):-1:1
                                        temp_within_perms(curr_sub,:,curr_perm) = orig_ord(curr_sub,randperm(out.num_rep_levs));
                                    end
                                end
                            else
                                temp_within_perms = [];
                            end
                            [out.thrAUC.rich_club_neg_bin_nodiscon.beta(con,curr_size,:), ...
                                out.thrAUC.rich_club_neg_bin_nodiscon.test_stat(con,curr_size,:), ...
                                out.thrAUC.rich_club_neg_bin_nodiscon.crit_val(con,curr_size,:), ...
                                out.thrAUC.rich_club_neg_bin_nodiscon.p_1tail(con,curr_size,:), ...
                                out.thrAUC.rich_club_neg_bin_nodiscon.p_2tail(con,curr_size,:), ...
                                ~, ...
                                out.thrAUC.rich_club_neg_bin_nodiscon.hadimag(con,curr_size,:), ...
                                out.thrAUC.rich_club_neg_bin_nodiscon.hadInf(con,curr_size,:), ...
                                permdis, ...
                                out.thrAUC.rich_club_neg_bin_nodiscon.nonperm_p(con,curr_size,:)] = ...
                                GTG_GLM(tempdata(notnan,:), ...
                                curr_full_desmat(notnan,:), ...
                                currcovars_desmat(notnan,:), ...
                                out.Contrast_or_F, ...
                                temp_perms, ...
                                out.HLA_reg_type, ...
                                use_parfor, ...
                                temp_within_perms, ...
                                out.alpha,out.MC_corr);
%                             if out.MC_corr==1
%                                 MC_permdis = min(MC_permdis,permdis);
%                             end
                        end
                    end
                end
            end
        end
        curr_props_calcd = curr_props_calcd+1;
        prog             = curr_props_calcd/tot_props_to_calc;
        progressbar(prog)
    end
    
    if out.MC_corr==1
        MC_permdis = sort(MC_permdis);
    end
    
    fprintf('Done running permutation analyses for contrast %s!\n\n',num2str(con))
    fprintf('Locating significant findings for contrast %s ...\n\n',num2str(con))
    
    if out.num_rep_levs>2
        lev_vals = out.num_rep_levs+1;
    else
        lev_vals = out.num_rep_levs;
    end
    for lev = 1:lev_vals
        if lev_vals>1
            if lev==1
                fprintf(sigeffects_fid,'Significant findings for Contrast #%s, between-subject effect (mean across repeated levels)\n\n',num2str(con));
            else
                fprintf(sigeffects_fid,'Significant findings for Contrast #%s, within-subject (repeated) contrast #%s\n\n',num2str(con),num2str(lev-1));
            end
        else
            fprintf(sigeffects_fid,'Significant findings for Contrast #%s\n\n',num2str(con));
        end
        
        % Network-wide properties
        % Full matrices
        for net_prop = 1:net_full_props_to_calc
            curr_short_name = net_full_short_names{net_prop};
            curr_long_name  = net_full_long_names{net_prop};
            if isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1']) && ~strcmp(curr_short_name,'clust_coef_signed_tot')
                if eval(['out.full.',curr_short_name,'_pos.p_2tail(con,lev)<=out.alpha && ~isnan(out.full.',curr_short_name,'_pos.beta(con,lev)) && out.full.',curr_short_name,'_pos.beta(con,lev) ~= 0'])
                    eval(['out.sig_find(con).',curr_short_name,'_pos_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    if out.MC_corr==1
                        eval(['corrected_p = sum(MC_permdis<=out.full.',curr_short_name,'_pos.nonperm_p(con,lev))/length(MC_permdis);']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_full{2,lev};num2cell([out.full.',curr_short_name,'_pos.beta(con,lev);out.full.',curr_short_name,'_pos.test_stat(con,lev);out.full.',curr_short_name,'_pos.crit_val(con,lev);out.full.',curr_short_name,'_pos.p_2tail(con,lev);corrected_p;out.full.',curr_short_name,'_pos.hadNaN(con);out.full.',curr_short_name,'_pos.hadimag(con);out.full.',curr_short_name,'_pos.hadInf(con)]'')];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, network-wide, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.full.',curr_short_name,'_pos.beta(con,lev),out.full.',curr_short_name,'_pos.test_stat(con,lev),out.full.',curr_short_name,'_pos.crit_val(con,lev),out.full.',curr_short_name,'_pos.p_2tail(con,lev),corrected_p,out.full.',curr_short_name,'_pos.hadNaN(con),out.full.',curr_short_name,'_pos.hadimag(con),out.full.',curr_short_name,'_pos.hadInf(con));']);
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_full{2,lev};num2cell([out.full.',curr_short_name,'_pos.beta(con,lev);out.full.',curr_short_name,'_pos.test_stat(con,lev);out.full.',curr_short_name,'_pos.crit_val(con,lev);out.full.',curr_short_name,'_pos.p_2tail(con,lev);out.full.',curr_short_name,'_pos.hadNaN(con);out.full.',curr_short_name,'_pos.hadimag(con);out.full.',curr_short_name,'_pos.hadInf(con)]'')];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, network-wide, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.full.',curr_short_name,'_pos.beta(con,lev),out.full.',curr_short_name,'_pos.test_stat(con,lev),out.full.',curr_short_name,'_pos.crit_val(con,lev),out.full.',curr_short_name,'_pos.p_2tail(con,lev),out.full.',curr_short_name,'_pos.hadNaN(con),out.full.',curr_short_name,'_pos.hadimag(con),out.full.',curr_short_name,'_pos.hadInf(con));']);
                    end
                end
                if strcmp(out.weightdirec,'Positive and Negative')
                    if eval(['out.full.',curr_short_name,'_neg.p_2tail(con,lev)<=out.alpha && ~isnan(out.full.',curr_short_name,'_neg.beta(con,lev)) && out.full.',curr_short_name,'_neg.beta(con,lev) ~= 0'])
                        eval(['out.sig_find(con).',curr_short_name,'_neg_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['corrected_p = sum(MC_permdis<=out.full.',curr_short_name,'_neg.nonperm_p(con,lev))/length(MC_permdis);']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_full{2,lev};num2cell([out.full.',curr_short_name,'_neg.beta(con,lev);out.full.',curr_short_name,'_neg.test_stat(con,lev);out.full.',curr_short_name,'_neg.crit_val(con,lev);out.full.',curr_short_name,'_neg.p_2tail(con,lev);corrected_p;out.full.',curr_short_name,'_neg.hadNaN(con);out.full.',curr_short_name,'_neg.hadimag(con);out.full.',curr_short_name,'_neg.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, network-wide, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.full.',curr_short_name,'_neg.beta(con,lev),out.full.',curr_short_name,'_neg.test_stat(con,lev),out.full.',curr_short_name,'_neg.crit_val(con,lev),out.full.',curr_short_name,'_neg.p_2tail(con,lev),corrected_p,out.full.',curr_short_name,'_neg.hadNaN(con),out.full.',curr_short_name,'_neg.hadimag(con),out.full.',curr_short_name,'_neg.hadInf(con));']);
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_full{2,lev};num2cell([out.full.',curr_short_name,'_neg.beta(con,lev);out.full.',curr_short_name,'_neg.test_stat(con,lev);out.full.',curr_short_name,'_neg.crit_val(con,lev);out.full.',curr_short_name,'_neg.p_2tail(con,lev);out.full.',curr_short_name,'_neg.hadNaN(con);out.full.',curr_short_name,'_neg.hadimag(con);out.full.',curr_short_name,'_neg.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, network-wide, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.full.',curr_short_name,'_neg.beta(con,lev),out.full.',curr_short_name,'_neg.test_stat(con,lev),out.full.',curr_short_name,'_neg.crit_val(con,lev),out.full.',curr_short_name,'_neg.p_2tail(con,lev),out.full.',curr_short_name,'_neg.hadNaN(con),out.full.',curr_short_name,'_neg.hadimag(con),out.full.',curr_short_name,'_neg.hadInf(con));']);
                        end
                    end
                end
            elseif isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1'])
                if eval(['out.full.',curr_short_name,'.p_2tail(con,lev)<=out.alpha && ~isnan(out.full.',curr_short_name,'.beta(con,lev)) && out.full.',curr_short_name,'.beta(con,lev) ~= 0'])
                    eval(['out.sig_find(con).',curr_short_name,'_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    if out.MC_corr==1
                        eval(['corrected_p = sum(MC_permdis<=out.full.',curr_short_name,'.nonperm_p(con,lev))/length(MC_permdis);']);
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = [out.sig_find(con).',curr_short_name,'_full{2,lev};num2cell([out.full.',curr_short_name,'.beta(con,lev);out.full.',curr_short_name,'.test_stat(con,lev);out.full.',curr_short_name,'.crit_val(con,lev);out.full.',curr_short_name,'.p_2tail(con,lev);corrected_p;out.full.',curr_short_name,'.hadNaN(con);out.full.',curr_short_name,'.hadimag(con);out.full.',curr_short_name,'.hadInf(con)]'')];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, network-wide):\n'');']);
                        fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.full.',curr_short_name,'.beta(con,lev),out.full.',curr_short_name,'.test_stat(con,lev),out.full.',curr_short_name,'.crit_val(con,lev),out.full.',curr_short_name,'.p_2tail(con,lev),corrected_p,out.full.',curr_short_name,'.hadNaN(con),out.full.',curr_short_name,'.hadimag(con),out.full.',curr_short_name,'.hadInf(con));']);
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = [out.sig_find(con).',curr_short_name,'_full{2,lev};num2cell([out.full.',curr_short_name,'.beta(con,lev);out.full.',curr_short_name,'.test_stat(con,lev);out.full.',curr_short_name,'.crit_val(con,lev);out.full.',curr_short_name,'.p_2tail(con,lev);out.full.',curr_short_name,'.hadNaN(con);out.full.',curr_short_name,'.hadimag(con);out.full.',curr_short_name,'.hadInf(con)]'')];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, network-wide):\n'');']);
                        fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.full.',curr_short_name,'.beta(con,lev),out.full.',curr_short_name,'.test_stat(con,lev),out.full.',curr_short_name,'.crit_val(con,lev),out.full.',curr_short_name,'.p_2tail(con,lev),out.full.',curr_short_name,'.hadNaN(con),out.full.',curr_short_name,'.hadimag(con),out.full.',curr_short_name,'.hadInf(con));']);
                    end
                end
            end
        end
        % Thresholded matrices
        for net_prop = 1:net_thr_props_to_calc
            curr_short_name = net_thr_short_names{net_prop};
            curr_long_name  = net_thr_long_names{net_prop};
            if isfield(out.test_props_thrmat,curr_short_name) && eval(['out.test_props_thrmat.',curr_short_name,'==1'])
                if eval(['out.thrAUC.',curr_short_name,'_pos.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_pos.beta(con,lev)) && out.thrAUC.',curr_short_name,'_pos.beta(con,lev) ~= 0'])
                    eval(['out.sig_find(con).',curr_short_name,'_pos_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    if out.MC_corr==1
                        eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos.nonperm_p(con,lev))/length(MC_permdis);']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos.beta(con,lev);out.thrAUC.',curr_short_name,'_pos.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos.hadNaN(con);out.thrAUC.',curr_short_name,'_pos.hadimag(con);out.thrAUC.',curr_short_name,'_pos.hadInf(con)]'')];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos.beta(con,lev),out.thrAUC.',curr_short_name,'_pos.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_pos.hadNaN(con),out.thrAUC.',curr_short_name,'_pos.hadimag(con),out.thrAUC.',curr_short_name,'_pos.hadInf(con));']);
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos.beta(con,lev);out.thrAUC.',curr_short_name,'_pos.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_pos.hadNaN(con);out.thrAUC.',curr_short_name,'_pos.hadimag(con);out.thrAUC.',curr_short_name,'_pos.hadInf(con)]'')];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos.beta(con,lev),out.thrAUC.',curr_short_name,'_pos.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_pos.hadNaN(con),out.thrAUC.',curr_short_name,'_pos.hadimag(con),out.thrAUC.',curr_short_name,'_pos.hadInf(con));']);
                    end
                end
                
                if isfield(out.test_props_thrmat,[curr_short_name,'_pos_bin'])
                    if eval(['out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_pos_bin.beta(con,lev)) && out.thrAUC.',curr_short_name,'_pos_bin.beta(con,lev) ~= 0'])
                        eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos_bin.nonperm_p(con,lev))/length(MC_permdis);']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos_bin.beta(con,lev);out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con);out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con);out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights, binarized matrices):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos_bin.beta(con,lev),out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con));']);
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos_bin.beta(con,lev);out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con);out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con);out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights, binarized matrices):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos_bin.beta(con,lev),out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con));']);
                        end
                    end
                end
                
                if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                    if eval(['out.thrAUC.',curr_short_name,'_neg.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_neg.beta(con,lev)) && out.thrAUC.',curr_short_name,'_neg.beta(con,lev) ~= 0'])
                        eval(['out.sig_find(con).',curr_short_name,'_neg_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg.nonperm_p(con,lev))/length(MC_permdis);']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg.beta(con,lev);out.thrAUC.',curr_short_name,'_neg.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg.hadNaN(con);out.thrAUC.',curr_short_name,'_neg.hadimag(con);out.thrAUC.',curr_short_name,'_neg.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg.beta(con,lev),out.thrAUC.',curr_short_name,'_neg.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_neg.hadNaN(con),out.thrAUC.',curr_short_name,'_neg.hadimag(con),out.thrAUC.',curr_short_name,'_neg.hadInf(con));']);
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg.beta(con,lev);out.thrAUC.',curr_short_name,'_neg.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_neg.hadNaN(con);out.thrAUC.',curr_short_name,'_neg.hadimag(con);out.thrAUC.',curr_short_name,'_neg.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg.beta(con,lev),out.thrAUC.',curr_short_name,'_neg.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_neg.hadNaN(con),out.thrAUC.',curr_short_name,'_neg.hadimag(con),out.thrAUC.',curr_short_name,'_neg.hadInf(con));']);
                        end
                    end
                    
                    if isfield(out.test_props_thrmat,[curr_short_name,'_neg_bin'])
                        if eval(['out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_neg_bin.beta(con,lev)) && out.thrAUC.',curr_short_name,'_neg_bin.beta(con,lev) ~= 0'])
                            eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            if out.MC_corr==1
                                eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg_bin.nonperm_p(con,lev))/length(MC_permdis);']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg_bin.beta(con,lev);out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con);out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con);out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con)]'')];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights, binarized matrices):\n'');']);
                                fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg_bin.beta(con,lev),out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con));']);
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg_bin.beta(con,lev);out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con);out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con);out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con))]'')];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights, binarized matrices):\n'');']);
                                fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg_bin.beta(con,lev),out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con));']);
                            end
                        end
                    end
                end
                if out.calcAUC_nodiscon==1
                    if eval(['out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,lev)) && out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,lev) ~= 0'])
                        eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos_nodiscon.nonperm_p(con,lev))/length(MC_permdis);']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights, excluding disconnected matrices in AUC):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con));']);
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con)]'')];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights, excluding disconnected matrices in AUC):\n'');']);
                            fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con));']);
                        end
                    end
                    
                    if isfield(out.test_props_thrmat,[curr_short_name,'_pos_bin_nodiscon'])
                        if eval(['out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,lev)) && out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,lev) ~= 0'])
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            if out.MC_corr==1
                                eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.nonperm_p(con,lev))/length(MC_permdis);']);
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con)]'')];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con));']);
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con)]'')];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, positive weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con));']);
                            end
                        end
                    end
                    
                    if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                        if eval(['out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,lev)) && out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,lev) ~= 0'])
                            eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            if out.MC_corr==1
                                eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg_nodiscon.nonperm_p(con,lev))/length(MC_permdis);']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con)]'')];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con));']);
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con)]'')];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con));']);
                            end
                        end
                        
                        if isfield(out.test_props_thrmat,[curr_short_name,'_neg_bin_nodiscon'])
                            if eval(['out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,lev)<=out.alpha && ~isnan(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,lev)) && out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,lev) ~= 0'])
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                                if out.MC_corr==1
                                    eval(['corrected_p = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.nonperm_p(con,lev))/length(MC_permdis);']);
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con)]'')];']);
                                    
                                    eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                    fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                    eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,lev),corrected_p,out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con));']);
                                else
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = {''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev};num2cell([out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con)]'')];']);
                                    
                                    eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, network-wide, negative weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                    fprintf(sigeffects_fid,'\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                    eval(['fprintf(sigeffects_fid,''\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n\n'',out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con));']);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % Node properties
        % Full matrices
        for node_prop = 1:node_full_props_to_calc
            curr_short_name = node_full_short_names{node_prop};
            curr_long_name  = node_full_long_names{node_prop};
            if isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1']) && ~strcmp(curr_short_name,'clust_coef_signed')
                eval(['ind = logical(out.full.',curr_short_name,'_pos.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.full.',curr_short_name,'_pos.beta(con,:,lev)) & out.full.',curr_short_name,'_pos.beta(con,:,lev) ~= 0);']);
                if any(ind)
                    eval(['out.sig_find(con).',curr_short_name,'_pos_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    f_ind = find(ind);
                    if out.MC_corr==1
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(f_ind)
                            eval(['corrected_p(ind_val) = sum(MC_permdis<=out.full.',curr_short_name,'_pos.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                            eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.full.',curr_short_name,'_pos.beta(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_pos.test_stat(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_pos.crit_val(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_pos.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.full.',curr_short_name,'_pos.hadNaN(con,f_ind(ind_val)),out.full.',curr_short_name,'_pos.hadimag(con,f_ind(ind_val)),out.full.',curr_short_name,'_pos.hadInf(con,f_ind(ind_val)));']);
                        end
                        
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_full{2,lev};[out.ROI_labels(ind),num2cell([out.full.',curr_short_name,'_pos.beta(con,ind,lev);out.full.',curr_short_name,'_pos.test_stat(con,ind,lev);out.full.',curr_short_name,'_pos.crit_val(con,ind,lev);out.full.',curr_short_name,'_pos.p_2tail(con,ind,lev);corrected_p;out.full.',curr_short_name,'_pos.hadNaN(con,ind);out.full.',curr_short_name,'_pos.hadimag(con,ind);out.full.',curr_short_name,'_pos.hadInf(con,ind)]'')]];']);
                        clear corrected_p
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_full{2,lev};[out.ROI_labels(ind),num2cell([out.full.',curr_short_name,'_pos.beta(con,ind,lev);out.full.',curr_short_name,'_pos.test_stat(con,ind,lev);out.full.',curr_short_name,'_pos.crit_val(con,ind,lev);out.full.',curr_short_name,'_pos.p_2tail(con,ind,lev);out.full.',curr_short_name,'_pos.hadNaN(con,ind);out.full.',curr_short_name,'_pos.hadimag(con,ind);out.full.',curr_short_name,'_pos.hadInf(con,ind)]'')]];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(f_ind)
                            eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.full.',curr_short_name,'_pos.beta(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_pos.test_stat(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_pos.crit_val(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_pos.p_2tail(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_pos.hadNaN(con,f_ind(ind_val)),out.full.',curr_short_name,'_pos.hadimag(con,f_ind(ind_val)),out.full.',curr_short_name,'_pos.hadInf(con,f_ind(ind_val)));']);
                        end
                    end
                    fprintf(sigeffects_fid,'\n');
                end
                if strcmp(out.weightdirec,'Positive and Negative')
                    eval(['ind = logical(out.full.',curr_short_name,'_neg.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.full.',curr_short_name,'_neg.beta(con,:,lev)) & out.full.',curr_short_name,'_neg.beta(con,:,lev) ~= 0);']);
                    if any(ind)
                        eval(['out.sig_find(con).',curr_short_name,'_neg_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        f_ind = find(ind);
                        if out.MC_corr==1
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['corrected_p(ind_val) = sum(MC_permdis<=out.full.',curr_short_name,'_neg.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.full.',curr_short_name,'_neg.beta(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_neg.test_stat(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_neg.crit_val(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_neg.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.full.',curr_short_name,'_neg.hadNaN(con,f_ind(ind_val)),out.full.',curr_short_name,'_neg.hadimag(con,f_ind(ind_val)),out.full.',curr_short_name,'_neg.hadInf(con,f_ind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_full{2,lev};[out.ROI_labels(ind),num2cell([out.full.',curr_short_name,'_neg.beta(con,ind,lev);out.full.',curr_short_name,'_neg.test_stat(con,ind,lev);out.full.',curr_short_name,'_neg.crit_val(con,ind,lev);out.full.',curr_short_name,'_neg.p_2tail(con,ind,lev);corrected_p;out.full.',curr_short_name,'_neg.hadNaN(con,ind);out.full.',curr_short_name,'_neg.hadimag(con,ind);out.full.',curr_short_name,'_neg.hadInf(con,ind)]'')]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_full{2,lev};[out.ROI_labels(ind),num2cell([out.full.',curr_short_name,'_neg.beta(con,ind,lev);out.full.',curr_short_name,'_neg.test_stat(con,ind,lev);out.full.',curr_short_name,'_neg.crit_val(con,ind,lev);out.full.',curr_short_name,'_neg.p_2tail(con,ind,lev);out.full.',curr_short_name,'_neg.hadNaN(con,ind);out.full.',curr_short_name,'_neg.hadimag(con,ind);out.full.',curr_short_name,'_neg.hadInf(con,ind)]'')]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.full.',curr_short_name,'_neg.beta(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_neg.test_stat(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_neg.crit_val(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_neg.p_2tail(con,f_ind(ind_val),lev),out.full.',curr_short_name,'_neg.hadNaN(con,f_ind(ind_val)),out.full.',curr_short_name,'_neg.hadimag(con,f_ind(ind_val)),out.full.',curr_short_name,'_neg.hadInf(con,f_ind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                end
            elseif isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1'])
                eval(['ind = logical(out.full.',curr_short_name,'.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.full.',curr_short_name,'.beta(con,:,lev)) & out.full.',curr_short_name,'.beta(con,:,lev) ~= 0);']);
                if any(ind)
                    eval(['out.sig_find(con).',curr_short_name,'_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    f_ind = find(ind);
                    if out.MC_corr==1
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices):\n'');']);
                        fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(f_ind)
                            eval(['corrected_p(ind_val) = sum(MC_permdis<=out.full.',curr_short_name,'.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                            eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.full.',curr_short_name,'.beta(con,f_ind(ind_val),lev),out.full.',curr_short_name,'.test_stat(con,f_ind(ind_val),lev),out.full.',curr_short_name,'.crit_val(con,f_ind(ind_val),lev),out.full.',curr_short_name,'.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.full.',curr_short_name,'.hadNaN(con,f_ind(ind_val)),out.full.',curr_short_name,'.hadimag(con,f_ind(ind_val)),out.full.',curr_short_name,'.hadInf(con,f_ind(ind_val)));']);
                        end
                        
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = [out.sig_find(con).',curr_short_name,'_full{2,lev};[out.ROI_labels(ind),num2cell([out.full.',curr_short_name,'.beta(con,ind,lev);out.full.',curr_short_name,'.test_stat(con,ind,lev);out.full.',curr_short_name,'.crit_val(con,ind,lev);out.full.',curr_short_name,'.p_2tail(con,ind,lev);corrected_p;out.full.',curr_short_name,'.hadNaN(con,ind);out.full.',curr_short_name,'.hadimag(con,ind);out.full.',curr_short_name,'.hadInf(con,ind)]'')]];']);
                        clear corrected_p
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_full{2,lev} = [out.sig_find(con).',curr_short_name,'_full{2,lev};[out.ROI_labels(ind),num2cell([out.full.',curr_short_name,'.beta(con,ind,lev);out.full.',curr_short_name,'.test_stat(con,ind,lev);out.full.',curr_short_name,'.crit_val(con,ind,lev);out.full.',curr_short_name,'.p_2tail(con,ind,lev);out.full.',curr_short_name,'.hadNaN(con,ind);out.full.',curr_short_name,'.hadimag(con,ind);out.full.',curr_short_name,'.hadInf(con,ind)]'')]];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices):\n'');']);
                        fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(f_ind)
                            eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.full.',curr_short_name,'.beta(con,f_ind(ind_val),lev),out.full.',curr_short_name,'.test_stat(con,f_ind(ind_val),lev),out.full.',curr_short_name,'.crit_val(con,f_ind(ind_val),lev),out.full.',curr_short_name,'.p_2tail(con,f_ind(ind_val),lev),out.full.',curr_short_name,'.hadNaN(con,f_ind(ind_val)),out.full.',curr_short_name,'.hadimag(con,f_ind(ind_val)),out.full.',curr_short_name,'.hadInf(con,f_ind(ind_val)));']);
                        end
                    end
                    fprintf(sigeffects_fid,'\n');
                end
            end
        end
        % Thresholded matrices
        for node_prop = 1:node_thr_props_to_calc
            curr_short_name = node_thr_short_names{node_prop};
            curr_long_name  = node_thr_long_names{node_prop};
            if isfield(out.test_props_thrmat,curr_short_name) && eval(['out.test_props_thrmat.',curr_short_name,'==1'])
                eval(['ind = logical(out.thrAUC.',curr_short_name,'_pos.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_pos.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_pos.beta(con,:,lev) ~= 0);']);
                if any(ind)
                    eval(['out.sig_find(con).',curr_short_name,'_pos_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    f_ind = find(ind);
                    if out.MC_corr==1
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(f_ind)
                            eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                            eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadInf(con,f_ind(ind_val)));']);
                        end
                        
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos.hadInf(con,ind)]'')]];']);
                        clear corrected_p
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_pos.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos.hadInf(con,ind)]'')]];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(f_ind)
                            eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadInf(con,f_ind(ind_val)));']);
                        end
                    end
                    fprintf(sigeffects_fid,'\n');
                end
                
                if isfield(out.test_props_thrmat,[curr_short_name,'_pos_bin'])
                    eval(['ind = logical(out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_pos_bin.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_pos_bin.beta(con,:,lev) ~= 0);']);
                    if any(ind)
                        eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        f_ind = find(ind);
                        if out.MC_corr==1
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized matrices):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos_bin.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos_bin.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,f_ind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos_bin.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,ind)]'')]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos_bin.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,ind)]'')]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized matrices):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos_bin.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,f_ind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                end
                
                if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                    eval(['ind = logical(out.thrAUC.',curr_short_name,'_neg.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_neg.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_neg.beta(con,:,lev) ~= 0);']);
                    if any(ind)
                        eval(['out.sig_find(con).',curr_short_name,'_neg_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        f_ind = find(ind);
                        if out.MC_corr==1
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadInf(con,f_ind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg.hadInf(con,ind)]'')]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_neg.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg.hadInf(con,ind)]'')]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadInf(con,f_ind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                    
                    if isfield(out.test_props_thrmat,[curr_short_name,'_neg_bin'])
                        eval(['ind = logical(out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_neg_bin.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_neg_bin.beta(con,:,lev) ~= 0);']);
                        if any(ind)
                            eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            f_ind = find(ind);
                            if out.MC_corr==1
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized matrices):\n'');']);
                                fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(f_ind)
                                    eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg_bin.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg_bin.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,f_ind(ind_val)));']);
                                end
                                
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg_bin.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,ind)]'')]];']);
                                clear corrected_p
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg_bin.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,ind)]'')]];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized matrices):\n'');']);
                                fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(f_ind)
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg_bin.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,f_ind(ind_val)));']);
                                end
                            end
                            fprintf(sigeffects_fid,'\n');
                        end
                    end
                end
                
                if out.calcAUC_nodiscon==1
                    eval(['ind = logical(out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,:,lev) ~= 0);']);
                    if any(ind)
                        eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        f_ind = find(ind);
                        if out.MC_corr==1
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, excluding disconnected matrices in AUC):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos_nodiscon.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,f_ind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,ind)]'')]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,ind)]'')]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, excluding disconnected matrices in AUC):\n'');']);
                            fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(f_ind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,f_ind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                    
                    if isfield(out.test_props_thrmat,[curr_short_name,'_pos_bin_nodiscon'])
                        eval(['ind = logical(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,:,lev) ~= 0);']);
                        if any(ind)
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            f_ind = find(ind);
                            if out.MC_corr==1
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(f_ind)
                                    eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,f_ind(ind_val)));']);
                                end
                                
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,ind)]'')]];']);
                                clear corrected_p
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,ind)]'')]];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(f_ind)
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,f_ind(ind_val)));']);
                                end
                            end
                            fprintf(sigeffects_fid,'\n');
                        end
                    end
                    
                    if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                        eval(['ind = logical(out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,:,lev) ~= 0);']);
                        if any(ind)
                            eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            f_ind = find(ind);
                            if out.MC_corr==1
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(f_ind)
                                    eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg_nodiscon.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,f_ind(ind_val)));']);
                                end
                                
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,ind)]'')]];']);
                                clear corrected_p
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,ind)]'')]];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(f_ind)
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,f_ind(ind_val)));']);
                                end
                            end
                            fprintf(sigeffects_fid,'\n');
                        end
                        
                        if isfield(out.test_props_thrmat,[curr_short_name,'_neg_bin_nodiscon'])
                            eval(['ind = logical(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,:,lev)) & out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,:,lev) ~= 0);']);
                            if any(ind)
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                                f_ind = find(ind);
                                if out.MC_corr==1
                                    eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                    fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                    for ind_val = 1:length(f_ind)
                                        eval(['corrected_p(ind_val) = sum(MC_permdis<=out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.nonperm_p(con,f_ind(ind_val),lev))/length(MC_permdis);']);
                                        eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,f_ind(ind_val),lev),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,f_ind(ind_val)));']);
                                    end
                                    
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,ind,lev);corrected_p;out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,ind)]'')]];']);
                                    clear corrected_p
                                else
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = {''node'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev};[out.ROI_labels(ind),num2cell([out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,ind,lev);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,ind);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,ind);out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,ind)]'')]];']);
                                    
                                    eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized matrices, excluding disconnected matrices in AUC):\n'');']);
                                    fprintf(sigeffects_fid,'\tnode\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                    for ind_val = 1:length(f_ind)
                                        eval(['fprintf(sigeffects_fid,''\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{f_ind(ind_val)},out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,f_ind(ind_val)));']);
                                    end
                                end
                                fprintf(sigeffects_fid,'\n');
                            end
                        end
                    end
                end
            end
        end
        
        % Edge properties
        % Full matrices
        for edge_prop = 1:edge_full_props_to_calc
            curr_short_name = edge_full_short_names{edge_prop};
            curr_long_name  = edge_full_long_names{edge_prop};
            if isfield(out.test_props_fullmat,curr_short_name) && eval(['out.test_props_fullmat.',curr_short_name,'==1'])
                eval(['ind = logical((squeeze(out.full.',curr_short_name,'_pos.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.full.',curr_short_name,'_pos.beta(con,:,:,lev))) & squeeze(out.full.',curr_short_name,'_pos.beta(con,:,:,lev) ~= 0));']);
                ind(logical(eye(size(ind)))) = 0; %#ok<*AGROW>
                [xind,yind]                  = find(ind==1); %#ok<*NASGU>
                if any(any(ind))
                    eval(['out.sig_find(con).',curr_short_name,'_pos_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    if out.MC_corr==1
                        eval(['tempbeta     = squeeze(out.full.',curr_short_name,'_pos.beta(con,:,:,lev));']);
                        eval(['tempteststat = squeeze(out.full.',curr_short_name,'_pos.test_stat(con,:,:,lev));']);
                        eval(['tempcritval  = squeeze(out.full.',curr_short_name,'_pos.crit_val(con,:,:,lev));']);
                        eval(['tempp        = squeeze(out.full.',curr_short_name,'_pos.p_2tail(con,:,:,lev));']);
                        eval(['tempnonpermp = squeeze(out.full.',curr_short_name,'_pos.nonperm_p(con,:,:,lev));']);
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(xind)
                            corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                            eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{xind(ind_val)},out.ROI_labels{yind(ind_val)},tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.full.',curr_short_name,'_pos.hadNaN(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_pos.hadimag(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_pos.hadInf(con,xind(ind_val),yind(ind_val)));']);
                        end
                        
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_full{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.full.',curr_short_name,'_pos.hadNaN(con,ind)'',out.full.',curr_short_name,'_pos.hadimag(con,ind)'',out.full.',curr_short_name,'_pos.hadInf(con,ind)''])]];']);
                        clear corrected_p
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['tempbeta                                              = squeeze(out.full.',curr_short_name,'_pos.beta(con,:,:,lev));']);
                        eval(['tempteststat                                          = squeeze(out.full.',curr_short_name,'_pos.test_stat(con,:,:,lev));']);
                        eval(['tempcritval                                           = squeeze(out.full.',curr_short_name,'_pos.crit_val(con,:,:,lev));']);
                        eval(['tempp                                                 = squeeze(out.full.',curr_short_name,'_pos.p_2tail(con,:,:,lev));']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_full{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_full{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.full.',curr_short_name,'_pos.hadNaN(con,ind)'',out.full.',curr_short_name,'_pos.hadimag(con,ind)'',out.full.',curr_short_name,'_pos.hadInf(con,ind)''])]];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(xind)
                            eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{xind(ind_val)},out.ROI_labels{yind(ind_val)},tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_pos.hadNaN(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_pos.hadimag(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_pos.hadInf(con,xind(ind_val),yind(ind_val)));']);
                        end
                    end
                    fprintf(sigeffects_fid,'\n');
                end
                if strcmp(out.weightdirec,'Positive and Negative')
                    eval(['ind = logical((squeeze(out.full.',curr_short_name,'_neg.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.full.',curr_short_name,'_neg.beta(con,:,:,lev))) & squeeze(out.full.',curr_short_name,'_neg.beta(con,:,:,lev) ~= 0));']);
                    ind(logical(eye(size(ind)))) = 0;
                    [xind,yind]                  = find(ind==1);
                    if any(any(ind))
                        eval(['out.sig_find(con).',curr_short_name,'_neg_full{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['tempbeta     = squeeze(out.full.',curr_short_name,'_neg.beta(con,:,:,lev));']);
                            eval(['tempteststat = squeeze(out.full.',curr_short_name,'_neg.test_stat(con,:,:,lev));']);
                            eval(['tempcritval  = squeeze(out.full.',curr_short_name,'_neg.crit_val(con,:,:,lev));']);
                            eval(['tempp        = squeeze(out.full.',curr_short_name,'_neg.p_2tail(con,:,:,lev));']);
                            eval(['tempnonpermp = squeeze(out.full.',curr_short_name,'_neg.nonperm_p(con,:,:,lev));']);
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{xind(ind_val)},out.ROI_labels{yind(ind_val)},tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.full.',curr_short_name,'_neg.hadNaN(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_neg.hadimag(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_neg.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_full{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.full.',curr_short_name,'_neg.hadNaN(con,ind)'',out.full.',curr_short_name,'_neg.hadimag(con,ind)'',out.full.',curr_short_name,'_neg.hadInf(con,ind)''])]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['tempbeta                                              = squeeze(out.full.',curr_short_name,'_neg.beta(con,:,:,lev));']);
                            eval(['tempteststat                                          = squeeze(out.full.',curr_short_name,'_neg.test_stat(con,:,:,lev));']);
                            eval(['tempcritval                                           = squeeze(out.full.',curr_short_name,'_neg.crit_val(con,:,:,lev));']);
                            eval(['tempp                                                 = squeeze(out.full.',curr_short_name,'_neg.p_2tail(con,:,:,lev));']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_full{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_full{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.full.',curr_short_name,'_neg.hadNaN(con,ind)'',out.full.',curr_short_name,'_neg.hadimag(con,ind)'',out.full.',curr_short_name,'_neg.hadInf(con,ind)''])]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (full matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels{xind(ind_val)},out.ROI_labels{yind(ind_val)},tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_neg.hadNaN(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_neg.hadimag(con,xind(ind_val),yind(ind_val)),out.full.',curr_short_name,'_neg.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                end
            end
        end
        % Thresholded matrices
        for edge_prop = 1:edge_thr_props_to_calc
            curr_short_name = edge_thr_short_names{edge_prop};
            curr_long_name = edge_thr_long_names{edge_prop};
            if isfield(out.test_props_thrmat,curr_short_name) && eval(['out.test_props_thrmat.',curr_short_name,'==1'])
                eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_pos.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_pos.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_pos.beta(con,:,:,lev) ~= 0));']);
                ind(logical(eye(size(ind)))) = 0;
                [xind,yind]                  = find(ind==1);
                if any(any(ind))
                    eval(['out.sig_find(con).',curr_short_name,'_pos_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                    if out.MC_corr==1
                        eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_pos.beta(con,:,:,lev));']);
                        eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_pos.test_stat(con,:,:,lev));']);
                        eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_pos.crit_val(con,:,:,lev));']);
                        eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_pos.p_2tail(con,:,:,lev));']);
                        eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_pos.nonperm_p(con,:,:,lev));']);
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(xind)
                            corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                            eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadInf(con,xind(ind_val),yind(ind_val)));']);
                        end
                        
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_pos.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos.hadInf(con,ind)''])]];']);
                        clear corrected_p
                    else
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                        eval(['tempbeta                                             = squeeze(out.thrAUC.',curr_short_name,'_pos.beta(con,:,:,lev));']);
                        eval(['tempteststat                                         = squeeze(out.thrAUC.',curr_short_name,'_pos.test_stat(con,:,:,lev));']);
                        eval(['tempcritval                                          = squeeze(out.thrAUC.',curr_short_name,'_pos.crit_val(con,:,:,lev));']);
                        eval(['tempp                                                = squeeze(out.thrAUC.',curr_short_name,'_pos.p_2tail(con,:,:,lev));']);
                        eval(['out.sig_find(con).',curr_short_name,'_pos_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos.hadInf(con,ind)''])]];']);
                        
                        eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights):\n'');']);
                        fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        for ind_val = 1:length(xind)
                            eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos.hadInf(con,xind(ind_val),yind(ind_val)));']);
                        end
                    end
                    fprintf(sigeffects_fid,'\n');
                end
                
                if isfield(out.test_props_thrmat,[curr_short_name,'_pos_bin'])
                    eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_pos_bin.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_pos_bin.beta(con,:,:,lev) ~= 0));']);
                    ind(logical(eye(size(ind)))) = 0;
                    [xind,yind]                  = find(ind==1);
                    if any(any(ind))
                        eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.beta(con,:,:,lev));']);
                            eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,:,:,lev));']);
                            eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,:,:,lev));']);
                            eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,:,:,lev));']);
                            eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.nonperm_p(con,:,:,lev));']);
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,ind)''])]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['tempbeta                                                 = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.beta(con,:,:,lev));']);
                            eval(['tempteststat                                             = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.test_stat(con,:,:,lev));']);
                            eval(['tempcritval                                              = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.crit_val(con,:,:,lev));']);
                            eval(['tempp                                                    = squeeze(out.thrAUC.',curr_short_name,'_pos_bin.p_2tail(con,:,:,lev));']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,ind)''])]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                end
                
                if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                    eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_neg.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_neg.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_neg.beta(con,:,:,lev) ~= 0));']);
                    ind(logical(eye(size(ind)))) = 0;
                    [xind,yind]                  = find(ind==1);
                    if any(any(ind))
                        eval(['out.sig_find(con).',curr_short_name,'_neg_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_neg.beta(con,:,:,lev));']);
                            eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_neg.test_stat(con,:,:,lev));']);
                            eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_neg.crit_val(con,:,:,lev));']);
                            eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_neg.p_2tail(con,:,:,lev));']);
                            eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_neg.nonperm_p(con,:,:,lev));']);
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_neg.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg.hadInf(con,ind)''])]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['tempbeta                                             = squeeze(out.thrAUC.',curr_short_name,'_neg.beta(con,:,:,lev));']);
                            eval(['tempteststat                                         = squeeze(out.thrAUC.',curr_short_name,'_neg.test_stat(con,:,:,lev));']);
                            eval(['tempcritval                                          = squeeze(out.thrAUC.',curr_short_name,'_neg.crit_val(con,:,:,lev));']);
                            eval(['tempp                                                = squeeze(out.thrAUC.',curr_short_name,'_neg.p_2tail(con,:,:,lev));']);
                            eval(['out.sig_find(con).',curr_short_name,'_neg_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg.hadInf(con,ind)''])]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                    
                    if isfield(out.test_props_thrmat,[curr_short_name,'_neg_bin'])
                        eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_neg_bin.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_neg_bin.beta(con,:,:,lev) ~= 0));']);
                        ind(logical(eye(size(ind)))) = 0;
                        [xind,yind]                  = find(ind==1);
                        if any(any(ind))
                            eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            if out.MC_corr==1
                                eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.beta(con,:,:,lev));']);
                                eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,:,:,lev));']);
                                eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,:,:,lev));']);
                                eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,:,:,lev));']);
                                eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.nonperm_p(con,:,:,lev));']);
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized weights):\n'');']);
                                fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(xind)
                                    corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                end
                                
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,ind)''])]];']);
                                clear corrected_p
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['tempbeta                                                 = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.beta(con,:,:,lev));']);
                                eval(['tempteststat                                             = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.test_stat(con,:,:,lev));']);
                                eval(['tempcritval                                              = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.crit_val(con,:,:,lev));']);
                                eval(['tempp                                                    = squeeze(out.thrAUC.',curr_short_name,'_neg_bin.p_2tail(con,:,:,lev));']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,ind)''])]];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized weights):\n'');']);
                                fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(xind)
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                end
                            end
                            fprintf(sigeffects_fid,'\n');
                        end
                    end
                end
                if out.calcAUC_nodiscon==1
                    eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,:,:,lev) ~= 0));']);
                    ind(logical(eye(size(ind)))) = 0;
                    [xind,yind]                  = find(ind==1);
                    if any(any(ind))
                        eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                        if out.MC_corr==1
                            eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,:,:,lev));']);
                            eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,:,:,lev));']);
                            eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,:,:,lev));']);
                            eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,:,:,lev));']);
                            eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.nonperm_p(con,:,:,lev));']);
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, excluding disconnected matrices in AUC):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                            
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,ind)''])]];']);
                            clear corrected_p
                        else
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                            eval(['tempbeta                                                      = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.beta(con,:,:,lev));']);
                            eval(['tempteststat                                                  = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.test_stat(con,:,:,lev));']);
                            eval(['tempcritval                                                   = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.crit_val(con,:,:,lev));']);
                            eval(['tempp                                                         = squeeze(out.thrAUC.',curr_short_name,'_pos_nodiscon.p_2tail(con,:,:,lev));']);
                            eval(['out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,ind)''])]];']);
                            
                            eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, excluding disconnected matrices in AUC):\n'');']);
                            fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            for ind_val = 1:length(xind)
                                eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                            end
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                    
                    if isfield(out.test_props_thrmat,[curr_short_name,'_pos_bin_nodiscon'])
                        eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,:,:,lev) ~= 0));']);
                        ind(logical(eye(size(ind)))) = 0;
                        [xind,yind]                  = find(ind==1);
                        if any(any(ind))
                            eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            if out.MC_corr==1
                                eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,:,:,lev));']);
                                eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,:,:,lev));']);
                                eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,:,:,lev));']);
                                eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,:,:,lev));']);
                                eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.nonperm_p(con,:,:,lev));']);
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(xind)
                                    corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                end
                                
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,ind)''])]];']);
                                clear corrected_p
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['tempbeta                                                          = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.beta(con,:,:,lev));']);
                                eval(['tempteststat                                                      = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.test_stat(con,:,:,lev));']);
                                eval(['tempcritval                                                       = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.crit_val(con,:,:,lev));']);
                                eval(['tempp                                                             = squeeze(out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.p_2tail(con,:,:,lev));']);
                                eval(['out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_pos_bin_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,ind)''])]];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, positive weights, binarized weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(xind)
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_pos_bin_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                end
                            end
                            fprintf(sigeffects_fid,'\n');
                        end
                    end
                    
                    if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                        eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,:,:,lev) ~= 0));']);
                        ind(logical(eye(size(ind)))) = 0;
                        [xind,yind]                  = find(ind==1);
                        if any(any(ind))
                            eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                            if out.MC_corr==1
                                eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,:,:,lev));']);
                                eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,:,:,lev));']);
                                eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,:,:,lev));']);
                                eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,:,:,lev));']);
                                eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.nonperm_p(con,:,:,lev));']);
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(xind)
                                    corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                end
                                
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,ind)''])]];']);
                                clear corrected_p
                            else
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                eval(['tempbeta                                                      = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.beta(con,:,:,lev));']);
                                eval(['tempteststat                                                  = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.test_stat(con,:,:,lev));']);
                                eval(['tempcritval                                                   = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.crit_val(con,:,:,lev));']);
                                eval(['tempp                                                         = squeeze(out.thrAUC.',curr_short_name,'_neg_nodiscon.p_2tail(con,:,:,lev));']);
                                eval(['out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,ind)''])]];']);
                                
                                eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, excluding disconnected matrices in AUC):\n'');']);
                                fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                for ind_val = 1:length(xind)
                                    eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                end
                            end
                            fprintf(sigeffects_fid,'\n');
                        end
                        
                        if isfield(out.test_props_thrmat,[curr_short_name,'_neg_bin_nodiscon'])
                            eval(['ind = logical((squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,:,:,lev))<=out.alpha) & ~isnan(squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,:,:,lev))) & squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,:,:,lev) ~= 0));']);
                            ind(logical(eye(size(ind)))) = 0;
                            [xind,yind]                  = find(ind==1);
                            if any(any(ind))
                                eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{1,lev} = [''Contrast/F-test #'',num2str(con)];']);
                                if out.MC_corr==1
                                    eval(['tempbeta     = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,:,:,lev));']);
                                    eval(['tempteststat = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,:,:,lev));']);
                                    eval(['tempcritval  = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,:,:,lev));']);
                                    eval(['tempp        = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,:,:,lev));']);
                                    eval(['tempnonpermp = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.nonperm_p(con,:,:,lev));']);
                                    eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized weights, excluding disconnected matrices in AUC):\n'');']);
                                    fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\tcorrected p\thad NaN\thad imag\thad Inf\n');
                                    for ind_val = 1:length(xind)
                                        corrected_p(ind_val) = sum(MC_permdis<=tempnonpermp(xind(ind_val),yind(ind_val)))/length(MC_permdis);
                                        eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),corrected_p(ind_val),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                    end
                                    
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''corrected_p'',''had NaN'',''had imag'',''had Inf''};']);
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),corrected_p,out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,ind)''])]];']);
                                    clear corrected_p
                                else
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = {''node_1'',''node_2'',''beta'',''stat'',''crit_val'',''p'',''had NaN'',''had imag'',''had Inf''};']);
                                    eval(['tempbeta                                                          = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.beta(con,:,:,lev));']);
                                    eval(['tempteststat                                                      = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.test_stat(con,:,:,lev));']);
                                    eval(['tempcritval                                                       = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.crit_val(con,:,:,lev));']);
                                    eval(['tempp                                                             = squeeze(out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.p_2tail(con,:,:,lev));']);
                                    eval(['out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev} = [out.sig_find(con).',curr_short_name,'_neg_bin_nodiscon_thr{2,lev};[out.ROI_labels(xind),out.ROI_labels(yind),num2cell([tempbeta(ind),tempteststat(ind),tempcritval(ind),tempp(ind),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,ind)'',out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,ind)''])]];']);
                                    
                                    eval(['fprintf(sigeffects_fid,''',curr_long_name,' (thresholded matrices, negative weights, binarized weights, excluding disconnected matrices in AUC):\n'');']);
                                    fprintf(sigeffects_fid,'\tnode_1\ttnode_2\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                                    for ind_val = 1:length(xind)
                                        eval(['fprintf(sigeffects_fid,''\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n'',out.ROI_labels(xind(ind_val)),out.ROI_labels(yind(ind_val)),tempbeta(xind(ind_val),yind(ind_val)),tempteststat(xind(ind_val),yind(ind_val)),tempcritval(xind(ind_val),yind(ind_val)),tempp(xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadNaN(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadimag(con,xind(ind_val),yind(ind_val)),out.thrAUC.',curr_short_name,'_neg_bin_nodiscon.hadInf(con,xind(ind_val),yind(ind_val)));']);
                                    end
                                end
                                fprintf(sigeffects_fid,'\n');
                            end
                        end
                    end
                end
            end
        end
        
        % Rich Club Networks
        if out.test_props_fullmat.rich_club==1
            ind = find(out.full.rich_club_pos.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.full.rich_club_pos.beta(con,:,lev)) & out.full.rich_club_pos.beta(con,:,lev) ~= 0);
            if any(ind)
                out.sig_find(con).rich_club_pos_full{1,lev} = ['Contrast/F-test #',num2str(con)];
                out.sig_find(con).rich_club_pos_full{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                out.sig_find(con).rich_club_pos_full{2,lev} = [out.sig_find(con).rich_club_pos_full{2,lev};num2cell([ind;out.full.rich_club_pos.beta(con,ind,lev);out.full.rich_club_pos.test_stat(con,ind,lev);out.full.rich_club_pos.crit_val(con,ind,lev);out.full.rich_club_pos.p_2tail(con,ind,lev);out.full.rich_club_pos.hadNaN(con,ind);out.full.rich_club_pos.hadimag(con,ind);out.full.rich_club_pos.hadInf(con,ind)]')];
                
                fprintf(sigeffects_fid,'Rich Club Networks (full matrices, positive weights):\n');
                fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                f_ind = find(ind);
                for ind_val = 1:length(f_ind)
                    fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.full.rich_club_pos.beta(con,f_ind(ind_val),lev),out.full.rich_club_pos.test_stat(con,f_ind(ind_val),lev),out.full.rich_club_pos.crit_val(con,f_ind(ind_val),lev),out.full.rich_club_pos.p_2tail(con,f_ind(ind_val),lev),out.full.rich_club_pos.hadNaN(con,f_ind(ind_val)),out.full.rich_club_pos.hadimag(con,f_ind(ind_val)),out.full.rich_club_pos.hadInf(con,f_ind(ind_val)));
                end
                fprintf(sigeffects_fid,'\n');
            end
            if strcmp(out.weightdirec,'Positive and Negative')
                ind = find(out.full.rich_club_neg.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.full.rich_club_neg.beta(con,:,lev)) & out.full.rich_club_neg.beta(con,:,lev) ~= 0);
                if any(ind)
                    out.sig_find(con).rich_club_neg_full{1,lev} = ['Contrast/F-test #',num2str(con)];
                    out.sig_find(con).rich_club_neg_full{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                    out.sig_find(con).rich_club_neg_full{2,lev} = [out.sig_find(con).rich_club_neg_full{2,lev};num2cell([ind;out.full.rich_club_neg.beta(con,ind,lev);out.full.rich_club_neg.test_stat(con,ind,lev);out.full.rich_club_neg.crit_val(con,ind,lev);out.full.rich_club_neg.p_2tail(con,ind,lev);out.full.rich_club_neg.hadNaN(con,ind);out.full.rich_club_neg.hadimag(con,ind);out.full.rich_club_neg.hadInf(con,ind)]')];
                    
                    fprintf(sigeffects_fid,'Rich Club Networks (full matrices, negative weights):\n');
                    fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                    f_ind = find(ind);
                    for ind_val = 1:length(f_ind)
                        fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.full.rich_club_neg.beta(con,f_ind(ind_val),lev),out.full.rich_club_neg.test_stat(con,f_ind(ind_val),lev),out.full.rich_club_neg.crit_val(con,f_ind(ind_val),lev),out.full.rich_club_neg.p_2tail(con,f_ind(ind_val),lev),out.full.rich_club_neg.hadNaN(con,f_ind(ind_val)),out.full.rich_club_neg.hadimag(con,f_ind(ind_val)),out.full.rich_club_neg.hadInf(con,f_ind(ind_val)));
                    end
                    fprintf(sigeffects_fid,'\n');
                end
            end
        end
        if out.test_props_thrmat.rich_club==1
            ind = find(out.thrAUC.rich_club_pos.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_pos.beta(con,:,lev)) & out.thrAUC.rich_club_pos.beta(con,:,lev) ~= 0);
            if any(ind)
                out.sig_find(con).rich_club_pos_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                out.sig_find(con).rich_club_pos_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                out.sig_find(con).rich_club_pos_thr{2,lev} = [out.sig_find(con).rich_club_pos_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_pos.beta(con,ind,lev);out.thrAUC.rich_club_pos.test_stat(con,ind,lev);out.thrAUC.rich_club_pos.crit_val(con,ind,lev);out.thrAUC.rich_club_pos.p_2tail(con,ind,lev);out.thrAUC.rich_club_pos.hadNaN(con,ind);out.thrAUC.rich_club_pos.hadimag(con,ind);out.thrAUC.rich_club_pos.hadInf(con,ind)]')];
                
                fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, positive weights):\n');
                fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                f_ind = find(ind);
                for ind_val = 1:length(f_ind)
                    fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_pos.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_pos.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_pos.hadInf(con,f_ind(ind_val)));
                end
                fprintf(sigeffects_fid,'\n');
            end
            
            if out.calcbinthresh==1
                ind = find(out.thrAUC.rich_club_pos_bin.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_pos_bin.beta(con,:,lev)) & out.thrAUC.rich_club_pos_bin.beta(con,:,lev) ~= 0);
                if any(ind)
                    out.sig_find(con).rich_club_pos_bin_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                    out.sig_find(con).rich_club_pos_bin_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                    out.sig_find(con).rich_club_pos_bin_thr{2,lev} = [out.sig_find(con).rich_club_pos_bin_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_pos_bin.beta(con,ind,lev);out.thrAUC.rich_club_pos_bin.test_stat(con,ind,lev);out.thrAUC.rich_club_pos_bin.crit_val(con,ind,lev);out.thrAUC.rich_club_pos_bin.p_2tail(con,ind,lev);out.thrAUC.rich_club_pos_bin.hadNaN(con,ind);out.thrAUC.rich_club_pos_bin.hadimag(con,ind);out.thrAUC.rich_club_pos_bin.hadInf(con,ind)]')];
                    
                    fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, binarized matrices, positive weights):\n');
                    fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                    f_ind = find(ind);
                    for ind_val = 1:length(f_ind)
                        fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_pos_bin.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_pos_bin.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_pos_bin.hadInf(con,f_ind(ind_val)));
                    end
                    fprintf(sigeffects_fid,'\n');
                end
            end
            
            if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                ind = find(out.thrAUC.rich_club_neg.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_neg.beta(con,:,lev)) & out.thrAUC.rich_club_neg.beta(con,:,lev) ~= 0);
                if any(ind)
                    out.sig_find(con).rich_club_neg_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                    out.sig_find(con).rich_club_neg_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                    out.sig_find(con).rich_club_neg_thr{2,lev} = [out.sig_find(con).rich_club_neg_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_neg.beta(con,ind,lev);out.thrAUC.rich_club_neg.test_stat(con,ind,lev);out.thrAUC.rich_club_neg.crit_val(con,ind,lev);out.thrAUC.rich_club_neg.p_2tail(con,ind,lev);out.thrAUC.rich_club_neg.hadNaN(con,ind);out.thrAUC.rich_club_neg.hadimag(con,ind);out.thrAUC.rich_club_neg.hadInf(con,ind)]')];
                    
                    fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, negative weights):\n');
                    fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                    f_ind = find(ind);
                    for ind_val = 1:length(f_ind)
                        fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_neg.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_neg.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_neg.hadInf(con,f_ind(ind_val)));
                    end
                    fprintf(sigeffects_fid,'\n');
                end
                
                if out.calcbinthresh==1
                    ind = find(out.thrAUC.rich_club_neg_bin.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_neg_bin.beta(con,:,lev)) & out.thrAUC.rich_club_neg_bin.beta(con,:,lev) ~= 0);
                    if any(ind)
                        out.sig_find(con).rich_club_neg_bin_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                        out.sig_find(con).rich_club_neg_bin_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                        out.sig_find(con).rich_club_neg_bin_thr{2,lev} = [out.sig_find(con).rich_club_neg_bin_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_neg_bin.beta(con,ind,lev);out.thrAUC.rich_club_neg_bin.test_stat(con,ind,lev);out.thrAUC.rich_club_neg_bin.crit_val(con,ind,lev);out.thrAUC.rich_club_neg_bin.p_2tail(con,ind,lev);out.thrAUC.rich_club_neg_bin.hadNaN(con,ind);out.thrAUC.rich_club_neg_bin.hadimag(con,ind);out.thrAUC.rich_club_neg_bin.hadInf(con,ind)]')];
                        
                        fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, binarized matrices, negative weights):\n');
                        fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        f_ind = find(ind);
                        for ind_val = 1:length(f_ind)
                            fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_neg_bin.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_neg_bin.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_neg_bin.hadInf(con,f_ind(ind_val)));
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                end
            end
            if out.calcAUC_nodiscon==1
                ind = find(out.thrAUC.rich_club_pos_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_pos_nodiscon.beta(con,:,lev)) & out.thrAUC.rich_club_pos_nodiscon.beta(con,:,lev) ~= 0);
                if any(ind)
                    out.sig_find(con).rich_club_pos_nodiscon_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                    out.sig_find(con).rich_club_pos_nodiscon_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                    out.sig_find(con).rich_club_pos_nodiscon_thr{2,lev} = [out.sig_find(con).rich_club_pos_nodiscon_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_pos_nodiscon.beta(con,ind,lev);out.thrAUC.rich_club_pos_nodiscon.test_stat(con,ind,lev);out.thrAUC.rich_club_pos_nodiscon.crit_val(con,ind,lev);out.thrAUC.rich_club_pos_nodiscon.p_2tail(con,ind,lev);out.thrAUC.rich_club_pos_nodiscon.hadNaN(con,ind);out.thrAUC.rich_club_pos_nodiscon.hadimag(con,ind);out.thrAUC.rich_club_pos_nodiscon.hadInf(con,ind)]')];
                    
                    fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, positive weights, excluding disconnected matrices in AUC):\n');
                    fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                    f_ind = find(ind);
                    for ind_val = 1:length(f_ind)
                        fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_pos_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_pos_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_pos_nodiscon.hadInf(con,f_ind(ind_val)));
                    end
                    fprintf(sigeffects_fid,'\n');
                end
                
                if out.calcbinthresh==1
                    ind = find(out.thrAUC.rich_club_pos_bin_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_pos_bin_nodiscon.beta(con,:,lev)) & out.thrAUC.rich_club_pos_bin_nodiscon.beta(con,:,lev) ~= 0);
                    if any(ind)
                        out.sig_find(con).rich_club_pos_bin_nodiscon_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                        out.sig_find(con).rich_club_pos_bin_nodiscon_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                        out.sig_find(con).rich_club_pos_bin_nodiscon_thr{2,lev} = [out.sig_find(con).rich_club_pos_bin_nodiscon_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_pos_bin_nodiscon.beta(con,ind,lev);out.thrAUC.rich_club_pos_bin_nodiscon.test_stat(con,ind,lev);out.thrAUC.rich_club_pos_bin_nodiscon.crit_val(con,ind,lev);out.thrAUC.rich_club_pos_bin_nodiscon.p_2tail(con,ind,lev);out.thrAUC.rich_club_pos_bin_nodiscon.hadNaN(con,ind);out.thrAUC.rich_club_pos_bin_nodiscon.hadimag(con,ind);out.thrAUC.rich_club_pos_bin_nodiscon.hadInf(con,ind)]')];
                        
                        fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, positive weights, binarized matrices, excluding disconnected matrices in AUC):\n');
                        fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        f_ind = find(ind);
                        for ind_val = 1:length(f_ind)
                            fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_pos_bin_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_pos_bin_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_pos_bin_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_pos_bin_nodiscon.hadInf(con,f_ind(ind_val)));
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                end
                
                if strcmp(out.weightdirec,'Positive and Negative') && out.neg_mindens_nan==0
                    ind = find(out.thrAUC.rich_club_neg_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_neg_nodiscon.beta(con,:,lev)) & out.thrAUC.rich_club_neg_nodiscon.beta(con,:,lev) ~= 0);
                    if any(ind)
                        out.sig_find(con).rich_club_neg_nodiscon_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                        out.sig_find(con).rich_club_neg_nodiscon_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                        out.sig_find(con).rich_club_neg_nodiscon_thr{2,lev} = [out.sig_find(con).rich_club_neg_nodiscon_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_neg_nodiscon.beta(con,ind,lev);out.thrAUC.rich_club_neg_nodiscon.test_stat(con,ind,lev);out.thrAUC.rich_club_neg_nodiscon.crit_val(con,ind,lev);out.thrAUC.rich_club_neg_nodiscon.p_2tail(con,ind,lev);out.thrAUC.rich_club_neg_nodiscon.hadNaN(con,ind);out.thrAUC.rich_club_neg_nodiscon.hadimag(con,ind);out.thrAUC.rich_club_neg_nodiscon.hadInf(con,ind)]')];
                        
                        fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, negative weights, excluding disconnected matrices in AUC):\n');
                        fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                        f_ind = find(ind);
                        for ind_val = 1:length(f_ind)
                            fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_neg_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_neg_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_neg_nodiscon.hadInf(con,f_ind(ind_val)));
                        end
                        fprintf(sigeffects_fid,'\n');
                    end
                    
                    if out.calcbinthresh==1
                        ind = find(out.thrAUC.rich_club_neg_bin_nodiscon.p_2tail(con,:,lev)<=out.alpha & ~isnan(out.thrAUC.rich_club_neg_bin_nodiscon.beta(con,:,lev)) & out.thrAUC.rich_club_neg_bin_nodiscon.beta(con,:,lev) ~= 0);
                        if any(ind)
                            out.sig_find(con).rich_club_neg_bin_nodiscon_thr{1,lev} = ['Contrast/F-test #',num2str(con)];
                            out.sig_find(con).rich_club_neg_bin_nodiscon_thr{2,lev} = {'node size','beta','stat','crit_val','p','had NaN','had imag','had Inf'};
                            out.sig_find(con).rich_club_neg_bin_nodiscon_thr{2,lev} = [out.sig_find(con).rich_club_neg_bin_nodiscon_thr{2,lev};num2cell([ind;out.thrAUC.rich_club_neg_bin_nodiscon.beta(con,ind,lev);out.thrAUC.rich_club_neg_bin_nodiscon.test_stat(con,ind,lev);out.thrAUC.rich_club_neg_bin_nodiscon.crit_val(con,ind,lev);out.thrAUC.rich_club_neg_bin_nodiscon.p_2tail(con,ind,lev);out.thrAUC.rich_club_neg_bin_nodiscon.hadNaN(con,ind);out.thrAUC.rich_club_neg_bin_nodiscon.hadimag(con,ind);out.thrAUC.rich_club_neg_bin_nodiscon.hadInf(con,ind)]')];
                            
                            fprintf(sigeffects_fid,'Rich Club Networks (thresholded matrices, negative weights, binarized matrices, excluding disconnected matrices in AUC):\n');
                            fprintf(sigeffects_fid,'\tnode size\tbeta\tstat\tcrit_val\tp\thad NaN\thad imag\thad Inf\n');
                            f_ind = find(ind);
                            for ind_val = 1:length(f_ind)
                                fprintf(sigeffects_fid,'\t%u\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%u\t%u\t%u\n',ind,out.thrAUC.rich_club_neg_bin_nodiscon.beta(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin_nodiscon.test_stat(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin_nodiscon.crit_val(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin_nodiscon.p_2tail(con,f_ind(ind_val),lev),out.thrAUC.rich_club_neg_bin_nodiscon.hadNaN(con,f_ind(ind_val)),out.thrAUC.rich_club_neg_bin_nodiscon.hadimag(con,f_ind(ind_val)),out.thrAUC.rich_club_neg_bin_nodiscon.hadInf(con,f_ind(ind_val)));
                            end
                            fprintf(sigeffects_fid,'\n');
                        end
                    end
                end
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

fclose(sigeffects_fid);
save(out.outname,'out')
fprintf('Done running permutation analyses!!!\n\n')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Calculate GLM & significance based on permutations
function [beta,test_stat,crit_val,p_1tail,p_2tail,hadNaN,hadimag,hadInf,p_permdis,nonperm_p] = GTG_GLM(DV,full_desmat,covars_desmat,Contrast_or_F,perms,reg_type,use_parfor,within_perms,alpha,MC_corr)

if size(DV,2)>size(DV,1)
    DV = DV';
end
if any(isnan(DV(:)))
    hadNaN = 1;
else
    hadNaN = 0;
end
if any(imag(DV(:))~=0)
    hadimag = 1;
else
    hadimag = 0;
end
if any(isinf(DV(:)))
    hadInf                         = 1;
    DV(isinf(DV) & (sign(DV)==1))  = max(DV(~isinf(DV)));
    DV(isinf(DV) & (sign(DV)==-1)) = min(DV(~isinf(DV)));
else
    hadInf = 0;
end

p_permdis = nan(size(perms,2),1);

if strcmp(Contrast_or_F,'Contrasts')
    if size(DV,2)==1
        if strcmp(reg_type,'OLS')
            temp                = regstats(DV,full_desmat,eye(size(full_desmat,2)),{'tstat'});
            beta                = temp.tstat.beta(1);
            test_stat           = temp.tstat.t(1);
            [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'OLS',use_parfor,MC_corr);
            nonperm_p           = temp.tstat.pval(1);
        elseif strcmp(reg_type,'Robust')
            [rob_b,rob_stats]   = robustfit_iterincrease(full_desmat,DV,'bisquare',4.685,'off');
            beta                = rob_b(1);
            test_stat           = rob_stats.t(1);
            [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'Robust',use_parfor,MC_corr);
            nonperm_p           = rob_stats.p(1);
        elseif strcmp(reg_type,'LTS')
            [rew_full,raw]      = ltsregres(full_desmat,DV,'plots',0,'intercept',0);
            beta                = raw.coefficients(1);
            rew_redu            = ltsregres(covars_desmat,DV,'plots',0,'intercept',0);
            test_stat           = ((rew_full.rsquared-rew_redu.rsquared)*(size(full_desmat,1)-size(full_desmat,2)))/(1-rew_full.rsquared);
            [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'LTS',use_parfor,MC_corr);
            nonperm_p           = fcdf(test_stat,1,(size(full_desmat,1)-size(full_desmat,2)),'upper');
        end
        if sign(test_stat)==1
            crit_val = permdis(round((1-alpha)*length(permdis)));
            p_1tail  = sum(permdis>=test_stat)/length(permdis);
            p_2tail  = sum(permdis<=-test_stat)/length(permdis)+p_1tail;
        else
            crit_val = permdis(round(alpha*length(permdis)));
            p_1tail  = sum(permdis<=test_stat)/length(permdis);
            p_2tail  = sum(permdis>=-test_stat)/length(permdis)+p_1tail;
        end
    elseif size(DV,3)==1
        K         = eye(size(full_desmat,2));
        M         = get_orthnorm_contmat(size(DV)); % Set M matrices, which specify what to do with the multiple DVs
        dof_error = size(full_desmat,1)-size(full_desmat,2)-1;                         % Degrees of freedom for the error
        
        % Estimate beta matrix
        B         = (full_desmat'*full_desmat)\full_desmat'*DV;
        beta      = M*B(1,:)';
        test_stat = zeros(size(M,2),1);
        crit_val  = zeros(size(M,2),1);
        p_1tail   = zeros(size(M,2),1);
        p_2tail   = zeros(size(M,2),1);
        
        for con = 1:size(M,2)
            H                        = M(:,con)'*(K(1,:)*B)'/(K(1,:)/(full_desmat'*full_desmat)*K(1,:)')*(K(1,:)*B)*M(:,con); % Estimate hypothesis matrix
            E                        = M(:,con)'*(DV'*DV-B'*(full_desmat'*full_desmat)*B)*M(:,con); % Estimate error matrix
            test_stat(con)           = H/(E/dof_error);                                         % Calculate F
            [permdis,temp_p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'~',use_parfor,MC_corr,within_perms,M(:,con));
            nonperm_p(con)           = fcdf(test_stat(con),1,dof_error,'upper');
            if MC_corr==1
                p_permdis = min(p_permdis,temp_p_permdis);
            end
            if sign(test_stat(con))==1
                crit_val(con) = permdis(round((1-alpha)*length(permdis)));
                p_1tail(con)  = sum(permdis>=test_stat(con))/length(permdis);
                p_2tail(con)  = sum(permdis<=-test_stat(con))/length(permdis)+p_1tail(con);
            else
                crit_val(con) = permdis(round(alpha*length(permdis)));
                p_1tail(con)  = sum(permdis<=test_stat(con))/length(permdis);
                p_2tail(con)  = sum(permdis>=-test_stat(con))/length(permdis)+p_1tail(con);
            end
        end
        if size(M,2)>2
            beta                     = [beta;sum(abs(beta(2:end)))];
            H                        = M(:,2:end)'*(K(1,:)*B)'/(K(1,:)/(full_desmat'*full_desmat)*K(1,:)')*(K(1,:)*B)*M(:,2:end); % Estimate hypothesis matrix
            E                        = M(:,2:end)'*(DV'*DV-B'*(full_desmat'*full_desmat)*B)*M(:,2:end); % Estimate error matrix
            test_stat(size(M,2)+1)   = det(E)/det(H+E);                                         % Calculate Wilks' Lambda
            [permdis,temp_p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'~',use_parfor,MC_corr,within_perms,M(:,2:end));
            
            if size(M(:,2:end),2)==2
                temp_F     = ((1-test_stat(size(M,2)+1))/test_stat(size(M,2)+1))*((numel(DV)-size(full_desmat,2)-1)/size(full_desmat,2));
                temp_dof_1 = size(full_desmat,2);
                temp_dof_2 = numel(DV)-size(full_desmat,2)-1;
            elseif size(M(:,2:end),2)==3
                temp_F     = ((1-sqrt(test_stat(size(M,2)+1)))/sqrt(test_stat(size(M,2)+1)))*((numel(DV)-size(full_desmat,2)-2)/size(full_desmat,2));
                temp_dof_1 = 2*size(full_desmat,2);
                temp_dof_2 = 2*(numel(DV)-size(full_desmat,2)-2);
            elseif size(full_desmat,2)==1
                temp_F     = ((1-test_stat(size(M,2)+1))/test_stat(size(M,2)+1))*((numel(DV)-size(M(:,2:end),2))/(size(M(:,2:end),2)-1));
                temp_dof_1 = size(M(:,2:end),2)-1;
                temp_dof_2 = numel(DV)-size(M(:,2:end),2);
            elseif size(full_desmat,2)==2
                temp_F     = ((1-sqrt(test_stat(size(M,2)+1)))/sqrt(test_stat(size(M,2)+1)))*((numel(DV)-size(M(:,2:end),2)-1)/(size(M(:,2:end),2)-1));
                temp_dof_1 = 2*(size(M(:,2:end),2)-1);
                temp_dof_2 = 2*(numel(DV)-size(M(:,2:end),2)-1);
            else
                pw = rank(H+E);
                qw = rank(K(1,:)/(full_desmat'*full_desmat)*K(1,:)');
                rw = dof_error-((pw-qw+1)/2);
                uw = ((pw*qw)-2)/4;
                if ((pw^2)+(qw^2)-5)>0
                    tw = sqrt((((pw^2)*(qw^2))-4)/((pw^2)+(qw^2)-5));
                else
                    tw = 1;
                end
                temp_F     = ((1-(test_stat(size(M,2)+1)^(1/tw)))/(test_stat(size(M,2)+1)^(1/tw)))*(((rw*tw)-(2*uw))/(pw*qw));
                temp_dof_1 = pw*qw;
                temp_dof_2 = (rw*tw)-(2*uw);
            end
            nonperm_p(size(M,2)+1) = fcdf(temp_F,temp_dof_1,temp_dof_2,'upper');
            if MC_corr==1
                p_permdis = min(p_permdis,temp_p_permdis);
            end
            if sign(test_stat(size(M,2)+1))==1
                crit_val(size(M,2)+1) = permdis(round((1-alpha)*length(permdis)));
                p_1tail(size(M,2)+1)  = sum(permdis>=test_stat(size(M,2)+1))/length(permdis);
                p_2tail(size(M,2)+1)  = sum(permdis<=-test_stat(size(M,2)+1))/length(permdis)+p_1tail(size(M,2)+1);
            else
                crit_val(size(M,2)+1) = permdis(round(alpha*length(permdis)));
                p_1tail(size(M,2)+1)  = sum(permdis<=test_stat(size(M,2)+1))/length(permdis);
                p_2tail(size(M,2)+1)  = sum(permdis>=-test_stat(size(M,2)+1))/length(permdis)+p_1tail(size(M,2)+1);
            end
        end
    end
else
    num_preds_in_F = size(full_desmat,2)-size(covars_desmat,2);
    if size(DV,2)==1
        if isempty(covars_desmat)
            if strcmp(reg_type,'OLS')
                temp                = regstats(DV,full_desmat,eye(size(full_desmat,2)),{'beta','fstat'});
                beta                = sum(abs(temp.beta(1:num_preds_in_F)));
                test_stat           = temp.fstat.f(1);
                [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'OLS',use_parfor,MC_corr);
                nonperm_p           = temp.fstat.pval(1);
            elseif strcmp(reg_type,'Robust')
                [rob_b,rob_stats]   = robustfit_iterincrease(full_desmat,DV,'bisquare',4.685,'off');
                beta                = sum(abs(rob_b(1:num_preds_in_F)));
                Rsq                 = 1-(corr(DV,rob_stats.resid))^2;
                test_stat           = ((Rsq)*(size(full_desmat,1)-num_preds_in_F-1))/((1-Rsq)*num_preds_in_F);
                [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'Robust',use_parfor,MC_corr);
                nonperm_p           = fcdf(test_stat,num_preds_in_F,(size(full_desmat,1)-num_preds_in_F-1),'upper');
            elseif strcmp(reg_type,'LTS')
                [rew,raw]           = ltsregres(full_desmat,DV,'plots',0,'intercept',0);
                beta                = sum(abs(raw.coefficients(1:num_preds_in_F)));
                test_stat           = (rew.rsquared*(size(full_desmat,1)-num_preds_in_F-1))/((1-rew.rsquared)*num_preds_in_F);
                [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'LTS',use_parfor,MC_corr);
                nonperm_p           = fcdf(test_stat,num_preds_in_F,(size(full_desmat,1)-num_preds_in_F-1),'upper');
            end
        else
            if strcmp(reg_type,'OLS')
                temp                = regstats(DV,full_desmat,eye(size(full_desmat,2)),{'beta','rsquare'});
                beta                = sum(abs(temp.beta(1:num_preds_in_F)));
                Rsq_full            = temp.rsquare(1);
                temp                = regstats(DV,covars_desmat,eye(size(covars_desmat,2)),{'rsquare'});
                Rsq_redu            = temp.rsquare(1);
                test_stat           = ((Rsq_full-Rsq_redu)*(size(full_desmat,1)-size(full_desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'OLS',use_parfor,MC_corr);
                nonperm_p           = fcdf(test_stat,num_preds_in_F,(size(full_desmat,1)-size(full_desmat,2)),'upper');
            elseif strcmp(reg_type,'Robust')
                [rob_b,rob_stats]   = robustfit_iterincrease(full_desmat,DV,'bisquare',4.685,'off');
                beta                = sum(abs(rob_b(1:num_preds_in_F)));
                Rsq_full            = 1-(corr(DV,rob_stats.resid))^2;
                [~,rob_stats]       = robustfit_iterincrease(covars_desmat,DV,'bisquare',4.685,'off');
                Rsq_redu            = 1-(corr(DV,rob_stats.resid))^2;
                test_stat           = ((Rsq_full-Rsq_redu)*(size(full_desmat,1)-size(full_desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'Robust',use_parfor,MC_corr);
                nonperm_p           = fcdf(test_stat,num_preds_in_F,(size(full_desmat,1)-size(full_desmat,2)),'upper');
            elseif strcmp(reg_type,'LTS')
                [rew_full,raw]      = ltsregres(full_desmat,DV,'plots',0,'intercept',0);
                beta                = sum(abs(raw.coefficients(1:num_preds_in_F)));
                rew_redu            = ltsregres(covars_desmat,DV,'plots',0,'intercept',0);
                test_stat           = ((rew_full.rsquared-rew_redu.rsquared)*(size(full_desmat,1)-size(full_desmat,2)))/((1-rew_full.rsquared)*num_preds_in_F);
                [permdis,p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'LTS',use_parfor,MC_corr);
                nonperm_p           = fcdf(test_stat,num_preds_in_F,(size(full_desmat,1)-size(full_desmat,2)),'upper');
            end
        end
        if sign(test_stat)==1
            crit_val = permdis(round((1-alpha)*length(permdis)));
            p_1tail  = sum(permdis>=test_stat)/length(permdis);
            p_2tail  = sum(permdis<=-test_stat)/length(permdis)+p_1tail;
        else
            crit_val = permdis(round(alpha*length(permdis)));
            p_1tail  = sum(permdis<=test_stat)/length(permdis);
            p_2tail  = sum(permdis>=-test_stat)/length(permdis)+p_1tail;
        end
    elseif size(DV,3)==1
        % Set M matrices, which specify what to do with the multiple DVs
        M = get_orthnorm_contmat(size(DV));
        if isempty(covars_desmat)
            % Estimate beta matrix
            B    = (full_desmat'*full_desmat)\full_desmat'*DV;
            beta = M*B';
            beta = sum(abs(beta(:,1:num_preds_in_F)));
            
            test_stat = zeros(size(M,2),1);
            crit_val  = zeros(size(M,2),1);
            p_1tail   = zeros(size(M,2),1);
            p_2tail   = zeros(size(M,2),1);
            
            for con = 1:size(M,2)
                Rsq                      = 1-(corr(DV*M(:,con),(DV*M(:,con)-full_desmat*B*M(:,con))))^2;
                test_stat(con)           = ((Rsq)*(size(full_desmat,1)-num_preds_in_F-1))/((1-Rsq)*num_preds_in_F);
                [permdis,temp_p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'~',use_parfor,MC_corr,within_perms,M(:,con));
                nonperm_p(con)           = fcdf(test_stat(con),num_preds_in_F,(size(full_desmat,1)-num_preds_in_F-1),'upper');
                if MC_corr==1
                    p_permdis = min(p_permdis,temp_p_permdis);
                end
                if sign(test_stat(con))==1
                    crit_val(con) = permdis(round((1-alpha)*length(permdis)));
                    p_1tail(con)  = sum(permdis>=test_stat(con))/length(permdis);
                    p_2tail(con)  = sum(permdis<=-test_stat(con))/length(permdis)+p_1tail(con);
                else
                    crit_val(con) = permdis(round(alpha*length(permdis)));
                    p_1tail(con)  = sum(permdis<=test_stat(con))/length(permdis);
                    p_2tail(con)  = sum(permdis>=-test_stat(con))/length(permdis)+p_1tail(con);
                end
            end
            if size(M,2)>2
                beta                     = [beta;sum(abs(beta(2:end)))];
                Rsq                      = 1-(corr2(DV*M(:,2:end),(DV*M(:,2:end)-full_desmat*B*M(:,2:end))))^2;
                test_stat(size(M,2)+1)   = ((Rsq)*(size(full_desmat,1)-num_preds_in_F-1))/((1-Rsq)*num_preds_in_F);
                [permdis,temp_p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'~',use_parfor,MC_corr,within_perms,M(:,2:end));
                nonperm_p(size(M,2)+1)   = fcdf(test_stat(size(M,2)+1),num_preds_in_F,(size(full_desmat,1)-num_preds_in_F-1),'upper');
                if MC_corr==1
                    p_permdis = min(p_permdis,temp_p_permdis);
                end
                if sign(test_stat(size(M,2)+1))==1
                    crit_val(size(M,2)+1) = permdis(round((1-alpha)*length(permdis)));
                    p_1tail(size(M,2)+1)  = sum(permdis>=test_stat(size(M,2)+1))/length(permdis);
                    p_2tail(size(M,2)+1)  = sum(permdis<=-test_stat(size(M,2)+1))/length(permdis)+p_1tail(size(M,2)+1);
                else
                    crit_val(size(M,2)+1) = permdis(round(alpha*length(permdis)));
                    p_1tail(size(M,2)+1)  = sum(permdis<=test_stat(size(M,2)+1))/length(permdis);
                    p_2tail(size(M,2)+1)  = sum(permdis>=-test_stat(size(M,2)+1))/length(permdis)+p_1tail(size(M,2)+1);
                end
            end
        else
            % Estimate beta matrix
            B_full = (full_desmat'*full_desmat)\full_desmat'*DV;
            B_redu = (covars_desmat'*covars_desmat)\covars_desmat'*DV;
            beta   = M*B_full';
            beta   = sum(abs(beta(:,1:num_preds_in_F)),2);
            
            test_stat = zeros(size(M,2),1);
            crit_val  = zeros(size(M,2),1);
            p_1tail   = zeros(size(M,2),1);
            p_2tail   = zeros(size(M,2),1);
            
            for con = 1:size(M,2)
                Rsq_full                 = 1-(corr(DV*M(:,con),(DV*M(:,con)-full_desmat*B_full*M(:,con))))^2;
                Rsq_redu                 = 1-(corr(DV*M(:,con),(DV*M(:,con)-covars_desmat*B_redu*M(:,con))))^2;
                test_stat(con)           = ((Rsq_full-Rsq_redu)*(size(full_desmat,1)-size(full_desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                [permdis,temp_p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'~',use_parfor,MC_corr,within_perms,M(:,con));
                nonperm_p(con)           = fcdf(test_stat(con),num_preds_in_F,(size(full_desmat,1)-size(fulldesmat,2)),'upper');
                if MC_corr==1
                    p_permdis = min(p_permdis,temp_p_permdis);
                end
                if sign(test_stat(con))==1
                    crit_val(con) = permdis(round((1-alpha)*length(permdis)));
                    p_1tail(con)  = sum(permdis>=test_stat(con))/length(permdis);
                    p_2tail(con)  = sum(permdis<=-test_stat(con))/length(permdis)+p_1tail(con);
                else
                    crit_val(con) = permdis(round(alpha*length(permdis)));
                    p_1tail(con)  = sum(permdis<=test_stat(con))/length(permdis);
                    p_2tail(con)  = sum(permdis>=-test_stat(con))/length(permdis)+p_1tail(con);
                end
            end
            if size(M,2)>2
                beta                     = [beta;sum(abs(beta(2:end)))];
                Rsq_full                 = 1-(corr2(DV*M(:,2:end),(DV*M(:,2:end)-full_desmat*B_full*M(:,2:end))))^2;
                Rsq_redu                 = 1-(corr2(DV*M(:,2:end),(DV*M(:,2:end)-covars_desmat*B_redu*M(:,2:end))))^2;
                test_stat(size(M,2)+1)   = ((Rsq_full-Rsq_redu)*(size(full_desmat,1)-size(full_desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                [permdis,temp_p_permdis] = create_perm_dist(DV,full_desmat,covars_desmat,Contrast_or_F,perms,'~',use_parfor,MC_corr,within_perms,M(:,2:end));
                nonperm_p(size(M,2)+1)   = fcdf(test_stat(size(M,2)+1),num_preds_in_F,(size(full_desmat,1)-size(full_desmat,2)),'upper');
                if MC_corr==1
                    p_permdis = min(p_permdis,temp_p_permdis);
                end
                if sign(test_stat(size(M,2)+1))==1
                    crit_val(size(M,2)+1) = permdis(round((1-alpha)*length(permdis)));
                    p_1tail(size(M,2)+1)  = sum(permdis>=test_stat(size(M,2)+1))/length(permdis);
                    p_2tail(size(M,2)+1)  = sum(permdis<=-test_stat(size(M,2)+1))/length(permdis)+p_1tail(size(M,2)+1);
                else
                    crit_val(size(M,2)+1) = permdis(round(alpha*length(permdis)));
                    p_1tail(size(M,2)+1)  = sum(permdis<=test_stat(size(M,2)+1))/length(permdis);
                    p_2tail(size(M,2)+1)  = sum(permdis>=-test_stat(size(M,2)+1))/length(permdis)+p_1tail(size(M,2)+1);
                end
            end
        end
    end
end

% if size(beta,1)>1
%     hadNaN(2:size(beta,1),1)  = hadNaN(1);
%     hadimag(2:size(beta,1),1) = hadimag(1);
%     hadInf(2:size(beta,1),1)  = hadInf(1);
% end



%%%% Obtain the orthonormal repeated contrast matrix
function M = get_orthnorm_contmat(size_DV)

if length(size_DV)==2
    switch size_DV(2)
        case 2                             % 2 levels
            M(:,1) = [1;1]./norm([1;1]);   % Orthonormal M Matrix for b/t effects
            M(:,2) = [-1;1]./norm([-1;1]); % Orthonormal M Matrix for w/in linear effects
        case 3                                 % 3 levels
            M(:,1) = [1;1;1]./norm([1;1;1]);   % Orthonormal M Matrix for b/t effects
            M(:,2) = [-1;0;1]./norm([-1;0;1]); % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [1;-2;1]./norm([1;-2;1]); % Orthonormal M Matrix for w/in quadratic effects
        case 4                                       % 4 levels
            M(:,1) = [1;1;1;1]./norm([1;1;1;1]);     % Orthonormal M Matrix for b/t effects
            M(:,2) = [-3;-1;1;3]./norm([-3;-1;1;3]); % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [1;-1;-1;1]./norm([1;-1;-1;1]); % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-1;3;-3;1]./norm([-1;3;-3;1]); % Orthonormal M Matrix for w/in cubic effects
        case 5                                             % 5 levels
            M(:,1) = [1;1;1;1;1]./norm([1;1;1;1;1]);       % Orthonormal M Matrix for b/t effects
            M(:,2) = [-2;-1;0;1;2]./norm([-2;-1;0;1;2]);   % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [2;-1;-2;-1;2]./norm([2;-1;-2;-1;2]); % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-1;2;0;-2;1]./norm([-1;2;0;-2;1]);   % Orthonormal M Matrix for w/in cubic effects
            M(:,5) = [1;-4;6;-4;1]./norm([1;-4;6;-4;1]);   % Orthonormal M Matrix for w/in quartic effects
        case 6                                                     % 6 levels
            M(:,1) = [1;1;1;1;1;1]./norm([1;1;1;1;1;1]);           % Orthonormal M Matrix for b/t effects
            M(:,2) = [-5;-3;-1;1;3;5]./norm([-5;-3;-1;1;3;5]);     % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [5;-1;-4;-4;-1;5]./norm([5;-1;-4;-4;-1;5]);   % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-5;7;4;-4;-7;5]./norm([-5;7;4;-4;-7;5]);     % Orthonormal M Matrix for w/in cubic effects
            M(:,5) = [1;-3;2;2;-3;1]./norm([1;-3;2;2;-3;1]);       % Orthonormal M Matrix for w/in quartic effects
            M(:,6) = [-1;5;-10;10;-5;1]./norm([-1;5;-10;10;-5;1]); % Orthonormal M Matrix for w/in quintic effects
        case 7                                                           % 7 levels
            M(:,1) = [1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1]);             % Orthonormal M Matrix for b/t effects
            M(:,2) = [-3;-2;-1;0;1;2;3]./norm([-3;-2;-1;0;1;2;3]);       % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [5;0;-3;-4;-3;0;5]./norm([5;0;-3;-4;-3;0;5]);       % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-1;1;1;0;-1;-1;1]./norm([-1;1;1;0;-1;-1;1]);       % Orthonormal M Matrix for w/in cubic effects
            M(:,5) = [3;-7;1;6;1;-7;3]./norm([3;-7;1;6;1;-7;3]);         % Orthonormal M Matrix for w/in quartic effects
            M(:,6) = [-1;4;-5;0;5;-4;1]./norm([-1;4;-5;0;5;-4;1]);       % Orthonormal M Matrix for w/in quintic effects
            M(:,7) = [1;-6;15;-20;15;-6;1]./norm([1;-6;15;-20;15;-6;1]); % Orthonormal M Matrix for w/in sextic effects
        case 8                                                                                           % 8 levels
            M(:,1) = [1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1]);                                           % Orthonormal M Matrix for b/t effects
            M(:,2) = [-7;-5;-3;-1;1;3;5;7]./norm([-7;-5;-3;-1;1;3;5;7]);                                 % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [7;1;-3;-5;-5;-3;1;7]./norm([7;1;-3;-5;-5;-3;1;7]);                                 % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-7;5;7;3;-3;-7;-5;7]./norm([-7;5;7;3;-3;-7;-5;7]);                                 % Orthonormal M Matrix for w/in cubic effects
            M(:,5) = [7;-13;-3;9;9;-3;-13;7]./norm([7;-13;-3;9;9;-3;-13;7]);                             % Orthonormal M Matrix for w/in quartic effects
            M(:,6) = [-150;492;-364;-321;321;364;-492;150]./norm([-150;492;-364;-321;321;364;-492;150]); % Orthonormal M Matrix for w/in quintic effects
            M(:,7) = [31;-154;277;-154;-154;277;-154;31]./norm([31;-154;277;-154;-154;277;-154;31]);     % Orthonormal M Matrix for w/in sextic effects
            M(:,8) = [-17;119;-358;597;-597;358;-119;17]./norm([-17;119;-358;597;-597;358;-119;17]);     % Orthonormal M Matrix for w/in septic effects
    end
elseif length(size_DV)==3
    switch size_DV(2)
        case 2                                       % 2 levels
            M(:,1) = [1;1;1;1]./norm([1;1;1;1]);     % Orthonormal M Matrix for b/t effects
            M(:,2) = [-1;1;-1;1]./norm([-1;1;-1;1]); % Orthonormal M Matrix for factor 1 w/in linear effects
            M(:,3) = [-1;-1;1;1]./norm([-1;-1;1;1]); % Orthonormal M Matrix for factor 2 w/in linear effects
            M(:,4) = [1;-1;-1;1]./norm([1;-1;-1;1]); % Orthonormal M Matrix for w/in linear x linear effects
        case 3                                                         % 3 levels
            switch size_DV(3)
                case 2
                    M(:,1) = [1;1;1;1;1;1]./norm([1;1;1;1;1;1]);       % Orthonormal M Matrix for b/t effects
                    M(:,2) = [-1;0;1;-1;0;1]./norm([-1;0;1;-1;0;1]);   % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3) = [1;-2;1;1;-2;1]./norm([1;-2;1;1;-2;1]);   % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4) = [-1;-1;-1;1;1;1]./norm([-1;-1;-1;1;1;1]); % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,5) = [1;0;-1;-1;0;1]./norm([1;0;-1;-1;0;1]);   % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,6) = [-1;2;-1;1;-2;1]./norm([-1;2;-1;1;-2;1]); % Orthonormal M Matrix for w/in quadratic x linear effects
                case 3
                    M(:,1) = [1;1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1;1;1]);         % Orthonormal M Matrix for b/t effects
                    M(:,2) = [-1;0;1;-1;0;1;-1;0;1]./norm([-1;0;1;-1;0;1;-1;0;1]);   % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3) = [1;-2;1;1;-2;1;1;-2;1]./norm([1;-2;1;1;-2;1;1;-2;1]);   % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4) = [-1;-1;-1;0;0;0;1;1;1]./norm([-1;-1;-1;0;0;0;1;1;1]);   % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,5) = [1;1;1;-2;-2;-2;1;1;1]./norm([1;1;1;-2;-2;-2;1;1;1]);   % Orthonormal M Matrix for factor 2 w/in quadratic effects
                    M(:,6) = [1;0;-1;0;0;0;-1;0;1]./norm([1;0;-1;0;0;0;-1;0;1]);     % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,7) = [-1;2;-1;0;0;0;1;-2;1]./norm([-1;2;-1;0;0;0;1;-2;1]);   % Orthonormal M Matrix for w/in quadratic x linear effects
                    M(:,8) = [-1;0;1;2;0;-2;-1;0;1]./norm([-1;0;1;2;0;-2;-1;0;1]);   % Orthonormal M Matrix for w/in linear x quadratic effects
                    M(:,9) = [1;-2;1;-2;4;-2;1;-2;1]./norm([1;-2;1;-2;4;-2;1;-2;1]); % Orthonormal M Matrix for w/in quadratic x quadratic effects
            end
        case 4                                                                   % 4 levels
            switch size_DV(3)
                case 2
                    M(:,1) = [1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1;1]);         % Orthonormal M Matrix for b/t effects
                    M(:,2) = [-3;-1;1;3;-3;-1;1;3]./norm([-3;-1;1;3;-3;-1;1;3]); % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3) = [1;-1;-1;1;1;-1;-1;1]./norm([1;-1;-1;1;1;-1;-1;1]); % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4) = [-1;3;-3;1;-1;3;-3;1]./norm([-1;3;-3;1;-1;3;-3;1]); % Orthonormal M Matrix for factor 1 w/in cubic effects
                    M(:,5) = [-1;-1;-1;-1;1;1;1;1]./norm([-1;-1;-1;-1;1;1;1;1]); % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,6) = [3;1;-1;-3;-3;-1;1;3]./norm([3;1;-1;-3;-3;-1;1;3]); % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,7) = [-1;1;1;-1;1;-1;-1;1]./norm([-1;1;1;-1;1;-1;-1;1]); % Orthonormal M Matrix for w/in quadratic x linear effects
                    M(:,8) = [1;-3;3;-1;-1;3;-3;1]./norm([1;-3;3;-1;-1;3;-3;1]); % Orthonormal M Matrix for w/in cubic x linear effects
                case 3
                    M(:,1)  = [1;1;1;1;1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1;1;1;1;1;1]);             % Orthonormal M Matrix for b/t effects
                    M(:,2)  = [-3;-1;1;3;-3;-1;1;3;-3;-1;1;3]./norm([-3;-1;1;3;-3;-1;1;3;-3;-1;1;3]); % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3)  = [1;-1;-1;1;1;-1;-1;1;1;-1;-1;1]./norm([1;-1;-1;1;1;-1;-1;1;1;-1;-1;1]); % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4)  = [-1;3;-3;1;-1;3;-3;1;-1;3;-3;1]./norm([-1;3;-3;1;-1;3;-3;1;-1;3;-3;1]); % Orthonormal M Matrix for factor 1 w/in cubic effects
                    M(:,5)  = [-1;-1;-1;-1;0;0;0;0;1;1;1;1]./norm([-1;-1;-1;-1;0;0;0;0;1;1;1;1]);     % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,6)  = [1;1;1;1;-2;-2;-2;-2;1;1;1;1]./norm([1;1;1;1;-2;-2;-2;-2;1;1;1;1]);     % Orthonormal M Matrix for factor 2 w/in quadratic effects
                    M(:,7)  = [3;1;-1;-3;0;0;0;0;-3;-1;1;3]./norm([3;1;-1;-3;0;0;0;0;-3;-1;1;3]);     % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,8)  = [-1;1;1;-1;0;0;0;0;1;-1;-1;1]./norm([-1;1;1;-1;0;0;0;0;1;-1;-1;1]);     % Orthonormal M Matrix for w/in quadratic x linear effects
                    M(:,9)  = [1;-3;3;-1;0;0;0;0;-1;3;-3;1]./norm([1;-3;3;-1;0;0;0;0;-1;3;-3;1]);     % Orthonormal M Matrix for w/in cubic x linear effects
                    M(:,10) = [-3;-1;1;3;6;2;-2;-6;-3;-1;1;3]./norm([-3;-1;1;3;6;2;-2;-6;-3;-1;1;3]); % Orthonormal M Matrix for w/in linear x quadratic effects
                    M(:,11) = [1;-1;-1;1;-2;2;2;-2;1;-1;-1;1]./norm([1;-1;-1;1;-2;2;2;-2;1;-1;-1;1]); % Orthonormal M Matrix for w/in quadratic x quadratic effects
                    M(:,12) = [-1;3;-3;1;2;-6;6;-2;-1;3;-3;1]./norm([-1;3;-3;1;2;-6;6;-2;-1;3;-3;1]); % Orthonormal M Matrix for w/in cubic x quadratic effects
                case 4
                    M(:,1)  = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]);                 % Orthonormal M Matrix for b/t effects
                    M(:,2)  = [-3;-1;1;3;-3;-1;1;3;-3;-1;1;3;-3;-1;1;3]./norm([-3;-1;1;3;-3;-1;1;3;-3;-1;1;3;-3;-1;1;3]); % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3)  = [1;-1;-1;1;1;-1;-1;1;1;-1;-1;1;1;-1;-1;1]./norm([1;-1;-1;1;1;-1;-1;1;1;-1;-1;1;1;-1;-1;1]); % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4)  = [-1;3;-3;1;-1;3;-3;1;-1;3;-3;1;-1;3;-3;1]./norm([-1;3;-3;1;-1;3;-3;1;-1;3;-3;1;-1;3;-3;1]); % Orthonormal M Matrix for factor 1 w/in cubic effects
                    M(:,5)  = [-3;-3;-3;-3;-1;-1;-1;-1;1;1;1;1;3;3;3;3]./norm([-3;-3;-3;-3;-1;-1;-1;-1;1;1;1;1;3;3;3;3]); % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,6)  = [1;1;1;1;-1;-1;-1;-1;-1;-1;-1;-1;1;1;1;1]./norm([1;1;1;1;-1;-1;-1;-1;-1;-1;-1;-1;1;1;1;1]); % Orthonormal M Matrix for factor 2 w/in quadratic effects
                    M(:,7)  = [-1;-1;-1;-1;3;3;3;3;-3;-3;-3;-3;1;1;1;1]./norm([-1;-1;-1;-1;3;3;3;3;-3;-3;-3;-3;1;1;1;1]); % Orthonormal M Matrix for factor 2 w/in cubic effects
                    M(:,8)  = [9;3;-3;-9;3;1;-1;-3;-3;-1;1;3;-9;-3;3;9]./norm([9;3;-3;-9;3;1;-1;-3;-3;-1;1;3;-9;-3;3;9]); % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,9)  = [-3;3;3;-3;-1;1;1;-1;1;-1;-1;1;3;-3;-3;3]./norm([-3;3;3;-3;-1;1;1;-1;1;-1;-1;1;3;-3;-3;3]); % Orthonormal M Matrix for w/in quadratic x linear effects
                    M(:,10) = [3;-9;9;-3;1;-3;3;-1;-1;3;-3;1;-3;9;-9;3]./norm([3;-9;9;-3;1;-3;3;-1;-1;3;-3;1;-3;9;-9;3]); % Orthonormal M Matrix for w/in cubic x linear effects
                    M(:,11) = [-3;-1;1;3;3;1;-1;-3;3;1;-1;-3;-3;-1;1;3]./norm([-3;-1;1;3;3;1;-1;-3;3;1;-1;-3;-3;-1;1;3]); % Orthonormal M Matrix for w/in linear x quadratic effects
                    M(:,12) = [1;-1;-1;1;-1;1;1;-1;-1;1;1;-1;1;-1;-1;1]./norm([1;-1;-1;1;-1;1;1;-1;-1;1;1;-1;1;-1;-1;1]); % Orthonormal M Matrix for w/in quadratic x quadratic effects
                    M(:,13) = [-1;3;-3;1;1;-3;3;-1;1;-3;3;-1;-1;3;-3;1]./norm([-1;3;-3;1;1;-3;3;-1;1;-3;3;-1;-1;3;-3;1]); % Orthonormal M Matrix for w/in cubic x quadratic effects
                    M(:,14) = [3;1;-1;-3;-9;-3;3;9;9;3;-3;-9;-3;-1;1;3]./norm([3;1;-1;-3;-9;-3;3;9;9;3;-3;-9;-3;-1;1;3]); % Orthonormal M Matrix for w/in linear x cubic effects
                    M(:,15) = [-1;1;1;-1;3;-3;-3;3;-3;3;3;-3;1;-1;-1;1]./norm([-1;1;1;-1;3;-3;-3;3;-3;3;3;-3;1;-1;-1;1]); % Orthonormal M Matrix for w/in quadratic x cubic effects
                    M(:,16) = [1;-3;3;-1;-3;9;-9;3;3;-9;9;-3;-1;3;-3;1]./norm([1;-3;3;-1;-3;9;-9;3;3;-9;9;-3;-1;3;-3;1]); % Orthonormal M Matrix for w/in cubic x cubic effects
            end
        case 5                                                                                % 5 levels
            switch size_DV(3)
                case 2
                    M(:,1)  = [1,1,1,1,1,1,1,1,1,1]./norm([1,1,1,1,1,1,1,1,1,1]);             % Orthonormal M Matrix for b/t effects
                    M(:,2)  = [-2;-1;0;1;2;-2;-1;0;1;2]./norm([-2;-1;0;1;2;-2;-1;0;1;2]);     % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3)  = [2;-1;-2;-1;2;2;-1;-2;-1;2]./norm([2;-1;-2;-1;2;2;-1;-2;-1;2]); % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4)  = [-1;2;0;-2;1;-1;2;0;-2;1]./norm([-1;2;0;-2;1;-1;2;0;-2;1]);     % Orthonormal M Matrix for factor 1 w/in cubic effects
                    M(:,5)  = [1;-4;6;-4;1;1;-4;6;-4;1]./norm([1;-4;6;-4;1;1;-4;6;-4;1]);     % Orthonormal M Matrix for factor 1 w/in quartic effects
                    M(:,6)  = [-1;-1;-1;-1-1;1;1;1;1;1]./norm([-1;-1;-1;-1-1;1;1;1;1;1]);     % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,7)  = [2;1;0;-1;-2;-2;-1;0;1;2]./norm([2;1;0;-1;-2;-2;-1;0;1;2]);     % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,8)  = [-2;1;2;1;-2;2;-1;-2;-1;2]./norm([-2;1;2;1;-2;2;-1;-2;-1;2]);   % Orthonormal M Matrix for w/in quadratic x linear effects
                    M(:,9)  = [1;-2;0;2;-1;-1;2;0;-2;1]./norm([1;-2;0;2;-1;-1;2;0;-2;1]);     % Orthonormal M Matrix for w/in cubic x linear effects
                    M(:,10) = [-1;4;-6;4;-1;1;-4;6;-4;1]./norm([-1;4;-6;4;-1;1;-4;6;-4;1]);   % Orthonormal M Matrix for w/in quartic x linear effects
                case 3
                    M(:,1)  = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]);                   % Orthonormal M Matrix for b/t effects
                    M(:,2)  = [-2;-1;0;1;2;-2;-1;0;1;2;-2;-1;0;1;2]./norm([-2;-1;0;1;2;-2;-1;0;1;2;-2;-1;0;1;2]);       % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3)  = [2;-1;-2;-1;2;2;-1;-2;-1;2;2;-1;-2;-1;2]./norm([2;-1;-2;-1;2;2;-1;-2;-1;2;2;-1;-2;-1;2]); % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4)  = [-1;2;0;-2;1;-1;2;0;-2;1;-1;2;0;-2;1]./norm([-1;2;0;-2;1;-1;2;0;-2;1;-1;2;0;-2;1]);       % Orthonormal M Matrix for factor 1 w/in cubic effects
                    M(:,5)  = [1;-4;6;-4;1;1;-4;6;-4;1;1;-4;6;-4;1]./norm([1;-4;6;-4;1;1;-4;6;-4;1;1;-4;6;-4;1]);       % Orthonormal M Matrix for factor 1 w/in quartic effects
                    M(:,6)  = [-1;-1;-1;-1;-1;0;0;0;0;0;1;1;1;1;1]./norm([-1;-1;-1;-1;-1;0;0;0;0;0;1;1;1;1;1]);         % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,7)  = [1;1;1;1;1;-2;-2;-2;-2;-2;1;1;1;1;1]./norm([1;1;1;1;1;-2;-2;-2;-2;-2;1;1;1;1;1]);         % Orthonormal M Matrix for factor 2 w/in quadratic effects
                    M(:,8)  = [2;1;0;-1;-2;0;0;0;0;0;-2;-1;0;1;2]./norm([2;1;0;-1;-2;0;0;0;0;0;-2;-1;0;1;2]);           % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,9)  = [-2;1;2;1;-2;0;0;0;0;0;2;-1;-2;-1;2]./norm([-2;1;2;1;-2;0;0;0;0;0;2;-1;-2;-1;2]);         % Orthonormal M Matrix for w/in quadratic x linear effects
                    M(:,10) = [1;-2;0;2;-1;0;0;0;0;0;-1;2;0;-2;1]./norm([1;-2;0;2;-1;0;0;0;0;0;-1;2;0;-2;1]);           % Orthonormal M Matrix for w/in cubic x linear effects
                    M(:,11) = [-1;4;-6;4;-1;0;0;0;0;0;1;-4;6;-4;1]./norm([-1;4;-6;4;-1;0;0;0;0;0;1;-4;6;-4;1]);         % Orthonormal M Matrix for w/in quartic x linear effects
                    M(:,12) = [-2;-1;0;1;2;4;2;0;-2;-4;-2;-1;0;1;2]./norm([-2;-1;0;1;2;4;2;0;-2;-4;-2;-1;0;1;2]);       % Orthonormal M Matrix for w/in linear x quadratic effects
                    M(:,13) = [2;-1;-2;-1;2;-4;2;4;2;-4;2;-1;-2;-1;2]./norm([2;-1;-2;-1;2;-4;2;4;2;-4;2;-1;-2;-1;2]);   % Orthonormal M Matrix for w/in quadratic x quadratic effects
                    M(:,14) = [-1;2;0;-2;1;2;-4;0;4;-2;-1;2;0;-2;1]./norm([-1;2;0;-2;1;2;-4;0;4;-2;-1;2;0;-2;1]);       % Orthonormal M Matrix for w/in cubic x quadratic effects
                    M(:,15) = [1;-4;6;-4;1;-2;8;-12;8;-2;1;-4;6;-4;1]./norm([1;-4;6;-4;1;-2;8;-12;8;-2;1;-4;6;-4;1]);   % Orthonormal M Matrix for w/in quartic x quadratic effects
                case 4
                    M(:,1)  = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]);                                 % Orthonormal M Matrix for b/t effects
                    M(:,2)  = [-2;-1;0;1;2;-2;-1;0;1;2;-2;-1;0;1;2;-2;-1;0;1;2]./norm([-2;-1;0;1;2;-2;-1;0;1;2;-2;-1;0;1;2;-2;-1;0;1;2]);                 % Orthonormal M Matrix for factor 1 w/in linear effects
                    M(:,3)  = [2;-1;-2;-1;2;2;-1;-2;-1;2;2;-1;-2;-1;2;2;-1;-2;-1;2]./norm([2;-1;-2;-1;2;2;-1;-2;-1;2;2;-1;-2;-1;2;2;-1;-2;-1;2]);         % Orthonormal M Matrix for factor 1 w/in quadratic effects
                    M(:,4)  = [-1;2;0;-2;1;-1;2;0;-2;1;-1;2;0;-2;1;-1;2;0;-2;1]./norm([-1;2;0;-2;1;-1;2;0;-2;1;-1;2;0;-2;1;-1;2;0;-2;1]);                 % Orthonormal M Matrix for factor 1 w/in cubic effects
                    M(:,5)  = [1;-4;6;-4;1;1;-4;6;-4;1;1;-4;6;-4;1;1;-4;6;-4;1]./norm([1;-4;6;-4;1;1;-4;6;-4;1;1;-4;6;-4;1;1;-4;6;-4;1]);                 % Orthonormal M Matrix for factor 1 w/in quartic effects
                    M(:,6)  = [-3;-3;-3;-3;-3;-1;-1;-1;-1;-1;1;1;1;1;1;3;3;3;3;3]./norm([-3;-3;-3;-3;-3;-1;-1;-1;-1;-1;1;1;1;1;1;3;3;3;3;3]);             % Orthonormal M Matrix for factor 2 w/in linear effects
                    M(:,7)  = [1;1;1;1;1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;1;1;1;1;1]./norm([1;1;1;1;1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;1;1;1;1;1]);             % Orthonormal M Matrix for factor 2 w/in quadratic effects
                    M(:,8)  = [-1;-1;-1;-1;-1;3;3;3;3;3;-3;-3;-3;-3;-3;1;1;1;1;1]./norm([-1;-1;-1;-1;-1;3;3;3;3;3;-3;-3;-3;-3;-3;1;1;1;1;1]);             % Orthonormal M Matrix for factor 2 w/in cubic effects
                    M(:,9)  = [6;3;0;-3;-6;2;1;0;-1;-2;-2;-1;0;1;2;-6;-3;0;3;6]./norm([6;3;0;-3;-6;2;1;0;-1;-2;-2;-1;0;1;2;-6;-3;0;3;6]);                 % Orthonormal M Matrix for w/in linear x linear effects
                    M(:,10) = [-6;3;6;3;-6;-2;1;2;1;-2;2;-1;-2;-1;2;6;-3;-6;-3;6]./norm([-6;3;6;3;-6;-2;1;2;1;-2;2;-1;-2;-1;2;6;-3;-6;-3;6]);             % Orthonormal M Matrix for w/in quadratic x linear effects
                    M(:,11) = [3;-6;0;6;-3;1;-2;0;2;-1;-1;2;0;-2;1;-3;6;0;-6;3]./norm([3;-6;0;6;-3;1;-2;0;2;-1;-1;2;0;-2;1;-3;6;0;-6;3]);                 % Orthonormal M Matrix for w/in cubic x linear effects
                    M(:,12) = [-3;12;-18;12;-3;-1;4;-6;4;-1;1;-4;6;-4;1;3;-12;18;-12;3]./norm([-3;12;-18;12;-3;-1;4;-6;4;-1;1;-4;6;-4;1;3;-12;18;-12;3]); % Orthonormal M Matrix for w/in quartic x linear effects
                    M(:,13) = [-2;-1;0;1;2;2;1;0;-1;-2;2;1;0;-1;-2;-2;-1;0;1;2]./norm([-2;-1;0;1;2;2;1;0;-1;-2;2;1;0;-1;-2;-2;-1;0;1;2]);                 % Orthonormal M Matrix for w/in linear x quadratic effects
                    M(:,14) = [2;-1;-2;-1;2;-2;1;2;1;-2;-2;1;2;1;-2;2;-1;-2;-1;2]./norm([2;-1;-2;-1;2;-2;1;2;1;-2;-2;1;2;1;-2;2;-1;-2;-1;2]);             % Orthonormal M Matrix for w/in quadratic x quadratic effects
                    M(:,15) = [-1;2;0;-2;1;1;-2;0;2;-1;1;-2;0;2;-1;-1;2;0;-2;1]./norm([-1;2;0;-2;1;1;-2;0;2;-1;1;-2;0;2;-1;-1;2;0;-2;1]);                 % Orthonormal M Matrix for w/in cubic x quadratic effects
                    M(:,16) = [1;-4;6;-4;1;-1;4;-6;4;-1;-1;4;-6;4;-1;1;-4;6;-4;1]./norm([1;-4;6;-4;1;-1;4;-6;4;-1;-1;4;-6;4;-1;1;-4;6;-4;1]);             % Orthonormal M Matrix for w/in quartic x quadratic effects
                    M(:,17) = [2;1;0;-1;-2;-6;-3;0;3;6;6;3;0;-3;-6;-2;-1;0;1;2]./norm([2;1;0;-1;-2;-6;-3;0;3;6;6;3;0;-3;-6;-2;-1;0;1;2]);                 % Orthonormal M Matrix for w/in linear x cubic effects
                    M(:,18) = [-2;1;2;1;-2;6;-3;-6;-3;6;-6;3;6;3;-6;2;-1;-2;-1;2]./norm([-2;1;2;1;-2;6;-3;-6;-3;6;-6;3;6;3;-6;2;-1;-2;-1;2]);             % Orthonormal M Matrix for w/in quadratic x cubic effects
                    M(:,19) = [1;-2;0;2;-1;-3;6;0;-6;3;3;-6;0;6;-3;-1;2;0;-2;1]./norm([1;-2;0;2;-1;-3;6;0;-6;3;3;-6;0;6;-3;-1;2;0;-2;1]);                 % Orthonormal M Matrix for w/in cubic x cubic effects
                    M(:,20) = [-1;4;-6;4;-1;3;-12;18;-12;3;-3;12;-18;12;-3;1;-4;6;-4;1]./norm([-1;4;-6;4;-1;3;-12;18;-12;3;-3;12;-18;12;-3;1;-4;6;-4;1]); % Orthonormal M Matrix for w/in quartic x cubic effects
            end
        case 6                                                     % 6 levels
            M(:,1) = [1;1;1;1;1;1]./norm([1;1;1;1;1;1]);           % Orthonormal M Matrix for b/t effects
            M(:,2) = [-5;-3;-1;1;3;5]./norm([-5;-3;-1;1;3;5]);     % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [5;-1;-4;-4;-1;5]./norm([5;-1;-4;-4;-1;5]);   % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-5;7;4;-4;-7;5]./norm([-5;7;4;-4;-7;5]);     % Orthonormal M Matrix for w/in cubic effects
            M(:,5) = [1;-3;2;2;-3;1]./norm([1;-3;2;2;-3;1]);       % Orthonormal M Matrix for w/in quartic effects
            M(:,6) = [-1;5;-10;10;-5;1]./norm([-1;5;-10;10;-5;1]); % Orthonormal M Matrix for w/in quintic effects
        case 7                                                           % 7 levels
            M(:,1) = [1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1]);             % Orthonormal M Matrix for b/t effects
            M(:,2) = [-3;-2;-1;0;1;2;3]./norm([-3;-2;-1;0;1;2;3]);       % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [5;0;-3;-4;-3;0;5]./norm([5;0;-3;-4;-3;0;5]);       % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-1;1;1;0;-1;-1;1]./norm([-1;1;1;0;-1;-1;1]);       % Orthonormal M Matrix for w/in cubic effects
            M(:,5) = [3;-7;1;6;1;-7;3]./norm([3;-7;1;6;1;-7;3]);         % Orthonormal M Matrix for w/in quartic effects
            M(:,6) = [-1;4;-5;0;5;-4;1]./norm([-1;4;-5;0;5;-4;1]);       % Orthonormal M Matrix for w/in quintic effects
            M(:,7) = [1;-6;15;-20;15;-6;1]./norm([1;-6;15;-20;15;-6;1]); % Orthonormal M Matrix for w/in sextic effects
        case 8                                                                                           % 8 levels
            M(:,1) = [1;1;1;1;1;1;1;1]./norm([1;1;1;1;1;1;1]);                                           % Orthonormal M Matrix for b/t effects
            M(:,2) = [-7;-5;-3;-1;1;3;5;7]./norm([-7;-5;-3;-1;1;3;5;7]);                                 % Orthonormal M Matrix for w/in linear effects
            M(:,3) = [7;1;-3;-5;-5;-3;1;7]./norm([7;1;-3;-5;-5;-3;1;7]);                                 % Orthonormal M Matrix for w/in quadratic effects
            M(:,4) = [-7;5;7;3;-3;-7;-5;7]./norm([-7;5;7;3;-3;-7;-5;7]);                                 % Orthonormal M Matrix for w/in cubic effects
            M(:,5) = [7;-13;-3;9;9;-3;-13;7]./norm([7;-13;-3;9;9;-3;-13;7]);                             % Orthonormal M Matrix for w/in quartic effects
            M(:,6) = [-150;492;-364;-321;321;364;-492;150]./norm([-150;492;-364;-321;321;364;-492;150]); % Orthonormal M Matrix for w/in quintic effects
            M(:,7) = [31;-154;277;-154;-154;277;-154;31]./norm([31;-154;277;-154;-154;277;-154;31]);     % Orthonormal M Matrix for w/in sextic effects
            M(:,8) = [-17;119;-358;597;-597;358;-119;17]./norm([-17;119;-358;597;-597;358;-119;17]);     % Orthonormal M Matrix for w/in septic effects
    end
end



%%%% Create a permutation distribution for a partial regression
%%%% coefficient
function [test_stats,p_stats] = create_perm_dist(DV,desmat,covariates,Contrast_or_F,perms,reg_type,use_parfor,MC_corr,varargin)

if ~isempty(varargin)
    within_perms = varargin{1};
    M            = varargin{2};
end

num_perms  = size(perms,2);
test_stats = nan(num_perms,1);
p_stats    = nan(num_perms,1);
dof_error  = size(desmat,1)-size(desmat,2)-1;                         % Degrees of freedom for the error

if strcmp(Contrast_or_F,'Contrasts')
    if use_parfor
        if size(DV,2)==1
            if strcmp(reg_type,'OLS')
                if isempty(covariates)
                    if unique(desmat)==1
                        parfor perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV.*curr_perm;
                            temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'tstat'});
                            test_stats(perm) = temp.tstat.t(1);
                            p_stats(perm)    = temp.tstat.pval(1);
                        end
                    else
                        parfor perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV(curr_perm); %#ok<*PFBNS>
                            temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'tstat'});
                            test_stats(perm) = temp.tstat.t(1);
                            p_stats(perm)    = temp.tstat.pval(1);
                        end
                    end
                else
                    DVvals = regstats(DV,covariates,eye(size(covariates,2)),{'yhat','r'});
                    resids = DVvals.r;
                    yhat   = DVvals.yhat;
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+yhat;
                        temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'tstat'});
                        test_stats(perm) = temp.tstat.t(1);
                        p_stats(perm)    = temp.tstat.pval(1);
                    end
                end
            elseif strcmp(reg_type,'Robust')
                if isempty(covariates)
                    if unique(desmat)==1
                        parfor perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV.*curr_perm;
                            [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                            test_stats(perm) = stats.t(1);
                            p_stats(perm)    = stats.p(1);
                        end
                    else
                        parfor perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV(curr_perm);
                            [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                            test_stats(perm) = stats.t(1);
                            p_stats(perm)    = stats.p(1);
                        end
                    end
                else
                    [~,orig_stats] = robustfit_iterincrease(covariates,DV,'bisquare',4.685,'off');
                    resids         = orig_stats.resid;
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+(DV-resids);
                        [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                        test_stats(perm) = stats.t(1);
                        p_stats(perm)    = stats.p(1);
                    end
                end
            elseif strcmp(reg_type,'LTS')
                if isempty(covariates)
                    if unique(desmat)==1
                        parfor perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV.*curr_perm;
                            rew              = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                            test_stats(perm) = (rew.rsquared*(size(desmat,1)-size(desmat,2)))/(1-rew.rsquared);
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    else
                        parfor perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV(curr_perm);
                            rew              = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                            test_stats(perm) = (rew.rsquared*(size(desmat,1)-size(desmat,2)))/(1-rew.rsquared);
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    end
                else
                    [~,raw] = ltsregres(covariates,DV,'plots',0,'intercept',0);
                    resids  = raw.res;
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+(DV-resids);
                        rew_full         = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                        rew_redu         = ltsregres(covariates,permed_DV,'plots',0,'intercept',0);
                        test_stats(perm) = ((rew_full.rsquared-rew_redu.rsquared)*(size(desmat,1)-size(desmat,2)))/(1-rew_full.rsquared);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            end
        else
            K = eye(size(desmat,2));
            if isempty(covariates)
                if unique(desmat)==1
                    parfor perm = 1:num_perms
                        permed_DV        = [];
                        curr_perm        = perms(:,perm);
                        curr_within_perm = within_perms(:,:,perm);
                        for lev = size(DV,2):-1:1
                            permed_DV(:,lev) = DV(:,lev).*curr_perm;
                        end
                        permed_DV = permed_DV(curr_within_perm);
                        B         = (desmat'*desmat)\desmat'*permed_DV;
                        H         = M'*(K(1,:)*B)'/(K(1,:)/(desmat'*desmat)*K(1,:)')*(K(1,:)*B)*M; % Estimate hypothesis matrix
                        E         = M'*(permed_DV'*permed_DV-B'*(desmat'*desmat)*B)*M; % Estimate error matrix
                        if size(M,2)>1
                            test_stats(perm) = det(E)/det(H+E);                                         % Calculate Wilks' Lambda
                            if MC_corr==1
                                if size(M,2)==2
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(desmat,2)-1)/size(desmat,2));
                                    temp_dof_1 = size(desmat,2);
                                    temp_dof_2 = numel(permed_DV)-size(desmat,2)-1;
                                elseif size(M,2)==3
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(desmat,2)-2)/size(desmat,2));
                                    temp_dof_1 = 2*size(desmat,2);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(desmat,2)-2);
                                elseif size(desmat,2)==1
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(M,2))/(size(M,2)-1));
                                    temp_dof_1 = size(M,2)-1;
                                    temp_dof_2 = numel(permed_DV)-size(M,2);
                                elseif size(desmat,2)==2
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(M,2)-1)/(size(M,2)-1));
                                    temp_dof_1 = 2*(size(M,2)-1);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(M,2)-1);
                                else
                                    pw = rank(H+E);
                                    qw = rank(K(1,:)/(desmat'*desmat)*K(1,:)');
                                    rw = dof_error-((pw-qw+1)/2);
                                    uw = ((pw*qw)-2)/4;
                                    if ((pw^2)+(qw^2)-5)>0
                                        tw = sqrt((((pw^2)*(qw^2))-4)/((pw^2)+(qw^2)-5));
                                    else
                                        tw = 1;
                                    end
                                    temp_F     = ((1-(test_stats(perm)^(1/tw)))/(test_stats(perm)^(1/tw)))*(((rw*tw)-(2*uw))/(pw*qw));
                                    temp_dof_1 = pw*qw;
                                    temp_dof_2 = (rw*tw)-(2*uw);
                                end
                                p_stats(perm) = fcdf(temp_F,temp_dof_1,temp_dof_2,'upper');
                            end
                        else
                            test_stats(perm) = H/(E/dof_error);                                         % Calculate F
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    end
                else
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        curr_within_perm = within_perms(:,:,perm);
                        permed_DV        = DV(curr_perm,:);
                        permed_DV        = permed_DV(curr_within_perm);
                        B                = (desmat'*desmat)\desmat'*permed_DV;
                        H                = M'*(K(1,:)*B)'/(K(1,:)/(desmat'*desmat)*K(1,:)')*(K(1,:)*B)*M; % Estimate hypothesis matrix
                        E                = M'*(permed_DV'*permed_DV-B'*(desmat'*desmat)*B)*M; % Estimate error matrix
                        if size(M,2)>1
                            test_stats(perm) = det(E)/det(H+E);                                         % Calculate Wilks' Lambda
                            if MC_corr==1
                                if size(M,2)==2
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(desmat,2)-1)/size(desmat,2));
                                    temp_dof_1 = size(desmat,2);
                                    temp_dof_2 = numel(permed_DV)-size(desmat,2)-1;
                                elseif size(M,2)==3
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(desmat,2)-2)/size(desmat,2));
                                    temp_dof_1 = 2*size(desmat,2);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(desmat,2)-2);
                                elseif size(desmat,2)==1
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(M,2))/(size(M,2)-1));
                                    temp_dof_1 = size(M,2)-1;
                                    temp_dof_2 = numel(permed_DV)-size(M,2);
                                elseif size(desmat,2)==2
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(M,2)-1)/(size(M,2)-1));
                                    temp_dof_1 = 2*(size(M,2)-1);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(M,2)-1);
                                else
                                    pw = rank(H+E);
                                    qw = rank(K(1,:)/(desmat'*desmat)*K(1,:)');
                                    rw = dof_error-((pw-qw+1)/2);
                                    uw = ((pw*qw)-2)/4;
                                    if ((pw^2)+(qw^2)-5)>0
                                        tw = sqrt((((pw^2)*(qw^2))-4)/((pw^2)+(qw^2)-5));
                                    else
                                        tw = 1;
                                    end
                                    temp_F     = ((1-(test_stats(perm)^(1/tw)))/(test_stats(perm)^(1/tw)))*(((rw*tw)-(2*uw))/(pw*qw));
                                    temp_dof_1 = pw*qw;
                                    temp_dof_2 = (rw*tw)-(2*uw);
                                end
                                p_stats(perm) = fcdf(temp_F,temp_dof_1,temp_dof_2,'upper');
                            end
                        else
                            test_stats(perm) = H/(E/dof_error);                                         % Calculate F
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    end
                end
            else
                yhat  = covariates/(covariates'*covariates)*covariates'*DV;
                resid = DV-covariates/(covariates'*covariates)*covariates'*DV;
                parfor perm = 1:num_perms
                    curr_perm        = perms(:,perm);
                    curr_within_perm = within_perms(:,:,perm);
                    permed_resid     = resid(curr_perm,:);
                    permed_resid     = permed_resid(curr_within_perm);
                    permed_DV        = permed_resid+yhat;
                    B                = (desmat'*desmat)\desmat'*permed_DV;
                    H                = M'*(K(1,:)*B)'/(K(1,:)/(desmat'*desmat)*K(1,:)')*(K(1,:)*B)*M; % Estimate hypothesis matrix
                    E                = M'*(permed_DV'*permed_DV-B'*(desmat'*desmat)*B)*M; % Estimate error matrix
                    if size(M,2)>1
                        test_stats(perm) = det(E)/det(H+E);                                         % Calculate Wilks' Lambda
                        if MC_corr==1
                            if size(M,2)==2
                                temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(desmat,2)-1)/size(desmat,2));
                                temp_dof_1 = size(desmat,2);
                                temp_dof_2 = numel(permed_DV)-size(desmat,2)-1;
                            elseif size(M,2)==3
                                temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(desmat,2)-2)/size(desmat,2));
                                temp_dof_1 = 2*size(desmat,2);
                                temp_dof_2 = 2*(numel(permed_DV)-size(desmat,2)-2);
                            elseif size(desmat,2)==1
                                temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(M,2))/(size(M,2)-1));
                                temp_dof_1 = size(M,2)-1;
                                temp_dof_2 = numel(permed_DV)-size(M,2);
                            elseif size(desmat,2)==2
                                temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(M,2)-1)/(size(M,2)-1));
                                temp_dof_1 = 2*(size(M,2)-1);
                                temp_dof_2 = 2*(numel(permed_DV)-size(M,2)-1);
                            else
                                pw = rank(H+E);
                                qw = rank(K(1,:)/(desmat'*desmat)*K(1,:)');
                                rw = dof_error-((pw-qw+1)/2);
                                uw = ((pw*qw)-2)/4;
                                if ((pw^2)+(qw^2)-5)>0
                                    tw = sqrt((((pw^2)*(qw^2))-4)/((pw^2)+(qw^2)-5));
                                else
                                    tw = 1;
                                end
                                temp_F     = ((1-(test_stats(perm)^(1/tw)))/(test_stats(perm)^(1/tw)))*(((rw*tw)-(2*uw))/(pw*qw));
                                temp_dof_1 = pw*qw;
                                temp_dof_2 = (rw*tw)-(2*uw);
                            end
                            p_stats(perm) = fcdf(temp_F,temp_dof_1,temp_dof_2,'upper');
                        end
                    else
                        test_stats(perm) = H/(E/dof_error);                                         % Calculate F
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            end
        end
    else
        if size(DV,2)==1
            if strcmp(reg_type,'OLS')
                if isempty(covariates)
                    if unique(desmat)==1
                        for perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV.*curr_perm;
                            temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'tstat'});
                            test_stats(perm) = temp.tstat.t(1);
                            p_stats(perm)    = temp.tstat.pval(1);
                        end
                    else
                        for perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV(curr_perm);
                            temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'tstat'});
                            test_stats(perm) = temp.tstat.t(1);
                            p_stats(perm)    = temp.tstat.pval(1);
                        end
                    end
                else
                    DVvals = regstats(DV,covariates,eye(size(covariates,2)),{'yhat','r'});
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DVvals.r(curr_perm)+DVvals.yhat;
                        temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'tstat'});
                        test_stats(perm) = temp.tstat.t(1);
                        p_stats(perm)    = temp.tstat.pval(1);
                    end
                end
            elseif strcmp(reg_type,'Robust')
                if isempty(covariates)
                    if unique(desmat)==1
                        for perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV.*curr_perm;
                            [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                            test_stats(perm) = stats.t(1);
                            p_stats(perm)    = stats.p(1);
                        end
                    else
                        for perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV(curr_perm);
                            [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                            test_stats(perm) = stats.t(1);
                            p_stats(perm)    = stats.p(1);
                        end
                    end
                else
                    [~,orig_stats] = robustfit_iterincrease(covariates,DV,'bisquare',4.685,'off');
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = orig_stats.resid(curr_perm)+(DV-orig_stats.resid);
                        [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                        test_stats(perm) = stats.t(1);
                        p_stats(perm)    = stats.p(1);
                    end
                end
            elseif strcmp(reg_type,'LTS')
                if isempty(covariates)
                    if unique(desmat)==1
                        for perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV.*curr_perm;
                            rew              = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                            test_stats(perm) = (rew.rsquared*(size(desmat,1)-size(desmat,2)))/(1-rew.rsquared);
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    else
                        for perm = 1:num_perms
                            curr_perm        = perms(:,perm);
                            permed_DV        = DV(curr_perm);
                            rew              = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                            test_stats(perm) = (rew.rsquared*(size(desmat,1)-size(desmat,2)))/(1-rew.rsquared);
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    end
                else
                    [~,raw] = ltsregres(covariates,DV,'plots',0,'intercept',0);
                    resids  = raw.res;
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+(DV-resids);
                        rew_full         = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                        rew_redu         = ltsregres(covariates,permed_DV,'plots',0,'intercept',0);
                        test_stats(perm) = ((rew_full.rsquared-rew_redu.rsquared)*(size(desmat,1)-size(desmat,2)))/(1-rew_full.rsquared);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            end
        else
            K = eye(size(desmat,2));
            if isempty(covariates)
                if unique(desmat)==1
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        curr_within_perm = within_perms(:,:,perm);
                        for lev = size(DV,2):-1:1
                            permed_DV(:,lev) = DV(:,lev).*curr_perm;
                        end
                        permed_DV = permed_DV(curr_within_perm);
                        B         = (desmat'*desmat)\desmat'*permed_DV;
                        H         = M'*(K(1,:)*B)'/(K(1,:)/(desmat'*desmat)*K(1,:)')*(K(1,:)*B)*M; % Estimate hypothesis matrix
                        E         = M'*(permed_DV'*permed_DV-B'*(desmat'*desmat)*B)*M; % Estimate error matrix
                        if size(M,2)>1
                            test_stats(perm) = det(E)/det(H+E);                                         % Calculate Wilks' Lambda
                            if MC_corr==1
                                if size(M,2)==2
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(desmat,2)-1)/size(desmat,2));
                                    temp_dof_1 = size(desmat,2);
                                    temp_dof_2 = numel(permed_DV)-size(desmat,2)-1;
                                elseif size(M,2)==3
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(desmat,2)-2)/size(desmat,2));
                                    temp_dof_1 = 2*size(desmat,2);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(desmat,2)-2);
                                elseif size(desmat,2)==1
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(M,2))/(size(M,2)-1));
                                    temp_dof_1 = size(M,2)-1;
                                    temp_dof_2 = numel(permed_DV)-size(M,2);
                                elseif size(desmat,2)==2
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(M,2)-1)/(size(M,2)-1));
                                    temp_dof_1 = 2*(size(M,2)-1);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(M,2)-1);
                                else
                                    pw = rank(H+E);
                                    qw = rank(K(1,:)/(desmat'*desmat)*K(1,:)');
                                    rw = dof_error-((pw-qw+1)/2);
                                    uw = ((pw*qw)-2)/4;
                                    if ((pw^2)+(qw^2)-5)>0
                                        tw = sqrt((((pw^2)*(qw^2))-4)/((pw^2)+(qw^2)-5));
                                    else
                                        tw = 1;
                                    end
                                    temp_F     = ((1-(test_stats(perm)^(1/tw)))/(test_stats(perm)^(1/tw)))*(((rw*tw)-(2*uw))/(pw*qw));
                                    temp_dof_1 = pw*qw;
                                    temp_dof_2 = (rw*tw)-(2*uw);
                                end
                                p_stats(perm) = fcdf(temp_F,temp_dof_1,temp_dof_2,'upper');
                            end
                        else
                            test_stats(perm) = H/(E/dof_error);                                         % Calculate F
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    end
                else
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        curr_within_perm = within_perms(:,:,perm);
                        permed_DV        = DV(curr_perm,:);
                        permed_DV        = permed_DV(curr_within_perm);
                        B                = (desmat'*desmat)\desmat'*permed_DV;
                        H                = M'*(K(1,:)*B)'/(K(1,:)/(desmat'*desmat)*K(1,:)')*(K(1,:)*B)*M; % Estimate hypothesis matrix
                        E                = M'*(permed_DV'*permed_DV-B'*(desmat'*desmat)*B)*M; % Estimate error matrix
                        if size(M,2)>1
                            test_stats(perm) = det(E)/det(H+E);                                         % Calculate Wilks' Lambda
                            if MC_corr==1
                                if size(M,2)==2
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(desmat,2)-1)/size(desmat,2));
                                    temp_dof_1 = size(desmat,2);
                                    temp_dof_2 = numel(permed_DV)-size(desmat,2)-1;
                                elseif size(M,2)==3
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(desmat,2)-2)/size(desmat,2));
                                    temp_dof_1 = 2*size(desmat,2);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(desmat,2)-2);
                                elseif size(desmat,2)==1
                                    temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(M,2))/(size(M,2)-1));
                                    temp_dof_1 = size(M,2)-1;
                                    temp_dof_2 = numel(permed_DV)-size(M,2);
                                elseif size(desmat,2)==2
                                    temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(M,2)-1)/(size(M,2)-1));
                                    temp_dof_1 = 2*(size(M,2)-1);
                                    temp_dof_2 = 2*(numel(permed_DV)-size(M,2)-1);
                                else
                                    pw = rank(H+E);
                                    qw = rank(K(1,:)/(desmat'*desmat)*K(1,:)');
                                    rw = dof_error-((pw-qw+1)/2);
                                    uw = ((pw*qw)-2)/4;
                                    if ((pw^2)+(qw^2)-5)>0
                                        tw = sqrt((((pw^2)*(qw^2))-4)/((pw^2)+(qw^2)-5));
                                    else
                                        tw = 1;
                                    end
                                    temp_F     = ((1-(test_stats(perm)^(1/tw)))/(test_stats(perm)^(1/tw)))*(((rw*tw)-(2*uw))/(pw*qw));
                                    temp_dof_1 = pw*qw;
                                    temp_dof_2 = (rw*tw)-(2*uw);
                                end
                                p_stats(perm) = fcdf(temp_F,temp_dof_1,temp_dof_2,'upper');
                            end
                        else
                            test_stats(perm) = H/(E/dof_error);                                         % Calculate F
                            if MC_corr==1
                                p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                            end
                        end
                    end
                end
            else
                yhat  = covariates/(covariates'*covariates)*covariates'*DV;
                resid = DV-covariates/(covariates'*covariates)*covariates'*DV;
                for perm = 1:num_perms
                    curr_perm        = perms(:,perm);
                    curr_within_perm = within_perms(:,:,perm);
                    permed_resid     = resid(curr_perm,:);
                    permed_resid     = permed_resid(curr_within_perm);
                    permed_DV        = permed_resid+yhat;
                    B                = (desmat'*desmat)\desmat'*permed_DV;
                    H                = M'*(K(1,:)*B)'/(K(1,:)/(desmat'*desmat)*K(1,:)')*(K(1,:)*B)*M; % Estimate hypothesis matrix
                    E                = M'*(permed_DV'*permed_DV-B'*(desmat'*desmat)*B)*M; % Estimate error matrix
                    if size(M,2)>1
                        test_stats(perm) = det(E)/det(H+E);                                         % Calculate Wilks' Lambda
                        if MC_corr==1
                            if size(M,2)==2
                                temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(desmat,2)-1)/size(desmat,2));
                                temp_dof_1 = size(desmat,2);
                                temp_dof_2 = numel(permed_DV)-size(desmat,2)-1;
                            elseif size(M,2)==3
                                temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(desmat,2)-2)/size(desmat,2));
                                temp_dof_1 = 2*size(desmat,2);
                                temp_dof_2 = 2*(numel(permed_DV)-size(desmat,2)-2);
                            elseif size(desmat,2)==1
                                temp_F     = ((1-test_stats(perm))/test_stats(perm))*((numel(permed_DV)-size(M,2))/(size(M,2)-1));
                                temp_dof_1 = size(M,2)-1;
                                temp_dof_2 = numel(permed_DV)-size(M,2);
                            elseif size(desmat,2)==2
                                temp_F     = ((1-sqrt(test_stats(perm)))/sqrt(test_stats(perm)))*((numel(permed_DV)-size(M,2)-1)/(size(M,2)-1));
                                temp_dof_1 = 2*(size(M,2)-1);
                                temp_dof_2 = 2*(numel(permed_DV)-size(M,2)-1);
                            else
                                pw = rank(H+E);
                                qw = rank(K(1,:)/(desmat'*desmat)*K(1,:)');
                                rw = dof_error-((pw-qw+1)/2);
                                uw = ((pw*qw)-2)/4;
                                if ((pw^2)+(qw^2)-5)>0
                                    tw = sqrt((((pw^2)*(qw^2))-4)/((pw^2)+(qw^2)-5));
                                else
                                    tw = 1;
                                end
                                temp_F     = ((1-(test_stats(perm)^(1/tw)))/(test_stats(perm)^(1/tw)))*(((rw*tw)-(2*uw))/(pw*qw));
                                temp_dof_1 = pw*qw;
                                temp_dof_2 = (rw*tw)-(2*uw);
                            end
                            p_stats(perm) = fcdf(temp_F,temp_dof_1,temp_dof_2,'upper');
                        end
                    else
                        test_stats(perm) = H/(E/dof_error);                                         % Calculate F
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),1,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            end
        end
    end
else
    num_preds_in_F = size(desmat,2)-size(covariates,2);
    if use_parfor
        if size(DV,2)==1
            if strcmp(reg_type,'OLS')
                if isempty(covariates)
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DV(curr_perm); %#ok<*PFBNS>
                        temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'fstat'});
                        test_stats(perm) = temp.fstat.f(1);
                        p_stats(perm)    = temp.fstat.pval(1);
                    end
                else
                    DVvals = regstats(DV,covariates,eye(size(covariates,2)),{'yhat','r'});
                    resids = DVvals.r;
                    yhat   = DVvals.yhat;
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+yhat;
                        temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'rsquare'});
                        Rsq_full         = temp.rsquare(1);
                        temp             = regstats(permed_DV,covariates,eye(size(covariates,2)),{'rsquare'});
                        Rsq_redu         = temp.rsquare(1);
                        test_stats(perm) = ((Rsq_full-Rsq_redu)*(size(desmat,1)-size(desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            elseif strcmp(reg_type,'Robust')
                if isempty(covariates)
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DV(curr_perm);
                        [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                        Rsq              = 1-(corr(permed_DV,stats.resid))^2;
                        test_stats(perm) = ((Rsq)*(size(desmat,1)-num_preds_in_F-1))/((1-Rsq)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-num_preds_in_F-1),'upper');
                        end
                    end
                else
                    [~,orig_stats] = robustfit_iterincrease(covariates,DV,'bisquare',4.685,'off');
                    resids         = orig_stats.resid;
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+(DV-resids);
                        [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                        Rsq_full         = 1-(corr(permed_DV,stats.resid))^2;
                        [~,stats]        = robustfit_iterincrease(covariates,permed_DV,'bisquare',4.685,'off');
                        Rsq_redu         = 1-(corr(permed_DV,stats.resid))^2;
                        test_stats(perm) = ((Rsq_full-Rsq_redu)*(size(desmat,1)-size(desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            elseif strcmp(reg_type,'LTS')
                if isempty(covariates)
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DV(curr_perm);
                        rew              = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                        test_stats(perm) = (rew.rsquared*(size(desmat,1)-num_preds_in_F-1))/((1-rew.rsquared)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-num_preds_in_F-1),'upper');
                        end
                    end
                else
                    [~,raw] = ltsregres(covariates,DV,'plots',0,'intercept',0);
                    resids  = raw.res;
                    parfor perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+(DV-resids);
                        rew_full         = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                        rew_redu         = ltsregres(covariates,permed_DV,'plots',0,'intercept',0);
                        test_stats(perm) = ((rew_full.rsquared-rew_redu.rsquared)*(size(desmat,1)-size(desmat,2)))/((1-rew_full.rsquared)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            end
        else
            if isempty(covariates)
                parfor perm = 1:num_perms
                    curr_perm        = perms(:,perm);
                    curr_within_perm = within_perms(:,:,perm);
                    permed_DV        = DV(curr_perm,:);
                    permed_DV        = permed_DV(curr_within_perm);
                    B                = (desmat'*desmat)\desmat'*permed_DV;
                    Rsq              = 1-(corr(permed_DV*M,(permed_DV*M-desmat*B*M)))^2;
                    test_stats(perm) = ((Rsq)*(size(desmat,1)-num_preds_in_F-1))/((1-Rsq)*num_preds_in_F);
                    if MC_corr==1
                        p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-num_preds_in_F-1),'upper');
                    end
                end
            else
                yhat  = covariates/(covariates'*covariates)*covariates'*DV;
                resid = DV-covariates/(covariates'*covariates)*covariates'*DV;
                parfor perm = 1:num_perms
                    curr_perm        = perms(:,perm);
                    curr_within_perm = within_perms(:,:,perm);
                    permed_resid     = resid(curr_perm,:);
                    permed_resid     = permed_resid(curr_within_perm);
                    permed_DV        = permed_resid+yhat;
                    B_full           = (desmat'*desmat)\desmat'*permed_DV;
                    B_redu           = (covariates'*covariates)\covariates'*permed_DV;
                    Rsq_full         = 1-(corr(permed_DV*M,(permed_DV*M-desmat*B_full*M)))^2;
                    Rsq_redu         = 1-(corr(permed_DV*M,(permed_DV*M-covariates*B_redu*M)))^2;
                    test_stats(perm) = ((Rsq_full-Rsq_redu)*(size(desmat,1)-size(desmat,2)))/((1-Rsq_full)*num_preds_in_F);                                         % Calculate F
                    if MC_corr==1
                        p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                    end
                end
            end
        end
    else
        if size(DV,2)==1
            if strcmp(reg_type,'OLS')
                if isempty(covariates)
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DV(curr_perm);
                        temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'fstat'});
                        test_stats(perm) = temp.fstat.f(1);
                        p_stats(perm)    = temp.fstat.pval(1);
                    end
                else
                    DVvals = regstats(DV,covariates,eye(size(covariates,2)),{'yhat','r'});
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DVvals.r(curr_perm)+DVvals.yhat;
                        temp             = regstats(permed_DV,desmat,eye(size(desmat,2)),{'rsquare'});
                        Rsq_full         = temp.rsquare(1);
                        temp             = regstats(permed_DV,covariates,eye(size(covariates,2)),{'rsquare'});
                        Rsq_redu         = temp.rsquare(1);
                        test_stats(perm) = ((Rsq_full-Rsq_redu)*(size(desmat,1)-size(desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            elseif strcmp(reg_type,'Robust')
                if isempty(covariates)
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DV(curr_perm);
                        [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                        Rsq              = 1-(corr(permed_DV,stats.resid))^2;
                        test_stats(perm) = ((Rsq)*(size(desmat,1)-num_preds_in_F-1))/((1-Rsq)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-num_preds_in_F-1),'upper');
                        end
                    end
                else
                    [~,orig_stats] = robustfit_iterincrease(covariates,DV,'bisquare',4.685,'off');
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = orig_stats.resid(curr_perm)+(DV-orig_stats.resid);
                        [~,stats]        = robustfit_iterincrease(desmat,permed_DV,'bisquare',4.685,'off');
                        Rsq_full         = 1-(corr(permed_DV,stats.resid))^2;
                        [~,stats]        = robustfit_iterincrease(covariates,permed_DV,'bisquare',4.685,'off');
                        Rsq_redu         = 1-(corr(permed_DV,stats.resid))^2;
                        test_stats(perm) = ((Rsq_full-Rsq_redu)*(size(desmat,1)-size(desmat,2)))/((1-Rsq_full)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            elseif strcmp(reg_type,'LTS')
                if isempty(covariates)
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = DV(curr_perm);
                        rew              = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                        test_stats(perm) = (rew.rsquared*(size(desmat,1)-num_preds_in_F-1))/((1-rew.rsquared)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-num_preds_in_F-1),'upper');
                        end
                    end
                else
                    [~,raw] = ltsregres(covariates,DV,'plots',0,'intercept',0);
                    resids  = raw.res;
                    for perm = 1:num_perms
                        curr_perm        = perms(:,perm);
                        permed_DV        = resids(curr_perm)+(DV-resids);
                        rew_full         = ltsregres(desmat,permed_DV,'plots',0,'intercept',0);
                        rew_redu         = ltsregres(covariates,permed_DV,'plots',0,'intercept',0);
                        test_stats(perm) = ((rew_full.rsquared-rew_redu.rsquared)*(size(desmat,1)-size(desmat,2)))/((1-rew_full.rsquared)*num_preds_in_F);
                        if MC_corr==1
                            p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                        end
                    end
                end
            end
        else
            if isempty(covariates)
                for perm = 1:num_perms
                    curr_perm        = perms(:,perm);
                    curr_within_perm = within_perms(:,:,perm);
                    permed_DV        = DV(curr_perm,:);
                    permed_DV        = permed_DV(curr_within_perm);
                    B                = (desmat'*desmat)\desmat'*permed_DV;
                    Rsq              = 1-(corr2(permed_DV*M,(permed_DV*M-desmat*B*M)))^2;
                    test_stats(perm) = ((Rsq)*(size(desmat,1)-num_preds_in_F-1))/((1-Rsq)*num_preds_in_F);
                    if MC_corr==1
                        p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-num_preds_in_F-1),'upper');
                    end
                end
            else
                yhat  = covariates/(covariates'*covariates)*covariates'*DV;
                resid = DV-covariates/(covariates'*covariates)*covariates'*DV;
                for perm = 1:num_perms
                    curr_perm        = perms(:,perm);
                    curr_within_perm = within_perms(:,:,perm);
                    permed_resid     = resid(curr_perm,:);
                    permed_resid     = permed_resid(curr_within_perm);
                    permed_DV        = permed_resid+yhat;
                    B_full           = (desmat'*desmat)\desmat'*permed_DV;
                    B_redu           = (covariates'*covariates)\covariates'*permed_DV;
                    Rsq_full         = 1-(corr2(permed_DV*M,(permed_DV*M-desmat*B_full*M)))^2;
                    Rsq_redu         = 1-(corr2(permed_DV*M,(permed_DV*M-covariates*B_redu*M)))^2;
                    test_stats(perm) = ((Rsq_full-Rsq_redu)*(size(desmat,1)-size(desmat,2)))/((1-Rsq_full)*num_preds_in_F);                                         % Calculate F
                    if MC_corr==1
                        p_stats(perm) = fcdf(test_stats(perm),num_preds_in_F,(size(desmat,1)-size(desmat,2)),'upper');
                    end
                end
            end
        end
    end
end
test_stats = sort(test_stats);
