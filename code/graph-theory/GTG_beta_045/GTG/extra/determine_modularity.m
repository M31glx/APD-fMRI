function [mod_grps,Q] = determine_modularity(conmats,num_mod_runs,gamma,B)

% Determine modular organization of matrix. Basicaly a wrapper for 
% community_louvain.m
% 
% conmats: one (usually mean) or more connectivity matrices. If multiple
%          matrices are given, the mean will be calculated and used to
%          determine modularity. Multiple matrices must be P x P x N, where
%          P = the number of nodes and N = the number of matrices
% 
% num_mod_runs: (optional, default = 5000) # of iterations (run w/highest
%               modularity chosen)
%
% gamma: (optional, default = 1) resolution parameter 
%        gamma>1,        detects smaller modules
%        0<=gamma<1,     detects larger modules
%        gamma=1,        classic modularity
% 
% B: (optional, default = 'negative_asym') objective-function type or 
%    custom objective matrix
%       'modularity'    = modularity
%       'potts'         = Potts-model Hamiltonian (for binary networks)
%       'negative_sym'  = symmetric treatment of negative weights
%       'negative_asym' = asymmetric treatment of negative weights
%       B               = custom objective-function matrix
% 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 03.09.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if ~exist('num_mod_runs','var')
    num_mod_runs = 5000;
end

if ~exist('gamma','var')
    gamma = 1;
end

if ~exist('B','var')
    B = 'negative_asym';
end

if size(conmats,3) > 1
    fish_z = atanh(conmats);
    full_mean_conmat = atan(mean(fish_z,3));
else
    full_mean_conmat = conmats;
end

for run = num_mod_runs:-1:1                                                         % Loop on run
    % Compute starting point partition
    M = community_louvain(full_mean_conmat,gamma,[],B);                      % Calculate the initial organization
    
    % Initialize modularity values
    Q0 = -1;
    Q1 = 0;
    
    % Iteratively refine partition
    while Q1-Q0>1e-5
        Q0     = Q1;
        [M,Q1] = community_louvain(full_mean_conmat,gamma,M,B);
    end
    temp_grps(:,run) = M;
    Q(run)           = Q1;
end

mod_grps = temp_grps(:,find(Q==max(Q),1,'first'));
Q        = max(Q);
