function S = HumphriesGurney_smallworldness_bu(inmat,varargin)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 06.27.14
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2015. Software can be modified and redistributed, but modifed, 
% redistributed versions must have the same rights

num_rand_nets = 1000;

if nargin > 1
    for arg = 1:nargin
        switch varargin{arg}
            case '-num_rand_nets'
                num_rand_nets = varargin{arg+1};
            case '-rand_data'
                mean_rand_clustcoef = varargin{arg+1};
                mean_rand_charpath = varargin{arg+2};
        end
    end
end

if length(unique(inmat)) > 2
    inmat = weight_conversion(inmat,'binarize');
end

in_clustcoef = mean(clustering_coef_bu(inmat));
in_charpath = charpath(distance_bin(inmat));

if ~exist('mean_rand_clustcoef','var')
    num_nodes = size(inmat,1);
    num_edges = sum(sum(triu(inmat)));
    
    rand_clustcoef = zeros(num_rand_nets,1);
    rand_charpath = zeros(num_rand_nets,1);
    
    for q = 1:num_rand_nets
        rand_net = makerandCIJ_und(num_nodes,num_edges);
        rand_clustcoef(q) = mean(clustering_coef_bu(rand_net));
        rand_charpath(q) = charpath(distance_bin(rand_net));
    end
    
    mean_rand_clustcoef = mean(rand_clustcoef);
    mean_rand_charpath = mean(rand_charpath);
end

S = (in_clustcoef/mean_rand_clustcoef)/(in_charpath/mean_rand_charpath);
