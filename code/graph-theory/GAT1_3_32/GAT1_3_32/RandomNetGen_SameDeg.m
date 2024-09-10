function [gam,lam,sig] = RandomNetGen_SameDeg(inp,iter)

%   gets a binary graph ("inp") and generates "iter" number of random graphs
%   with the same total degree and degree distribution. Then it obtain the
%   network measures for each of the "iter" random graphs and average them.


%% Hadi Hosseini July 2011
%  
%   Ref. van der Heuvel, J Neuorsci 2009


%%

s=size(inp,1);
inp(1:s+1:end)=0;
inp(inp<0)=0;
inp(inp>0)=1;


%% BCT random graphs
clust=[];pathl=[];sworld=[];

for o=1:iter
    
    R = randomizer_bin_und(inp, 20);
    clust=[clust;mean(clustering_coef_bu(R))];
    pathl=[pathl;charpath(distance_bin(R))];

end

%% calculating gamma, lambda, and sigma
gam = mean(clustering_coef_bu(inp))/mean(clust);
lam = charpath(distance_bin(inp))/mean(pathl);
sig = gam/lam;
  

