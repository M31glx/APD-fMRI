function [MeanSizeT,ixT]=TargetAttack_AllNetMeasures(g,InMes,OutMes)
%This function do targtted attack analysis on the input binary graph
%The inputs are:
%G(binary graph),
%InMes: the measure based on which the nodes are removed: betw, deg; dist; clust, leff
%OutMes: the measure for which you want to see the effect of attack: comp, dist, deg, assort, dens, clust, trans, geff, leff, mod, modl, charpath, node_betw, edge_betw, lambda, gamma, sigma  

%% Hadi Hosseini, June 2011

switch InMes
    case 'betw'
        %%%order by node betwenness
        Mes=betweenness_bin(g);%node-betw
        [ordT,ixT]=sort(Mes,'descend');
    
    case 'deg'
        Mes=degrees_und(g);%degree
        [ordT,ixT]=sort(Mes,'descend');
    
    case 'dist'
         Mes=distance_bin(g);%distance
         Mes=mean(Mes);
         [ordT,ixT]=sort(Mes,'ascend');
    
    case 'clust'
        Mes=clustering_coef_bu(g);%clustering
        [ordT,ixT]=sort(Mes,'descend');
        ixT=ixT';ordT=ordT';
    
    case 'leff'
        Mes=efficiency(g,1);%local eff
        [ordT,ixT]=sort(Mes,'descend');
        ixT=ixT';ordT=ordT';
end


        
%%%%
zT=size(g,1);
%%%%
switch OutMes
    case 'comp' %largest components
        sz_fullT=size(largest_component(g),2);%largest connected component
    case 'dist'
        sz_fullT=mean(mean(distance_bin(g)));
    case 'deg'
        sz_fullT=mean(degrees_und(g));
    case 'assort'
        sz_fullT=assortativity(g,0);
    case 'dens'
        sz_fullT=density_und(g);
    case 'clust'
        sz_fullT=mean(clustering_coef_bu(g));
    case 'trans'
        sz_fullT=transitivity_bu(g);
    case 'geff'
        sz_fullT=efficiency(g);
    case 'leff'
        sz_fullT=mean(efficiency(g,1));
    case 'mod'
        sz_fullT=modularity_und(g);
    case 'modl'
        sz_fullT=modularity_louvain_und(g);
    case 'charpath'
        sz_fullT=charpath(distance_bin(g));
    case 'node_betw'
        sz_fullT=mean(betweenness_bin(g));
    %case 'edge_betw'
    %    sz_fullT=mean(edge_betweenness_bin(g));
    case 'lambda'
        [PPP,PAPA,sz_fullT] = conn_network_mindist(g);
    case 'gamma'
        [CCC,CACA,sz_fullT] = conn_network_clustering(g);
    case 'sigma'
        [PPP,PAPA,Lam] = conn_network_mindist(g);
        [CCC,CACA,Gam] = conn_network_clustering(g);
        sz_fullT =  Gam/Lam;
end

szT=[];
ii=0;
for i=ixT%in the order of betweenness
    ii=ii+1;
        gg=g;%a copy of the original graph
        InxT=ixT(1:ii);
        gg(:,InxT)=0;gg(InxT,:)=0;

        %calculate the size of the largest remained connected component
        
        switch OutMes
            case 'comp' %largest components
                szT=[szT;size(largest_component(gg),2)];%largest connected component
            case 'dist'
                szT=[szT;mean(mean(distance_bin(gg)))];
            case 'deg'
                szT=[szT;mean(degrees_und(gg))];
            case 'assort'
                szT=[szT;assortativity(gg,0)];
            case 'dens'
                szT=[szT;density_und(gg)];
            case 'clust'
                szT=[szT;mean(clustering_coef_bu(gg))];
            case 'trans'
                szT=[szT;transitivity_bu(gg)];
            case 'geff'
                szT=[szT;efficiency(gg)];
            case 'leff'
                szT=[szT;mean(efficiency(gg,1))];
            case 'mod'
                szT=[szT;modularity_und(gg)];
            case 'modl'
                szT=[szT;modularity_louvain_und(gg)];
            case 'charpath'
                szT=[szT;charpath(distance_bin(gg))];
            case 'node_betw'
                szT=[szT;mean(betweenness_bin(gg))];
            case 'edge_betw'
                szT=[szT;mean(edge_betweenness_bin(gg))];
            case 'lambda'
                [PPP,PAPA,LAM] = conn_network_mindist(gg);
                LAM(isnan(LAM))=inf;
                szT=[szT;LAM];
                
            case 'gamma'
                [CCC,CACA,GAM] = conn_network_clustering(gg);
                GAM(isnan(GAM))=inf;
                szT=[szT;GAM];
            case 'sigma'
                [PPP,PAPA,LAM] = conn_network_mindist(gg);
                [CCC,CACA,GAM] = conn_network_clustering(gg);
                SIG = GAM./LAM; SIG(isnan(SIG))=inf;
                szT=[szT;SIG];
        end
               
end

MeanSizeT=szT./sz_fullT;
