% This code tests the significance of the difference in area under the
% curve (auc) between two groups net measures. the auc is calculated based
% on several density criteria: 1. auc between MinDens and MaxDens; 2. auce
% between MinDens of full connectivity to Max Dens; 3. askes for the
% density interval that you are interested in

%%   
data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Net Measures)');
ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;
Group2 = ff.mat4GAT.g2;

if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha') 
        
    Alpha = ff.mat4GAT.alpha;
    Tail = ff.mat4GAT.tail;

else

    Alpha = .05;
    Tail = 2;

end
    

NetMesPerDens1 = spm_select(1,['NetMesPerDens_OrigNet_' Group1],'Select the "NetMesPerDens_OrigNet_Group1.mat" (output from "Compare Graphs Across Densities")');
NetMesPerDens2 = spm_select(1,['NetMesPerDens_OrigNet_' Group2],'Select the "NetMesPerDens_OrigNet_Group2.mat" (output from "Compare Graphs Across Densities")');
NetMesPerDens1_rand = spm_select(1,['NetMesPerDens_' Group1],'Select the "NetMesPerDens_Group1.mat" (output from "Compare Graphs Across Densities")');
NetMesPerDens2_rand = spm_select(1,['NetMesPerDens_' Group2],'Select the "NetMesPerDens_Group2.mat" (output from "Compare Graphs Across Densities")');
DensInt = spm_select(1,'DensityInterval*','Select the "DensityInterval.mat" (output from "Compare Graphs Across Densities")');


%% 

MinDens = input('type the min density of interest: ');
MaxDens = input('type the max density of interest: ');
DensStep = input('type the density interval (step): ');
iXax=[MinDens:DensStep:MaxDens];

%% 
xx = load(DensInt,'Xax');Xax = xx.Xax;

f1 = load(NetMesPerDens1,'MClust_1');MClust_1 = f1.MClust_1;
f1 = load(NetMesPerDens1,'Trans_1');Trans_1 = f1.Trans_1;
f1 = load(NetMesPerDens1,'Assort_1');Assort_1 = f1.Assort_1;
f1 = load(NetMesPerDens1,'GEff_1');GEff_1 = f1.GEff_1;
f1 = load(NetMesPerDens1,'MLocEff_1');MLocEff_1 = f1.MLocEff_1;
f1 = load(NetMesPerDens1,'Mod_1');Mod_1 = f1.Mod_1;
f1 = load(NetMesPerDens1,'ModL_1');ModL_1 = f1.ModL_1;
f1 = load(NetMesPerDens1,'PathL_1');PathL_1 = f1.PathL_1;
f1 = load(NetMesPerDens1,'MNodeBetw_1');MNodeBetw_1 = f1.MNodeBetw_1;
f1 = load(NetMesPerDens1,'MEdgeBetw_1');MEdgeBetw_1 = f1.MEdgeBetw_1;
f1 = load(NetMesPerDens1,'Lambda_1');Lambda_1 = f1.Lambda_1;
f1 = load(NetMesPerDens1,'Gamma_1');Gamma_1 = f1.Gamma_1;
f1 = load(NetMesPerDens1,'Sigma_1');Sigma_1 = f1.Sigma_1;

f2 = load(NetMesPerDens2,'MClust_2');MClust_2 = f2.MClust_2;
f2 = load(NetMesPerDens2,'Trans_2');Trans_2 = f2.Trans_2;
f2 = load(NetMesPerDens2,'Assort_2');Assort_2 = f2.Assort_2;
f2 = load(NetMesPerDens2,'GEff_2');GEff_2 = f2.GEff_2;
f2 = load(NetMesPerDens2,'MLocEff_2');MLocEff_2 = f2.MLocEff_2;
f2 = load(NetMesPerDens2,'Mod_2');Mod_2 = f2.Mod_2;
f2 = load(NetMesPerDens2,'ModL_2');ModL_2 = f2.ModL_2;
f2 = load(NetMesPerDens2,'PathL_2');PathL_2 = f2.PathL_2;
f2 = load(NetMesPerDens2,'MNodeBetw_2');MNodeBetw_2 = f2.MNodeBetw_2;
f2 = load(NetMesPerDens2,'MEdgeBetw_2');MEdgeBetw_2 = f2.MEdgeBetw_2;
f2 = load(NetMesPerDens2,'Lambda_2');Lambda_2 = f2.Lambda_2;
f2 = load(NetMesPerDens2,'Gamma_2');Gamma_2 = f2.Gamma_2;
f2 = load(NetMesPerDens2,'Sigma_2');Sigma_2 = f2.Sigma_2;

f1_rand = load(NetMesPerDens1_rand,'MClust_1_rand');MClust_1_rand = f1_rand.MClust_1_rand;
f1_rand = load(NetMesPerDens1_rand,'Trans_1_rand');Trans_1_rand = f1_rand.Trans_1_rand;
f1_rand = load(NetMesPerDens1_rand,'Assort_1_rand');Assort_1_rand = f1_rand.Assort_1_rand;
f1_rand = load(NetMesPerDens1_rand,'GEff_1_rand');GEff_1_rand = f1_rand.GEff_1_rand;
f1_rand = load(NetMesPerDens1_rand,'MLocEff_1_rand');MLocEff_1_rand = f1_rand.MLocEff_1_rand;
f1_rand = load(NetMesPerDens1_rand,'Mod_1_rand');Mod_1_rand = f1_rand.Mod_1_rand;
f1_rand = load(NetMesPerDens1_rand,'ModL_1_rand');ModL_1_rand = f1_rand.ModL_1_rand;
f1_rand = load(NetMesPerDens1_rand,'PathL_1_rand');PathL_1_rand = f1_rand.PathL_1_rand;
f1_rand = load(NetMesPerDens1_rand,'MNodeBetw_1_rand');MNodeBetw_1_rand = f1_rand.MNodeBetw_1_rand;
f1_rand = load(NetMesPerDens1_rand,'MEdgeBetw_1_rand');MEdgeBetw_1_rand = f1_rand.MEdgeBetw_1_rand;
f1_rand = load(NetMesPerDens1_rand,'Lambda_1_rand');Lambda_1_rand = f1_rand.Lambda_1_rand;
f1_rand = load(NetMesPerDens1_rand,'Gamma_1_rand');Gamma_1_rand = f1_rand.Gamma_1_rand;
f1_rand = load(NetMesPerDens1_rand,'Sigma_1_rand');Sigma_1_rand = f1_rand.Sigma_1_rand;

f2_rand = load(NetMesPerDens2_rand,'MClust_2_rand');MClust_2_rand = f2_rand.MClust_2_rand;
f2_rand = load(NetMesPerDens2_rand,'Trans_2_rand');Trans_2_rand = f2_rand.Trans_2_rand;
f2_rand = load(NetMesPerDens2_rand,'Assort_2_rand');Assort_2_rand = f2_rand.Assort_2_rand;
f2_rand = load(NetMesPerDens2_rand,'GEff_2_rand');GEff_2_rand = f2_rand.GEff_2_rand;
f2_rand = load(NetMesPerDens2_rand,'MLocEff_2_rand');MLocEff_2_rand = f2_rand.MLocEff_2_rand;
f2_rand = load(NetMesPerDens2_rand,'Mod_2_rand');Mod_2_rand = f2_rand.Mod_2_rand;
f2_rand = load(NetMesPerDens2_rand,'ModL_2_rand');ModL_2_rand = f2_rand.ModL_2_rand;
f2_rand = load(NetMesPerDens2_rand,'PathL_2_rand');PathL_2_rand = f2_rand.PathL_2_rand;
f2_rand = load(NetMesPerDens2_rand,'MNodeBetw_2_rand');MNodeBetw_2_rand = f2_rand.MNodeBetw_2_rand;
f2_rand = load(NetMesPerDens2_rand,'MEdgeBetw_2_rand');MEdgeBetw_2_rand = f2_rand.MEdgeBetw_2_rand;
f2_rand = load(NetMesPerDens2_rand,'Lambda_2_rand');Lambda_2_rand = f2_rand.Lambda_2_rand;
f2_rand = load(NetMesPerDens2_rand,'Gamma_2_rand');Gamma_2_rand = f2_rand.Gamma_2_rand;
f2_rand = load(NetMesPerDens2_rand,'Sigma_2_rand');Sigma_2_rand = f2_rand.Sigma_2_rand;


%% 
if ~isequal(MinDens,Xax(1))
   
    Precision = .1;Ind = find(Xax >= MinDens-Precision & Xax <= MinDens+Precision );
    Dind = Xax(Ind)-MinDens;
    Dind = Xax(Ind)-MinDens;
    Ind_Ind = find(Dind == min(abs(Dind)) | Dind == -min(abs(Dind)));
    Ind_min = Ind(min(Ind_Ind));
    
else
    
    Ind_min=1;
    
end

if ~isequal(MaxDens,Xax(end))
    Precision = .1;Ind = find(Xax >= MaxDens-Precision & Xax <= MaxDens+Precision );
    Dind = Xax(Ind)-MaxDens;
    Dind = Xax(Ind)-MaxDens;
    Ind_Ind = find(Dind == min(abs(Dind)) | Dind == -min(abs(Dind)));
    Ind_max = Ind(min(Ind_Ind));
    
else
    
    Ind_max=size(Xax,2);
    
end

%% 

auc_MClust_1 = trapz(iXax,MClust_1(Ind_min:Ind_max));
auc_Trans_1 = trapz(iXax,Trans_1(Ind_min:Ind_max));
auc_Assort_1 = trapz(iXax,Assort_1(Ind_min:Ind_max));
auc_GEff_1 = trapz(iXax,GEff_1(Ind_min:Ind_max));
auc_MLocEff_1 = trapz(iXax,MLocEff_1(Ind_min:Ind_max));
auc_Mod_1 = trapz(iXax,Mod_1(Ind_min:Ind_max));
auc_ModL_1 = trapz(iXax,ModL_1(Ind_min:Ind_max));
auc_PathL_1 = trapz(iXax,PathL_1(Ind_min:Ind_max));
auc_MNodeBetw_1 = trapz(iXax,MNodeBetw_1(Ind_min:Ind_max));
auc_MEdgeBetw_1 = trapz(iXax,MEdgeBetw_1(Ind_min:Ind_max));
auc_Lambda_1 = trapz(iXax,Lambda_1(Ind_min:Ind_max));
auc_Gamma_1 = trapz(iXax,Gamma_1(Ind_min:Ind_max));
auc_Sigma_1 = trapz(iXax,Sigma_1(Ind_min:Ind_max));

auc_MClust_2 = trapz(iXax,MClust_2(Ind_min:Ind_max));
auc_Trans_2 = trapz(iXax,Trans_2(Ind_min:Ind_max));
auc_Assort_2 = trapz(iXax,Assort_2(Ind_min:Ind_max));
auc_GEff_2 = trapz(iXax,GEff_2(Ind_min:Ind_max));
auc_MLocEff_2 = trapz(iXax,MLocEff_2(Ind_min:Ind_max));
auc_Mod_2 = trapz(iXax,Mod_2(Ind_min:Ind_max));
auc_ModL_2 = trapz(iXax,ModL_2(Ind_min:Ind_max));
auc_PathL_2 = trapz(iXax,PathL_2(Ind_min:Ind_max));
auc_MNodeBetw_2 = trapz(iXax,MNodeBetw_2(Ind_min:Ind_max));
auc_MEdgeBetw_2 = trapz(iXax,MEdgeBetw_2(Ind_min:Ind_max));
auc_Lambda_2 = trapz(iXax,Lambda_2(Ind_min:Ind_max));
auc_Gamma_2 = trapz(iXax,Gamma_2(Ind_min:Ind_max));
auc_Sigma_2 = trapz(iXax,Sigma_2(Ind_min:Ind_max));

auc_MClust_1_rand = trapz(iXax,MClust_1_rand(Ind_min:Ind_max,:));
auc_Trans_1_rand = trapz(iXax,Trans_1_rand(Ind_min:Ind_max,:));
auc_Assort_1_rand = trapz(iXax,Assort_1_rand(Ind_min:Ind_max,:));
auc_GEff_1_rand = trapz(iXax,GEff_1_rand(Ind_min:Ind_max,:));
auc_MLocEff_1_rand = trapz(iXax,MLocEff_1_rand(Ind_min:Ind_max,:));
auc_Mod_1_rand = trapz(iXax,Mod_1_rand(Ind_min:Ind_max,:));
auc_ModL_1_rand = trapz(iXax,ModL_1_rand(Ind_min:Ind_max,:));
auc_PathL_1_rand = trapz(iXax,PathL_1_rand(Ind_min:Ind_max,:));
auc_MNodeBetw_1_rand = trapz(iXax,MNodeBetw_1_rand(Ind_min:Ind_max,:));
auc_MEdgeBetw_1_rand = trapz(iXax,MEdgeBetw_1_rand(Ind_min:Ind_max,:));
auc_Lambda_1_rand = trapz(iXax,Lambda_1_rand(Ind_min:Ind_max,:));
auc_Gamma_1_rand = trapz(iXax,Gamma_1_rand(Ind_min:Ind_max,:));
auc_Sigma_1_rand = trapz(iXax,Sigma_1_rand(Ind_min:Ind_max,:));

auc_MClust_2_rand = trapz(iXax,MClust_2_rand(Ind_min:Ind_max,:));
auc_Trans_2_rand = trapz(iXax,Trans_2_rand(Ind_min:Ind_max,:));
auc_Assort_2_rand = trapz(iXax,Assort_2_rand(Ind_min:Ind_max,:));
auc_GEff_2_rand = trapz(iXax,GEff_2_rand(Ind_min:Ind_max,:));
auc_MLocEff_2_rand = trapz(iXax,MLocEff_2_rand(Ind_min:Ind_max,:));
auc_Mod_2_rand = trapz(iXax,Mod_2_rand(Ind_min:Ind_max,:));
auc_ModL_2_rand = trapz(iXax,ModL_2_rand(Ind_min:Ind_max,:));
auc_PathL_2_rand = trapz(iXax,PathL_2_rand(Ind_min:Ind_max,:));
auc_MNodeBetw_2_rand = trapz(iXax,MNodeBetw_2_rand(Ind_min:Ind_max,:));
auc_MEdgeBetw_2_rand = trapz(iXax,MEdgeBetw_2_rand(Ind_min:Ind_max,:));
auc_Lambda_2_rand = trapz(iXax,Lambda_2_rand(Ind_min:Ind_max,:));
auc_Gamma_2_rand = trapz(iXax,Gamma_2_rand(Ind_min:Ind_max,:));
auc_Sigma_2_rand = trapz(iXax,Sigma_2_rand(Ind_min:Ind_max,:));

%% test the significance of AUC between groups
p_auc_MClust = CL_Pval(auc_MClust_2_rand - auc_MClust_1_rand, auc_MClust_2 - auc_MClust_1,'AUC_MClust', Tail)
p_auc_Trans = CL_Pval(auc_Trans_2_rand - auc_Trans_1_rand, auc_Trans_2 - auc_Trans_1,'AUC_Trans', Tail)
p_auc_Assort = CL_Pval(auc_Assort_2_rand - auc_Assort_1_rand, auc_Assort_2 - auc_Assort_1,'AUC_Assort', Tail)
p_auc_GEff = CL_Pval(auc_GEff_2_rand - auc_GEff_1_rand, auc_GEff_2 - auc_GEff_1,'AUC_GEff', Tail)
p_auc_MLocEff = CL_Pval(auc_MLocEff_2_rand - auc_MLocEff_1_rand, auc_MLocEff_2 - auc_MLocEff_1,'AUC_MLocEff', Tail)
p_auc_Mod = CL_Pval(auc_Mod_2_rand - auc_Mod_1_rand, auc_Mod_2 - auc_Mod_1,'AUC_Mod', Tail)
p_auc_ModL = CL_Pval(auc_ModL_2_rand - auc_ModL_1_rand, auc_ModL_2 - auc_ModL_1,'AUC_ModL', Tail)
p_auc_PathL = CL_Pval(auc_PathL_2_rand - auc_PathL_1_rand, auc_PathL_2 - auc_PathL_1,'AUC_PathL', Tail)
p_auc_MNodeBetw = CL_Pval(auc_MNodeBetw_2_rand - auc_MNodeBetw_1_rand, auc_MNodeBetw_2 - auc_MNodeBetw_1,'AUC_MNodeBetw', Tail)
p_auc_MEdgeBetw = CL_Pval(auc_MEdgeBetw_2_rand - auc_MEdgeBetw_1_rand, auc_MEdgeBetw_2 - auc_MEdgeBetw_1,'AUC_MEdgeBetw', Tail)
p_auc_Lambda = CL_Pval(auc_Lambda_2_rand - auc_Lambda_1_rand, auc_Lambda_2 - auc_Lambda_1,'AUC_Lambda', Tail)
p_auc_Gamma = CL_Pval(auc_Gamma_2_rand - auc_Gamma_1_rand, auc_Gamma_2 - auc_Gamma_1,'AUC_Gamma', Tail)
p_auc_Sigma = CL_Pval(auc_Sigma_2_rand - auc_Sigma_1_rand, auc_Sigma_2 - auc_Sigma_1,'AUC_Sigma', Tail)


