% This code utilizes functional data analysis to test the significance of 
% the difference in network masures curve between two groups.

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

fda_MClust_1 = MClust_1(Ind_min:Ind_max);
fda_Trans_1 = Trans_1(Ind_min:Ind_max);
fda_Assort_1 = Assort_1(Ind_min:Ind_max);
fda_GEff_1 = GEff_1(Ind_min:Ind_max);
fda_MLocEff_1 = MLocEff_1(Ind_min:Ind_max);
fda_Mod_1 = Mod_1(Ind_min:Ind_max);
fda_ModL_1 = ModL_1(Ind_min:Ind_max);
fda_PathL_1 = PathL_1(Ind_min:Ind_max);
fda_MNodeBetw_1 = MNodeBetw_1(Ind_min:Ind_max);
fda_MEdgeBetw_1 = MEdgeBetw_1(Ind_min:Ind_max);
fda_Lambda_1 = Lambda_1(Ind_min:Ind_max);
fda_Gamma_1 = Gamma_1(Ind_min:Ind_max);
fda_Sigma_1 = Sigma_1(Ind_min:Ind_max);

fda_MClust_2 = MClust_2(Ind_min:Ind_max);
fda_Trans_2 = Trans_2(Ind_min:Ind_max);
fda_Assort_2 = Assort_2(Ind_min:Ind_max);
fda_GEff_2 = GEff_2(Ind_min:Ind_max);
fda_MLocEff_2 = MLocEff_2(Ind_min:Ind_max);
fda_Mod_2 = Mod_2(Ind_min:Ind_max);
fda_ModL_2 = ModL_2(Ind_min:Ind_max);
fda_PathL_2 = PathL_2(Ind_min:Ind_max);
fda_MNodeBetw_2 = MNodeBetw_2(Ind_min:Ind_max);
fda_MEdgeBetw_2 = MEdgeBetw_2(Ind_min:Ind_max);
fda_Lambda_2 = Lambda_2(Ind_min:Ind_max);
fda_Gamma_2 = Gamma_2(Ind_min:Ind_max);
fda_Sigma_2 = Sigma_2(Ind_min:Ind_max);

fda_MClust_1_rand = MClust_1_rand(Ind_min:Ind_max,:);
fda_Trans_1_rand = Trans_1_rand(Ind_min:Ind_max,:);
fda_Assort_1_rand = Assort_1_rand(Ind_min:Ind_max,:);
fda_GEff_1_rand = GEff_1_rand(Ind_min:Ind_max,:);
fda_MLocEff_1_rand = MLocEff_1_rand(Ind_min:Ind_max,:);
fda_Mod_1_rand = Mod_1_rand(Ind_min:Ind_max,:);
fda_ModL_1_rand = ModL_1_rand(Ind_min:Ind_max,:);
fda_PathL_1_rand = PathL_1_rand(Ind_min:Ind_max,:);
fda_MNodeBetw_1_rand = MNodeBetw_1_rand(Ind_min:Ind_max,:);
fda_MEdgeBetw_1_rand = MEdgeBetw_1_rand(Ind_min:Ind_max,:);
fda_Lambda_1_rand = Lambda_1_rand(Ind_min:Ind_max,:);
fda_Gamma_1_rand = Gamma_1_rand(Ind_min:Ind_max,:);
fda_Sigma_1_rand = Sigma_1_rand(Ind_min:Ind_max,:);

fda_MClust_2_rand = MClust_2_rand(Ind_min:Ind_max,:);
fda_Trans_2_rand = Trans_2_rand(Ind_min:Ind_max,:);
fda_Assort_2_rand = Assort_2_rand(Ind_min:Ind_max,:);
fda_GEff_2_rand = GEff_2_rand(Ind_min:Ind_max,:);
fda_MLocEff_2_rand = MLocEff_2_rand(Ind_min:Ind_max,:);
fda_Mod_2_rand = Mod_2_rand(Ind_min:Ind_max,:);
fda_ModL_2_rand = ModL_2_rand(Ind_min:Ind_max,:);
fda_PathL_2_rand = PathL_2_rand(Ind_min:Ind_max,:);
fda_MNodeBetw_2_rand = MNodeBetw_2_rand(Ind_min:Ind_max,:);
fda_MEdgeBetw_2_rand = MEdgeBetw_2_rand(Ind_min:Ind_max,:);
fda_Lambda_2_rand = Lambda_2_rand(Ind_min:Ind_max,:);
fda_Gamma_2_rand = Gamma_2_rand(Ind_min:Ind_max,:);
fda_Sigma_2_rand = Sigma_2_rand(Ind_min:Ind_max,:);

%% 
p_fda_MClust = CL_Pval(sum(fda_MClust_2_rand - fda_MClust_1_rand), sum(fda_MClust_2 - fda_MClust_1),'fda_MClust',Tail)
p_fda_Trans = CL_Pval(sum(fda_Trans_2_rand - fda_Trans_1_rand), sum(fda_Trans_2 - fda_Trans_1),'fda_Trans',Tail)
p_fda_Assort = CL_Pval(sum(fda_Assort_2_rand - fda_Assort_1_rand), sum(fda_Assort_2 - fda_Assort_1),'fda_Assort',Tail)
p_fda_GEff = CL_Pval(sum(fda_GEff_2_rand - fda_GEff_1_rand), sum(fda_GEff_2 - fda_GEff_1),'fda_GEff',Tail)
p_fda_MLocEff = CL_Pval(sum(fda_MLocEff_2_rand - fda_MLocEff_1_rand), sum(fda_MLocEff_2 - fda_MLocEff_1),'fda_MLocEff',Tail)
p_fda_Mod = CL_Pval(sum(fda_Mod_2_rand - fda_Mod_1_rand), sum(fda_Mod_2 - fda_Mod_1),'fda_Mod',Tail)
p_fda_ModL = CL_Pval(sum(fda_ModL_2_rand - fda_ModL_1_rand), sum(fda_ModL_2 - fda_ModL_1),'fda_ModL',Tail)
p_fda_PathL = CL_Pval(sum(fda_PathL_2_rand - fda_PathL_1_rand), sum(fda_PathL_2 - fda_PathL_1),'fda_PathL',Tail)
p_fda_MNodeBetw = CL_Pval(sum(fda_MNodeBetw_2_rand - fda_MNodeBetw_1_rand), sum(fda_MNodeBetw_2 - fda_MNodeBetw_1),'fda_MNodeBetw',Tail)
p_fda_MEdgeBetw = CL_Pval(sum(fda_MEdgeBetw_2_rand - fda_MEdgeBetw_1_rand), sum(fda_MEdgeBetw_2 - fda_MEdgeBetw_1),'fda_MEdgeBetw',Tail)
p_fda_Lambda = CL_Pval(sum(fda_Lambda_2_rand - fda_Lambda_1_rand), sum(fda_Lambda_2 - fda_Lambda_1),'fda_Lambda',Tail)
p_fda_Gamma = CL_Pval(sum(fda_Gamma_2_rand - fda_Gamma_1_rand), sum(fda_Gamma_2 - fda_Gamma_1),'fda_Gamma',Tail)
p_fda_Sigma = CL_Pval(sum(fda_Sigma_2_rand - fda_Sigma_1_rand), sum(fda_Sigma_2 - fda_Sigma_1),'fda_Sigma',Tail)


