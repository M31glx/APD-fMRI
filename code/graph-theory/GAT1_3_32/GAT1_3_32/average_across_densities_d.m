function average_across_densities_d%(NetM,MinMax,Group1,Group2,varargin)
%averagig net measures across densities and then compare the mean between groups (either
%permutation or t-test)

%NetM : network measures for each group and permutations (all the measures written in the NetMeasures_perDensity folder)
%MinMax: mininmum, maximum thresholds and the thresholded step

%varargin: optional input arguments are those (covariates)for which you need to correct the data:
%Standard format:
%'var1': an array with the same number of rows as data1 whose columns
%       contin different covariates 
%'var2': an array with the same number of rows as data1 whose columns
%       contin different covariates

%output
%%
%-- Hadi Hosseini, Mar 21, 2011
%-- updated on Apr 20,2011 for GUI
%-- updated on August 21,2011 for functionals


%%

GUIinput=0;

if nargin<1
    
    data_log = spm_select(1,'mat4GATd*','Select mat4GATd .mat (output from "Network Measures")');

    ff = load(data_log,'mat4GATd');
    Group1 = ff.mat4GATd.g1;Group2 = ff.mat4GATd.g2;
    
    if isfield(ff.mat4GATd,'tail') && isfield(ff.mat4GATd,'alpha') 
        
        Alpha = ff.mat4GATd.alpha;
        Tail = ff.mat4GATd.tail;
        
    else
        
        Alpha = .05;
        Tail = 2;
    
    end
    
    
    NetM = spm_select(1,'NetMesPerDens_D.mat','Select NetMesPerDens_D.mat file ');
    OutputFName='Average_Across_Densities_Results_D.mat';

    fprintf('%-4s\n',['loading inputs...']);
    load(NetM);
    
    MinMesPlot=ff.mat4GATd.MinThr;
    MaxMesPlot=ff.mat4GATd.MaxThr;
    MesStepPlot=ff.mat4GATd.MesStep;
    
    GUIinput=1;
end

%%

if isequal(GUIinput,0)
    if nargin <3

        load(NetMeasures);
        Group1='Group1';
        Group2='Group2';

    elseif nargin <2
    
        'error: please provide one data set for each group)' 
    
    end
    
end


%% 
xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinThr=input('select your final decision on minimum threshold:');MinIdx=find(single(xxx)==single(MinThr));
MaxThr=input('select your final decision on maximum threshold:');MaxIdx=find(single(xxx)==single(MaxThr));

if MaxIdx > size(Dens_1_rand,2)
    
    MaxIdx = size(Dens_1_rand,2);
end

Dens_1_rand=Dens_1_rand(:,MinIdx:MaxIdx,:);
MClust_1_rand=MClust_1_rand(:,MinIdx:MaxIdx,:);
MDeg_1_rand=MDeg_1_rand(:,MinIdx:MaxIdx,:);
Trans_1_rand=Trans_1_rand(:,MinIdx:MaxIdx,:);
Assort_1_rand=Assort_1_rand(:,MinIdx:MaxIdx,:);
GEff_1_rand=GEff_1_rand(:,MinIdx:MaxIdx,:);
MLocEff_1_rand=MLocEff_1_rand(:,MinIdx:MaxIdx,:); 
Mod_1_rand=Mod_1_rand(:,MinIdx:MaxIdx,:); 
ModL_1_rand=ModL_1_rand(:,MinIdx:MaxIdx,:);     
PathL_1_rand=PathL_1_rand(:,MinIdx:MaxIdx,:); 
MNodeBetw_1_rand=MNodeBetw_1_rand(:,MinIdx:MaxIdx,:); 
MEdgeBetw_1_rand=MEdgeBetw_1_rand(:,MinIdx:MaxIdx,:); 
Lambda_1_rand=Lambda_1_rand(:,MinIdx:MaxIdx,:); 
Gamma_1_rand=Gamma_1_rand(:,MinIdx:MaxIdx,:);     
Sigma_1_rand=Sigma_1_rand(:,MinIdx:MaxIdx,:);

Dens_2_rand=Dens_2_rand(:,MinIdx:MaxIdx,:); 
MClust_2_rand=MClust_2_rand(:,MinIdx:MaxIdx,:); 
MDeg_2_rand=MDeg_2_rand(:,MinIdx:MaxIdx,:); 
Trans_2_rand=Trans_2_rand(:,MinIdx:MaxIdx,:);      
Assort_2_rand=Assort_2_rand(:,MinIdx:MaxIdx,:); 
GEff_2_rand=GEff_2_rand(:,MinIdx:MaxIdx,:); 
MLocEff_2_rand=MLocEff_2_rand(:,MinIdx:MaxIdx,:); 
Mod_2_rand=Mod_2_rand(:,MinIdx:MaxIdx,:); 
ModL_2_rand=ModL_2_rand(:,MinIdx:MaxIdx,:); 
PathL_2_rand=PathL_2_rand(:,MinIdx:MaxIdx,:); 
MNodeBetw_2_rand=MNodeBetw_2_rand(:,MinIdx:MaxIdx,:); 
MEdgeBetw_2_rand=MEdgeBetw_2_rand(:,MinIdx:MaxIdx,:); 
Lambda_2_rand=Lambda_2_rand(:,MinIdx:MaxIdx,:); 
Gamma_2_rand=Gamma_2_rand(:,MinIdx:MaxIdx,:);     
Sigma_2_rand=Sigma_2_rand(:,MinIdx:MaxIdx,:); 

Dens_1=Dens_1(:,MinIdx:MaxIdx); 
MClust_1=MClust_1(:,MinIdx:MaxIdx); 
MDeg_1=MDeg_1(:,MinIdx:MaxIdx); 
Trans_1=Trans_1(:,MinIdx:MaxIdx);
Assort_1=Assort_1(:,MinIdx:MaxIdx);
GEff_1=GEff_1(:,MinIdx:MaxIdx);
MLocEff_1=MLocEff_1(:,MinIdx:MaxIdx);
Mod_1=Mod_1(:,MinIdx:MaxIdx);
ModL_1=ModL_1(:,MinIdx:MaxIdx);
PathL_1=PathL_1(:,MinIdx:MaxIdx);
MNodeBetw_1=MNodeBetw_1(:,MinIdx:MaxIdx);
MEdgeBetw_1=MEdgeBetw_1(:,MinIdx:MaxIdx);
Lambda_1=Lambda_1(:,MinIdx:MaxIdx);
Gamma_1=Gamma_1(:,MinIdx:MaxIdx);
Sigma_1=Sigma_1(:,MinIdx:MaxIdx);

Dens_2=Dens_2(:,MinIdx:MaxIdx);
MClust_2=MClust_2(:,MinIdx:MaxIdx);
MDeg_2=MDeg_2(:,MinIdx:MaxIdx);
Trans_2=Trans_2(:,MinIdx:MaxIdx);
Assort_2=Assort_2(:,MinIdx:MaxIdx);
GEff_2=GEff_2(:,MinIdx:MaxIdx);
MLocEff_2=MLocEff_2(:,MinIdx:MaxIdx);
Mod_2=Mod_2(:,MinIdx:MaxIdx);
ModL_2=ModL_2(:,MinIdx:MaxIdx);
PathL_2=PathL_2(:,MinIdx:MaxIdx);
MNodeBetw_2=MNodeBetw_2(:,MinIdx:MaxIdx);
MEdgeBetw_2=MEdgeBetw_2(:,MinIdx:MaxIdx);
Lambda_2=Lambda_2(:,MinIdx:MaxIdx);
Gamma_2=Gamma_2(:,MinIdx:MaxIdx);
Sigma_2=Sigma_2(:,MinIdx:MaxIdx);


fprintf('%-4s\n',['saving net measures across new threshold range...']);
dd=pwd;     
mkdir('AverageAcrossDens_Results');
cd([dd '/AverageAcrossDens_Results']);

save NetMesPerDens_AveAcrossDens_D Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
     Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
     PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
     Sigma_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
     Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
     PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
     Sigma_2_rand Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
     Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
     Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ...
     ModL_2 PathL_2 MNodeBetw_2 MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 


MeanAcrossDens_MClust_1_rand(1,:)=mean(mean(MClust_1_rand,2));
MeanAcrossDens_MDeg_1_rand(1,:)=mean(mean(MDeg_1_rand,2));
MeanAcrossDens_Trans_1_rand(1,:)=mean(mean(Trans_1_rand,2));
MeanAcrossDens_Assort_1_rand(1,:)=mean(mean(Assort_1_rand,2));
MeanAcrossDens_GEff_1_rand(1,:)=mean(mean(GEff_1_rand,2));
MeanAcrossDens_MLocEff_1_rand(1,:)=mean(mean(MLocEff_1_rand,2));
MeanAcrossDens_Mod_1_rand(1,:)=mean(mean(Mod_1_rand,2));
MeanAcrossDens_ModL_1_rand(1,:)=mean(mean(ModL_1_rand,2));
MeanAcrossDens_PathL_1_rand(1,:)=mean(mean(PathL_1_rand,2));
MeanAcrossDens_MNodeBetw_1_rand(1,:)=mean(mean(MNodeBetw_1_rand,2));
MeanAcrossDens_MEdgeBetw_1_rand(1,:)=mean(mean(MEdgeBetw_1_rand,2));
MeanAcrossDens_Lambda_1_rand(1,:)=mean(mean(Lambda_1_rand,2));
MeanAcrossDens_Gamma_1_rand(1,:)=mean(mean(Gamma_1_rand,2));
MeanAcrossDens_Sigma_1_rand(1,:)=mean(mean(Sigma_1_rand,2));


MeanAcrossDens_MClust_2_rand(1,:)=mean(mean(MClust_2_rand,2));
MeanAcrossDens_MDeg_2_rand(1,:)=mean(mean(MDeg_2_rand,2));
MeanAcrossDens_Trans_2_rand(1,:)=mean(mean(Trans_2_rand,2));
MeanAcrossDens_Assort_2_rand(1,:)=mean(mean(Assort_2_rand,2));
MeanAcrossDens_GEff_2_rand(1,:)=mean(mean(GEff_2_rand,2));
MeanAcrossDens_MLocEff_2_rand(1,:)=mean(mean(MLocEff_2_rand,2));
MeanAcrossDens_Mod_2_rand(1,:)=mean(mean(Mod_2_rand,2));
MeanAcrossDens_ModL_2_rand(1,:)=mean(mean(ModL_2_rand,2));
MeanAcrossDens_PathL_2_rand(1,:)=mean(mean(PathL_2_rand,2));
MeanAcrossDens_MNodeBetw_2_rand(1,:)=mean(mean(MNodeBetw_2_rand,2));
MeanAcrossDens_MEdgeBetw_2_rand(1,:)=mean(mean(MEdgeBetw_2_rand,2));
MeanAcrossDens_Lambda_2_rand(1,:)=mean(mean(Lambda_2_rand,2));
MeanAcrossDens_Gamma_2_rand(1,:)=mean(mean(Gamma_2_rand,2));
MeanAcrossDens_Sigma_2_rand(1,:)=mean(mean(Sigma_2_rand,2));


MeanAcrossDens_MClust_1=mean(mean(MClust_1,2));
MeanAcrossDens_MDeg_1=mean(mean(MDeg_1,2));
MeanAcrossDens_Trans_1=mean(mean(Trans_1,2));
MeanAcrossDens_Assort_1=mean(mean(Assort_1,2));
MeanAcrossDens_GEff_1=mean(mean(GEff_1,2));
MeanAcrossDens_MLocEff_1=mean(mean(MLocEff_1,2));
MeanAcrossDens_Mod_1=mean(mean(Mod_1,2));
MeanAcrossDens_ModL_1=mean(mean(ModL_1,2));
MeanAcrossDens_PathL_1=mean(mean(PathL_1,2));
MeanAcrossDens_MNodeBetw_1=mean(mean(MNodeBetw_1,2));
MeanAcrossDens_MEdgeBetw_1=mean(mean(MEdgeBetw_1,2));
MeanAcrossDens_Lambda_1=mean(mean(Lambda_1,2));
MeanAcrossDens_Gamma_1=mean(mean(Gamma_1,2));
MeanAcrossDens_Sigma_1=mean(mean(Sigma_1,2));


MeanAcrossDens_MClust_2=mean(mean(MClust_2,2));
MeanAcrossDens_MDeg_2=mean(mean(MDeg_2,2));
MeanAcrossDens_Trans_2=mean(mean(Trans_2,2));
MeanAcrossDens_Assort_2=mean(mean(Assort_2,2));
MeanAcrossDens_GEff_2=mean(mean(GEff_2,2));
MeanAcrossDens_MLocEff_2=mean(mean(MLocEff_2,2));
MeanAcrossDens_Mod_2=mean(mean(Mod_2,2));
MeanAcrossDens_ModL_2=mean(mean(ModL_2,2));
MeanAcrossDens_PathL_2=mean(mean(PathL_2,2));
MeanAcrossDens_MNodeBetw_2=mean(mean(MNodeBetw_2,2));
MeanAcrossDens_MEdgeBetw_2=mean(mean(MEdgeBetw_2,2));
MeanAcrossDens_Lambda_2=mean(mean(Lambda_2,2));
MeanAcrossDens_Gamma_2=mean(mean(Gamma_2,2));
MeanAcrossDens_Sigma_2=mean(mean(Sigma_2,2));


%% 

[mu_MClust1_rand,sigma_MClust1_rand,muCi_MClust1_rand,sigmaCi_MClust1_rand]=normfit(MeanAcrossDens_MClust_1_rand',Alpha);
[mu_MDeg1_rand,sigma_MDeg1_rand,muCi_MDeg1_rand,sigmaCi_MDeg1_rand]=normfit(MeanAcrossDens_MDeg_1_rand',Alpha);
[mu_Trans1_rand,sigma_Trans1_rand,muCi_Trans1_rand,sigmaCi_Trans1_rand]=normfit(MeanAcrossDens_Trans_1_rand',Alpha);
[mu_Assort1_rand,sigma_Assort1_rand,muCi_Assort1_rand,sigmaCi_Assort1_rand]=normfit(MeanAcrossDens_Assort_1_rand',Alpha);
[mu_GEff1_rand,sigma_GEff1_rand,muCi_GEff1_rand,sigmaCi_GEff1_rand]=normfit(MeanAcrossDens_GEff_1_rand',Alpha);
[mu_MLocEff1_rand,sigma_MLocEff1_rand,muCi_MLocEff1_rand,sigmaCi_MLocEff1_rand]=normfit(MeanAcrossDens_MLocEff_1_rand',Alpha);
[mu_Mod1_rand,sigma_Mod1_rand,muCi_Mod1_rand,sigmaCi_Mod1_rand]=normfit(MeanAcrossDens_Mod_1_rand',Alpha);
[mu_ModL1_rand,sigma_ModL1_rand,muCi_ModL1_rand,sigmaCi_ModL1_rand]=normfit(MeanAcrossDens_ModL_1_rand',Alpha);
[mu_PathL1_rand,sigma_PathL1_rand,muCi_PathL1_rand,sigmaCi_PathL1_rand]=normfit(MeanAcrossDens_PathL_1_rand',Alpha);
[mu_MNodeBetw1_rand,sigma_MNodeBetw1_rand,muCi_MNodeBetw1_rand,sigmaCi_MNodeBetw1_rand]=normfit(MeanAcrossDens_MNodeBetw_1_rand',Alpha);
[mu_MEdgeBetw1_rand,sigma_MEdgeBetw1_rand,muCi_MEdgeBetw1_rand,sigmaCi_MEdgeBetw1_rand]=normfit(MeanAcrossDens_MEdgeBetw_1_rand',Alpha);
[mu_Lambda1_rand,sigma_Lambda1_rand,muCi_Lambda1_rand,sigmaCi_Lambda1_rand]=normfit(MeanAcrossDens_Lambda_1_rand',Alpha);
[mu_Gamma1_rand,sigma_Gamma1_rand,muCi_Gamma1_rand,sigmaCi_Gamma1_rand]=normfit(MeanAcrossDens_Gamma_1_rand',Alpha);
[mu_Sigma1_rand,sigma_Sigma1_rand,muCi_Sigma1_rand,sigmaCi_Sigma1_rand]=normfit(MeanAcrossDens_Sigma_1_rand',Alpha);

[mu_MClust2_rand,sigma_MClust2_rand,muCi_MClust2_rand,sigmaCi_MClust2_rand]=normfit(MeanAcrossDens_MClust_2_rand',Alpha);
 [mu_MDeg2_rand,sigma_MDeg2_rand,muCi_MDeg2_rand,sigmaCi_MDeg2_rand]=normfit(MeanAcrossDens_MDeg_2_rand',Alpha);
[mu_Trans2_rand,sigma_Trans2_rand,muCi_Trans2_rand,sigmaCi_Trans2_rand]=normfit(MeanAcrossDens_Trans_2_rand',Alpha);
[mu_Assort2_rand,sigma_Assort2_rand,muCi_Assort2_rand,sigmaCi_Assort2_rand]=normfit(MeanAcrossDens_Assort_2_rand',Alpha);
[mu_GEff2_rand,sigma_GEff2_rand,muCi_GEff2_rand,sigmaCi_GEff2_rand]=normfit(MeanAcrossDens_GEff_2_rand',Alpha);
[mu_MLocEff2_rand,sigma_MLocEff2_rand,muCi_MLocEff2_rand,sigmaCi_MLocEff2_rand]=normfit(MeanAcrossDens_MLocEff_2_rand',Alpha);
[mu_Mod2_rand,sigma_Mod2_rand,muCi_Mod2_rand,sigmaCi_Mod2_rand]=normfit(MeanAcrossDens_Mod_2_rand',Alpha);
[mu_ModL2_rand,sigma_ModL2_rand,muCi_ModL2_rand,sigmaCi_ModL2_rand]=normfit(MeanAcrossDens_ModL_2_rand',Alpha);
[mu_PathL2_rand,sigma_PathL2_rand,muCi_PathL2_rand,sigmaCi_PathL2_rand]=normfit(MeanAcrossDens_PathL_2_rand',Alpha);
[mu_MNodeBetw2_rand,sigma_MNodeBetw2_rand,muCi_MNodeBetw2_rand,sigmaCi_MNodeBetw2_rand]=normfit(MeanAcrossDens_MNodeBetw_2_rand',Alpha);
[mu_MEdgeBetw2_rand,sigma_MEdgeBetw2_rand,muCi_MEdgeBetw2_rand,sigmaCi_MEdgeBetw2_rand]=normfit(MeanAcrossDens_MEdgeBetw_2_rand',Alpha);
[mu_Lambda2_rand,sigma_Lambda2_rand,muCi_Lambda2_rand,sigmaCi_Lambda2_rand]=normfit(MeanAcrossDens_Lambda_2_rand',Alpha);
[mu_Gamma2_rand,sigma_Gamma2_rand,muCi_Gamma2_rand,sigmaCi_Gamma2_rand]=normfit(MeanAcrossDens_Gamma_2_rand',Alpha);
[mu_Sigma2_rand,sigma_Sigma2_rand,muCi_Sigma2_rand,sigmaCi_Sigma2_rand]=normfit(MeanAcrossDens_Sigma_2_rand',Alpha);


%%

N_rand=size(MeanAcrossDens_MClust_1_rand,2);

Pvalue = Alpha;

Ci_MClust=CL_per(MeanAcrossDens_MClust_2_rand'-MeanAcrossDens_MClust_1_rand',Pvalue);
Ci_MDeg=CL_per(MeanAcrossDens_MDeg_2_rand'-MeanAcrossDens_MDeg_1_rand',Pvalue);
Ci_Trans=CL_per(MeanAcrossDens_Trans_2_rand'-MeanAcrossDens_Trans_1_rand',Pvalue);
Ci_Assort=CL_per(MeanAcrossDens_Assort_2_rand'-MeanAcrossDens_Assort_1_rand',Pvalue);
Ci_GEff=CL_per(MeanAcrossDens_GEff_2_rand'-MeanAcrossDens_GEff_1_rand',Pvalue);
Ci_MLocEff=CL_per(MeanAcrossDens_MLocEff_2_rand'-MeanAcrossDens_MLocEff_1_rand',Pvalue);
Ci_Mod=CL_per(MeanAcrossDens_Mod_2_rand'-MeanAcrossDens_Mod_1_rand',Pvalue);
Ci_ModL=CL_per(MeanAcrossDens_ModL_2_rand'-MeanAcrossDens_ModL_1_rand',Pvalue);
Ci_PathL=CL_per(MeanAcrossDens_PathL_2_rand'-MeanAcrossDens_PathL_1_rand',Pvalue);
Ci_MNodeBetw=CL_per(MeanAcrossDens_MNodeBetw_2_rand'-MeanAcrossDens_MNodeBetw_1_rand',Pvalue);
Ci_MEdgeBetw=CL_per(MeanAcrossDens_MEdgeBetw_2_rand'-MeanAcrossDens_MEdgeBetw_1_rand',Pvalue);
Ci_Lambda=CL_per(MeanAcrossDens_Lambda_2_rand'-MeanAcrossDens_Lambda_1_rand',Pvalue);
Ci_Gamma=CL_per(MeanAcrossDens_Gamma_2_rand'-MeanAcrossDens_Gamma_1_rand',Pvalue);
Ci_Sigma=CL_per(MeanAcrossDens_Sigma_2_rand'-MeanAcrossDens_Sigma_1_rand',Pvalue);



%% 
   
p_Assort=CL_Pval(MeanAcrossDens_Assort_2_rand-MeanAcrossDens_Assort_1_rand,MeanAcrossDens_Assort_2-MeanAcrossDens_Assort_1,'Assort',Tail);
p_GEff=CL_Pval(MeanAcrossDens_GEff_2_rand-MeanAcrossDens_GEff_1_rand,MeanAcrossDens_GEff_2-MeanAcrossDens_GEff_1,'GEff',Tail);
p_Gamma=CL_Pval(MeanAcrossDens_Gamma_2_rand-MeanAcrossDens_Gamma_1_rand,MeanAcrossDens_Gamma_2-MeanAcrossDens_Gamma_1,'Gamma',Tail);
p_Lambda=CL_Pval(MeanAcrossDens_Lambda_2_rand-MeanAcrossDens_Lambda_1_rand,MeanAcrossDens_Lambda_2-MeanAcrossDens_Lambda_1,'Lambda',Tail);
p_MClust=CL_Pval(MeanAcrossDens_MClust_2_rand-MeanAcrossDens_MClust_1_rand,MeanAcrossDens_MClust_2-MeanAcrossDens_MClust_1,'MClust',Tail);
p_MEdgeBetw=CL_Pval(MeanAcrossDens_MEdgeBetw_2_rand-MeanAcrossDens_MEdgeBetw_1_rand,MeanAcrossDens_MEdgeBetw_2-MeanAcrossDens_MEdgeBetw_1,'MEdgeBetw',Tail);
p_MLocEff=CL_Pval(MeanAcrossDens_MLocEff_2_rand-MeanAcrossDens_MLocEff_1_rand,MeanAcrossDens_MLocEff_2-MeanAcrossDens_MLocEff_1,'MLocEff',Tail);
p_MNodeBetw=CL_Pval(MeanAcrossDens_MNodeBetw_2_rand-MeanAcrossDens_MNodeBetw_1_rand,MeanAcrossDens_MNodeBetw_2-MeanAcrossDens_MNodeBetw_1,'MNodeBetw',Tail);
p_ModL=CL_Pval(MeanAcrossDens_ModL_2_rand-MeanAcrossDens_ModL_1_rand,MeanAcrossDens_ModL_2-MeanAcrossDens_ModL_1,'ModL',Tail);
p_Mod=CL_Pval(MeanAcrossDens_Mod_2_rand-MeanAcrossDens_Mod_1_rand,MeanAcrossDens_Mod_2-MeanAcrossDens_Mod_1,'Mod',Tail);
p_PathL=CL_Pval(MeanAcrossDens_PathL_2_rand-MeanAcrossDens_PathL_1_rand,MeanAcrossDens_PathL_2-MeanAcrossDens_PathL_1,'PathL',Tail);
p_Sigma=CL_Pval(MeanAcrossDens_Sigma_2_rand-MeanAcrossDens_Sigma_1_rand,MeanAcrossDens_Sigma_2-MeanAcrossDens_Sigma_1,'Sigma',Tail);
p_Trans=CL_Pval(MeanAcrossDens_Trans_2_rand-MeanAcrossDens_Trans_1_rand,MeanAcrossDens_Trans_2-MeanAcrossDens_Trans_1,'Trans',Tail);
p_MDeg=CL_Pval(MeanAcrossDens_MDeg_2_rand-MeanAcrossDens_MDeg_1_rand,MeanAcrossDens_MDeg_2-MeanAcrossDens_MDeg_1,'MDeg',Tail);


%% 

MeanAcrossDens_MClust_1_ttest=mean(MClust_1,2);
MeanAcrossDens_MDeg_1_ttest=mean(MDeg_1,2);
MeanAcrossDens_Trans_1_ttest=mean(Trans_1,2);
MeanAcrossDens_Assort_1_ttest=mean(Assort_1,2);
MeanAcrossDens_GEff_1_ttest=mean(GEff_1,2);
MeanAcrossDens_MLocEff_1_ttest=mean(MLocEff_1,2);
MeanAcrossDens_Mod_1_ttest=mean(Mod_1,2);
MeanAcrossDens_ModL_1_ttest=mean(ModL_1,2);
MeanAcrossDens_PathL_1_ttest=mean(PathL_1,2);
MeanAcrossDens_MNodeBetw_1_ttest=mean(MNodeBetw_1,2);
MeanAcrossDens_MEdgeBetw_1_ttest=mean(MEdgeBetw_1,2);
MeanAcrossDens_Lambda_1_ttest=mean(Lambda_1,2);
MeanAcrossDens_Gamma_1_ttest=mean(Gamma_1,2);
MeanAcrossDens_Sigma_1_ttest=mean(Sigma_1,2);

MeanAcrossDens_MClust_2_ttest=mean(MClust_2,2);
MeanAcrossDens_MDeg_2_ttest=mean(MDeg_2,2);
%MeanAcrossDens_MStr_2_ttest=mean(MStr_2,2);
MeanAcrossDens_Trans_2_ttest=mean(Trans_2,2);
MeanAcrossDens_Assort_2_ttest=mean(Assort_2,2);
MeanAcrossDens_GEff_2_ttest=mean(GEff_2,2);
MeanAcrossDens_MLocEff_2_ttest=mean(MLocEff_2,2);
MeanAcrossDens_Mod_2_ttest=mean(Mod_2,2);
MeanAcrossDens_ModL_2_ttest=mean(ModL_2,2);
MeanAcrossDens_PathL_2_ttest=mean(PathL_2,2);
MeanAcrossDens_MNodeBetw_2_ttest=mean(MNodeBetw_2,2);
MeanAcrossDens_MEdgeBetw_2_ttest=mean(MEdgeBetw_2,2);
MeanAcrossDens_Lambda_2_ttest=mean(Lambda_2,2);
MeanAcrossDens_Gamma_2_ttest=mean(Gamma_2,2);
MeanAcrossDens_Sigma_2_ttest=mean(Sigma_2,2);


save MeanAcrossDens_NetMeasVals MeanAcrossDens_MClust_1_ttest MeanAcrossDens_MDeg_1_ttest...
MeanAcrossDens_Trans_1_ttest MeanAcrossDens_Assort_1_ttest ...
MeanAcrossDens_GEff_1_ttest MeanAcrossDens_MLocEff_1_ttest...
MeanAcrossDens_Mod_1_ttest MeanAcrossDens_ModL_1_ttest ...
MeanAcrossDens_PathL_1_ttest MeanAcrossDens_MNodeBetw_1_ttest...
MeanAcrossDens_MEdgeBetw_1_ttest MeanAcrossDens_Lambda_1_ttest...
MeanAcrossDens_Gamma_1_ttest MeanAcrossDens_Sigma_1_ttest...
MeanAcrossDens_MClust_2_ttest MeanAcrossDens_MDeg_2_ttest...
MeanAcrossDens_Trans_2_ttest MeanAcrossDens_Assort_2_ttest...
MeanAcrossDens_GEff_2_ttest MeanAcrossDens_MLocEff_2_ttest...
MeanAcrossDens_Mod_2_ttest MeanAcrossDens_ModL_2_ttest...
MeanAcrossDens_PathL_2_ttest MeanAcrossDens_MNodeBetw_2_ttest...
MeanAcrossDens_MEdgeBetw_2_ttest MeanAcrossDens_Lambda_2_ttest...
MeanAcrossDens_Gamma_2_ttest MeanAcrossDens_Sigma_2_ttest ...

p = Alpha;

if isequal(Tail,1)
    
    Model = 'left';

elseif isequal(Tail,2)
    
    Model = 'both';

end

[h_Trans,pTTest_Trans,ci_Trans,stats_Trans]=ttest2(MeanAcrossDens_Trans_1_ttest,MeanAcrossDens_Trans_2_ttest,p,Model);
[h_Assort,pTTest_Assort,ci_Assort,stats_Assort]=ttest2(MeanAcrossDens_Assort_1_ttest,MeanAcrossDens_Assort_2_ttest,p,Model);
[h_LocEff,pTTest_LocEff,ci_LocEff,stats_LocEff]=ttest2(MeanAcrossDens_MLocEff_1_ttest,MeanAcrossDens_MLocEff_2_ttest,p,Model);
[h_ModL,pTTest_ModL,ci_ModL,stats_ModL]=ttest2(MeanAcrossDens_ModL_1_ttest,MeanAcrossDens_ModL_2_ttest,p,Model);
[h_MNodeBetw,pTTest_MNodeBetw,ci_MNodeBetw,stats_MNodeBetw]=ttest2(MeanAcrossDens_MNodeBetw_1_ttest,MeanAcrossDens_MNodeBetw_2_ttest,p,Model);
[h_MEdgeBetw,pTTest_MEdgeBetw,ci_MEdgeBetw,stats_MEdgeBetw]=ttest2(MeanAcrossDens_MEdgeBetw_1_ttest,MeanAcrossDens_MEdgeBetw_2_ttest,p,Model);
[h_Lambda,pTTest_Lambda,ci_Lambda,stats_Lambda]=ttest2(MeanAcrossDens_Lambda_1_ttest,MeanAcrossDens_Lambda_2_ttest,p,Model);
[h_Gamma,pTTest_Gamma,ci_Gamma,stats_Gamma]=ttest2(MeanAcrossDens_Gamma_1_ttest,MeanAcrossDens_Gamma_2_ttest,p,Model);
[h_Sigma,pTTest_Sigma,ci_Sigma,stats_Sigma]=ttest2(MeanAcrossDens_Sigma_1_ttest,MeanAcrossDens_Sigma_2_ttest,p,Model);


save(['pVal_TTest'],'pTTest_Trans','pTTest_Assort','pTTest_LocEff','pTTest_ModL','pTTest_MNodeBetw', ...
     'pTTest_MEdgeBetw','pTTest_Lambda','pTTest_Gamma','pTTest_Sigma');
