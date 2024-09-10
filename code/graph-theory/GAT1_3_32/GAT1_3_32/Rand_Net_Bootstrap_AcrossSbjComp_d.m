function [N1,N2,N_G1,N_G2]=Rand_Net_Bootstrap_AcrossSbjComp_d(data1,data2,MinMesPlot,MesStepPlot,MaxMesPlot,Nperm,Group1,Group2,varargin)

%data1/2 : group 1/2 functional network measures (output NetMesBin_* from NetMesvsDensity..._Functionals)
%MinMax: mninmum, maximum thresholds and the thresholded step
%Nperm: number of permutation (sampling) (default 100)

%varargin: optional input arguments are those (covariates)for which you need to correct the data:
%Standard format:
%'var1': an array with the same number of rows as data1 whose columns
%       contin different covariates 
%'var2': an array with the same number of rows as data1 whose columns
%       contin different covariates

%output
%'N1': Original unthresholded correlation matrix for group1
%'N2': Original unthresholded correlation matrix for group2
%'N_G1': Cell array containing Randomly generated unthresholded correlation matrix for group1
%'N_G2': Cell array containing Randomly generated unthresholded correlation matrix for group2

%%

%-- Hadi Hosseini, Mar 21, 2011
%-- updated on Apr 20,2011 for GUI
%-- updated on August 18,2011 for functionals


%% 

GUIinput=0;

if nargin<1

    Nperm=input('please type the number of random networks to be generated (e.g. 1000) and press enter  ');
    
    if isempty(Nperm)
    
        Nperm=100;
    
    end
    
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
    
    OutputFName='RandNetMesBootstrap_Results_D.mat';
    
    cov=ff.mat4GATd.flagCov;
    
    if isequal(cov,1)
    
        Data1 = spm_select(1,['NetMesBin_D_Adjusted_' Group1 '_FinalThrRange.mat'],'Select 1st group residuals (NetMesBin_D_adjusted.mat .mat file)');
        Data2 = spm_select(1,['NetMesBin_D_Adjusted_' Group2 '_FinalThrRange.mat'],'Select 2nd Group residuals (NetMesBin_D_adjusted.mat .mat file)');
        
    elseif isequal(cov,0)
        
        Data1 = spm_select(1,['NetMesBin_D_' Group1 '_FinalThrRange.mat'],'Select 1st group NetMesBin_D.mat file ');
        Data2 = spm_select(1,['NetMesBin_D_' Group2 '_FinalThrRange.mat'],'Select 2nd group NetMesBin_D.mat file ');
        
    end
    
    fprintf('%-4s\n',['loading inputs...']);
    f1=load(Data1,'NetMes_Bin');data1=f1.NetMes_Bin;
    f2=load(Data2,'NetMes_Bin');data2=f2.NetMes_Bin;
    
    MinMesPlot=ff.mat4GATd.MinThr;
    MaxMesPlot=ff.mat4GATd.MaxThr;
    MesStepPlot=ff.mat4GATd.MesStep;
    
    GUIinput=1;
    
end


%%

if isequal(GUIinput,0)
    
    if nargin <6
    
        Nperm=100;
        Group1='Group1';
        Group2='Group2';
        
    elseif nargin <5

        'error: please provide the mat files including net measures and the density interval of interest' 
        
    end

end
 
N1=data1;N2=data2;

%% 

Sz1=size(data1,1);Sz2=size(data2,1);

data=data1;data(Sz1+1:Sz1+Sz2,:)=data2;
RandIndex=randperm(Sz1+Sz2);
Randata(1:Sz1+Sz2,:)=data(RandIndex(1:Sz1+Sz2),:);

N_G1=cell(1,Nperm);N_G2=cell(1,Nperm);

for i=1:Nperm
    
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    Samp1=randsample(Sz1+Sz2,Sz1,'true');
    Samp2=randsample(Sz1+Sz2,Sz2,'true');
    N_G1{i}=Randata(Samp1,:);
    N_G2{i}=Randata(Samp2,:);

end


%% 

fprintf('%-4s\n',' extracting network measures...');

Dens_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));%[];
MClust_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MDeg_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MStr_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Trans_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Assort_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
GEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MLocEff_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Mod_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));ModL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
PathL_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MNodeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));MEdgeBetw_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));
Lambda_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Gamma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));Sigma_1_rand=zeros(size(N1,1),size(N1,2),size(N_G1,2));

Dens_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
MClust_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MDeg_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MStr_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Trans_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Assort_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
GEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MLocEff_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Mod_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));ModL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
PathL_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MNodeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));MEdgeBetw_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));
Lambda_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Gamma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));Sigma_2_rand=zeros(size(N2,1),size(N2,2),size(N_G2,2));


Dens_1=zeros(size(N1));
MClust_1=zeros(size(N1));MDeg_1=zeros(size(N1));MStr_1=zeros(size(N1));Trans_1=zeros(size(N1));Assort_1=zeros(size(N1));
GEff_1=zeros(size(N1));MLocEff_1=zeros(size(N1));Mod_1=zeros(size(N1));ModL_1=zeros(size(N1));
PathL_1=zeros(size(N1));MNodeBetw_1=zeros(size(N1));MEdgeBetw_1=zeros(size(N1));
Lambda_1=zeros(size(N1));Gamma_1=zeros(size(N1));Sigma_1=zeros(size(N1));

Dens_2=zeros(size(N2));
MClust_2=zeros(size(N2));MDeg_2=zeros(size(N2));MStr_2=zeros(size(N2));Trans_2=zeros(size(N2));Assort_2=zeros(size(N2));
GEff_2=zeros(size(N2));MLocEff_2=zeros(size(N2));Mod_2=zeros(size(N2));ModL_2=zeros(size(N2));
PathL_2=zeros(size(N2));MNodeBetw_2=zeros(size(N2));MEdgeBetw_2=zeros(size(N2));
Lambda_2=zeros(size(N2));Gamma_2=zeros(size(N2));Sigma_2=zeros(size(N2));


for i=1:size(N_G1,2)
    
    for j=1:size(N1,1)
    
        for k=1:size(N1,2)
        
            Dens_1_rand(j,k,i)=N_G1{i}{j,k}{6};
            MClust_1_rand(j,k,i)=N_G1{i}{j,k}{8};
            MDeg_1_rand(j,k,i)=N_G1{i}{j,k}{2};
            Trans_1_rand(j,k,i)=N_G1{i}{j,k}{9};
            Assort_1_rand(j,k,i)=N_G1{i}{j,k}{5};
            GEff_1_rand(j,k,i)=N_G1{i}{j,k}{10};
            MLocEff_1_rand(j,k,i)=N_G1{i}{j,k}{12};
            Mod_1_rand(j,k,i)=N_G1{i}{j,k}{13};
            ModL_1_rand(j,k,i)=N_G1{i}{j,k}{14};
            PathL_1_rand(j,k,i)=N_G1{i}{j,k}{15};
            MNodeBetw_1_rand(j,k,i)=N_G1{i}{j,k}{17};
            MEdgeBetw_1_rand(j,k,i)=N_G1{i}{j,k}{19};
            Lambda_1_rand(j,k,i)=N_G1{i}{j,k}{20};
            Gamma_1_rand(j,k,i)=N_G1{i}{j,k}{21};
            Sigma_1_rand(j,k,i)=N_G1{i}{j,k}{22};
            
        end
        
    end
    
    for j=1:size(N2,1)
    
        for k=1:size(N2,2)
            
            Dens_2_rand(j,k,i)=N_G2{i}{j,k}{6};
            MClust_2_rand(j,k,i)=N_G2{i}{j,k}{8};
            MDeg_2_rand(j,k,i)=N_G2{i}{j,k}{2};
            Trans_2_rand(j,k,i)=N_G2{i}{j,k}{9};
            Assort_2_rand(j,k,i)=N_G2{i}{j,k}{5};
            GEff_2_rand(j,k,i)=N_G2{i}{j,k}{10};
            MLocEff_2_rand(j,k,i)=N_G2{i}{j,k}{12};
            Mod_2_rand(j,k,i)=N_G2{i}{j,k}{13};
            ModL_2_rand(j,k,i)=N_G2{i}{j,k}{14};
            PathL_2_rand(j,k,i)=N_G2{i}{j,k}{15};
            MNodeBetw_2_rand(j,k,i)=N_G2{i}{j,k}{17};
            MEdgeBetw_2_rand(j,k,i)=N_G2{i}{j,k}{19};
            Lambda_2_rand(j,k,i)=N_G2{i}{j,k}{20};
            Gamma_2_rand(j,k,i)=N_G2{i}{j,k}{21};
            Sigma_2_rand(j,k,i)=N_G2{i}{j,k}{22};       
            
        end
        
    end
    
end

for j=1:size(N1,1)
    
    for k=1:size(N1,2)
        
        Dens_1(j,k)=N1{j,k}{6};
        MClust_1(j,k)=N1{j,k}{8};
        MDeg_1(j,k)=N1{j,k}{2};
        Trans_1(j,k)=N1{j,k}{9};
        Assort_1(j,k)=N1{j,k}{5};
        GEff_1(j,k)=N1{j,k}{10};
        MLocEff_1(j,k)=N1{j,k}{12};
        Mod_1(j,k)=N1{j,k}{13};
        ModL_1(j,k)=N1{j,k}{14};
        PathL_1(j,k)=N1{j,k}{15};
        MNodeBetw_1(j,k)=N1{j,k}{17};
        MEdgeBetw_1(j,k)=N1{j,k}{19};
        Lambda_1(j,k)=N1{j,k}{20};
        Gamma_1(j,k)=N1{j,k}{21};
        Sigma_1(j,k)=N1{j,k}{22};
        
    end
    
end

for j=1:size(N2,1)

    for k=1:size(N2,2)
        
        Dens_2(j,k)=N2{j,k}{6};
        MClust_2(j,k)=N2{j,k}{8};
        MDeg_2(j,k)=N2{j,k}{2};
        Trans_2(j,k)=N2{j,k}{9};
        Assort_2(j,k)=N2{j,k}{5};
        GEff_2(j,k)=N2{j,k}{10};
        MLocEff_2(j,k)=N2{j,k}{12};
        Mod_2(j,k)=N2{j,k}{13};
        ModL_2(j,k)=N2{j,k}{14};
        PathL_2(j,k)=N2{j,k}{15};
        MNodeBetw_2(j,k)=N2{j,k}{17};
        MEdgeBetw_2(j,k)=N2{j,k}{19};
        Lambda_2(j,k)=N2{j,k}{20};
        Gamma_2(j,k)=N2{j,k}{21};
        Sigma_2(j,k)=N2{j,k}{22};
        
    end
    
end

fprintf('%-4s\n',' saving output...');


save NetMesPerDens_D Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
     Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
     PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
     Sigma_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand Assort_2_rand GEff_2_rand ...
     MLocEff_2_rand Mod_2_rand ModL_2_rand PathL_2_rand MNodeBetw_2_rand ...
     MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand Sigma_2_rand ...
     Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
     Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
     Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ...
     ModL_2 PathL_2 MNodeBetw_2 MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2


%%

fprintf('%-4s\n',' calculating confidence intervals ...');

%% Method1: (Bruno, Hosseini, Kesler, 2012 Neurobiology of Disease)

MeanAcrossSbj_MClust_1_rand(:,:)=mean(MClust_1_rand,1);
MeanAcrossSbj_MDeg_1_rand(:,:)=mean(MDeg_1_rand,1);
MeanAcrossSbj_Trans_1_rand(:,:)=mean(Trans_1_rand,1);
MeanAcrossSbj_Assort_1_rand(:,:)=mean(Assort_1_rand,1);
MeanAcrossSbj_GEff_1_rand(:,:)=mean(GEff_1_rand,1);
MeanAcrossSbj_MLocEff_1_rand(:,:)=mean(MLocEff_1_rand,1);
MeanAcrossSbj_Mod_1_rand(:,:)=mean(Mod_1_rand,1);
MeanAcrossSbj_ModL_1_rand(:,:)=mean(ModL_1_rand,1);
MeanAcrossSbj_PathL_1_rand(:,:)=mean(PathL_1_rand,1);
MeanAcrossSbj_MNodeBetw_1_rand(:,:)=mean(MNodeBetw_1_rand,1);
MeanAcrossSbj_MEdgeBetw_1_rand(:,:)=mean(MEdgeBetw_1_rand,1);
MeanAcrossSbj_Lambda_1_rand(:,:)=mean(Lambda_1_rand,1);
MeanAcrossSbj_Gamma_1_rand(:,:)=mean(Gamma_1_rand,1);
MeanAcrossSbj_Sigma_1_rand(:,:)=mean(Sigma_1_rand,1);

MeanAcrossSbj_MClust_2_rand(:,:)=mean(MClust_2_rand,1);
MeanAcrossSbj_MDeg_2_rand(:,:)=mean(MDeg_2_rand,1);
MeanAcrossSbj_Trans_2_rand(:,:)=mean(Trans_2_rand,1);
MeanAcrossSbj_Assort_2_rand(:,:)=mean(Assort_2_rand,1);
MeanAcrossSbj_GEff_2_rand(:,:)=mean(GEff_2_rand,1);
MeanAcrossSbj_MLocEff_2_rand(:,:)=mean(MLocEff_2_rand,1);
MeanAcrossSbj_Mod_2_rand(:,:)=mean(Mod_2_rand,1);
MeanAcrossSbj_ModL_2_rand(:,:)=mean(ModL_2_rand,1);
MeanAcrossSbj_PathL_2_rand(:,:)=mean(PathL_2_rand,1);
MeanAcrossSbj_MNodeBetw_2_rand(:,:)=mean(MNodeBetw_2_rand,1);
MeanAcrossSbj_MEdgeBetw_2_rand(:,:)=mean(MEdgeBetw_2_rand,1);
MeanAcrossSbj_Lambda_2_rand(:,:)=mean(Lambda_2_rand,1);
MeanAcrossSbj_Gamma_2_rand(:,:)=mean(Gamma_2_rand,1);
MeanAcrossSbj_Sigma_2_rand(:,:)=mean(Sigma_2_rand,1);

MeanAcrossSbj_MClust_1=mean(MClust_1);
MeanAcrossSbj_MDeg_1=mean(MDeg_1);
MeanAcrossSbj_Trans_1=mean(Trans_1);
MeanAcrossSbj_Assort_1=mean(Assort_1);
MeanAcrossSbj_GEff_1=mean(GEff_1);
MeanAcrossSbj_MLocEff_1=mean(MLocEff_1);
MeanAcrossSbj_Mod_1=mean(Mod_1);
MeanAcrossSbj_ModL_1=mean(ModL_1);
MeanAcrossSbj_PathL_1=mean(PathL_1);
MeanAcrossSbj_MNodeBetw_1=mean(MNodeBetw_1);
MeanAcrossSbj_MEdgeBetw_1=mean(MEdgeBetw_1);
MeanAcrossSbj_Lambda_1=mean(Lambda_1);
MeanAcrossSbj_Gamma_1=mean(Gamma_1);
MeanAcrossSbj_Sigma_1=mean(Sigma_1);

MeanAcrossSbj_MClust_2=mean(MClust_2);
MeanAcrossSbj_MDeg_2=mean(MDeg_2);
MeanAcrossSbj_Trans_2=mean(Trans_2);
MeanAcrossSbj_Assort_2=mean(Assort_2);
MeanAcrossSbj_GEff_2=mean(GEff_2);
MeanAcrossSbj_MLocEff_2=mean(MLocEff_2);
MeanAcrossSbj_Mod_2=mean(Mod_2);
MeanAcrossSbj_ModL_2=mean(ModL_2);
MeanAcrossSbj_PathL_2=mean(PathL_2);
MeanAcrossSbj_MNodeBetw_2=mean(MNodeBetw_2);
MeanAcrossSbj_MEdgeBetw_2=mean(MEdgeBetw_2);
MeanAcrossSbj_Lambda_2=mean(Lambda_2);
MeanAcrossSbj_Gamma_2=mean(Gamma_2);
MeanAcrossSbj_Sigma_2=mean(Sigma_2);

%% 

[mu_MClust1_rand,sigma_MClust1_rand,muCi_MClust1_rand,sigmaCi_MClust1_rand]=normfit(MeanAcrossSbj_MClust_1_rand',Alpha);
[mu_MDeg1_rand,sigma_MDeg1_rand,muCi_MDeg1_rand,sigmaCi_MDeg1_rand]=normfit(MeanAcrossSbj_MDeg_1_rand',Alpha);
[mu_Trans1_rand,sigma_Trans1_rand,muCi_Trans1_rand,sigmaCi_Trans1_rand]=normfit(MeanAcrossSbj_Trans_1_rand',Alpha);
[mu_Assort1_rand,sigma_Assort1_rand,muCi_Assort1_rand,sigmaCi_Assort1_rand]=normfit(MeanAcrossSbj_Assort_1_rand',Alpha);
[mu_GEff1_rand,sigma_GEff1_rand,muCi_GEff1_rand,sigmaCi_GEff1_rand]=normfit(MeanAcrossSbj_GEff_1_rand',Alpha);
[mu_MLocEff1_rand,sigma_MLocEff1_rand,muCi_MLocEff1_rand,sigmaCi_MLocEff1_rand]=normfit(MeanAcrossSbj_MLocEff_1_rand',Alpha);
[mu_Mod1_rand,sigma_Mod1_rand,muCi_Mod1_rand,sigmaCi_Mod1_rand]=normfit(MeanAcrossSbj_Mod_1_rand',Alpha);
[mu_ModL1_rand,sigma_ModL1_rand,muCi_ModL1_rand,sigmaCi_ModL1_rand]=normfit(MeanAcrossSbj_ModL_1_rand',Alpha);
[mu_PathL1_rand,sigma_PathL1_rand,muCi_PathL1_rand,sigmaCi_PathL1_rand]=normfit(MeanAcrossSbj_PathL_1_rand',Alpha);
[mu_MNodeBetw1_rand,sigma_MNodeBetw1_rand,muCi_MNodeBetw1_rand,sigmaCi_MNodeBetw1_rand]=normfit(MeanAcrossSbj_MNodeBetw_1_rand',Alpha);
[mu_MEdgeBetw1_rand,sigma_MEdgeBetw1_rand,muCi_MEdgeBetw1_rand,sigmaCi_MEdgeBetw1_rand]=normfit(MeanAcrossSbj_MEdgeBetw_1_rand',Alpha);
[mu_Lambda1_rand,sigma_Lambda1_rand,muCi_Lambda1_rand,sigmaCi_Lambda1_rand]=normfit(MeanAcrossSbj_Lambda_1_rand',Alpha);
[mu_Gamma1_rand,sigma_Gamma1_rand,muCi_Gamma1_rand,sigmaCi_Gamma1_rand]=normfit(MeanAcrossSbj_Gamma_1_rand',Alpha);
[mu_Sigma1_rand,sigma_Sigma1_rand,muCi_Sigma1_rand,sigmaCi_Sigma1_rand]=normfit(MeanAcrossSbj_Sigma_1_rand',Alpha);

[mu_MClust2_rand,sigma_MClust2_rand,muCi_MClust2_rand,sigmaCi_MClust2_rand]=normfit(MeanAcrossSbj_MClust_2_rand',Alpha);
 [mu_MDeg2_rand,sigma_MDeg2_rand,muCi_MDeg2_rand,sigmaCi_MDeg2_rand]=normfit(MeanAcrossSbj_MDeg_2_rand',Alpha);
[mu_Trans2_rand,sigma_Trans2_rand,muCi_Trans2_rand,sigmaCi_Trans2_rand]=normfit(MeanAcrossSbj_Trans_2_rand',Alpha);
[mu_Assort2_rand,sigma_Assort2_rand,muCi_Assort2_rand,sigmaCi_Assort2_rand]=normfit(MeanAcrossSbj_Assort_2_rand',Alpha);
[mu_GEff2_rand,sigma_GEff2_rand,muCi_GEff2_rand,sigmaCi_GEff2_rand]=normfit(MeanAcrossSbj_GEff_2_rand',Alpha);
[mu_MLocEff2_rand,sigma_MLocEff2_rand,muCi_MLocEff2_rand,sigmaCi_MLocEff2_rand]=normfit(MeanAcrossSbj_MLocEff_2_rand',Alpha);
[mu_Mod2_rand,sigma_Mod2_rand,muCi_Mod2_rand,sigmaCi_Mod2_rand]=normfit(MeanAcrossSbj_Mod_2_rand',Alpha);
[mu_ModL2_rand,sigma_ModL2_rand,muCi_ModL2_rand,sigmaCi_ModL2_rand]=normfit(MeanAcrossSbj_ModL_2_rand',Alpha);
[mu_PathL2_rand,sigma_PathL2_rand,muCi_PathL2_rand,sigmaCi_PathL2_rand]=normfit(MeanAcrossSbj_PathL_2_rand',Alpha);
[mu_MNodeBetw2_rand,sigma_MNodeBetw2_rand,muCi_MNodeBetw2_rand,sigmaCi_MNodeBetw2_rand]=normfit(MeanAcrossSbj_MNodeBetw_2_rand',Alpha);
[mu_MEdgeBetw2_rand,sigma_MEdgeBetw2_rand,muCi_MEdgeBetw2_rand,sigmaCi_MEdgeBetw2_rand]=normfit(MeanAcrossSbj_MEdgeBetw_2_rand',Alpha);
[mu_Lambda2_rand,sigma_Lambda2_rand,muCi_Lambda2_rand,sigmaCi_Lambda2_rand]=normfit(MeanAcrossSbj_Lambda_2_rand',Alpha);
[mu_Gamma2_rand,sigma_Gamma2_rand,muCi_Gamma2_rand,sigmaCi_Gamma2_rand]=normfit(MeanAcrossSbj_Gamma_2_rand',Alpha);
[mu_Sigma2_rand,sigma_Sigma2_rand,muCi_Sigma2_rand,sigmaCi_Sigma2_rand]=normfit(MeanAcrossSbj_Sigma_2_rand',Alpha);

%%

N_rand=size(N_G1,2);
Pvalue = Alpha;

Ci_MClust=CL_per(MeanAcrossSbj_MClust_2_rand'-MeanAcrossSbj_MClust_1_rand',Pvalue);
Ci_Trans=CL_per(MeanAcrossSbj_Trans_2_rand'-MeanAcrossSbj_Trans_1_rand',Pvalue);
Ci_Assort=CL_per(MeanAcrossSbj_Assort_2_rand'-MeanAcrossSbj_Assort_1_rand',Pvalue);
Ci_GEff=CL_per(MeanAcrossSbj_GEff_2_rand'-MeanAcrossSbj_GEff_1_rand',Pvalue);
Ci_MLocEff=CL_per(MeanAcrossSbj_MLocEff_2_rand'-MeanAcrossSbj_MLocEff_1_rand',Pvalue);
Ci_Mod=CL_per(MeanAcrossSbj_Mod_2_rand'-MeanAcrossSbj_Mod_1_rand',Pvalue);
Ci_ModL=CL_per(MeanAcrossSbj_ModL_2_rand'-MeanAcrossSbj_ModL_1_rand',Pvalue);
Ci_PathL=CL_per(MeanAcrossSbj_PathL_2_rand'-MeanAcrossSbj_PathL_1_rand',Pvalue);
Ci_MNodeBetw=CL_per(MeanAcrossSbj_MNodeBetw_2_rand'-MeanAcrossSbj_MNodeBetw_1_rand',Pvalue);
Ci_MEdgeBetw=CL_per(MeanAcrossSbj_MEdgeBetw_2_rand'-MeanAcrossSbj_MEdgeBetw_1_rand',Pvalue);
Ci_Lambda=CL_per(MeanAcrossSbj_Lambda_2_rand'-MeanAcrossSbj_Lambda_1_rand',Pvalue);
Ci_Gamma=CL_per(MeanAcrossSbj_Gamma_2_rand'-MeanAcrossSbj_Gamma_1_rand',Pvalue);
Ci_Sigma=CL_per(MeanAcrossSbj_Sigma_2_rand'-MeanAcrossSbj_Sigma_1_rand',Pvalue);


%%

Xax=[MinMesPlot:MesStepPlot:MaxMesPlot];

if ~isequal(size(Xax,2),size(MeanAcrossSbj_MClust_1,2))

        Xax = Xax(1:size(MeanAcrossSbj_MClust_1,2));
end


figure
plot(Xax,MeanAcrossSbj_MClust_2'-MeanAcrossSbj_MClust_1','r*');hold 
plot(Xax,mu_MClust2_rand-mu_MClust1_rand,'bx');
plot(Xax,Ci_MClust(1,:)','b--');
plot(Xax,Ci_MClust(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in clustering coefficient','fontsize',12,'fontweight','b')
title('Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['Clustering_vs_Null_' num2str(Nperm) '.fig'])



figure
plot(Xax,MeanAcrossSbj_Trans_2'-MeanAcrossSbj_Trans_1','r*');hold 
plot(Xax,mu_Trans2_rand-mu_Trans1_rand,'bx');
plot(Xax,Ci_Trans(1,:)','b--');
plot(Xax,Ci_Trans(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in transitivity coefficient','fontsize',12,'fontweight','b')
title('transitivity','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['Transitivity_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_PathL_2'-MeanAcrossSbj_PathL_1','r*');hold 
plot(Xax,mu_PathL2_rand-mu_PathL1_rand,'bx');
plot(Xax,Ci_PathL(1,:)','b--');
plot(Xax,Ci_PathL(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in path length','fontsize',12,'fontweight','b')
title('Char Path length','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['CharPath_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Assort_2'-MeanAcrossSbj_Assort_1','r*');hold 
plot(Xax,mu_Assort2_rand-mu_Assort1_rand,'bx');
plot(Xax,Ci_Assort(1,:)','b--');
plot(Xax,Ci_Assort(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in assortativity','fontsize',12,'fontweight','b')
title('Mean assortativity','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['MeanAssortativity_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_GEff_2'-MeanAcrossSbj_GEff_1','r*');hold 
plot(Xax,mu_GEff2_rand-mu_GEff1_rand,'bx');
plot(Xax,Ci_GEff(1,:)','b--');
plot(Xax,Ci_GEff(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in global efficiency','fontsize',12,'fontweight','b')
title('Global efficiency','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['GlobalEfficiecny_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_MLocEff_2'-MeanAcrossSbj_MLocEff_1','r*');hold 
plot(Xax,mu_MLocEff2_rand-mu_MLocEff1_rand,'bx');
plot(Xax,Ci_MLocEff(1,:)','b--');
plot(Xax,Ci_MLocEff(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in mean local efficiency','fontsize',12,'fontweight','b')
title('Mean local efficiency','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['MeanLocalEfficiecny_vs_Null_' num2str(Nperm) '.fig'])



figure
plot(Xax,MeanAcrossSbj_Mod_2'-MeanAcrossSbj_Mod_1','r*');hold 
plot(Xax,mu_Mod2_rand-mu_Mod1_rand,'bx');
plot(Xax,Ci_Mod(1,:)','b--');
plot(Xax,Ci_Mod(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in modularity','fontsize',12,'fontweight','b')
title('Modularity','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['Modularity_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_ModL_2'-MeanAcrossSbj_ModL_1','r*');hold 
plot(Xax,mu_ModL2_rand-mu_ModL1_rand,'bx');
plot(Xax,Ci_ModL(1,:)','b--');
plot(Xax,Ci_ModL(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in modularity_L','fontsize',12,'fontweight','b')
title('Modularity_L','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['ModularityL_vs_Null_' num2str(Nperm) '.fig'])



figure
plot(Xax,MeanAcrossSbj_MNodeBetw_2'-MeanAcrossSbj_MNodeBetw_1','r*');hold 
plot(Xax,mu_MNodeBetw2_rand-mu_MNodeBetw1_rand,'bx');
plot(Xax,Ci_MNodeBetw(1,:)','b--');
plot(Xax,Ci_MNodeBetw(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in Mean Node Betweenness','fontsize',12,'fontweight','b')
title('Mean Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['MeanNodeBetweenness_vs_Null_' num2str(Nperm) '.fig'])



figure
plot(Xax,MeanAcrossSbj_MEdgeBetw_2'-MeanAcrossSbj_MEdgeBetw_1','r*');hold 
plot(Xax,mu_MEdgeBetw2_rand-mu_MEdgeBetw1_rand,'bx');
plot(Xax,Ci_MEdgeBetw(1,:)','b--');
plot(Xax,Ci_MEdgeBetw(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in Mean Edge Betweenness','fontsize',12,'fontweight','b')
title('Mean Edge Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['MeanEdgeBetweenness_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Lambda_2'-MeanAcrossSbj_Lambda_1','r*');hold 
plot(Xax,mu_Lambda2_rand-mu_Lambda1_rand,'bx');
plot(Xax,Ci_Lambda(1,:)','b--');
plot(Xax,Ci_Lambda(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in Lambda','fontsize',12,'fontweight','b')
title('Lambda','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['Lambda_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Gamma_2'-MeanAcrossSbj_Gamma_1','r*');hold 
plot(Xax,mu_Gamma2_rand-mu_Gamma1_rand,'bx');
plot(Xax,Ci_Gamma(1,:)','b--');
plot(Xax,Ci_Gamma(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in Gamma','fontsize',12,'fontweight','b')
title('Gamma','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['Gamma_vs_Null_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Sigma_2'-MeanAcrossSbj_Sigma_1','r*');hold 
plot(Xax,mu_Sigma2_rand-mu_Sigma1_rand,'bx');
plot(Xax,Ci_Sigma(1,:)','b--');
plot(Xax,Ci_Sigma(2,:)','b--');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Difference in Sigma','fontsize',12,'fontweight','b')
title('Sigma','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave(['Sigma_vs_Null_' num2str(Nperm) '.fig'])


%% 

figure
plot(Xax,MeanAcrossSbj_MClust_1','b*');hold 
plot(Xax,MeanAcrossSbj_MClust_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Mean clustering coefficient','fontsize',12,'fontweight','b')
title('Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['Clustering_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Trans_1','b*');hold 
plot(Xax,MeanAcrossSbj_Trans_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Transitivity','fontsize',12,'fontweight','b')
title('Transitivity','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['Transitivity_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Assort_1','b*');hold 
plot(Xax,MeanAcrossSbj_Assort_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Assortativity','fontsize',12,'fontweight','b')
title('Assortativity','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['Assortativity_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_GEff_1','b*');hold 
plot(Xax,MeanAcrossSbj_GEff_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Global Efficiency','fontsize',12,'fontweight','b')
title('Global Efficicency','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['GEfficiency_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_MLocEff_1','b*');hold 
plot(Xax,MeanAcrossSbj_MLocEff_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Local Efficiency','fontsize',12,'fontweight','b')
title('Local Efficicency','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['LocalEfficiency_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Mod_1','b*');hold 
plot(Xax,MeanAcrossSbj_Mod_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Modularity','fontsize',12,'fontweight','b')
title('Modularity','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['Modularity_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_ModL_1','b*');hold 
plot(Xax,MeanAcrossSbj_ModL_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Modularity_L','fontsize',12,'fontweight','b')
title('Modularity_L','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['ModularityL_Value_' num2str(Nperm) '.fig'])

figure
plot(Xax,MeanAcrossSbj_PathL_1','b*');hold 
plot(Xax,MeanAcrossSbj_PathL_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Characteristic path length','fontsize',12,'fontweight','b')
title('Characteristic path length','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['CharPathLength_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_MNodeBetw_1','b*');hold 
plot(Xax,MeanAcrossSbj_MNodeBetw_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Mean Node Betweenness','fontsize',12,'fontweight','b')
title('Mean Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['MNodeBetweeness_Value_' num2str(Nperm) '.fig'])

figure
plot(Xax,MeanAcrossSbj_MEdgeBetw_1','b*');hold 
plot(Xax,MeanAcrossSbj_MEdgeBetw_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Mean Edge Betweenness','fontsize',12,'fontweight','b')
title('Mean Edge Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['MEdgeBetweeness_Value_' num2str(Nperm) '.fig'])

figure
plot(Xax,MeanAcrossSbj_Lambda_1','b*');hold 
plot(Xax,MeanAcrossSbj_Lambda_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Lambda','fontsize',12,'fontweight','b')
title('Lambda','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['Lambda_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Gamma_1','b*');hold 
plot(Xax,MeanAcrossSbj_Gamma_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Gamma','fontsize',12,'fontweight','b')
title('Gamma','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['Gamma_Value_' num2str(Nperm) '.fig'])


figure
plot(Xax,MeanAcrossSbj_Sigma_1','b*');hold 
plot(Xax,MeanAcrossSbj_Sigma_2','r*');
xlabel('Density','fontsize',12,'fontweight','b')
ylabel('Sigma','fontsize',12,'fontweight','b')
title('Sigma','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1,Group2)
grid on
hgsave(['Sigma_Value_' num2str(Nperm) '.fig'])

%%

p_assort=CL_Pval(MeanAcrossSbj_Assort_2_rand-MeanAcrossSbj_Assort_1_rand,MeanAcrossSbj_Assort_2-MeanAcrossSbj_Assort_1,'Assort',Tail);
p_GEff=CL_Pval(MeanAcrossSbj_GEff_2_rand-MeanAcrossSbj_GEff_1_rand,MeanAcrossSbj_GEff_2-MeanAcrossSbj_GEff_1,'GEff',Tail);
p_Gamma=CL_Pval(MeanAcrossSbj_Gamma_2_rand-MeanAcrossSbj_Gamma_1_rand,MeanAcrossSbj_Gamma_2-MeanAcrossSbj_Gamma_1,'Gamma',Tail);
p_Lambda=CL_Pval(MeanAcrossSbj_Lambda_2_rand-MeanAcrossSbj_Lambda_1_rand,MeanAcrossSbj_Lambda_2-MeanAcrossSbj_Lambda_1,'Lambda',Tail);
p_MClust=CL_Pval(MeanAcrossSbj_MClust_2_rand-MeanAcrossSbj_MClust_1_rand,MeanAcrossSbj_MClust_2-MeanAcrossSbj_MClust_1,'MClust',Tail);
p_MEdgeBetw=CL_Pval(MeanAcrossSbj_MEdgeBetw_2_rand-MeanAcrossSbj_MEdgeBetw_1_rand,MeanAcrossSbj_MEdgeBetw_2-MeanAcrossSbj_MEdgeBetw_1,'MEdgeBetw',Tail);
p_MLocEff=CL_Pval(MeanAcrossSbj_MLocEff_2_rand-MeanAcrossSbj_MLocEff_1_rand,MeanAcrossSbj_MLocEff_2-MeanAcrossSbj_MLocEff_1,'MLocEff',Tail);
p_MNodeBetw=CL_Pval(MeanAcrossSbj_MNodeBetw_2_rand-MeanAcrossSbj_MNodeBetw_1_rand,MeanAcrossSbj_MNodeBetw_2-MeanAcrossSbj_MNodeBetw_1,'MNodeBetw',Tail);
p_ModL=CL_Pval(MeanAcrossSbj_ModL_2_rand-MeanAcrossSbj_ModL_1_rand,MeanAcrossSbj_ModL_2-MeanAcrossSbj_ModL_1,'ModL',Tail);
p_Mod=CL_Pval(MeanAcrossSbj_Mod_2_rand-MeanAcrossSbj_Mod_1_rand,MeanAcrossSbj_Mod_2-MeanAcrossSbj_Mod_1,'Mod',Tail);
p_PathL=CL_Pval(MeanAcrossSbj_PathL_2_rand-MeanAcrossSbj_PathL_1_rand,MeanAcrossSbj_PathL_2-MeanAcrossSbj_PathL_1,'PathL',Tail);
p_Sigma=CL_Pval(MeanAcrossSbj_Sigma_2_rand-MeanAcrossSbj_Sigma_1_rand,MeanAcrossSbj_Sigma_2-MeanAcrossSbj_Sigma_1,'Sigma',Tail);
p_Trans=CL_Pval(MeanAcrossSbj_Trans_2_rand-MeanAcrossSbj_Trans_1_rand,MeanAcrossSbj_Trans_2-MeanAcrossSbj_Trans_1,'Trans',Tail);
p_MDeg=CL_Pval(MeanAcrossSbj_MDeg_2_rand-MeanAcrossSbj_MDeg_1_rand,MeanAcrossSbj_MDeg_2-MeanAcrossSbj_MDeg_1,'MDeg',Tail);




%% 2-SAMPLE T-TEST

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinThr=input('select your final decision on minimum threshold:');MinIdx=find(single(xxx)==single(MinThr));
MaxThr=input('select your final decision on maximum threshold:');MaxIdx=find(single(xxx)==single(MaxThr));

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


p = Alpha;
if isequal(Tail,1)
    
    Model = 'left';

elseif isequal(Tail,2)
    
    Model = 'both';

end


[h_Trans,pTTest_Trans,ci_Trans,stats_Trans]=ttest2(Trans_1,Trans_2,p,Model);
[h_Assort,pTTest_Assort,ci_Assort,stats_Assort]=ttest2(Assort_1,Assort_2,p,Model);
[h_LocEff,pTTest_LocEff,ci_LocEff,stats_LocEff]=ttest2(MLocEff_1,MLocEff_2,p,Model);
[h_ModL,pTTest_ModL,ci_ModL,stats_ModL]=ttest2(ModL_1,ModL_2,p,Model);
[h_MNodeBetw,pTTest_MNodeBetw,ci_MNodeBetw,stats_MNodeBetw]=ttest2(MNodeBetw_1,MNodeBetw_2,p,Model);
[h_MEdgeBetw,pTTest_MEdgeBetw,ci_MEdgeBetw,stats_MEdgeBetw]=ttest2(MEdgeBetw_1,MEdgeBetw_2,p,Model);
[h_Lambda,pTTest_Lambda,ci_Lambda,stats_Lambda]=ttest2(Lambda_1,Lambda_2,p,Model);
[h_Gamma,pTTest_Gamma,ci_Gamma,stats_Gamma]=ttest2(Gamma_1,Gamma_2,p,Model);
[h_Sigma,pTTest_Sigma,ci_Sigma,stats_Sigma]=ttest2(Sigma_1,Sigma_2,p,Model);


MeanAcrossSbj_MClust_1_ttest=mean(MClust_1,1);
MeanAcrossSbj_MDeg_1_ttest=mean(MDeg_1,1);
MeanAcrossSbj_Trans_1_ttest=mean(Trans_1,1);
MeanAcrossSbj_Assort_1_ttest=mean(Assort_1,1);
MeanAcrossSbj_GEff_1_ttest=mean(GEff_1,1);
MeanAcrossSbj_MLocEff_1_ttest=mean(MLocEff_1,1);
MeanAcrossSbj_Mod_1_ttest=mean(Mod_1,1);
MeanAcrossSbj_ModL_1_ttest=mean(ModL_1,1);
MeanAcrossSbj_PathL_1_ttest=mean(PathL_1,1);
MeanAcrossSbj_MNodeBetw_1_ttest=mean(MNodeBetw_1,1);
MeanAcrossSbj_MEdgeBetw_1_ttest=mean(MEdgeBetw_1,1);
MeanAcrossSbj_Lambda_1_ttest=mean(Lambda_1,1);
MeanAcrossSbj_Gamma_1_ttest=mean(Gamma_1,1);
MeanAcrossSbj_Sigma_1_ttest=mean(Sigma_1,1);


MeanAcrossSbj_MClust_2_ttest=mean(MClust_2,1);
MeanAcrossSbj_MDeg_2_ttest=mean(MDeg_2,1);
%MeanAcrossSbj_MStr_2_ttest=mean(MStr_2,1);
MeanAcrossSbj_Trans_2_ttest=mean(Trans_2,1);
MeanAcrossSbj_Assort_2_ttest=mean(Assort_2,1);
MeanAcrossSbj_GEff_2_ttest=mean(GEff_2,1);
MeanAcrossSbj_MLocEff_2_ttest=mean(MLocEff_2,1);
MeanAcrossSbj_Mod_2_ttest=mean(Mod_2,1);
MeanAcrossSbj_ModL_2_ttest=mean(ModL_2,1);
MeanAcrossSbj_PathL_2_ttest=mean(PathL_2,1);
MeanAcrossSbj_MNodeBetw_2_ttest=mean(MNodeBetw_2,1);
MeanAcrossSbj_MEdgeBetw_2_ttest=mean(MEdgeBetw_2,1);
MeanAcrossSbj_Lambda_2_ttest=mean(Lambda_2,1);
MeanAcrossSbj_Gamma_2_ttest=mean(Gamma_2,1);
MeanAcrossSbj_Sigma_2_ttest=mean(Sigma_2,1);


dd=pwd;     
mkdir('TTest_Results');
cd([dd '/TTest_Results']);

save(['pVal_TTest'],'pTTest_Trans','pTTest_Assort','pTTest_LocEff','pTTest_ModL','pTTest_MNodeBetw', ...
     'pTTest_MEdgeBetw','pTTest_Lambda','pTTest_Gamma','pTTest_Sigma');


fprintf('%-4s\n','......Done ...');

