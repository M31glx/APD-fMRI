function AUC_Analysis_NetMesf


GUIinput=0;

if nargin<1

    data_log = spm_select(1,'mat4GATf*','Select mat4GATf .mat (output from "Network Measures")');
    
    ff = load(data_log,'mat4GATf');
    Group1 = ff.mat4GATf.g1;Group2 = ff.mat4GATf.g2;
    
    if isfield(ff.mat4GATf,'tail') && isfield(ff.mat4GATf,'alpha') 
        
        Alpha = ff.mat4GATf.alpha;
        Tail = ff.mat4GATf.tail;
        
    else
        
        Alpha = .05;
        Tail = 2;
    
    end
    
    NetM = spm_select(1,'NetMesPerDens_f.mat','Select NetMesPerDens_f.mat file ');
    OutputFName='Average_Across_Densities_Results_f.mat';
    
    fprintf('%-4s\n',['loading inputs...']);
    
    load(NetM);
    MinMesPlot=ff.mat4GATf.MinThr;
    MaxMesPlot=ff.mat4GATf.MaxThr;
    MesStepPlot=ff.mat4GATf.MesStep;
    
    GUIinput=1;
end


xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinThr=input('select your final decision on minimum threshold:');MinIdx=find(single(xxx)==single(MinThr));
MaxThr=input('select your final decision on maximum threshold:');MaxIdx=find(single(xxx)==single(MaxThr));

if MaxIdx > size(Dens_1_rand,2)
    
    MaxIdx = size(Dens_1_rand,2);
    
end

%% 

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

%%

fprintf('%-4s\n',['saving net measures across new threshold range...']);

dd=pwd;     
mkdir('AUC_Results');
cd([dd '/AUC_Results']);

save NetMesPerDens_AveAcrossDens_f Dens_1_rand MClust_1_rand MDeg_1_rand Trans_1_rand ...
     Assort_1_rand GEff_1_rand MLocEff_1_rand Mod_1_rand ModL_1_rand ...
     PathL_1_rand MNodeBetw_1_rand MEdgeBetw_1_rand Lambda_1_rand Gamma_1_rand ...
     Sigma_1_rand Dens_2_rand MClust_2_rand MDeg_2_rand Trans_2_rand ...
     Assort_2_rand GEff_2_rand MLocEff_2_rand Mod_2_rand ModL_2_rand ...
     PathL_2_rand MNodeBetw_2_rand MEdgeBetw_2_rand Lambda_2_rand Gamma_2_rand ...
     Sigma_2_rand Dens_1 MClust_1 MDeg_1 Trans_1 Assort_1 GEff_1 MLocEff_1 ...
     Mod_1 ModL_1 PathL_1 MNodeBetw_1 MEdgeBetw_1 Lambda_1 Gamma_1 Sigma_1 ...
     Dens_2 MClust_2 MDeg_2 Trans_2 Assort_2 GEff_2 MLocEff_2 Mod_2 ModL_2 PathL_2 MNodeBetw_2 ...
     MEdgeBetw_2 Lambda_2 Gamma_2 Sigma_2 


%% 

Xax=MinMesPlot:MesStepPlot:MaxMesPlot;
iXax=Xax(MinIdx:MaxIdx);


auc_MClust_1 = trapz(iXax,MClust_1,2);
auc_Trans_1 = trapz(iXax,Trans_1,2);
auc_Assort_1 = trapz(iXax,Assort_1,2);
auc_GEff_1 = trapz(iXax,GEff_1,2);
auc_MLocEff_1 = trapz(iXax,MLocEff_1,2);
auc_Mod_1 = trapz(iXax,Mod_1,2);
auc_ModL_1 = trapz(iXax,ModL_1,2);
auc_PathL_1 = trapz(iXax,PathL_1,2);
auc_MNodeBetw_1 = trapz(iXax,MNodeBetw_1,2);
auc_MEdgeBetw_1 = trapz(iXax,MEdgeBetw_1,2);
auc_Lambda_1 = trapz(iXax,Lambda_1,2);
auc_Gamma_1 = trapz(iXax,Gamma_1,2);
auc_Sigma_1 = trapz(iXax,Sigma_1,2);


auc_MClust_2 = trapz(iXax,MClust_2,2);
auc_Trans_2 = trapz(iXax,Trans_2,2);
auc_Assort_2 = trapz(iXax,Assort_2,2);
auc_GEff_2 = trapz(iXax,GEff_2,2);
auc_MLocEff_2 = trapz(iXax,MLocEff_2,2);
auc_Mod_2 = trapz(iXax,Mod_2,2);
auc_ModL_2 = trapz(iXax,ModL_2,2);
auc_PathL_2 = trapz(iXax,PathL_2,2);
auc_MNodeBetw_2 = trapz(iXax,MNodeBetw_2,2);
auc_MEdgeBetw_2 = trapz(iXax,MEdgeBetw_2,2);
auc_Lambda_2 = trapz(iXax,Lambda_2,2);
auc_Gamma_2 = trapz(iXax,Gamma_2,2);
auc_Sigma_2 = trapz(iXax,Sigma_2,2);


auc_MClust_1_rand = trapz(iXax,MClust_1_rand,2);
auc_Trans_1_rand = trapz(iXax,Trans_1_rand,2);
auc_Assort_1_rand = trapz(iXax,Assort_1_rand,2);
auc_GEff_1_rand = trapz(iXax,GEff_1_rand,2);
auc_MLocEff_1_rand = trapz(iXax,MLocEff_1_rand,2);
auc_Mod_1_rand = trapz(iXax,Mod_1_rand,2);
auc_ModL_1_rand = trapz(iXax,ModL_1_rand,2);
auc_PathL_1_rand = trapz(iXax,PathL_1_rand,2);
auc_MNodeBetw_1_rand = trapz(iXax,MNodeBetw_1_rand,2);
auc_MEdgeBetw_1_rand = trapz(iXax,MEdgeBetw_1_rand,2);
auc_Lambda_1_rand = trapz(iXax,Lambda_1_rand,2);
auc_Gamma_1_rand = trapz(iXax,Gamma_1_rand,2);
auc_Sigma_1_rand = trapz(iXax,Sigma_1_rand,2);


auc_MClust_2_rand = trapz(iXax,MClust_2_rand,2);
auc_Trans_2_rand = trapz(iXax,Trans_2_rand,2);
auc_Assort_2_rand = trapz(iXax,Assort_2_rand,2);
auc_GEff_2_rand = trapz(iXax,GEff_2_rand,2);
auc_MLocEff_2_rand = trapz(iXax,MLocEff_2_rand,2);
auc_Mod_2_rand = trapz(iXax,Mod_2_rand,2);
auc_ModL_2_rand = trapz(iXax,ModL_2_rand,2);
auc_PathL_2_rand = trapz(iXax,PathL_2_rand,2);
auc_MNodeBetw_2_rand = trapz(iXax,MNodeBetw_2_rand,2);
auc_MEdgeBetw_2_rand = trapz(iXax,MEdgeBetw_2_rand,2);
auc_Lambda_2_rand = trapz(iXax,Lambda_2_rand,2);
auc_Gamma_2_rand = trapz(iXax,Gamma_2_rand,2);
auc_Sigma_2_rand = trapz(iXax,Sigma_2_rand,2);


save AUC_NetMesPerDens_AveAcrossDens_f  auc_MClust_1_rand auc_Trans_1_rand ...
     auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
     auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
     auc_Sigma_1_rand auc_MClust_2_rand auc_Trans_2_rand ...
     auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
     auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
     auc_Sigma_2_rand auc_MClust_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
     auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 ...
     auc_MClust_2 auc_Trans_2 auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
     auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2


%% 

nRand=size(auc_MClust_1_rand,3);

auc_MClust_1 = mean(auc_MClust_1);
auc_Trans_1 = mean(auc_Trans_1);
auc_Assort_1 = mean(auc_Assort_1);
auc_GEff_1 = mean(auc_GEff_1);
auc_MLocEff_1 = mean(auc_MLocEff_1);
auc_Mod_1 = mean(auc_Mod_1);
auc_ModL_1 = mean(auc_ModL_1);
auc_PathL_1 = mean(auc_PathL_1);
auc_MNodeBetw_1 = mean(auc_MNodeBetw_1);
auc_MEdgeBetw_1 = mean(auc_MEdgeBetw_1);
auc_Lambda_1 = mean(auc_Lambda_1);
auc_Gamma_1 = mean(auc_Gamma_1);
auc_Sigma_1 = mean(auc_Sigma_1);


auc_MClust_2 = mean(auc_MClust_2);
auc_Trans_2 = mean(auc_Trans_2);
auc_Assort_2 = mean(auc_Assort_2);
auc_GEff_2 = mean(auc_GEff_2);
auc_MLocEff_2 = mean(auc_MLocEff_2);
auc_Mod_2 = mean(auc_Mod_2);
auc_ModL_2 = mean(auc_ModL_2);
auc_PathL_2 = mean(auc_PathL_2);
auc_MNodeBetw_2 = mean(auc_MNodeBetw_2);
auc_MEdgeBetw_2 = mean(auc_MEdgeBetw_2);
auc_Lambda_2 = mean(auc_Lambda_2);
auc_Gamma_2 = mean(auc_Gamma_2);
auc_Sigma_2 = mean(auc_Sigma_2);


auc_MClust_1_rand = mean(auc_MClust_1_rand,1);auc_MClust_1_rand = reshape(auc_MClust_1_rand,1,nRand);
auc_Trans_1_rand = mean(auc_Trans_1_rand,1);auc_Trans_1_rand = reshape(auc_Trans_1_rand ,1,nRand);
auc_Assort_1_rand = mean(auc_Assort_1_rand,1);auc_Assort_1_rand = reshape(auc_Assort_1_rand ,1,nRand);
auc_GEff_1_rand = mean(auc_GEff_1_rand,1);auc_GEff_1_rand = reshape(auc_GEff_1_rand ,1,nRand);
auc_MLocEff_1_rand = mean(auc_MLocEff_1_rand,1);auc_MLocEff_1_rand = reshape(auc_MLocEff_1_rand ,1,nRand);
auc_Mod_1_rand = mean(auc_Mod_1_rand,1);auc_Mod_1_rand = reshape(auc_Mod_1_rand ,1,nRand);
auc_ModL_1_rand = mean(auc_ModL_1_rand,1);auc_ModL_1_rand = reshape(auc_ModL_1_rand ,1,nRand);
auc_PathL_1_rand = mean(auc_PathL_1_rand,1);auc_PathL_1_rand = reshape(auc_PathL_1_rand ,1,nRand);
auc_MNodeBetw_1_rand = mean(auc_MNodeBetw_1_rand,1);auc_MNodeBetw_1_rand = reshape(auc_MNodeBetw_1_rand ,1,nRand);
auc_MEdgeBetw_1_rand = mean(auc_MEdgeBetw_1_rand,1);auc_MEdgeBetw_1_rand = reshape(auc_MEdgeBetw_1_rand ,1,nRand);
auc_Lambda_1_rand = mean(auc_Lambda_1_rand,1);auc_Lambda_1_rand = reshape( auc_Lambda_1_rand,1,nRand);
auc_Gamma_1_rand = mean(auc_Gamma_1_rand,1);auc_Gamma_1_rand = reshape(auc_Gamma_1_rand ,1,nRand);
auc_Sigma_1_rand = mean(auc_Sigma_1_rand,1);auc_Sigma_1_rand = reshape(auc_Sigma_1_rand ,1,nRand);


auc_MClust_2_rand = mean(auc_MClust_2_rand,1);auc_MClust_2_rand = reshape(auc_MClust_2_rand,1,nRand);
auc_Trans_2_rand = mean(auc_Trans_2_rand,1);auc_Trans_2_rand = reshape(auc_Trans_2_rand ,1,nRand);
auc_Assort_2_rand = mean(auc_Assort_2_rand,1);auc_Assort_2_rand = reshape(auc_Assort_2_rand ,1,nRand);
auc_GEff_2_rand = mean(auc_GEff_2_rand,1);auc_GEff_2_rand = reshape(auc_GEff_2_rand ,1,nRand);
auc_MLocEff_2_rand = mean(auc_MLocEff_2_rand,1);auc_MLocEff_2_rand = reshape(auc_MLocEff_2_rand ,1,nRand);
auc_Mod_2_rand = mean(auc_Mod_2_rand,1);auc_Mod_2_rand = reshape(auc_Mod_2_rand ,1,nRand);
auc_ModL_2_rand = mean(auc_ModL_2_rand,1);auc_ModL_2_rand = reshape(auc_ModL_2_rand ,1,nRand);
auc_PathL_2_rand = mean(auc_PathL_2_rand,1);auc_PathL_2_rand = reshape(auc_PathL_2_rand ,1,nRand);
auc_MNodeBetw_2_rand = mean(auc_MNodeBetw_2_rand,1);auc_MNodeBetw_2_rand = reshape(auc_MNodeBetw_2_rand ,1,nRand);
auc_MEdgeBetw_2_rand = mean(auc_MEdgeBetw_2_rand,1);auc_MEdgeBetw_2_rand = reshape(auc_MEdgeBetw_2_rand ,1,nRand);
auc_Lambda_2_rand = mean(auc_Lambda_2_rand,1);auc_Lambda_2_rand = reshape( auc_Lambda_2_rand,1,nRand);
auc_Gamma_2_rand = mean(auc_Gamma_2_rand,1);auc_Gamma_2_rand = reshape(auc_Gamma_2_rand ,1,nRand);
auc_Sigma_2_rand = mean(auc_Sigma_2_rand,1);auc_Sigma_2_rand = reshape(auc_Sigma_2_rand ,1,nRand);

save AUC_Mean auc_MClust_1_rand auc_Trans_1_rand ...
     auc_Assort_1_rand auc_GEff_1_rand auc_MLocEff_1_rand auc_Mod_1_rand auc_ModL_1_rand ...
     auc_PathL_1_rand auc_MNodeBetw_1_rand auc_MEdgeBetw_1_rand auc_Lambda_1_rand auc_Gamma_1_rand ...
     auc_Sigma_1_rand auc_MClust_2_rand auc_Trans_2_rand ...
     auc_Assort_2_rand auc_GEff_2_rand auc_MLocEff_2_rand auc_Mod_2_rand auc_ModL_2_rand ...
     auc_PathL_2_rand auc_MNodeBetw_2_rand auc_MEdgeBetw_2_rand auc_Lambda_2_rand auc_Gamma_2_rand ...
     auc_Sigma_2_rand auc_MClust_1 auc_Trans_1 auc_Assort_1 auc_GEff_1 auc_MLocEff_1 ...
     auc_Mod_1 auc_ModL_1 auc_PathL_1 auc_MNodeBetw_1 auc_MEdgeBetw_1 auc_Lambda_1 auc_Gamma_1 auc_Sigma_1 ...
     auc_MClust_2 auc_Trans_2 ...
     auc_Assort_2 auc_GEff_2 auc_MLocEff_2 auc_Mod_2 auc_ModL_2 auc_PathL_2 auc_MNodeBetw_2 ...
     auc_MEdgeBetw_2 auc_Lambda_2 auc_Gamma_2 auc_Sigma_2 

%%

p_auc_MClust = CL_Pval(auc_MClust_2_rand - auc_MClust_1_rand, auc_MClust_2 - auc_MClust_1,'AUC_MClust',Tail)
p_auc_Trans = CL_Pval(auc_Trans_2_rand - auc_Trans_1_rand, auc_Trans_2 - auc_Trans_1,'AUC_Trans',Tail)
p_auc_Assort = CL_Pval(auc_Assort_2_rand - auc_Assort_1_rand, auc_Assort_2 - auc_Assort_1,'AUC_Assort',Tail)
p_auc_GEff = CL_Pval(auc_GEff_2_rand - auc_GEff_1_rand, auc_GEff_2 - auc_GEff_1,'AUC_GEff',Tail)
p_auc_MLocEff = CL_Pval(auc_MLocEff_2_rand - auc_MLocEff_1_rand, auc_MLocEff_2 - auc_MLocEff_1,'AUC_MLocEff',Tail)
p_auc_Mod = CL_Pval(auc_Mod_2_rand - auc_Mod_1_rand, auc_Mod_2 - auc_Mod_1,'AUC_Mod',Tail)
p_auc_ModL = CL_Pval(auc_ModL_2_rand - auc_ModL_1_rand, auc_ModL_2 - auc_ModL_1,'AUC_ModL',Tail)
p_auc_PathL = CL_Pval(auc_PathL_2_rand - auc_PathL_1_rand, auc_PathL_2 - auc_PathL_1,'AUC_PathL',Tail)
p_auc_MNodeBetw = CL_Pval(auc_MNodeBetw_2_rand - auc_MNodeBetw_1_rand, auc_MNodeBetw_2 - auc_MNodeBetw_1,'AUC_MNodeBetw',Tail)
p_auc_MEdgeBetw = CL_Pval(auc_MEdgeBetw_2_rand - auc_MEdgeBetw_1_rand, auc_MEdgeBetw_2 - auc_MEdgeBetw_1,'AUC_MEdgeBetw',Tail)
p_auc_Lambda = CL_Pval(auc_Lambda_2_rand - auc_Lambda_1_rand, auc_Lambda_2 - auc_Lambda_1,'AUC_Lambda',Tail)
p_auc_Gamma = CL_Pval(auc_Gamma_2_rand - auc_Gamma_1_rand, auc_Gamma_2 - auc_Gamma_1,'AUC_Gamma',Tail)
p_auc_Sigma = CL_Pval(auc_Sigma_2_rand - auc_Sigma_1_rand, auc_Sigma_2 - auc_Sigma_1,'AUC_Sigma',Tail)


p = Alpha;


CL_auc_Lambda = CL_per(auc_Lambda_2_rand' - auc_Lambda_1_rand', p);
figure; errorbar(1, mean(auc_Lambda_2_rand - auc_Lambda_1_rand), CL_auc_Lambda(1),CL_auc_Lambda(2),'O','LineWidth',1,'MarkerSize',10);
hold;
plot(1,auc_Lambda_2 - auc_Lambda_1,'rv','LineWidth',1.2,'MarkerSize',8);
%set(gca,'YLim',[-0.005 0.005]);

CL_auc_Gamma = CL_per(auc_Gamma_2_rand' - auc_Gamma_1_rand', p);
errorbar(2, mean(auc_Gamma_2_rand - auc_Gamma_1_rand), CL_auc_Gamma(1),CL_auc_Gamma(2),'O','LineWidth',1,'MarkerSize',10);
plot(2,auc_Gamma_2 - auc_Gamma_1,'rv','LineWidth',1.2,'MarkerSize',8);

CL_auc_Sigma = CL_per(auc_Sigma_2_rand' - auc_Sigma_1_rand', p);
errorbar(3, mean(auc_Sigma_2_rand - auc_Sigma_1_rand), CL_auc_Sigma(1),CL_auc_Sigma(2),'O','LineWidth',1,'MarkerSize',10);
plot(3,auc_Sigma_2 - auc_Sigma_1,'rv','LineWidth',1.2,'MarkerSize',8);

set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'path length','clustering','small-worldness'},'FontSize',12)


ylabel('difference in AUC','FontSize',12)
title('Difference in AUC of Network Measures','fontsize',14,'Color','red')
legend('Null mean and 95% CI',[Group2 ' vs. ' Group1])
grid on

% set(gcf,'WindowStyle','normal')
% set(gcf,'PaperPositionMode', 'auto');
% scrsz = get(0,'ScreenSize');
FIG='AUC_Diff_vs_Null';

hgsave(FIG) 
print('-depsc',FIG)
print('-djpeg','-r300',FIG)
print('-dtiff','-r600',FIG)
print('-dpng','-r300',FIG)




