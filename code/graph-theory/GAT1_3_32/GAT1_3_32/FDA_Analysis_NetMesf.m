function FDA_Analysis_NetMesf


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



fprintf('%-4s\n',['saving net measures across new threshold range...']);
dd=pwd;     
mkdir('FDA_Results');
cd([dd '/FDA_Results']);

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

fda_MClust_1 = sum(MClust_1,2);
fda_Trans_1 = sum(Trans_1,2);
fda_Assort_1 = sum(Assort_1,2);
fda_GEff_1 = sum(GEff_1,2);
fda_MLocEff_1 = sum(MLocEff_1,2);
fda_Mod_1 = sum(Mod_1,2);
fda_ModL_1 = sum(ModL_1,2);
fda_PathL_1 = sum(PathL_1,2);
fda_MNodeBetw_1 = sum(MNodeBetw_1,2);
fda_MEdgeBetw_1 = sum(MEdgeBetw_1,2);
fda_Lambda_1 = sum(Lambda_1,2);
fda_Gamma_1 = sum(Gamma_1,2);
fda_Sigma_1 = sum(Sigma_1,2);


fda_MClust_2 = sum(MClust_2,2);
fda_Trans_2 = sum(Trans_2,2);
fda_Assort_2 = sum(Assort_2,2);
fda_GEff_2 = sum(GEff_2,2);
fda_MLocEff_2 = sum(MLocEff_2,2);
fda_Mod_2 = sum(Mod_2,2);
fda_ModL_2 = sum(ModL_2,2);
fda_PathL_2 = sum(PathL_2,2);
fda_MNodeBetw_2 = sum(MNodeBetw_2,2);
fda_MEdgeBetw_2 = sum(MEdgeBetw_2,2);
fda_Lambda_2 = sum(Lambda_2,2);
fda_Gamma_2 = sum(Gamma_2,2);
fda_Sigma_2 = sum(Sigma_2,2);


fda_MClust_1_rand = sum(MClust_1_rand,2);
fda_Trans_1_rand = sum(Trans_1_rand,2);
fda_Assort_1_rand = sum(Assort_1_rand,2);
fda_GEff_1_rand = sum(GEff_1_rand,2);
fda_MLocEff_1_rand = sum(MLocEff_1_rand,2);
fda_Mod_1_rand = sum(Mod_1_rand,2);
fda_ModL_1_rand = sum(ModL_1_rand,2);
fda_PathL_1_rand = sum(PathL_1_rand,2);
fda_MNodeBetw_1_rand = sum(MNodeBetw_1_rand,2);
fda_MEdgeBetw_1_rand = sum(MEdgeBetw_1_rand,2);
fda_Lambda_1_rand = sum(Lambda_1_rand,2);
fda_Gamma_1_rand = sum(Gamma_1_rand,2);
fda_Sigma_1_rand = sum(Sigma_1_rand,2);


fda_MClust_2_rand = sum(MClust_2_rand,2);
fda_Trans_2_rand = sum(Trans_2_rand,2);
fda_Assort_2_rand = sum(Assort_2_rand,2);
fda_GEff_2_rand = sum(GEff_2_rand,2);
fda_MLocEff_2_rand = sum(MLocEff_2_rand,2);
fda_Mod_2_rand = sum(Mod_2_rand,2);
fda_ModL_2_rand = sum(ModL_2_rand,2);
fda_PathL_2_rand = sum(PathL_2_rand,2);
fda_MNodeBetw_2_rand = sum(MNodeBetw_2_rand,2);
fda_MEdgeBetw_2_rand = sum(MEdgeBetw_2_rand,2);
fda_Lambda_2_rand = sum(Lambda_2_rand,2);
fda_Gamma_2_rand = sum(Gamma_2_rand,2);
fda_Sigma_2_rand = sum(Sigma_2_rand,2);


save fda_NetMesPerDens_AveAcrossDens_f  fda_MClust_1_rand fda_Trans_1_rand ...
     fda_Assort_1_rand fda_GEff_1_rand fda_MLocEff_1_rand fda_Mod_1_rand fda_ModL_1_rand ...
     fda_PathL_1_rand fda_MNodeBetw_1_rand fda_MEdgeBetw_1_rand fda_Lambda_1_rand fda_Gamma_1_rand ...
     fda_Sigma_1_rand fda_MClust_2_rand fda_Trans_2_rand ...
     fda_Assort_2_rand fda_GEff_2_rand fda_MLocEff_2_rand fda_Mod_2_rand fda_ModL_2_rand ...
     fda_PathL_2_rand fda_MNodeBetw_2_rand fda_MEdgeBetw_2_rand fda_Lambda_2_rand fda_Gamma_2_rand ...
     fda_Sigma_2_rand fda_MClust_1 fda_Trans_1 fda_Assort_1 fda_GEff_1 fda_MLocEff_1 ...
     fda_Mod_1 fda_ModL_1 fda_PathL_1 fda_MNodeBetw_1 fda_MEdgeBetw_1 fda_Lambda_1 fda_Gamma_1 fda_Sigma_1 ...
     fda_MClust_2 fda_Trans_2 fda_Assort_2 fda_GEff_2 fda_MLocEff_2 fda_Mod_2 fda_ModL_2 fda_PathL_2 fda_MNodeBetw_2 ...
     fda_MEdgeBetw_2 fda_Lambda_2 fda_Gamma_2 fda_Sigma_2


%% 

nRand=size(fda_MClust_1_rand,3);

fda_MClust_1 = mean(fda_MClust_1);
fda_Trans_1 = mean(fda_Trans_1);
fda_Assort_1 = mean(fda_Assort_1);
fda_GEff_1 = mean(fda_GEff_1);
fda_MLocEff_1 = mean(fda_MLocEff_1);
fda_Mod_1 = mean(fda_Mod_1);
fda_ModL_1 = mean(fda_ModL_1);
fda_PathL_1 = mean(fda_PathL_1);
fda_MNodeBetw_1 = mean(fda_MNodeBetw_1);
fda_MEdgeBetw_1 = mean(fda_MEdgeBetw_1);
fda_Lambda_1 = mean(fda_Lambda_1);
fda_Gamma_1 = mean(fda_Gamma_1);
fda_Sigma_1 = mean(fda_Sigma_1);


fda_MClust_2 = mean(fda_MClust_2);
fda_Trans_2 = mean(fda_Trans_2);
fda_Assort_2 = mean(fda_Assort_2);
fda_GEff_2 = mean(fda_GEff_2);
fda_MLocEff_2 = mean(fda_MLocEff_2);
fda_Mod_2 = mean(fda_Mod_2);
fda_ModL_2 = mean(fda_ModL_2);
fda_PathL_2 = mean(fda_PathL_2);
fda_MNodeBetw_2 = mean(fda_MNodeBetw_2);
fda_MEdgeBetw_2 = mean(fda_MEdgeBetw_2);
fda_Lambda_2 = mean(fda_Lambda_2);
fda_Gamma_2 = mean(fda_Gamma_2);
fda_Sigma_2 = mean(fda_Sigma_2);


fda_MClust_1_rand = mean(fda_MClust_1_rand,1);fda_MClust_1_rand = reshape(fda_MClust_1_rand,1,nRand);
fda_Trans_1_rand = mean(fda_Trans_1_rand,1);fda_Trans_1_rand = reshape(fda_Trans_1_rand ,1,nRand);
fda_Assort_1_rand = mean(fda_Assort_1_rand,1);fda_Assort_1_rand = reshape(fda_Assort_1_rand ,1,nRand);
fda_GEff_1_rand = mean(fda_GEff_1_rand,1);fda_GEff_1_rand = reshape(fda_GEff_1_rand ,1,nRand);
fda_MLocEff_1_rand = mean(fda_MLocEff_1_rand,1);fda_MLocEff_1_rand = reshape(fda_MLocEff_1_rand ,1,nRand);
fda_Mod_1_rand = mean(fda_Mod_1_rand,1);fda_Mod_1_rand = reshape(fda_Mod_1_rand ,1,nRand);
fda_ModL_1_rand = mean(fda_ModL_1_rand,1);fda_ModL_1_rand = reshape(fda_ModL_1_rand ,1,nRand);
fda_PathL_1_rand = mean(fda_PathL_1_rand,1);fda_PathL_1_rand = reshape(fda_PathL_1_rand ,1,nRand);
fda_MNodeBetw_1_rand = mean(fda_MNodeBetw_1_rand,1);fda_MNodeBetw_1_rand = reshape(fda_MNodeBetw_1_rand ,1,nRand);
fda_MEdgeBetw_1_rand = mean(fda_MEdgeBetw_1_rand,1);fda_MEdgeBetw_1_rand = reshape(fda_MEdgeBetw_1_rand ,1,nRand);
fda_Lambda_1_rand = mean(fda_Lambda_1_rand,1);fda_Lambda_1_rand = reshape( fda_Lambda_1_rand,1,nRand);
fda_Gamma_1_rand = mean(fda_Gamma_1_rand,1);fda_Gamma_1_rand = reshape(fda_Gamma_1_rand ,1,nRand);
fda_Sigma_1_rand = mean(fda_Sigma_1_rand,1);fda_Sigma_1_rand = reshape(fda_Sigma_1_rand ,1,nRand);


fda_MClust_2_rand = mean(fda_MClust_2_rand,1);fda_MClust_2_rand = reshape(fda_MClust_2_rand,1,nRand);
fda_Trans_2_rand = mean(fda_Trans_2_rand,1);fda_Trans_2_rand = reshape(fda_Trans_2_rand ,1,nRand);
fda_Assort_2_rand = mean(fda_Assort_2_rand,1);fda_Assort_2_rand = reshape(fda_Assort_2_rand ,1,nRand);
fda_GEff_2_rand = mean(fda_GEff_2_rand,1);fda_GEff_2_rand = reshape(fda_GEff_2_rand ,1,nRand);
fda_MLocEff_2_rand = mean(fda_MLocEff_2_rand,1);fda_MLocEff_2_rand = reshape(fda_MLocEff_2_rand ,1,nRand);
fda_Mod_2_rand = mean(fda_Mod_2_rand,1);fda_Mod_2_rand = reshape(fda_Mod_2_rand ,1,nRand);
fda_ModL_2_rand = mean(fda_ModL_2_rand,1);fda_ModL_2_rand = reshape(fda_ModL_2_rand ,1,nRand);
fda_PathL_2_rand = mean(fda_PathL_2_rand,1);fda_PathL_2_rand = reshape(fda_PathL_2_rand ,1,nRand);
fda_MNodeBetw_2_rand = mean(fda_MNodeBetw_2_rand,1);fda_MNodeBetw_2_rand = reshape(fda_MNodeBetw_2_rand ,1,nRand);
fda_MEdgeBetw_2_rand = mean(fda_MEdgeBetw_2_rand,1);fda_MEdgeBetw_2_rand = reshape(fda_MEdgeBetw_2_rand ,1,nRand);
fda_Lambda_2_rand = mean(fda_Lambda_2_rand,1);fda_Lambda_2_rand = reshape( fda_Lambda_2_rand,1,nRand);
fda_Gamma_2_rand = mean(fda_Gamma_2_rand,1);fda_Gamma_2_rand = reshape(fda_Gamma_2_rand ,1,nRand);
fda_Sigma_2_rand = mean(fda_Sigma_2_rand,1);fda_Sigma_2_rand = reshape(fda_Sigma_2_rand ,1,nRand);

save fda_Mean fda_MClust_1_rand fda_Trans_1_rand ...
     fda_Assort_1_rand fda_GEff_1_rand fda_MLocEff_1_rand fda_Mod_1_rand fda_ModL_1_rand ...
     fda_PathL_1_rand fda_MNodeBetw_1_rand fda_MEdgeBetw_1_rand fda_Lambda_1_rand fda_Gamma_1_rand ...
     fda_Sigma_1_rand fda_MClust_2_rand fda_Trans_2_rand ...
     fda_Assort_2_rand fda_GEff_2_rand fda_MLocEff_2_rand fda_Mod_2_rand fda_ModL_2_rand ...
     fda_PathL_2_rand fda_MNodeBetw_2_rand fda_MEdgeBetw_2_rand fda_Lambda_2_rand fda_Gamma_2_rand ...
     fda_Sigma_2_rand fda_MClust_1 fda_Trans_1 fda_Assort_1 fda_GEff_1 fda_MLocEff_1 ...
     fda_Mod_1 fda_ModL_1 fda_PathL_1 fda_MNodeBetw_1 fda_MEdgeBetw_1 fda_Lambda_1 fda_Gamma_1 fda_Sigma_1 ...
     fda_MClust_2 fda_Trans_2 ...
     fda_Assort_2 fda_GEff_2 fda_MLocEff_2 fda_Mod_2 fda_ModL_2 fda_PathL_2 fda_MNodeBetw_2 ...
     fda_MEdgeBetw_2 fda_Lambda_2 fda_Gamma_2 fda_Sigma_2 

%% 

p_fda_MClust = CL_Pval(fda_MClust_2_rand - fda_MClust_1_rand, fda_MClust_2 - fda_MClust_1,'fda_MClust',Tail)
p_fda_Trans = CL_Pval(fda_Trans_2_rand - fda_Trans_1_rand, fda_Trans_2 - fda_Trans_1,'fda_Trans',Tail)
p_fda_Assort = CL_Pval(fda_Assort_2_rand - fda_Assort_1_rand, fda_Assort_2 - fda_Assort_1,'fda_Assort',Tail)
p_fda_GEff = CL_Pval(fda_GEff_2_rand - fda_GEff_1_rand, fda_GEff_2 - fda_GEff_1,'fda_GEff',Tail)
p_fda_MLocEff = CL_Pval(fda_MLocEff_2_rand - fda_MLocEff_1_rand, fda_MLocEff_2 - fda_MLocEff_1,'fda_MLocEff',Tail)
p_fda_Mod = CL_Pval(fda_Mod_2_rand - fda_Mod_1_rand, fda_Mod_2 - fda_Mod_1,'fda_Mod',Tail)
p_fda_ModL = CL_Pval(fda_ModL_2_rand - fda_ModL_1_rand, fda_ModL_2 - fda_ModL_1,'fda_ModL',Tail)
p_fda_PathL = CL_Pval(fda_PathL_2_rand - fda_PathL_1_rand, fda_PathL_2 - fda_PathL_1,'fda_PathL',Tail)
p_fda_MNodeBetw = CL_Pval(fda_MNodeBetw_2_rand - fda_MNodeBetw_1_rand, fda_MNodeBetw_2 - fda_MNodeBetw_1,'fda_MNodeBetw',Tail)
p_fda_MEdgeBetw = CL_Pval(fda_MEdgeBetw_2_rand - fda_MEdgeBetw_1_rand, fda_MEdgeBetw_2 - fda_MEdgeBetw_1,'fda_MEdgeBetw',Tail)
p_fda_Lambda = CL_Pval(fda_Lambda_2_rand - fda_Lambda_1_rand, fda_Lambda_2 - fda_Lambda_1,'fda_Lambda',Tail)
p_fda_Gamma = CL_Pval(fda_Gamma_2_rand - fda_Gamma_1_rand, fda_Gamma_2 - fda_Gamma_1,'fda_Gamma',Tail)
p_fda_Sigma = CL_Pval(fda_Sigma_2_rand - fda_Sigma_1_rand, fda_Sigma_2 - fda_Sigma_1,'fda_Sigma',Tail)


%%
p = Alpha;

CL_fda_Lambda = CL_per(fda_Lambda_2_rand' - fda_Lambda_1_rand', p);
figure; errorbar(1, mean(fda_Lambda_2_rand - fda_Lambda_1_rand), CL_fda_Lambda(1),CL_fda_Lambda(2),'O','LineWidth',1,'MarkerSize',10);
hold;
plot(1,fda_Lambda_2 - fda_Lambda_1,'rv','LineWidth',1.2,'MarkerSize',8);
%set(gca,'YLim',[-0.005 0.005]);

CL_fda_Gamma = CL_per(fda_Gamma_2_rand' - fda_Gamma_1_rand', p);
errorbar(2, mean(fda_Gamma_2_rand - fda_Gamma_1_rand), CL_fda_Gamma(1),CL_fda_Gamma(2),'O','LineWidth',1,'MarkerSize',10);
plot(2,fda_Gamma_2 - fda_Gamma_1,'rv','LineWidth',1.2,'MarkerSize',8);

CL_fda_Sigma = CL_per(fda_Sigma_2_rand' - fda_Sigma_1_rand', p);
errorbar(3, mean(fda_Sigma_2_rand - fda_Sigma_1_rand), CL_fda_Sigma(1),CL_fda_Sigma(2),'O','LineWidth',1,'MarkerSize',10);
plot(3,fda_Sigma_2 - fda_Sigma_1,'rv','LineWidth',1.2,'MarkerSize',8);

set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'path length','clustering','small-worldness'},'FontSize',12)

%plot properties
ylabel('difference in fda','FontSize',12)
title('Difference in fda of Network Measures','fontsize',14,'Color','red')
legend('Null mean and 95% CI',[Group2 ' vs. ' Group1])
grid on

% set(gcf,'WindowStyle','normal')
% set(gcf,'PaperPositionMode', 'auto');
% scrsz = get(0,'ScreenSize');
FIG='fda_Diff_vs_Null';

hgsave(FIG) 
print('-depsc',FIG)
print('-djpeg','-r300',FIG)
print('-dtiff','-r600',FIG)
print('-dpng','-r300',FIG)




