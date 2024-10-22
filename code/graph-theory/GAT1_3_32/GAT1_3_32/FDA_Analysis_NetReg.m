function FDA_Analysis_NetReg(Group1,Group2,NetMes1,NetMes2,NetMes1_rand,NetMes2_rand,Xax,MinDens,MaxDens,DensStep,ROI)


%% 

GUIinput=0;

if nargin<1
    MinDens = input('type the min density of interest: ');
    MaxDens = input('type the max density of interest: ');
    DensStep = input('type the density interval (step): ');
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
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
    
    
    ROI = ff.mat4GAT.roi1;
    Data1 = spm_select(1,['NetMesBin_' Group1 '.mat'],'Select the "NetMes_Bin_Group1.mat" file ("output from Compare Graphs Across Densities")');
    Data2 = spm_select(1,['NetMesBin_' Group2 '.mat'],'Select the "NetMes_Bin_Group2.mat" file ("output from Compare Graphs Across Densities")');
    Data1r = spm_select(1,['NetMesBin_rand_' Group1 '.mat'],'Select the "NetMes_Bin_rand_Group1.mat" file (output from "Compare Graphs Across Densities")');
    Data2r = spm_select(1,['NetMesBin_rand_' Group2 '.mat'],'Select the "NetMes_Bin_rand_Group2.mat" file (output from "Compare Graphs Across Densities")');
    DensInt = spm_select(1,'DensityInterval*','Select the "DensityInterval.mat" (output from "Compare Graphs Across Densities")');
    fprintf('%-4s\n',' loading input data....');
    f1 = load(Data1,'NetMes_B1');NetMes1 = f1.NetMes_B1;
    f2 = load(Data2,'NetMes_B2');NetMes2 = f2.NetMes_B2;
    e1 = load(Data1r,'NetMes_B1_rand');NetMes1_rand = e1.NetMes_B1_rand;
    e2 = load(Data2r,'NetMes_B2_rand');NetMes2_rand = e2.NetMes_B2_rand;
    xx = load(DensInt,'Xax');Xax = xx.Xax;
    
    GUIinput=1;
end
fprintf('%-4s\n',' extracting regional network measures....');


%% 
iXax=[MinDens:DensStep:MaxDens];
nRand = size(NetMes1_rand,2);

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
    
    Ind_max = size(Xax,2);
    
end

%%

MClust1rn = []; MDeg1rn = []; MNodeBetw1rn=[]; MLocEff1rn = [];

MClust2rn = []; MDeg2rn = []; MNodeBetw2rn=[];

cn=0;

for i=Ind_min:Ind_max
    cn=cn+1;
    NM1=NetMes1{i};
    NM2=NetMes2{i};
    
    
    MClust1o(cn,:)=NM1{7,1}';
    MDeg1o(cn,:)=NM1{1,1};
    MNodeBetw1o(cn,:)=NM1{16,1};
    MLocEff1o(cn,:)=NM1{11,1}';
    
    MClust1n(cn,:)=NM1{7,1}'/mean(NM1{7,1});
    MDeg1n(cn,:)=NM1{1,1}/mean(NM1{1,1});
    MNodeBetw1n(cn,:)=NM1{16,1}/mean(NM1{16,1});
    MLocEff1n(cn,:)=NM1{11,1}'/mean(NM1{11,1});
    
    MClust2o(cn,:)=NM2{7,1}';
    MDeg2o(cn,:)=NM2{1,1};
    MNodeBetw2o(cn,:)=NM2{16,1};
    MLocEff2o(cn,:)=NM2{11,1}';
    
    MClust2n(cn,:)=NM2{7,1}'/mean(NM2{7,1});
    MDeg2n(cn,:)=NM2{1,1}/mean(NM2{1,1});
    MNodeBetw2n(cn,:)=NM2{16,1}/mean(NM2{16,1});
    MLocEff2n(cn,:)=NM2{11,1}'/mean(NM2{11,1});
    
    for j=1:size(NetMes1_rand,2)
        
        NM1r=NetMes1_rand{j}{i};
        MClust1rn(cn,:,j)=NM1r{7,1}'/mean(NM1r{7,1});
        MDeg1rn(cn,:,j)=NM1r{1,1}'/mean(NM1r{1,1});
        MNodeBetw1rn(cn,:,j)=NM1r{16,1}'/mean(NM1r{16,1});
        MLocEff1rn(cn,:,j)=NM1r{11,1}'/mean(NM1r{11,1});
        
    end
    
    for j=1:size(NetMes2_rand,2)
    
        NM2r=NetMes2_rand{j}{i};
        MClust2rn(cn,:,j)=NM2r{7,1}'/mean(NM2r{7,1});
        MDeg2rn(cn,:,j)=NM2r{1,1}'/mean(NM2r{1,1});
        MNodeBetw2rn(cn,:,j)=NM2r{16,1}'/mean(NM2r{16,1});
        MLocEff2rn(cn,:,j)=NM2r{11,1}'/mean(NM2r{11,1});
        
    end
    
end

%% 

save(['NetMesReg_Rand_Norm_' Group1],'MClust1rn','MDeg1rn','MNodeBetw1rn','MLocEff1rn');
save(['NetMesReg_Rand_Norm_' Group2],'MClust2rn','MDeg2rn','MNodeBetw2rn','MLocEff2rn');

save(['NetMesReg_Norm' Group1],'MClust1n','MDeg1n','MNodeBetw1n','MLocEff1n');
save(['NetMesReg_Norm' Group2],'MClust2n','MDeg2n','MNodeBetw2n','MLocEff2n');

save(['NetMesReg_' Group1],'MClust1o','MDeg1o','MNodeBetw1o','MLocEff1o');
save(['NetMesReg_' Group2],'MClust2o','MDeg2o','MNodeBetw2o','MLocEff2o');

nROI = size(MClust1n,2);

%% 

fda_MClust1o = sum(MClust1o);
fda_MDeg1o = sum(MDeg1o);
fda_MNodeBetw1o = sum(MNodeBetw1o);
fda_MLocEff1o = sum(MLocEff1o);

fda_MClust1n = sum(MClust1n);
fda_MDeg1n = sum(MDeg1n);
fda_MNodeBetw1n = sum(MNodeBetw1n);
fda_MLocEff1n = sum(MLocEff1n);

fda_MClust2o = sum(MClust2o);
fda_MDeg2o = sum(MDeg2o);
fda_MNodeBetw2o = sum(MNodeBetw2o);
fda_MLocEff2o = sum(MLocEff2o);

fda_MClust2n = sum(MClust2n);
fda_MDeg2n = sum(MDeg2n);
fda_MNodeBetw2n = sum(MNodeBetw2n);
fda_MLocEff2n = sum(MLocEff2n);

fda_MClust1rn = sum(MClust1rn);
fda_MClust1rn = reshape(fda_MClust1rn,nROI,nRand);
fda_MDeg1rn = sum(MDeg1rn);
fda_MDeg1rn = reshape(fda_MDeg1rn,nROI,nRand);
fda_MNodeBetw1rn = sum(MNodeBetw1rn);
fda_MNodeBetw1rn = reshape(fda_MNodeBetw1rn,nROI,nRand);
fda_MLocEff1rn = sum(MLocEff1rn);
fda_MLocEff1rn = reshape(fda_MLocEff1rn,nROI,nRand);

fda_MClust2rn = sum(MClust2rn);
fda_MClust2rn = reshape(fda_MClust2rn,nROI,nRand);
fda_MDeg2rn = sum(MDeg2rn);
fda_MDeg2rn = reshape(fda_MDeg2rn,nROI,nRand);
fda_MNodeBetw2rn = sum(MNodeBetw2rn);
fda_MNodeBetw2rn = reshape(fda_MNodeBetw2rn,nROI,nRand);
fda_MLocEff2rn = sum(MLocEff2rn);
fda_MLocEff2rn = reshape(fda_MLocEff2rn,nROI,nRand);


%
save(['fda_NetMesReg_Norm' Group1],'fda_MClust1n','fda_MDeg1n','fda_MNodeBetw1n','fda_MLocEff1n');
save(['fda_NetMesReg_Norm' Group2],'fda_MClust2n','fda_MDeg2n','fda_MNodeBetw2n','fda_MLocEff2n');

%
save(['fda_NetMesReg_Norm_Rand' Group1],'fda_MClust1rn','fda_MDeg1rn','fda_MNodeBetw1rn','fda_MLocEff1rn');
save(['fda_NetMesReg_Norm_Rand' Group2],'fda_MClust2rn','fda_MDeg2rn','fda_MNodeBetw2rn','fda_MLocEff2rn');


%
save(['fda_NetMesReg_Orig' Group1],'fda_MClust1o','fda_MDeg1o','fda_MNodeBetw1o','fda_MLocEff1o');
save(['fda_NetMesReg_Orig' Group2],'fda_MClust2o','fda_MDeg2o','fda_MNodeBetw2o','fda_MLocEff2o');

%% 

p_fda_MClust = CL_Pval(fda_MClust2rn - fda_MClust1rn, fda_MClust2n' - fda_MClust1n','fda_MClustReg',Tail)
p_fda_MDeg = CL_Pval(fda_MDeg2rn - fda_MDeg1rn, fda_MDeg2n' - fda_MDeg1n','fda_MDegReg',Tail)
p_fda_MNodeBetw = CL_Pval(fda_MNodeBetw2rn - fda_MNodeBetw1rn, fda_MNodeBetw2n' - fda_MNodeBetw1n','fda_MNodeBetwReg',Tail)
p_fda_MLocEff = CL_Pval(fda_MLocEff2rn - fda_MLocEff1rn, fda_MLocEff2n' - fda_MLocEff1n','fda_MLocEffReg',Tail)

%% 

fdr_cor(p_fda_MClust,ROI,'fda_MClust_fdr_pval');
fdr_cor(p_fda_MDeg,ROI,'fda_MDeg_fdr_pval');
fdr_cor(p_fda_MNodeBetw,ROI,'fda_MNodeBetw_fdr_pval');
fdr_cor(p_fda_MLocEff,ROI,'fda_MLocEff_fdr_pval');


%% 

[mu_MClust1_rand,sigma_MClust1_rand,muCi_MClust1_rand,sigmaCi_MClust1_rand]=normfit(fda_MClust1rn',Alpha);
[mu_MDeg1_rand,sigma_MDeg1_rand,muCi_MDeg1_rand,sigmaCi_MDeg1_rand]=normfit(fda_MDeg1rn',Alpha);
[mu_MNodeBetw1_rand,sigma_MNodeBetw1_rand,muCi_MNodeBetw1_rand,sigmaCi_MNodeBetw1_rand]=normfit(fda_MNodeBetw1rn',Alpha);
[mu_MLocEff1_rand,sigma_MLocEff1_rand,muCi_MLocEff1_rand,sigmaCi_MLocEff1_rand]=normfit(fda_MLocEff1rn',Alpha);

[mu_MClust2_rand,sigma_MClust2_rand,muCi_MClust2_rand,sigmaCi_MClust2_rand]=normfit(fda_MClust2rn',Alpha);
[mu_MDeg2_rand,sigma_MDeg2_rand,muCi_MDeg2_rand,sigmaCi_MDeg2_rand]=normfit(fda_MDeg2rn',Alpha);
[mu_MNodeBetw2_rand,sigma_MNodeBetw2_rand,muCi_MNodeBetw2_rand,sigmaCi_MNodeBetw2_rand]=normfit(fda_MNodeBetw2rn',Alpha);
[mu_MLocEff2_rand,sigma_MLocEff2_rand,muCi_MLocEff2_rand,sigmaCi_MLocEff2_rand]=normfit(fda_MLocEff2rn',Alpha);

nROI=size(MClust1n,2);
N_rand=size(NetMes1_rand,2);
Pvalue = Alpha;

Ci_MClust=CL_per(fda_MClust2rn'-fda_MClust1rn',Pvalue);
Ci_MDeg=CL_per(fda_MDeg2rn'-fda_MDeg1rn',Pvalue);
Ci_MNodeBetw=CL_per(fda_MNodeBetw2rn'-fda_MNodeBetw1rn',Pvalue);
Ci_MLocEff=CL_per(fda_MLocEff2rn'-fda_MLocEff1rn',Pvalue);


figure
plot([1:nROI],fda_MClust2n-fda_MClust1n,'r*');hold 
plot([1:nROI],mu_MClust2_rand-mu_MClust1_rand,'bx');
plot([1:nROI],Ci_MClust(1,:)','b--');
plot([1:nROI],Ci_MClust(2,:)','b--');

xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of clustering','fontsize',12,'fontweight','b')
title(' fda of Normalized Clustering ','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('fda_Reg_NormClustering_vs_Null_Nperm.fig')


figure
plot([1:nROI],fda_MDeg2n-fda_MDeg1n,'r*');hold 
plot([1:nROI],mu_MDeg2_rand-mu_MDeg1_rand,'bx');
plot([1:nROI],Ci_MDeg(1,:)','b--');
plot([1:nROI],Ci_MDeg(2,:)','b--');

xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of degree','fontsize',12,'fontweight','b')
title(' fda of Normalized Degree ','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('fda_Reg_NormDegree_vs_Null_Nperm.fig')



figure
plot([1:nROI],fda_MNodeBetw2n-fda_MNodeBetw1n,'r*');hold 
plot([1:nROI],mu_MNodeBetw2_rand-mu_MNodeBetw1_rand,'bx');
plot([1:nROI],Ci_MNodeBetw(1,:)','b--');
plot([1:nROI],Ci_MNodeBetw(2,:)','b--');

xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of Betweenness','fontsize',12,'fontweight','b')
title(' fda of Normalized Betweenness ','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('fda_Reg_NormBetweenness_vs_Null_Nperm.fig')


figure
plot([1:nROI],fda_MLocEff2n-fda_MLocEff1n,'r*');hold 
plot([1:nROI],mu_MLocEff2_rand-mu_MLocEff1_rand,'bx');
plot([1:nROI],Ci_MLocEff(1,:)','b--');
plot([1:nROI],Ci_MLocEff(2,:)','b--');

xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of Local Efficiency','fontsize',12,'fontweight','b')
title(' fda of Normalized Local Efficiency ','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('fda_Reg_NormLocalEfficiency_vs_Null_Nperm.fig')


