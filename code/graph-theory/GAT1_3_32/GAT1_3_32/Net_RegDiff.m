function Net_RegDiff(NetMes1,NetMes2,NetMes1_rand,NetMes2_rand,DensInd,ROI)

% NetMes1/2: the output (NetMes_B1/B2) network measures from NetMesvsDensity... 
% NetMes1_rand/NetMes2_rand: the output (NetMes_B1_rand/NetMes_B2_rand) network measures from NetMesvsDensity...
%DensInd: The index corresponding to minimum density for fully connected graph (output from KAN_FixedDensity) 

%%
%-- Hadi Hosseini (Created Apr 12)
%-- Hadi Hosseini (Updated Apr 21 for GUI)


%%
GUIinput=0;

if nargin<1
  
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;
    Group2 = ff.mat4GAT.g2;
    ROI = ff.mat4GAT.roi1;
    
    if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha')
        
        Alpha = ff.mat4GAT.alpha;
        Tail = ff.mat4GAT.tail;
        
    else
        
        Alpha = .05;
        Tail = 2;
        
    end
    
    
    Data1 = spm_select(1,['NetMesBin_' Group1 '.mat'],'Select the "NetMes_Bin_Group1.mat" file ("output from Compare Graphs Across Densities")');
    Data2 = spm_select(1,['NetMesBin_' Group2 '.mat'],'Select the "NetMes_Bin_Group2.mat" file ("output from Compare Graphs Across Densities")');
    Data1r = spm_select(1,['NetMesBin_rand_' Group1 '.mat'],'Select the "NetMes_Bin_rand_Group1.mat" file (output from "Compare Graphs Across Densities")');
    Data2r = spm_select(1,['NetMesBin_rand_' Group2 '.mat'],'Select the "NetMes_Bin_rand_Group2.mat" file (output from "Compare Graphs Across Densities")');
    DensInp= spm_select(1,'KAN_FixedDens*','Select the KAN_FixedDens_Results.mat file (output from "Compare Graphs at Minimum Density")');
    
    fprintf('%-4s\n',' loading input data....');
    
    f1=load(Data1,'NetMes_B1');NetMes1=f1.NetMes_B1;
    f2=load(Data2,'NetMes_B2');NetMes2=f2.NetMes_B2;
    e1=load(Data1r,'NetMes_B1_rand');NetMes1_rand=e1.NetMes_B1_rand;
    e2=load(Data2r,'NetMes_B2_rand');NetMes2_rand=e2.NetMes_B2_rand;
    e3=load(DensInp,'NetMes_Bin1');DensMin=e3.NetMes_Bin1;
    
    MinDens=DensMin{6,1};
    
    for i=1:size(NetMes1,2)
       
        densit(i)=NetMes1{i}{6,1};
        Err(i)=abs(densit(i)-MinDens);
        
    end
    
    [Min,IMin]=min(Err);DensInd=IMin;
    
    GUIinput=1;
    
end

%%

fprintf('%-4s\n',' calculating regional network measures....');

NetMes_1=NetMes1{DensInd};
NetMes_2=NetMes2{DensInd};

for i=1:size(NetMes1_rand,2)

    NetMes_1_rand{i}=NetMes1_rand{i}{DensInd};
    
end

for i=1:size(NetMes2_rand,2)
    
    NetMes_2_rand{i}=NetMes2_rand{i}{DensInd};
    
end

nROI=size(NetMes_1{1},2);

%% 

MClust_1_rand=[];MDeg_1_rand=[];MNodeBetw_1_rand=[];
MClust_1_randnorm=[];MDeg_1_randnorm=[];MNodeBetw_1_randnorm=[];

MClust_2_rand=[];MDeg_2_rand=[];MNodeBetw_2_rand=[];
MClust_2_randnorm=[];MDeg_2_randnorm=[];MNodeBetw_2_randnorm=[];

MClust_1=NetMes_1{7}';MDeg_1=NetMes_1{1};
MNodeBetw_1=NetMes_1{16};

MClust_2=NetMes_2{7}';MDeg_2=NetMes_2{1};
MNodeBetw_2=NetMes_2{16};

MClust_1norm=MClust_1/mean(MClust_1);
MDeg_1norm=MDeg_1/mean(MDeg_1);
MNodeBetw_1norm=MNodeBetw_1/mean(MNodeBetw_1);

MClust_2norm=MClust_2/mean(MClust_2);
MDeg_2norm=MDeg_2/mean(MDeg_2);
MNodeBetw_2norm=MNodeBetw_2/mean(MNodeBetw_2);


%%

for i=1:size(NetMes_1_rand,2)

    MClust_1_rand=[MClust_1_rand;NetMes_1_rand{i}{7}'];
    MDeg_1_rand=[MDeg_1_rand;NetMes_1_rand{i}{1}];
    MNodeBetw_1_rand=[MNodeBetw_1_rand;NetMes_1_rand{i}{16}];
    
    MClust_2_rand=[MClust_2_rand;NetMes_2_rand{i}{7}'];
    MDeg_2_rand=[MDeg_2_rand;NetMes_2_rand{i}{1}];
    MNodeBetw_2_rand=[MNodeBetw_2_rand;NetMes_2_rand{i}{16}];
    
    MClust_1_randnorm=[MClust_1_randnorm;MClust_1_rand(i,:)/mean(MClust_1_rand(i,:))];
    MDeg_1_randnorm=[MDeg_1_randnorm;MDeg_1_rand(i,:)/mean(MDeg_1_rand(i,:))];
    MNodeBetw_1_randnorm=[MNodeBetw_1_randnorm;MNodeBetw_1_rand(i,:)/mean(MNodeBetw_1_rand(i,:))];    
    
    MClust_2_randnorm=[MClust_2_randnorm;MClust_2_rand(i,:)/mean(MClust_2_rand(i,:))];
    MDeg_2_randnorm=[MDeg_2_randnorm;MDeg_2_rand(i,:)/mean(MDeg_2_rand(i,:))];
    MNodeBetw_2_randnorm=[MNodeBetw_2_randnorm;MNodeBetw_2_rand(i,:)/mean(MNodeBetw_2_rand(i,:))];    
       
end


save(['NetMesReg_rand_' Group1],'MClust_1_rand','MDeg_1_rand','MNodeBetw_1_rand','MClust_1_randnorm','MDeg_1_randnorm','MNodeBetw_1_randnorm');
save(['NetMesReg_rand_' Group2],'MClust_2_rand','MDeg_2_rand','MNodeBetw_2_rand','MClust_2_randnorm','MDeg_2_randnorm','MNodeBetw_2_randnorm');

save(['NetMesReg_' Group1],'MClust_1','MDeg_1','MNodeBetw_1','MClust_1norm','MDeg_1norm','MNodeBetw_1norm');
save(['NetMesReg_' Group2],'MClust_2','MDeg_2','MNodeBetw_2','MClust_2norm','MDeg_2norm','MNodeBetw_2norm');

%%

[mu_MClust1_rand,sigma_MClust1_rand,muCi_MClust1_rand,sigmaCi_MClust1_rand]=normfit(MClust_1_rand,0.05);
[mu_MDeg1_rand,sigma_MDeg1_rand,muCi_MDeg1_rand,sigmaCi_MDeg1_rand]=normfit(MDeg_1_rand,0.05);
[mu_MNodeBetw1_rand,sigma_MNodeBetw1_rand,muCi_MNodeBetw1_rand,sigmaCi_MNodeBetw1_rand]=normfit(MNodeBetw_1_rand,0.05);
%
[mu_MClust2_rand,sigma_MClust2_rand,muCi_MClust2_rand,sigmaCi_MClust2_rand]=normfit(MClust_2_rand,0.05);
[mu_MDeg2_rand,sigma_MDeg2_rand,muCi_MDeg2_rand,sigmaCi_MDeg2_rand]=normfit(MDeg_2_rand,0.05);
[mu_MNodeBetw2_rand,sigma_MNodeBetw2_rand,muCi_MNodeBetw2_rand,sigmaCi_MNodeBetw2_rand]=normfit(MNodeBetw_2_rand,0.05);

[mu_MClust1_randnorm,sigma_MClust1_randnorm,muCi_MClust1_randnorm,sigmaCi_MClust1_randnorm]=normfit(MClust_1_randnorm,0.05);
[mu_MDeg1_randnorm,sigma_MDeg1_randnorm,muCi_MDeg1_randnorm,sigmaCi_MDeg1_randnorm]=normfit(MDeg_1_randnorm,0.05);
[mu_MNodeBetw1_randnorm,sigma_MNodeBetw1_randnorm,muCi_MNodeBetw1_randnorm,sigmaCi_MNodeBetw1_randnorm]=normfit(MNodeBetw_1_randnorm,0.05);

[mu_MClust2_randnorm,sigma_MClust2_randnorm,muCi_MClust2_randnorm,sigmaCi_MClust2_randnorm]=normfit(MClust_2_randnorm,0.05);
[mu_MDeg2_randnorm,sigma_MDeg2_randnorm,muCi_MDeg2_randnorm,sigmaCi_MDeg2_randnorm]=normfit(MDeg_2_randnorm,0.05);
[mu_MNodeBetw2_randnorm,sigma_MNodeBetw2_randnorm,muCi_MNodeBetw2_randnorm,sigmaCi_MNodeBetw2_randnorm]=normfit(MNodeBetw_2_randnorm,0.05);


N_rand=size(NetMes1_rand,2);
Pvalue = Alpha;


Ci_MClust=CL_per(MClust_2_rand-MClust_1_rand,Pvalue);
Ci_MDeg=CL_per(MDeg_2_rand-MDeg_1_rand,Pvalue);
Ci_MNodeBetw=CL_per(MNodeBetw_2_rand-MNodeBetw_1_rand,Pvalue);

Ci_MClustnorm=CL_per(MClust_2_randnorm-MClust_1_randnorm,Pvalue);
Ci_MDegnorm=CL_per(MDeg_2_randnorm-MDeg_1_randnorm,Pvalue);
Ci_MNodeBetwnorm=CL_per(MNodeBetw_2_randnorm-MNodeBetw_1_randnorm,Pvalue);


%%

p_RegDeg = CL_Pval((MDeg_2_rand-MDeg_1_rand)',(MDeg_2'-MDeg_1'),'RegDeg',Tail);
p_RegNodeBetw = CL_Pval((MNodeBetw_2_rand-MNodeBetw_1_rand)',(MNodeBetw_2'-MNodeBetw_1'),'RegNodeBetw',Tail);
p_RegClust = CL_Pval((MClust_2_rand-MClust_1_rand)',(MClust_2'-MClust_1'),'RegClust',Tail);

p_RegDeg_norm = CL_Pval((MDeg_2_randnorm-MDeg_1_randnorm)',(MDeg_2norm'-MDeg_1norm'),'RegDegNorm',Tail);
p_RegNodeBetw_norm = CL_Pval((MNodeBetw_2_randnorm-MNodeBetw_1_randnorm)',(MNodeBetw_2norm'-MNodeBetw_1norm'),'RegNodeBetwNorm',Tail);
p_RegClust_norm = CL_Pval((MClust_2_randnorm-MClust_1_randnorm)',(MClust_2norm'-MClust_1norm'),'RegClustNorm',Tail);


fdr_cor(p_RegDeg,ROI,'RegDeg');
fdr_cor(p_RegNodeBetw,ROI,'RegNodeBetw');
fdr_cor(p_RegClust,ROI,'RegClust');

fdr_cor(p_RegDeg_norm,ROI,'RegDegNorm');
fdr_cor(p_RegNodeBetw_norm,ROI,'RegNodeBetwNorm');
fdr_cor(p_RegClust_norm,ROI,'RegClustNorm');



%%

%clustering
figure
plot([1:nROI],MClust_2norm-MClust_1norm,'r*');hold 
plot([1:nROI],mu_MClust2_randnorm-mu_MClust1_randnorm,'bx');%mean random difference
plot([1:nROI],Ci_MClustnorm(1,:)','b--');%mean randnormom difference
plot([1:nROI],Ci_MClustnorm(2,:)','b--');%mean randnormom difference
%plot properties
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in clustering coefficient','fontsize',12,'fontweight','b')
title('Normalized Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_NormClustering_vs_Null_Nperm.fig')


%Degree
figure
plot([1:nROI],MDeg_2norm'-MDeg_1norm','r*');hold 
plot([1:nROI],mu_MDeg2_randnorm-mu_MDeg1_randnorm,'bx');%mean randnormom difference
plot([1:nROI],Ci_MDegnorm(1,:)','b--');%mean randnormom difference
plot([1:nROI],Ci_MDegnorm(2,:)','b--');%mean randnormom difference
%plot properties
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in mean degree','fontsize',12,'fontweight','b')
title('Normalized degree','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_NormDegree_vs_Null_Nperm.fig')


%Mean Node Betweenness
figure
plot([1:nROI],MNodeBetw_2norm'-MNodeBetw_1norm','r*');hold 
plot([1:nROI],mu_MNodeBetw2_randnorm-mu_MNodeBetw1_randnorm,'bx');%mean randnormom difference
plot([1:nROI],Ci_MNodeBetwnorm(1,:)','b--');%mean randnormom difference
plot([1:nROI],Ci_MNodeBetwnorm(2,:)','b--');%mean randnormom difference
%plot properties
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in Mean Node Betweenness','fontsize',12,'fontweight','b')
title('Nomrmalized Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_NormNodeBetweenness_vs_Null_Nperm.fig')



fprintf('%-4s\n','.... done ....');
