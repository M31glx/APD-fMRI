function Net_RegDiff_f_withBootstrap_d(data1,data2,MinMesPlot,MesStepPlot,MaxMesPlot,Nperm,Group1,Group2,varargin)

% data1/2 : group 1/2 functional network measures (output NetMesBin_* from NetMesvsDensity..._Functionals)
% MinMax: mninmum, maximum thresholds and the thresholded step
% Nperm: number of permutation (sampling) (default 100)

% varargin: optional input arguments are those (covariates)for which you need to correct the data:
% Standard format:
% 'var1': an array with the same number of rows as data1 whose columns
%       contin different covariates 
% 'var2': an array with the same number of rows as data1 whose columns
%       contin different covariates


%%

%-- Hadi Hosseini (Created Apr 12,2011)
%-- Hadi Hosseini (Updated Apr 21,2011 for GUI)
%-- Hadi Hosseini (Updated Sept 14,2011 for functionals)

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
    ROI = ff.mat4GATd.roi1;
    
    if isfield(ff.mat4GATd,'tail') && isfield(ff.mat4GATd,'alpha')
        
        Alpha = ff.mat4GATd.alpha;
        Tail = ff.mat4GATd.tail;
            
    else
        
        Alpha = .05;
        Tail = 2;
        
    end
    
    MinMesPlot=ff.mat4GATd.MinThr;
    MaxMesPlot=ff.mat4GATd.MaxThr;
    MesStepPlot=ff.mat4GATd.MesStep;
    
    OutputFName='RegionalNetMesBootstrap_Results_D.mat';
    
    cov=ff.mat4GATd.flagCov;
    
    if isequal(cov,1)
        
        Data1 = spm_select(1,['NetMesBin_D_Adjusted_' Group1 '_FinalThrRange.mat'],'Select 1st group residuals (NetMesBin_D_adjusted.mat .mat file)');
        Data2 = spm_select(1,['NetMesBin_D_Adjusted_' Group2 '_FinalThrRange.mat'],'Select 2nd Group residuals (NetMesBin_D_adjusted.mat .mat file)');
        
    elseif isequal(cov,0)
        
        Data1 = spm_select(1,['NetMesBin_D_' Group1 '_FinalThrRange.mat'],'Select 1st group NetMesBin_D.mat file ');
        Data2 = spm_select(1,['NetMesBin_D_' Group2 '_FinalThrRange.mat'],'Select 2nd group NetMesBin_D.mat file ');
        
    end
    
    fprintf('%-4s\n',['loading inputs...']);
    
    f1=load(Data1,'NetMes_Bin');NetMes1=f1.NetMes_Bin;
    f2=load(Data2,'NetMes_Bin');NetMes2=f2.NetMes_Bin;
    
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



%% 

Sz1=size(NetMes1,1);Sz2=size(NetMes2,1);
data=NetMes1;data(Sz1+1:Sz1+Sz2,:)=NetMes2;
RandIndex=randperm(Sz1+Sz2);
Randata(1:Sz1+Sz2,:)=data(RandIndex(1:Sz1+Sz2),:);

%% 

NetMes1_rand=cell(1,Nperm);NetMes2_rand=cell(1,Nperm);

for i=1:Nperm
    
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    
    Samp1=randsample(Sz1+Sz2,Sz1,'true');
    Samp2=randsample(Sz1+Sz2,Sz2,'true');
    NetMes1_rand{i}=Randata(Samp1,:);
    NetMes2_rand{i}=Randata(Samp2,:);
    
end

%% 

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinThr=input('select your final decision on minimum threshold:');MinIdx=find(single(xxx)==single(MinThr));
MaxThr=input('select your final decision on maximum threshold:');MaxIdx=find(single(xxx)==single(MaxThr));
Xax = [MinThr:MesStepPlot:MaxThr];

%%

dd=pwd;     
mkdir('Regional');
cd([dd '/Regional']);


fprintf('%-4s\n',' calculating regional network measures....');

NetMes1=NetMes1(:,MinIdx:MaxIdx);
NetMes2=NetMes2(:,MinIdx:MaxIdx);

MClust_1norm=[];MDeg_1norm=[];MNodeBetw_1norm=[];MLocEff_1norm=[];
MClust_2norm=[];MDeg_2norm=[];MNodeBetw_2norm=[];MLocEff_2norm=[];
AUC_MClust_1norm=[];AUC_MDeg_1norm=[];AUC_MNodeBetw_1norm=[];AUC_MLocEff_1norm=[];
AUC_MClust_2norm=[];AUC_MDeg_2norm=[];AUC_MNodeBetw_2norm=[];AUC_MLocEff_2norm=[];

fda_MClust_1norm=[];fda_MDeg_1norm=[];fda_MNodeBetw_1norm=[];fda_MLocEff_1norm=[];
fda_MClust_2norm=[];fda_MDeg_2norm=[];fda_MNodeBetw_2norm=[];fda_MLocEff_2norm=[];


for i=1:size(NetMes1,1)
    fprintf('%-4s\n',['calculating group1 subject ' num2str(i) ' regional network measures....']);
    
    temp_clust1=[];temp_deg1=[];temp_nodeb1=[];temp_leff1=[];
    for j=1:size(NetMes1,2)
        
        
        temp_clust1=[temp_clust1;NetMes1{i,j}{7,3}'];
        temp_deg1=[temp_deg1;NetMes1{i,j}{1,3}];
        temp_nodeb1=[temp_nodeb1;NetMes1{i,j}{16,3}];
        temp_leff1=[temp_leff1;NetMes1{i,j}{11,3}'];

    end
    
    MClust_1norm=[MClust_1norm;mean(temp_clust1)];
    MDeg_1norm=[MDeg_1norm;mean(temp_deg1)];
    MNodeBetw_1norm=[MNodeBetw_1norm;mean(temp_nodeb1)];
    MLocEff_1norm=[MLocEff_1norm;mean(temp_leff1)];
    
    
    AUC_MClust_1norm=[AUC_MClust_1norm;trapz(Xax,temp_clust1)];
    AUC_MDeg_1norm=[AUC_MDeg_1norm;trapz(Xax,temp_deg1)];
    AUC_MNodeBetw_1norm=[AUC_MNodeBetw_1norm;trapz(Xax,temp_nodeb1)];
    AUC_MLocEff_1norm=[AUC_MLocEff_1norm;trapz(Xax,temp_leff1)];
    
    
    fda_MClust_1norm=[fda_MClust_1norm;sum(temp_clust1)];
    fda_MDeg_1norm=[fda_MDeg_1norm;sum(temp_deg1)];
    fda_MNodeBetw_1norm=[fda_MNodeBetw_1norm;sum(temp_nodeb1)];
    fda_MLocEff_1norm=[fda_MLocEff_1norm;sum(temp_leff1)];

end


save(['Indiv_NetMesReg_' Group1],'MClust_1norm','MDeg_1norm','MNodeBetw_1norm','MLocEff_1norm');
save(['Indiv_AUC_NetMesReg_' Group1],'AUC_MClust_1norm','AUC_MDeg_1norm','AUC_MNodeBetw_1norm','AUC_MLocEff_1norm');
save(['Indiv_fda_NetMesReg_' Group1],'fda_MClust_1norm','fda_MDeg_1norm','fda_MNodeBetw_1norm','fda_MLocEff_1norm');

MClust_1norm=mean(MClust_1norm);
MDeg_1norm=mean(MDeg_1norm);
MNodeBetw_1norm=mean(MNodeBetw_1norm);
MLocEff_1norm=mean(MLocEff_1norm);

AUC_MClust_1norm=mean(AUC_MClust_1norm);
AUC_MDeg_1norm=mean(AUC_MDeg_1norm);
AUC_MNodeBetw_1norm=mean(AUC_MNodeBetw_1norm);
AUC_MLocEff_1norm=mean(AUC_MLocEff_1norm);

fda_MClust_1norm=mean(fda_MClust_1norm);
fda_MDeg_1norm=mean(fda_MDeg_1norm);
fda_MNodeBetw_1norm=mean(fda_MNodeBetw_1norm);
fda_MLocEff_1norm=mean(fda_MLocEff_1norm);


for i=1:size(NetMes2,1)
    
    fprintf('%-4s\n',['calculating group2 subject ' num2str(i) ' regional network measures....']);
    
    temp_clust2=[];temp_deg2=[];temp_nodeb2=[];temp_leff2=[];
    for j=1:size(NetMes2,2)
        
        
        temp_clust2=[temp_clust2;NetMes2{i,j}{7,3}'];
        temp_deg2=[temp_deg2;NetMes2{i,j}{1,3}];
        temp_nodeb2=[temp_nodeb2;NetMes2{i,j}{16,3}];
        temp_leff2=[temp_leff2;NetMes2{i,j}{11,3}'];
        
    end

    
    MClust_2norm=[MClust_2norm;mean(temp_clust2)];
    MDeg_2norm=[MDeg_2norm;mean(temp_deg2)];
    MNodeBetw_2norm=[MNodeBetw_2norm;mean(temp_nodeb2)];
    MLocEff_2norm=[MLocEff_2norm;mean(temp_leff2)];
    
    
    AUC_MClust_2norm=[AUC_MClust_2norm;trapz(Xax,temp_clust2)];
    AUC_MDeg_2norm=[AUC_MDeg_2norm;trapz(Xax,temp_deg2)];
    AUC_MNodeBetw_2norm=[AUC_MNodeBetw_2norm;trapz(Xax,temp_nodeb2)];
    AUC_MLocEff_2norm=[AUC_MLocEff_2norm;trapz(Xax,temp_leff2)];
    
    
    fda_MClust_2norm=[fda_MClust_2norm;sum(temp_clust2)];
    fda_MDeg_2norm=[fda_MDeg_2norm;sum(temp_deg2)];
    fda_MNodeBetw_2norm=[fda_MNodeBetw_2norm;sum(temp_nodeb2)];
    fda_MLocEff_2norm=[fda_MLocEff_2norm;sum(temp_leff2)];
    
end

save(['Indiv_NetMesReg_' Group2],'MClust_2norm','MDeg_2norm','MNodeBetw_2norm','MLocEff_2norm');
save(['Indiv_AUC_NetMesReg_' Group2],'AUC_MClust_2norm','AUC_MDeg_2norm','AUC_MNodeBetw_2norm','AUC_MLocEff_2norm');
save(['Indiv_fda_NetMesReg_' Group2],'fda_MClust_2norm','fda_MDeg_2norm','fda_MNodeBetw_2norm','fda_MLocEff_2norm');

MClust_2norm=mean(MClust_2norm);
MDeg_2norm=mean(MDeg_2norm);
MNodeBetw_2norm=mean(MNodeBetw_2norm);
MLocEff_2norm=mean(MLocEff_2norm);

AUC_MClust_2norm=mean(AUC_MClust_2norm);
AUC_MDeg_2norm=mean(AUC_MDeg_2norm);
AUC_MNodeBetw_2norm=mean(AUC_MNodeBetw_2norm);
AUC_MLocEff_2norm=mean(AUC_MLocEff_2norm);

fda_MClust_2norm=mean(fda_MClust_2norm);
fda_MDeg_2norm=mean(fda_MDeg_2norm);
fda_MNodeBetw_2norm=mean(fda_MNodeBetw_2norm);
fda_MLocEff_2norm=mean(fda_MLocEff_2norm);


%%

MClust_1_randnorm=[];MDeg_1_randnorm=[];MNodeBetw_1_randnorm=[];MLocEff_1_randnorm=[];
MClust_2_randnorm=[];MDeg_2_randnorm=[];MNodeBetw_2_randnorm=[];MLocEff_2_randnorm=[];

AUC_MClust_1_randnorm=[];AUC_MDeg_1_randnorm=[];AUC_MNodeBetw_1_randnorm=[];AUC_MLocEff_1_randnorm=[];
AUC_MClust_2_randnorm=[];AUC_MDeg_2_randnorm=[];AUC_MNodeBetw_2_randnorm=[];AUC_MLocEff_2_randnorm=[];

fda_MClust_1_randnorm=[];fda_MDeg_1_randnorm=[];fda_MNodeBetw_1_randnorm=[];fda_MLocEff_1_randnorm=[];
fda_MClust_2_randnorm=[];fda_MDeg_2_randnorm=[];fda_MNodeBetw_2_randnorm=[];fda_MLocEff_2_randnorm=[];

for k=1:size(NetMes1_rand,2)

    fprintf('%-4s\n',['calculating ranodm net ' num2str(k) ' regional network measures....']);
    temp_rand1=NetMes1_rand{1,k};
    temp_rand2=NetMes2_rand{1,k};
    
    temp_rand1=temp_rand1(:,MinIdx:MaxIdx);
    temp_rand2=temp_rand2(:,MinIdx:MaxIdx);
    
    MClust_1norm_rand=[];MDeg_1norm_rand=[];MNodeBetw_1norm_rand=[];MLocEff_1norm_rand=[];
    MClust_2norm_rand=[];MDeg_2norm_rand=[];MNodeBetw_2norm_rand=[];MLocEff_2norm_rand=[];
    AUC_MClust_1norm_rand=[];AUC_MDeg_1norm_rand=[];AUC_MNodeBetw_1norm_rand=[];AUC_MLocEff_1norm_rand=[];
    AUC_MClust_2norm_rand=[];AUC_MDeg_2norm_rand=[];AUC_MNodeBetw_2norm_rand=[];AUC_MLocEff_2norm_rand=[];
    fda_MClust_1norm_rand=[];fda_MDeg_1norm_rand=[];fda_MNodeBetw_1norm_rand=[];fda_MLocEff_1norm_rand=[];
    fda_MClust_2norm_rand=[];fda_MDeg_2norm_rand=[];fda_MNodeBetw_2norm_rand=[];fda_MLocEff_2norm_rand=[];
    
    
    for i=1:size(temp_rand1,1)
        
        temp_rand_clust1=[];temp_rand_deg1=[];temp_rand_nodeb1=[];temp_rand_leff1=[];
        
        for j=1:size(temp_rand1,2)
            
            
            temp_rand_clust1=[temp_rand_clust1;temp_rand1{i,j}{7,3}'];
            temp_rand_deg1=[temp_rand_deg1;temp_rand1{i,j}{1,3}];
            temp_rand_nodeb1=[temp_rand_nodeb1;temp_rand1{i,j}{16,3}];
            temp_rand_leff1=[temp_rand_leff1;temp_rand1{i,j}{11,3}'];

        end
        
        MClust_1norm_rand=[MClust_1norm_rand;mean(temp_rand_clust1)];
        MDeg_1norm_rand=[MDeg_1norm_rand;mean(temp_rand_deg1)];
        MNodeBetw_1norm_rand=[MNodeBetw_1norm_rand;mean(temp_rand_nodeb1)];
        MLocEff_1norm_rand=[MLocEff_1norm_rand;mean(temp_rand_leff1)];
        
        
        AUC_MClust_1norm_rand=[AUC_MClust_1norm_rand;trapz(Xax,temp_rand_clust1)];
        AUC_MDeg_1norm_rand=[AUC_MDeg_1norm_rand;trapz(Xax,temp_rand_deg1)];
        AUC_MNodeBetw_1norm_rand=[AUC_MNodeBetw_1norm_rand;trapz(Xax,temp_rand_nodeb1)];
        AUC_MLocEff_1norm_rand=[AUC_MLocEff_1norm_rand;trapz(Xax,temp_rand_leff1)];
        
        
        fda_MClust_1norm_rand=[fda_MClust_1norm_rand;sum(temp_rand_clust1)];
        fda_MDeg_1norm_rand=[fda_MDeg_1norm_rand;sum(temp_rand_deg1)];
        fda_MNodeBetw_1norm_rand=[fda_MNodeBetw_1norm_rand;sum(temp_rand_nodeb1)];
        fda_MLocEff_1norm_rand=[fda_MLocEff_1norm_rand;sum(temp_rand_leff1)];

    end
    
    MClust_1_randnorm=[MClust_1_randnorm;mean(MClust_1norm_rand)];
    MDeg_1_randnorm=[MDeg_1_randnorm;mean(MDeg_1norm_rand)];
    MNodeBetw_1_randnorm=[MNodeBetw_1_randnorm;mean(MNodeBetw_1norm_rand)];
    MLocEff_1_randnorm=[MLocEff_1_randnorm;mean(MLocEff_1norm_rand)];
    
    AUC_MClust_1_randnorm=[AUC_MClust_1_randnorm;mean(AUC_MClust_1norm_rand)];
    AUC_MDeg_1_randnorm=[AUC_MDeg_1_randnorm;mean(AUC_MDeg_1norm_rand)];
    AUC_MNodeBetw_1_randnorm=[AUC_MNodeBetw_1_randnorm;mean(AUC_MNodeBetw_1norm_rand)];
    AUC_MLocEff_1_randnorm=[AUC_MLocEff_1_randnorm;mean(AUC_MLocEff_1norm_rand)];
    
    fda_MClust_1_randnorm=[fda_MClust_1_randnorm;mean(fda_MClust_1norm_rand)];
    fda_MDeg_1_randnorm=[fda_MDeg_1_randnorm;mean(fda_MDeg_1norm_rand)];
    fda_MNodeBetw_1_randnorm=[fda_MNodeBetw_1_randnorm;mean(fda_MNodeBetw_1norm_rand)];
    fda_MLocEff_1_randnorm=[fda_MLocEff_1_randnorm;mean(fda_MLocEff_1norm_rand)];
    
    
    for i=1:size(temp_rand2,1)
        temp_rand_clust2=[];temp_rand_deg2=[];temp_rand_nodeb2=[];temp_rand_leff2=[];
        for j=1:size(temp_rand2,2)
            
            temp_rand_clust2=[temp_rand_clust2;temp_rand2{i,j}{7,3}'];
            temp_rand_deg2=[temp_rand_deg2;temp_rand2{i,j}{1,3}];
            temp_rand_nodeb2=[temp_rand_nodeb2;temp_rand2{i,j}{16,3}];
            temp_rand_leff2=[temp_rand_leff2;temp_rand2{i,j}{11,3}'];
        end

        
        MClust_2norm_rand=[MClust_2norm_rand;mean(temp_rand_clust2)];
        MDeg_2norm_rand=[MDeg_2norm_rand;mean(temp_rand_deg2)];
        MNodeBetw_2norm_rand=[MNodeBetw_2norm_rand;mean(temp_rand_nodeb2)];
        MLocEff_2norm_rand=[MLocEff_2norm_rand;mean(temp_rand_leff2)];
        
        AUC_MClust_2norm_rand=[AUC_MClust_2norm_rand;trapz(Xax,temp_rand_clust2)];
        AUC_MDeg_2norm_rand=[AUC_MDeg_2norm_rand;trapz(Xax,temp_rand_deg2)];
        AUC_MNodeBetw_2norm_rand=[AUC_MNodeBetw_2norm_rand;trapz(Xax,temp_rand_nodeb2)];
        AUC_MLocEff_2norm_rand=[AUC_MLocEff_2norm_rand;trapz(Xax,temp_rand_leff2)];
        
        fda_MClust_2norm_rand=[fda_MClust_2norm_rand;sum(temp_rand_clust2)];
        fda_MDeg_2norm_rand=[fda_MDeg_2norm_rand;sum(temp_rand_deg2)];
        fda_MNodeBetw_2norm_rand=[fda_MNodeBetw_2norm_rand;sum(temp_rand_nodeb2)];
        fda_MLocEff_2norm_rand=[fda_MLocEff_2norm_rand;sum(temp_rand_leff2)];
        
    end
    
    MClust_2_randnorm=[MClust_2_randnorm;mean(MClust_2norm_rand)];
    MDeg_2_randnorm=[MDeg_2_randnorm;mean(MDeg_2norm_rand)];
    MNodeBetw_2_randnorm=[MNodeBetw_2_randnorm;mean(MNodeBetw_2norm_rand)];
    MLocEff_2_randnorm=[MLocEff_2_randnorm;mean(MLocEff_2norm_rand)];
    
    AUC_MClust_2_randnorm=[AUC_MClust_2_randnorm;mean(AUC_MClust_2norm_rand)];
    AUC_MDeg_2_randnorm=[AUC_MDeg_2_randnorm;mean(AUC_MDeg_2norm_rand)];
    AUC_MNodeBetw_2_randnorm=[AUC_MNodeBetw_2_randnorm;mean(AUC_MNodeBetw_2norm_rand)];
    AUC_MLocEff_2_randnorm=[AUC_MLocEff_2_randnorm;mean(AUC_MLocEff_2norm_rand)];
    
    fda_MClust_2_randnorm=[fda_MClust_2_randnorm;mean(fda_MClust_2norm_rand)];
    fda_MDeg_2_randnorm=[fda_MDeg_2_randnorm;mean(fda_MDeg_2norm_rand)];
    fda_MNodeBetw_2_randnorm=[fda_MNodeBetw_2_randnorm;mean(fda_MNodeBetw_2norm_rand)];
    fda_MLocEff_2_randnorm=[fda_MLocEff_2_randnorm;mean(fda_MLocEff_2norm_rand)];
    
end

nROI=size(MClust_1norm,2);

%% 

save(['NetMesReg_rand_' Group1],'MClust_1_randnorm','MDeg_1_randnorm','MNodeBetw_1_randnorm','MLocEff_1_randnorm');
save(['NetMesReg_rand_' Group2],'MClust_2_randnorm','MDeg_2_randnorm','MNodeBetw_2_randnorm','MLocEff_2_randnorm');
save(['AUC_NetMesReg_rand_' Group1],'AUC_MClust_1_randnorm','AUC_MDeg_1_randnorm','AUC_MNodeBetw_1_randnorm','AUC_MLocEff_1_randnorm');
save(['AUC_NetMesReg_rand_' Group2],'AUC_MClust_2_randnorm','AUC_MDeg_2_randnorm','AUC_MNodeBetw_2_randnorm','AUC_MLocEff_2_randnorm');
save(['fda_NetMesReg_rand_' Group1],'fda_MClust_1_randnorm','fda_MDeg_1_randnorm','fda_MNodeBetw_1_randnorm','fda_MLocEff_1_randnorm');
save(['fda_NetMesReg_rand_' Group2],'fda_MClust_2_randnorm','fda_MDeg_2_randnorm','fda_MNodeBetw_2_randnorm','fda_MLocEff_2_randnorm');

save(['NetMesReg_' Group1],'MClust_1norm','MDeg_1norm','MNodeBetw_1norm','MLocEff_1norm');
save(['NetMesReg_' Group2],'MClust_2norm','MDeg_2norm','MNodeBetw_2norm','MLocEff_2norm');
save(['AUC_NetMesReg_' Group1],'AUC_MClust_1norm','AUC_MDeg_1norm','AUC_MNodeBetw_1norm','AUC_MLocEff_1norm');
save(['AUC_NetMesReg_' Group2],'AUC_MClust_2norm','AUC_MDeg_2norm','AUC_MNodeBetw_2norm','AUC_MLocEff_2norm');
save(['fda_NetMesReg_' Group1],'fda_MClust_1norm','fda_MDeg_1norm','fda_MNodeBetw_1norm','fda_MLocEff_1norm');
save(['fda_NetMesReg_' Group2],'fda_MClust_2norm','fda_MDeg_2norm','fda_MNodeBetw_2norm','fda_MLocEff_2norm');

[mu_MClust1_randnorm,sigma_MClust1_randnorm,muCi_MClust1_randnorm,sigmaCi_MClust1_randnorm]=normfit(MClust_1_randnorm,0.05);
[mu_MDeg1_randnorm,sigma_MDeg1_randnorm,muCi_MDeg1_randnorm,sigmaCi_MDeg1_randnorm]=normfit(MDeg_1_randnorm,0.05);
[mu_MNodeBetw1_randnorm,sigma_MNodeBetw1_randnorm,muCi_MNodeBetw1_randnorm,sigmaCi_MNodeBetw1_randnorm]=normfit(MNodeBetw_1_randnorm,0.05);
[mu_MLocEff1_randnorm,sigma_MLocEff1_randnorm,muCi_MLocEff1_randnorm,sigmaCi_MLocEff1_randnorm]=normfit(MLocEff_1_randnorm,0.05);

[mu_MClust2_randnorm,sigma_MClust2_randnorm,muCi_MClust2_randnorm,sigmaCi_MClust2_randnorm]=normfit(MClust_2_randnorm,0.05);
[mu_MDeg2_randnorm,sigma_MDeg2_randnorm,muCi_MDeg2_randnorm,sigmaCi_MDeg2_randnorm]=normfit(MDeg_2_randnorm,0.05);
[mu_MNodeBetw2_randnorm,sigma_MNodeBetw2_randnorm,muCi_MNodeBetw2_randnorm,sigmaCi_MNodeBetw2_randnorm]=normfit(MNodeBetw_2_randnorm,0.05);
[mu_MLocEff2_randnorm,sigma_MLocEff2_randnorm,muCi_MLocEff2_randnorm,sigmaCi_MLocEff2_randnorm]=normfit(MLocEff_2_randnorm,0.05);


[AUC_mu_MClust1_randnorm,sigma_MClust1_randnorm,muCi_MClust1_randnorm,sigmaCi_MClust1_randnorm]=normfit(AUC_MClust_1_randnorm,0.05);
[AUC_mu_MDeg1_randnorm,sigma_MDeg1_randnorm,muCi_MDeg1_randnorm,sigmaCi_MDeg1_randnorm]=normfit(AUC_MDeg_1_randnorm,0.05);
[AUC_mu_MNodeBetw1_randnorm,sigma_MNodeBetw1_randnorm,muCi_MNodeBetw1_randnorm,sigmaCi_MNodeBetw1_randnorm]=normfit(AUC_MNodeBetw_1_randnorm,0.05);
[AUC_mu_MLocEff1_randnorm,sigma_MLocEff1_randnorm,muCi_MLocEff1_randnorm,sigmaCi_MLocEff1_randnorm]=normfit(AUC_MLocEff_1_randnorm,0.05);

[AUC_mu_MClust2_randnorm,sigma_MClust2_randnorm,muCi_MClust2_randnorm,sigmaCi_MClust2_randnorm]=normfit(AUC_MClust_2_randnorm,0.05);
[AUC_mu_MDeg2_randnorm,sigma_MDeg2_randnorm,muCi_MDeg2_randnorm,sigmaCi_MDeg2_randnorm]=normfit(AUC_MDeg_2_randnorm,0.05);
[AUC_mu_MNodeBetw2_randnorm,sigma_MNodeBetw2_randnorm,muCi_MNodeBetw2_randnorm,sigmaCi_MNodeBetw2_randnorm]=normfit(AUC_MNodeBetw_2_randnorm,0.05);
[AUC_mu_MLocEff2_randnorm,sigma_MLocEff2_randnorm,muCi_MLocEff2_randnorm,sigmaCi_MLocEff2_randnorm]=normfit(AUC_MLocEff_2_randnorm,0.05);


[fda_mu_MClust1_randnorm,sigma_MClust1_randnorm,muCi_MClust1_randnorm,sigmaCi_MClust1_randnorm]=normfit(fda_MClust_1_randnorm,0.05);
[fda_mu_MDeg1_randnorm,sigma_MDeg1_randnorm,muCi_MDeg1_randnorm,sigmaCi_MDeg1_randnorm]=normfit(fda_MDeg_1_randnorm,0.05);
[fda_mu_MNodeBetw1_randnorm,sigma_MNodeBetw1_randnorm,muCi_MNodeBetw1_randnorm,sigmaCi_MNodeBetw1_randnorm]=normfit(fda_MNodeBetw_1_randnorm,0.05);
[fda_mu_MLocEff1_randnorm,sigma_MLocEff1_randnorm,muCi_MLocEff1_randnorm,sigmaCi_MLocEff1_randnorm]=normfit(fda_MLocEff_1_randnorm,0.05);

[fda_mu_MClust2_randnorm,sigma_MClust2_randnorm,muCi_MClust2_randnorm,sigmaCi_MClust2_randnorm]=normfit(fda_MClust_2_randnorm,0.05);
[fda_mu_MDeg2_randnorm,sigma_MDeg2_randnorm,muCi_MDeg2_randnorm,sigmaCi_MDeg2_randnorm]=normfit(fda_MDeg_2_randnorm,0.05);
[fda_mu_MNodeBetw2_randnorm,sigma_MNodeBetw2_randnorm,muCi_MNodeBetw2_randnorm,sigmaCi_MNodeBetw2_randnorm]=normfit(fda_MNodeBetw_2_randnorm,0.05);
[fda_mu_MLocEff2_randnorm,sigma_MLocEff2_randnorm,muCi_MLocEff2_randnorm,sigmaCi_MLocEff2_randnorm]=normfit(fda_MLocEff_2_randnorm,0.05);


%%

N_rand=size(NetMes1_rand,2);

Pvalue = Alpha;

Ci_MClustnorm=CL_per(MClust_2_randnorm-MClust_1_randnorm,Pvalue);
Ci_MDegnorm=CL_per(MDeg_2_randnorm-MDeg_1_randnorm,Pvalue);
Ci_MNodeBetwnorm=CL_per(MNodeBetw_2_randnorm-MNodeBetw_1_randnorm,Pvalue);
Ci_MLocEffnorm=CL_per(MLocEff_2_randnorm-MLocEff_1_randnorm,Pvalue);


AUC_Ci_MClustnorm=CL_per(AUC_MClust_2_randnorm-AUC_MClust_1_randnorm,Pvalue);
AUC_Ci_MDegnorm=CL_per(AUC_MDeg_2_randnorm-AUC_MDeg_1_randnorm,Pvalue);
AUC_Ci_MNodeBetwnorm=CL_per(AUC_MNodeBetw_2_randnorm-AUC_MNodeBetw_1_randnorm,Pvalue);
AUC_Ci_MLocEffnorm=CL_per(AUC_MLocEff_2_randnorm-AUC_MLocEff_1_randnorm,Pvalue);

fda_Ci_MClustnorm=CL_per(fda_MClust_2_randnorm-fda_MClust_1_randnorm,Pvalue);
fda_Ci_MDegnorm=CL_per(fda_MDeg_2_randnorm-fda_MDeg_1_randnorm,Pvalue);
fda_Ci_MNodeBetwnorm=CL_per(fda_MNodeBetw_2_randnorm-fda_MNodeBetw_1_randnorm,Pvalue);
fda_Ci_MLocEffnorm=CL_per(fda_MLocEff_2_randnorm-fda_MLocEff_1_randnorm,Pvalue);

%%

p_RegClust_norm = CL_Pval((MClust_2_randnorm-MClust_1_randnorm)',(MClust_2norm'-MClust_1norm'),'RegClustNorm',Tail);
p_RegDeg_norm = CL_Pval((MDeg_2_randnorm-MDeg_1_randnorm)',(MDeg_2norm'-MDeg_1norm'),'RegDegNorm',Tail);
p_RegNodeBetw_norm = CL_Pval((MNodeBetw_2_randnorm-MNodeBetw_1_randnorm)',(MNodeBetw_2norm'-MNodeBetw_1norm'),'RegNodeBetwNorm',Tail);
p_RegLocEff_norm = CL_Pval((MLocEff_2_randnorm-MLocEff_1_randnorm)',(MLocEff_2norm'-MLocEff_1norm'),'RegLocEffNorm',Tail);


p_AUC_RegClust_norm = CL_Pval((AUC_MClust_2_randnorm-AUC_MClust_1_randnorm)',(AUC_MClust_2norm'-AUC_MClust_1norm'),'AUC_RegClustNorm',Tail);
p_AUC_RegDeg_norm = CL_Pval((AUC_MDeg_2_randnorm-AUC_MDeg_1_randnorm)',(AUC_MDeg_2norm'-AUC_MDeg_1norm'),'AUC_RegDegNorm',Tail);
p_AUC_RegNodeBetw_norm = CL_Pval((AUC_MNodeBetw_2_randnorm-AUC_MNodeBetw_1_randnorm)',(AUC_MNodeBetw_2norm'-AUC_MNodeBetw_1norm'),'AUC_RegNodeBetwNorm',Tail);
p_AUC_RegLocEff_norm = CL_Pval((AUC_MLocEff_2_randnorm-AUC_MLocEff_1_randnorm)',(AUC_MLocEff_2norm'-AUC_MLocEff_1norm'),'AUC_RegLocEffNorm',Tail);


p_fda_RegClust_norm = CL_Pval((fda_MClust_2_randnorm-fda_MClust_1_randnorm)',(fda_MClust_2norm'-fda_MClust_1norm'),'fda_RegClustNorm',Tail);
p_fda_RegDeg_norm = CL_Pval((fda_MDeg_2_randnorm-fda_MDeg_1_randnorm)',(fda_MDeg_2norm'-fda_MDeg_1norm'),'fda_RegDegNorm',Tail);
p_fda_RegNodeBetw_norm = CL_Pval((fda_MNodeBetw_2_randnorm-fda_MNodeBetw_1_randnorm)',(fda_MNodeBetw_2norm'-fda_MNodeBetw_1norm'),'fda_RegNodeBetwNorm',Tail);
p_fda_RegLocEff_norm = CL_Pval((fda_MLocEff_2_randnorm-fda_MLocEff_1_randnorm)',(fda_MLocEff_2norm'-fda_MLocEff_1norm'),'fda_RegLocEffNorm',Tail);


fdr_cor(p_RegClust_norm,ROI,'RegClustNorm');
fdr_cor(p_RegDeg_norm,ROI,'RegDegNorm');
fdr_cor(p_RegNodeBetw_norm,ROI,'RegNodeBetwNorm');
fdr_cor(p_RegLocEff_norm,ROI,'RegLocEffNorm');

fdr_cor(p_AUC_RegClust_norm,ROI,'AUC_RegClustNorm');
fdr_cor(p_AUC_RegDeg_norm,ROI,'AUC_RegDegNorm');
fdr_cor(p_AUC_RegNodeBetw_norm,ROI,'AUC_RegNodeBetwNorm');
fdr_cor(p_AUC_RegLocEff_norm,ROI,'AUC_RegLocEffNorm');

fdr_cor(p_fda_RegClust_norm,ROI,'fda_RegClustNorm');
fdr_cor(p_fda_RegDeg_norm,ROI,'fda_RegDegNorm');
fdr_cor(p_fda_RegNodeBetw_norm,ROI,'fda_RegNodeBetwNorm');
fdr_cor(p_fda_RegLocEff_norm,ROI,'fda_RegLocEffNorm');


%%

figure
plot([1:nROI],MClust_2norm-MClust_1norm,'r*');hold 
plot([1:nROI],mu_MClust2_randnorm-mu_MClust1_randnorm,'bx');
plot([1:nROI],Ci_MClustnorm(1,:)','b--');
plot([1:nROI],Ci_MClustnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in clustering coefficient','fontsize',12,'fontweight','b')
title('Normalized Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_NormClustering_vs_Null_Nperm.fig')
close(gcf)


figure
plot([1:nROI],MDeg_2norm'-MDeg_1norm','r*');hold 
plot([1:nROI],mu_MDeg2_randnorm-mu_MDeg1_randnorm,'bx');
plot([1:nROI],Ci_MDegnorm(1,:)','b--');
plot([1:nROI],Ci_MDegnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in mean degree','fontsize',12,'fontweight','b')
title('Normalized degree','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_NormDegree_vs_Null_Nperm.fig')
close(gcf)


figure
plot([1:nROI],MNodeBetw_2norm'-MNodeBetw_1norm','r*');hold 
plot([1:nROI],mu_MNodeBetw2_randnorm-mu_MNodeBetw1_randnorm,'bx');
plot([1:nROI],Ci_MNodeBetwnorm(1,:)','b--');
plot([1:nROI],Ci_MNodeBetwnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in Mean Node Betweenness','fontsize',12,'fontweight','b')
title('Nomrmalized Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_NormNodeBetweenness_vs_Null_Nperm.fig')
close(gcf)


figure
plot([1:nROI],MLocEff_2norm'-MLocEff_1norm','r*');hold 
plot([1:nROI],mu_MLocEff2_randnorm-mu_MLocEff1_randnorm,'bx');
plot([1:nROI],Ci_MLocEffnorm(1,:)','b--');
plot([1:nROI],Ci_MLocEffnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in mean local efficiency','fontsize',12,'fontweight','b')
title('Normalized local efficiency','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_NormlocEff_vs_Null_Nperm.fig')
close(gcf)


%% 

figure
plot([1:nROI],AUC_MClust_2norm-AUC_MClust_1norm,'r*');hold 
plot([1:nROI],AUC_mu_MClust2_randnorm-AUC_mu_MClust1_randnorm,'bx');
plot([1:nROI],AUC_Ci_MClustnorm(1,:)','b--');
plot([1:nROI],AUC_Ci_MClustnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in AUC of clustering coefficient','fontsize',12,'fontweight','b')
title('AUC of Normalized Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_AUC_NormClustering_vs_Null_Nperm.fig')



figure
plot([1:nROI],AUC_MDeg_2norm'-AUC_MDeg_1norm','r*');hold 
plot([1:nROI],AUC_mu_MDeg2_randnorm-AUC_mu_MDeg1_randnorm,'bx');
plot([1:nROI],AUC_Ci_MDegnorm(1,:)','b--');
plot([1:nROI],AUC_Ci_MDegnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in AUC of mean degree','fontsize',12,'fontweight','b')
title('AUC of Normalized degree','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_AUC_NormDegree_vs_Null_Nperm.fig')



figure
plot([1:nROI],AUC_MNodeBetw_2norm'-AUC_MNodeBetw_1norm','r*');hold 
plot([1:nROI],AUC_mu_MNodeBetw2_randnorm-AUC_mu_MNodeBetw1_randnorm,'bx');
plot([1:nROI],AUC_Ci_MNodeBetwnorm(1,:)','b--');
plot([1:nROI],AUC_Ci_MNodeBetwnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in AUC of Mean Node Betweenness','fontsize',12,'fontweight','b')
title('AUC of Nomrmalized Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_AUC_NormNodeBetweenness_vs_Null_Nperm.fig')


figure
plot([1:nROI],AUC_MLocEff_2norm'-AUC_MLocEff_1norm','r*');hold 
plot([1:nROI],AUC_mu_MLocEff2_randnorm-AUC_mu_MLocEff1_randnorm,'bx');
plot([1:nROI],AUC_Ci_MLocEffnorm(1,:)','b--');
plot([1:nROI],AUC_Ci_MLocEffnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in AUC of mean local efficiency','fontsize',12,'fontweight','b')
title('AUC of Normalized local efficiency','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_AUC_NormlocEff_vs_Null_Nperm.fig')


%%

figure
plot([1:nROI],fda_MClust_2norm-fda_MClust_1norm,'r*');hold 
plot([1:nROI],fda_mu_MClust2_randnorm-fda_mu_MClust1_randnorm,'bx');
plot([1:nROI],fda_Ci_MClustnorm(1,:)','b--');
plot([1:nROI],fda_Ci_MClustnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of clustering coefficient','fontsize',12,'fontweight','b')
title('fda of Normalized Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_fda_NormClustering_vs_Null_Nperm.fig')



figure
plot([1:nROI],fda_MDeg_2norm'-fda_MDeg_1norm','r*');hold 
plot([1:nROI],fda_mu_MDeg2_randnorm-fda_mu_MDeg1_randnorm,'bx');
plot([1:nROI],fda_Ci_MDegnorm(1,:)','b--');
plot([1:nROI],fda_Ci_MDegnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of mean degree','fontsize',12,'fontweight','b')
title('fda of Normalized degree','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_fda_NormDegree_vs_Null_Nperm.fig')


figure
plot([1:nROI],fda_MNodeBetw_2norm'-fda_MNodeBetw_1norm','r*');hold 
plot([1:nROI],fda_mu_MNodeBetw2_randnorm-fda_mu_MNodeBetw1_randnorm,'bx');
plot([1:nROI],fda_Ci_MNodeBetwnorm(1,:)','b--');
plot([1:nROI],fda_Ci_MNodeBetwnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of Mean Node Betweenness','fontsize',12,'fontweight','b')
title('fda of Nomrmalized Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_fda_NormNodeBetweenness_vs_Null_Nperm.fig')


figure
plot([1:nROI],fda_MLocEff_2norm'-fda_MLocEff_1norm','r*');hold 
plot([1:nROI],fda_mu_MLocEff2_randnorm-fda_mu_MLocEff1_randnorm,'bx');
plot([1:nROI],fda_Ci_MLocEffnorm(1,:)','b--');
plot([1:nROI],fda_Ci_MLocEffnorm(2,:)','b--');
xlabel('ROI','fontsize',12,'fontweight','b')
ylabel('Difference in fda of mean local efficiency','fontsize',12,'fontweight','b')
title('fda of Normalized local efficiency','fontsize',14,'fontweight','b','fontangle','italic')
legend([Group2 'vs.' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
grid on
hgsave('Reg_fda_NormlocEff_vs_Null_Nperm.fig')



fprintf('%-4s\n','.... done ....');

