function NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_D%(Input1,Input2,rand_Input1,rand_Input2,BiggerIsBetter,Measure,MinMes,MaxMes) 

% The program accepts the Results_ROI (output from conn) of both groups
% combined (note that the groups should be in order e.g. 1:34 BC; 35:61 CON
%'MinMes': The minimum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.05 for density)
%'MaxMes': The maximum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.6 for density)


%%
%-- Hadi Hosseini March 22, 2011
%-- Hadi Hosseini August 18, 2011 (updated for functional connectivity)
%-- Hadi Hosseini August 18, 2011 (updated for covariates of no interest)


%% 

clear all;
Precision=10;
Step=0.0001;

GroupName=input('please type name of the Data set in single quotation:  ');    
MinMes=input('please type the lowest density of interest (e.g. 0.01):  ');
MaxMes=input('please type the highest density of interest (e.g. 0.6) :  ');
MesStep=input('please type the density interval (setps) of interest (e.g. 0.01):  ');

%%

cov=input('any covariates: 1 for yes , 2 for no: ');

if isequal(cov,1)
       
    CovModel = input('modeling covariates with(type 1) or wo(type 2) interactions: ');
        
    if isequal(CovModel,1)

        mat4GATd.CovModel = 'interaction';
        CovModel = 'interaction';

    else

        mat4GATd.CovModel = 'linear';
        CovModel = 'linear';

    end
        
    CovType = input('XLS (type 1) or Mat (type 2) file? ');
    
    if isequal(CovType,1)
    
        [CovName,CovPathName] = uigetfile('*.xls',['Select the Covariates XLS file for Group1']);
        Covar1 = xlsread([CovPathName CovName]);
        [CovName,CovPathName] = uigetfile('*.xls',['Select the Covariates XLS file for Group2']);
        Covar2 = xlsread([CovPathName CovName]);
        
        
    elseif isequal(CovType,2)
        
        Covar1 = spm_select(1,'.mat','Select 1st Group covariate .mat file (nSubjct*nCovar)');
        Covar2 = spm_select(1,'.mat','Select 2nd Group covariate .mat file (nSubjct*nCovar)');
        Covar1=load(Covar1,'covar');Covar1=Covar1.covar;
        Covar2=load(Covar2,'covar');Covar2=Covar2.covar;
        
    end
    
    Covar=[Covar1;Covar2];
    mat4GATd.cov1=Covar1;
    mat4GATd.cov2=Covar2;
    mat4GATd.flagCov=1;
    
else
    
    mat4GATd.flagCov=0;
    
end

%% 

pathname=pwd;
[filename, pathname] = uigetfile('results_diff*.mat', 'Select results_diff*.mat file');
filename=fullfile(pathname,filename);
fprintf('%-4s\n',' loading the input data....');

%% 
load(filename,'Z','names')
[nROI,nROI2]=size(Z);

[rois,ok]=listdlg('PromptString','Select ROIs:',...
        'SelectionMode','multiple',...
        'ListString',char(names),...
        'initialvalue',1:nROI);
if ~ok,return;end

names={names{rois}};Z=Z(rois,rois,:);nROI=length(rois);

%% 

%null_flag = input('type 2 (same deg dist), 3 (H-Q-S), 4 (RandCorr): '); 
null_flag = 2;%at this moment, only upports same deg dist
mat4GATd.null = null_flag;
NullType = null_flag;

Null_Iter = input('number of null networks (e.g. 20):  ');
mat4GATd.nullIter = Null_Iter;

Tail = input('one (type 1) or two (type 2) -tailed test:  ');
mat4GATd.tail = Tail;

Alpha = input('significance threshold level (e.g. 0.05):  ');
mat4GATd.alpha = Alpha;



%% 

MinMesPlot=MinMes;
MaxMesPlot=MaxMes;
MesStepPlot=MesStep;

%%

N = size(Z,1);
MinEdge=MinMes*((N^2-N)/2);MinMes=floor(MinEdge);
MaxEdge=MaxMes*((N^2-N)/2);MaxMes=ceil(MaxEdge);
EdgeStep=MesStep*((N^2-N)/2);MesStep=round(EdgeStep);


%% 

ZZ=Z;


%% 
Sz=size(ZZ,1);

for i=1:size(ZZ,3)

    z=ZZ(:,:,i);
    z(1:Sz+1:end)=0;z(z<0)=0;
    %normalize by mean for thresholding
    %z = z ./ max(max(z));
    
    ZZ(:,:,i)=z;
    
end

save(['UnThresh_CorrMatrices_D_' GroupName],'ZZ');

%% 

Step = max(max(max(ZZ)))./10000;


fprintf('%-4s\n',' calculating threshold level. it might take tens of minutes depending on the number of networks....');

Thresh_Density=zeros(size([Step:Step:max(max(max(ZZ)))],2),2,size(ZZ,3));

for n=1:size(ZZ,3)

    C=ZZ(:,:,n);
    cn=0;

    for i=Step:Step:max(max(C))

        cn=cn+1;
        C(C<=i)=0;
        Thresh_Density(cn,:,n)=[i nnz(triu(C))];

    end
    
end

save(['ThreshDens_D_' GroupName] ,'Thresh_Density');

count=0;

for u=MinMes:MesStep:MaxMes
    
    count=count+1;
    
    fprintf('%-4s\n',[' calculating network measures at each density: ' num2str(count) ' out of ' num2str(ceil((MaxMes-MinMes)./MesStep)) '...']);
    
    for n=1:size(Thresh_Density,3)
        
        T_Density=Thresh_Density(:,:,n);
        
        Ind=find(T_Density(:,2) >= u-Precision & T_Density(:,2) <= u+Precision );
        Dind=T_Density(Ind,2)-u;
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
        if isempty(Ind)
            
            fprintf('%-4s\n',[' No same density within the default precision: decreasing the precision for subject#' num2str(n)]);
            Ind=find(T_Density(:,2) >= u-100*Precision & T_Density(:,2) <= u+100*Precision );
            Dind=Thresh_Density(Ind,2)-u;
            Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
            Ind=Ind(min(Ind_Ind));
            
        end
        
        Thresh=T_Density(Ind,1);
        R=ZZ(:,:,n);R(R<=Thresh)=0;
        
        %if isequal(NullType,3) || isequal(NullType,4)
                
            %NetMes_Bin = NetMeasures_Binary(R,NullType,Null_Iter,d,Inp);

        %else

            NetMes = NetMeasures_Binary(R,NullType,Null_Iter);

        %end 
        
        
        NetMes{23,1}=Thresh;
        NetMes{24,1}=R;
        
        NetMes_Bin{n,count}=NetMes;
        
    end
    
end

save(['NetMesBin_D_' GroupName],'NetMes_Bin');

MinThr = MinMesPlot;MinIdx = 1;
MaxThr = MaxMesPlot;MaxIdx =  size(NetMes_Bin,2);


Group1=input('type name of the 1st group in single quotations: ');
Group2=input('type name of the 2nd group in single quotations: ');
n1=input('type number of subjects in the 1st group: ');
n2=input('type number of subjects in the 2nd group: ');
NetMes_Bin=NetMes_Bin(1:n1+n2,MinIdx:MaxIdx);
save MinMaxThr MinThr MaxThr MesStepPlot;
save(['NetMesBin_D_' GroupName '_FinalThrRange'],'NetMes_Bin','-v7.3');

mat4GATd.g1 = Group1;mat4GATd.g2 = Group2;
mat4GATd.n1 = n1;mat4GATd.n2 = n2;
mat4GATd.MinThr = MinThr;
mat4GATd.MaxThr = MaxThr;
mat4GATd.MesStep = MesStepPlot;
mat4GATd.roi1 = names;
mat4GATd.roi2 = names;
save(['mat4GATd_' Group1 '_' Group2],'mat4GATd');


Spare=NetMes_Bin;
NetMes_Bin=Spare(1:n1,:);save(['NetMesBin_D_' Group1 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin=Spare(n1+1:n1+n2,:);save(['NetMesBin_D_' Group2 '_FinalThrRange'],'NetMes_Bin','-v7.3');
NetMes_Bin=Spare;

    
%% 

N=NetMes_Bin(1:n1+n2,:);

if isequal(cov,1) && ~isempty(Covar) 
    
    ExcludeNans=find(isnan(Covar)==1);
    Included=setdiff(1:size(N,1),ExcludeNans);
    N=N(Included,:);
    Covar=Covar(Included,:);
    
    for thresh=1:size(N,2)
        
        deg=[];degN=[];degM=[];str=[];strN=[];strM=[];assort=[];dens=[];
        clust=[];clustN=[];clustM=[];trans=[];geff=[];leff=[];leffN=[];leffM=[];
        mod=[];modl=[];charp=[];nodebetw=[];nodebetwN=[];nodebetwM=[];
        edgebetw=[];edgebetwN=[];edgebetwM=[];lam=[];gam=[];sig=[];
        
        for subj=1:size(N,1)
        
            deg=[deg;N{subj,thresh}{1,1}];
            degN=[degN;N{subj,thresh}{1,3}];
            degM=[degM;N{subj,thresh}{2,1}];
            str=[str;N{subj,thresh}{3,1}];
            strN=[strN;N{subj,thresh}{3,3}];
            strM=[strM;N{subj,thresh}{4,1}];
            assort=[assort;N{subj,thresh}{5,1}];
            dens=[dens;N{subj,thresh}{6,1}];
            clust=[clust;transpose(N{subj,thresh}{7,1})];
            clustN=[clustN;transpose(N{subj,thresh}{7,3})];
            clustM=[clustM;N{subj,thresh}{8,1}];
            trans=[trans;N{subj,thresh}{9,1}];
            geff=[geff;N{subj,thresh}{10,1}];
            leff=[leff;transpose(N{subj,thresh}{11,1})];
            leffN=[leffN;transpose(N{subj,thresh}{11,3})];
            leffM=[leffM;N{subj,thresh}{12,1}];
            mod=[mod;N{subj,thresh}{13,1}];
            modl=[modl;N{subj,thresh}{14,1}];
            charp=[charp;N{subj,thresh}{15,1}];
            nodebetw=[nodebetw;N{subj,thresh}{16,1}];
            nodebetwN=[nodebetwN;N{subj,thresh}{16,3}];
            nodebetwM=[nodebetwM;N{subj,thresh}{17,1}];
            temp=N{subj,thresh}{18,1};edgebetw=[edgebetw;reshape(temp,1,size(temp,1).^2)];
            temp=N{subj,thresh}{18,3};edgebetwN=[edgebetwN;reshape(temp,1,size(temp,1).^2)];
            edgebetwM=[edgebetwM;N{subj,thresh}{19,1}];
            lam=[lam;N{subj,thresh}{20,1}];
            gam=[gam;N{subj,thresh}{21,1}];
            sig=[sig;N{subj,thresh}{22,1}];
            
        end
        
        for nroi=1:size(deg,2)
        
            stats=regstats(deg(:,nroi),Covar,CovModel);deg_res(:,nroi)=mean(deg(:,nroi))+stats.r;
            stats=regstats(degN(:,nroi),Covar,CovModel);degN_res(:,nroi)=mean(degN(:,nroi))+stats.r;
            stats=regstats(str(:,nroi),Covar,CovModel);str_res(:,nroi)=mean(str(:,nroi))+stats.r;
            stats=regstats(strN(:,nroi),Covar,CovModel);strN_res(:,nroi)=mean(strN(:,nroi))+stats.r;
            stats=regstats(clust(:,nroi),Covar,CovModel);clust_res(:,nroi)=mean(clust(:,nroi))+stats.r;
            stats=regstats(clustN(:,nroi),Covar,CovModel);clustN_res(:,nroi)=mean(clustN(:,nroi))+stats.r;
            stats=regstats(leff(:,nroi),Covar,CovModel);leff_res(:,nroi)=mean(leff(:,nroi))+stats.r;
            stats=regstats(leffN(:,nroi),Covar,CovModel);leffN_res(:,nroi)=mean(leffN(:,nroi))+stats.r;
            stats=regstats(nodebetw(:,nroi),Covar,CovModel);nodebetw_res(:,nroi)=mean(nodebetw(:,nroi))+stats.r;
            stats=regstats(nodebetwN(:,nroi),Covar,CovModel);nodebetwN_res(:,nroi)=mean(nodebetwN(:,nroi))+stats.r;
            
            if isequal(nroi,1)
                
                stats=regstats(degM,Covar,CovModel);degM_res(:,nroi)=mean(degM)+stats.r;
                stats=regstats(strM,Covar,CovModel);strM_res(:,nroi)=mean(strM)+stats.r;
                stats=regstats(assort,Covar,CovModel);assort_res(:,nroi)=mean(assort)+stats.r;
                stats=regstats(dens,Covar,CovModel);dens_res(:,nroi)=mean(dens)+stats.r;
                stats=regstats(clustM,Covar,CovModel);clustM_res(:,nroi)=mean(clustM)+stats.r;
                stats=regstats(trans,Covar,CovModel);trans_res(:,nroi)=mean(trans)+stats.r;
                stats=regstats(geff,Covar,CovModel);geff_res(:,nroi)=mean(geff)+stats.r;
                stats=regstats(leffM,Covar,CovModel);leffM_res(:,nroi)=mean(leffM)+stats.r;
                stats=regstats(mod,Covar,CovModel);mod_res(:,nroi)=mean(mod)+stats.r;
                stats=regstats(modl,Covar,CovModel);modl_res(:,nroi)=mean(modl)+stats.r;
                stats=regstats(charp,Covar,CovModel);charp_res(:,nroi)=mean(charp)+stats.r;
                stats=regstats(nodebetwM,Covar,CovModel);nodebetwM_res(:,nroi)=mean(nodebetwM)+stats.r;
                stats=regstats(edgebetwM,Covar,CovModel);edgebetwM_res(:,nroi)=mean(edgebetwM)+stats.r;
                stats=regstats(lam,Covar,CovModel);lam_res(:,nroi)=mean(lam)+stats.r;
                stats=regstats(gam,Covar,CovModel);gam_res(:,nroi)=mean(gam)+stats.r;
                stats=regstats(sig,Covar,CovModel);sig_res(:,nroi)=mean(sig)+stats.r;
                
            end
            
        end
        
        for nroi2=1:size(deg,2)^2
        
            stats=regstats(edgebetw(:,nroi2),Covar,CovModel);edgebetw_res(:,nroi2)=edgebetw(:,nroi2)+stats.r;
            stats=regstats(edgebetwN(:,nroi2),Covar,CovModel);edgebetwN_res(:,nroi2)=edgebetwN(:,nroi2)+stats.r;            
        
        end
        
        for subj=1:size(N,1)
            N{subj,thresh}{1,1}=deg_res(subj,:);
            N{subj,thresh}{1,3}=degN_res(subj,:);
            N{subj,thresh}{2,1}=degM_res(subj,:);
            N{subj,thresh}{3,1}=str_res(subj,:);
            N{subj,thresh}{3,3}=strN_res(subj,:);
            N{subj,thresh}{4,1}=strM_res(subj,:);
            N{subj,thresh}{5,1}=assort_res(subj,:);
            N{subj,thresh}{6,1}=dens_res(subj,:);
            N{subj,thresh}{7,1}=clust_res(subj,:)';
            N{subj,thresh}{7,3}=clustN_res(subj,:)';
            N{subj,thresh}{8,1}=clustM_res(subj,:);
            N{subj,thresh}{9,1}=trans_res(subj,:);
            N{subj,thresh}{10,1}=geff_res(subj,:);
            N{subj,thresh}{11,1}=leff_res(subj,:)';
            N{subj,thresh}{11,3}=leffN_res(subj,:)';
            N{subj,thresh}{12,1}=leffM_res(subj,:);
            N{subj,thresh}{13,1}=mod_res(subj,:);
            N{subj,thresh}{14,1}=modl_res(subj,:);
            N{subj,thresh}{15,1}=charp_res(subj,:);
            N{subj,thresh}{16,1}=nodebetw_res(subj,:);
            N{subj,thresh}{16,3}=nodebetwN_res(subj,:);
            N{subj,thresh}{17,1}=nodebetwM_res(subj,:);
            temp=edgebetw_res(subj,:);N{subj,thresh}{18,1}=reshape(temp,sqrt(size(temp,2)),sqrt(size(temp,2)));
            temp=edgebetwN_res(subj,:);N{subj,thresh}{18,3}=reshape(temp,sqrt(size(temp,2)),sqrt(size(temp,2)));
            N{subj,thresh}{19,1}=edgebetwM_res(subj,:);
            N{subj,thresh}{20,1}=lam_res(subj,:);
            N{subj,thresh}{21,1}=gam_res(subj,:);
            N{subj,thresh}{22,1}=sig_res(subj,:);
            
        end
        
    end
    
    
    NetMes_Bin=N;
    
    save(['NetMesBin_D_Adjusted_' GroupName '_FinalThrRange'],'NetMes_Bin','-v7.3');
    
    Inc1=setdiff(1:n1,ExcludeNans);
    Inc2=setdiff(n1+1:n1+n2,ExcludeNans);
    NetMes_Bin=N(1:size(Inc1,2),:);save(['NetMesBin_D_Adjusted_' Group1 '_FinalThrRange'],'NetMes_Bin','-v7.3');
    NetMes_Bin=N(size(Inc1,2)+1:size(Inc1,2)+size(Inc2,2),:);save(['NetMesBin_D_Adjusted_' Group2 '_FinalThrRange'],'NetMes_Bin','-v7.3');
    
    NetMes_Bin=Spare(Inc1,:);save(['NetMesBin_D_ExcludeNanUnadjusted_' Group1 '_FinalThrRange'],'NetMes_Bin','-v7.3');
    NetMes_Bin=Spare(Inc2,:);save(['NetMesBin_D_ExcludeNanUnadjusted_' Group2 '_FinalThrRange'],'NetMes_Bin','-v7.3');
    
end



