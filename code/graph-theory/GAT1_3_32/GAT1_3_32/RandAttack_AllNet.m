function [G1msz,G2msz,Grand1msz,Grand2msz]=RandAttack_AllNet(Group1, Group2, GG1, GG2, Gr1, Gr2, nsimul) 

%% 
%  This function do random failure analysis on the input graphs

%%

%-- Hadi Hosseini: June, 2011


%%

GUIinput=0;

if nargin<1
    
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;
    Group2 = ff.mat4GAT.g2;
    
    if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha') %&& isfield(mat4GAT,'dir') 
        
        Alpha = ff.mat4GAT.alpha;
        Tail = ff.mat4GAT.tail;
        
    else
        
        Alpha = .05;
        Tail = 2;
    
    end
    
    Range = spm_select(1,'DensityInterval*','Select DensityInterval .mat');
    Xr = load(Range,'Xax');Xrange=Xr.Xax;
    
    D1 = spm_select(1,['CorrMat_' Group1 '*'],['Select ' Group1 ' "CorrMat*.mat" ']);
    D2 = spm_select(1,['CorrMat_' Group2 '*'],['Select ' Group2 ' "CorrMat*.mat" ']);
    D1r = spm_select(1,['CorrMat_rand_' Group1 '*'],'Select random set 1 "CorrMat_rand_*.mat" ');
    D2r = spm_select(1,['CorrMat_rand_' Group2 '*'],'Select random set 2 "CorrMat_rand_*.mat" ');
    
    nsimul=input('please type number of simulations ');
    Dens = input('type the network density of interest: ');
    
    fprintf('%-4s\n',' loading the input graphs....');
    
    ff1=load(D1,'CorrMat1');GG1=ff1.CorrMat1;%group1 correlation matrix (nROI*nROI*nDensity)
    ff2=load(D2,'CorrMat2');GG2=ff2.CorrMat2;
    ff1r=load(D1r,'CorrMat1_rand');Gr1=ff1r.CorrMat1_rand;%random set 1 correlation matrices Gr1{nPerm}(nROI*nROI*nDensity)
    ff2r=load(D2r,'CorrMat2_rand');Gr2=ff2r.CorrMat2_rand;
    
    GUIinput=1;
end

%%

if isequal(GUIinput,0)%command-based run

    if nargin<4
    
        'error...not enough number of input arguments'
    
    elseif isequal(nargin,5)
    
        'default number of simulations was applied'
        nsimul=100;
    
    elseif nargin>5
    
        'error...too many input arguments'
    
    end
    
end

%%


[dd, idx] = find(Xrange == Dens);

G1 = GG1(:,:,idx);
G2 = GG2(:,:,idx);
for rr = 1:size(Gr1,2)%number of random networks
    
    Grand1{rr} = Gr1{rr}(:,:,idx);%cells are used just to be compatible with previous version
    Grand2{rr} = Gr2{rr}(:,:,idx);

end


%% 

Sz1=size(G1,1);Sz2=size(G2,1);
G1(1:Sz1+1:end)=0;G2(1:Sz2+1:end)=0;

Grand1=cellfun(@zero_diag,Grand1,'UniformOutput', false);
Grand2=cellfun(@zero_diag,Grand2,'UniformOutput', false);             

G1(G1<0)=0;G2(G2<0)=0;

Grand1=cellfun(@neg_zero,Grand1,'UniformOutput', false);
Grand2=cellfun(@neg_zero,Grand2,'UniformOutput', false);          

G1(G1>0)=1;G2(G2>0)=1;

Grand1=cellfun(@pos_one,Grand1,'UniformOutput', false);
Grand2=cellfun(@pos_one,Grand2,'UniformOutput', false);          


%% 

fprintf('%-4s\n',' random attack analysis: group 1 network....');
G1msz=RandomAttack(G1,nsimul);

fprintf('%-4s\n',' random attack analysis: group 2 network....');
G2msz=RandomAttack(G2,nsimul);

fprintf('%-4s\n',' random attack analysis: set 1&2 random networks....');

Grand1msz=[];Grand2msz=[];

for k=1:size(Grand1,2)
 
    fprintf('%-4s\n',[' random network ' num2str(k)]);
    
    Gtemp1=cell2mat(Grand1(k));
    Grand1msz=[Grand1msz RandomAttack(Gtemp1,nsimul)];
    
    Gtemp2=cell2mat(Grand2(k));
    Grand2msz=[Grand2msz RandomAttack(Gtemp2,nsimul)];
    
end

save(['MSizeRandAttack_' Group1 '_' num2str(100*Dens)],'G1msz');
save(['MSizeRandAttack_' Group2 '_' num2str(100*Dens)],'G2msz');
save(['MSizeRandAttack_Rand_' Group1 '_' num2str(100*Dens)],'Grand1msz','-v7.3');
save(['MSizeRandAttack_Rand_' Group2 '_' num2str(100*Dens)],'Grand2msz','-v7.3');


%% 

Pvalue = Alpha;

p_msz=CL_Pval(Grand2msz-Grand1msz,G2msz-G1msz,['MeanSizeTrg' '_' num2str(100*Dens)],Tail);

save(['PvalRandAttack_' Group2 'vs' Group1 '_' num2str(100*Dens)],'p_msz');


%%

figure
plot([1:Sz1]/Sz1,G1msz,'blackx');hold
plot([1:Sz1]/Sz1,G2msz,'redo');
G2vsG1=find(p_msz<=Pvalue & p_msz>=0);t1=ones(1,size(G2vsG1,1));t1(:)=1-Pvalue;
plot(G2vsG1/Sz1,t1,'red*');
G1vsG2=find(p_msz<0 & p_msz>=-Pvalue);t2=ones(1,size(G1vsG2,1));t2(:)=1-Pvalue;
plot(G1vsG2/Sz1,t2,'black*');
xlabel('fraction of removed nodes');ylabel('relative size of the remaining giant component');
legend(Group1,Group2);
hgsave(['MeanSizeDiff_RandAtt_' Group1 '_' Group2 '_' num2str(100*Dens) '.fig'])

