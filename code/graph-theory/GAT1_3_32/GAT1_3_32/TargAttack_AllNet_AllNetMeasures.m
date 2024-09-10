function [G1mszT,G2mszT,Grand1mszT,Grand2mszT]=TargAttack_AllNet_AllNetMeasures(GG1,GG2,Gr1,Gr2,InMeas,OutMeas) 

%%
%   This function do targtted attack analysis on the input graphs

%%

%-- Hadi Hosseini: June, 2011


%%   

GUIinput=0;

if nargin<1

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
    
    Range = spm_select(1,'DensityInterval*','Select DensityInterval .mat');
    Xr = load(Range,'Xax');Xrange=Xr.Xax;
    
    D1 = spm_select(1,['CorrMat_' Group1 '*'],['Select ' Group1 ' "CorrMat*.mat" ']);
    D2 = spm_select(1,['CorrMat_' Group2 '*'],['Select ' Group2 ' "CorrMat*.mat" ']);
    D1r = spm_select(1,['CorrMat_rand_' Group1 '*'],'Select random set 1 "CorrMat_rand_*.mat" ');
    D2r = spm_select(1,['CorrMat_rand_' Group2 '*'],'Select random set 2 "CorrMat_rand_*.mat" ');
    
        
    InMeas=input('please type the input attack param in single quotation mark (deg,betw,...)  ');
    OutMeas=input('please type the output attack param in single quotation mark (comp,dist,geff,...)  ');    
    
    Dens = input('type the network density of interest: ');
    
    fprintf('%-4s\n',' loading the input graphs....');
    
    ff1=load(D1,'CorrMat1');GG1=ff1.CorrMat1;%group1 correlation matrix (nROI*nROI*nDensity)
    ff2=load(D2,'CorrMat2');GG2=ff2.CorrMat2;
    ff1r=load(D1r,'CorrMat1_rand');Gr1=ff1r.CorrMat1_rand;%random set 1 correlation matrices Grand{nPerm}(nROI*nROI*nDensity)
    ff2r=load(D2r,'CorrMat2_rand');Gr2=ff2r.CorrMat2_rand;
    
    GUIinput=1;
    
end
%%

if isequal(GUIinput,0)

    if nargin<6
    
        'error...not enough number of input arguments'
    
    elseif nargin>6
    
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

fprintf('%-4s\n',' targetted attack analysis: group 1 network....');
[G1mszT,IndexG1]=TargetAttack_AllNetMeasures(G1,InMeas,OutMeas);

fprintf('%-4s\n',' targetted attack analysis: group 2 network....');
[G2mszT,IndexG2]=TargetAttack_AllNetMeasures(G2,InMeas,OutMeas);

fprintf('%-4s\n',' targetted attack analysis: set 1&2 random networks....');

Grand1mszT=[];Grand2mszT=[];

for k=1:size(Grand1,2)

    fprintf('%-4s\n',[' random network ' num2str(k)]);
    
    Gtemp1=cell2mat(Grand1(k));
    Grand1mszT=[Grand1mszT TargetAttack_AllNetMeasures(Gtemp1,InMeas,OutMeas)];
    
    Gtemp2=cell2mat(Grand2(k));
    Grand2mszT=[Grand2mszT TargetAttack_AllNetMeasures(Gtemp2,InMeas,OutMeas)];
    
end

save(['MSizeTargAttack_' Group1 '_' num2str(100*Dens)],'G1mszT','IndexG1');
save(['MSizeTargAttack_' Group2 '_' num2str(100*Dens)],'G2mszT','IndexG2');
save(['MSizeTargAttack_Rand_' Group1 '_' num2str(100*Dens)],'Grand1mszT','-v7.3');
save(['MSizeTargAttack_Rand_' Group2 '_' num2str(100*Dens)],'Grand2mszT','-v7.3');   


%% 

Pvalue = Alpha;

p_mszT=CL_Pval(Grand2mszT-Grand1mszT,G2mszT-G1mszT,['MeanSizeTrg' '_' num2str(100*Dens)],Tail);

save(['PvalTargAttack_' Group2 'vs' Group1 '_' num2str(100*Dens)],'p_mszT');


%%

figure
plot([1:Sz1]/Sz1,G1mszT,'blackx');hold
plot([1:Sz1]/Sz1,G2mszT,'redo');
G2vsG1T=find(p_mszT<=Pvalue & p_mszT>=0);t1T=ones(1,size(G2vsG1T,1));t1T(:)=1-Pvalue;
plot(G2vsG1T/Sz1,t1T,'red*');
G1vsG2T=find(p_mszT<0 & p_mszT>=-Pvalue);t2T=ones(1,size(G1vsG2T,1));t2T(:)=1-Pvalue;
plot(G1vsG2T/Sz1,t2T,'black*');
xlabel('fraction of removed nodes');ylabel('relative size of the remaining giant component');
legend(Group1,Group2);
hgsave(['MeanSizeDiff_TargAtt_' Group1 Group2 '_' num2str(100*Dens) '.fig'])


