function [G1mszT,G2mszT,Grand1mszT,Grand2mszT]=TargAttack_AllNet_AllNetMeasures_AllDensity(Group1,Group2,G1,G2,Grand1,Grand2,InMeas,OutMeas) 

%% Hadi Hosseini June 2011

%This function do targtted attack analysis on the input graphs

%% Note
% the input graphs (original and random) should be either binary or thresholded weighted graphs
%% Inputs
%The input graphs are:
%G1(first group binary graph),output R1 from 'Rand_Net_Bootstrap.m' (N*N array)
%G2(second group binary graph),output R2 from 'Rand_Net_Bootstrap.m' (N*N array)
%Grand1(set of randomly generated graphs corresponding to G1);output R_G1 from 'Rand_Net_Bootstrap.m' whihc is a 1*N cell array containing N random graphs
%Gradn2(set of randomly generated graphs corresponding to G2);output R_G2 from 'Rand_Net_Bootstrap.m' whihc is a 1*N cell array containing N random graphs
%nsimul is the number of simualtions (default 100)

%% outputs
%G1msz/T: N*1 vector of mean size of the remaining vector after consecutive random/targetted deletion of G1 nodes
%G2msz/T: N*1 vector of mean size of the remaining vector after consecutive random/targetted deletion of G2 nodes
%Grand1msz/T: N*M array of mean size of the remaining vector after consecutive random/targetted deletion nodes of each of the M nrandom networks of Grand1 (M is size(Grand1,2)
%Grand1msz/T: N*M array of mean size of the remaining vector after consecutive random/targetted deletion nodes of each of the M nrandom networks of Grand2 (M is size(Grand2,2)

%% Input dialog  
GUIinput=0;
if nargin<1
    %Data = spm_select(1,'RandNets_ThreshAtMin*','Select the "RandNet_ThreshAtMin.mat" (output from "compare graphs at min density")');
    %KAN = spm_select(1,'KAN_FixedDens*','Select "KAN_FixedDens_Results" (output from "graph measures at min density")');%to find the minimum edg(density) for thresholding
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;
    Group2 = ff.mat4GAT.g2;
    
    if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha') %&& isfield(mat4GAT,'dir') 
        
        Alpha = ff.mat4GAT.alpha;
        Tail = ff.mat4GAT.tail;
        %Dir = ff.mat4GAT.dir;
        
    else
        
        Alpha = .05;
        Tail = 1;
        %Dir = 1;
    
    end
    
    
    
    Range = spm_select(1,'DensityInterval*','Select DensityInterval .mat');
    Xrange = load(Range,'Xax');
    
    D1 = spm_select(1,['CorrMat_' Group1 '*'],['Select ' Group1 ' "CorrMat*.mat" ']);
    D2 = spm_select(1,['CorrMat_' Group2 '*'],['Select ' Group2 ' "CorrMat*.mat" ']);
    D1r = spm_select(1,['CorrMat_rand_' Group1 '*'],'Select random set 1 "CorrMat_rand_*.mat" ');
    D2r = spm_select(1,['CorrMat_rand' Group2 '*'],'Select random set 2 "CorrMat_rand_*.mat" ');
    
    InMeas=input('type the input attack param in single quotation mark (deg,betw,...)  ');
    OutMeas=input('type the output attack param in single quotation mark (comp,dist,geff,...)  ');    
    
    Dens = input('type the network density of interest: ');
    
    fprintf('%-4s\n',' loading the input graphs....');
    ff1=load(D1,'CorrMat1');GG1=ff1.CorrMat1;%group1 correlation matrix (nROI*nROI*nDensity)
    ff2=load(D2,'CorrMat2');GG2=ff2.CorrMat2;%group2
    ff1r=load(D1r,'CorrMat1_rand');Gr1=ff1r.CorrMat1_rand;%random set 1 correlation matrices
    ff2r=load(D2r,'CorrMat2_rand');Gr2=ff2r.CorrMat2_rand;%random set 2 correlation matrices Grand{nPerm}(nROI*nROI*nDensity)
    
    %%%
    GUIinput=1;%GUI based
end
%%
if isequal(GUIinput,0)%command-based run
    if nargin<6
        'error...not enough number of input arguments'
    elseif nargin>4%5
        'error...too many input arguments'
    end
end

%extract the networks in the specified density
[dd, idx] = find(Xrange == Dens);

G1 = GG1(:,:,idx);
G2 = GG2(:,:,idx);
for kk = 1:size(Gr1,2)%number of random networks
    
    Grand1{kk} = Gr1{kk}(:,:,idx);%cells are used just to be compatible with previous version
    Grand2{kk} = Gr2{kk}(:,:,idx);

end

%% converting to binary

%making diagonal zero 
% on original matrices
Sz1=size(G1,1);Sz2=size(G2,1);
G1(1:Sz1+1:end)=0;G2(1:Sz2+1:end)=0;

%on random generated matrices
Grand1=cellfun(@zero_diag,Grand1,'UniformOutput', false);%apply zero_diag function on each cell of the cell array
Grand2=cellfun(@zero_diag,Grand2,'UniformOutput', false);             

%change negative corrrelations to zero
%original graphs
G1(G1<0)=0;G2(G2<0)=0;

%random graphs
Grand1=cellfun(@neg_zero,Grand1,'UniformOutput', false);%apply zero_diag function on each cell of the cell array
Grand2=cellfun(@neg_zero,Grand2,'UniformOutput', false);          

%convert weights to one
G1(G1>0)=1;G2(G2>0)=1;

%to random graphs
Grand1=cellfun(@pos_one,Grand1,'UniformOutput', false);%apply pos_one function on each cell of the cell array
Grand2=cellfun(@pos_one,Grand2,'UniformOutput', false);          

%% attack analysis
%G1
fprintf('%-4s\n',' targetted attack analysis: group 1 network....');
[G1mszT,IndexG1]=TargetAttack_AllNetMeasures(G1,InMeas,OutMeas);%,nsimul);
%G2
fprintf('%-4s\n',' targetted attack analysis: group 2 network....');
[G2mszT,IndexG2]=TargetAttack_AllNetMeasures(G2,InMeas,OutMeas);%,nsimul);
%Grand1
fprintf('%-4s\n',' targetted attack analysis: set 1&2 random networks....');

Grand1mszT=[];Grand2mszT=[];
for k=1:size(Grand1,2)
    fprintf('%-4s\n',[' random network ' num2str(k)]);
    %Grand1
    Gtemp1=cell2mat(Grand1(k));
    Grand1mszT=[Grand1mszT TargetAttack_AllNetMeasures(Gtemp1,InMeas,OutMeas)];%,nsimul)];
    %Grand2
    Gtemp2=cell2mat(Grand2(k));
    Grand2mszT=[Grand2mszT TargetAttack_AllNetMeasures(Gtemp2,InMeas,OutMeas)];%,nsimul)];
end
save(['MSizeTargAttack_' Group1 '_' num2str(100*idx)],'G1mszT','IndexG1');
save(['MSizeTargAttack_' Group2 '_' num2str(100*idx)],'G2mszT','IndexG2');
save(['MSizeTargAttack_Rand_' Group1 '_' num2str(100*idx)],'Grand1mszT','-v7.3');
save(['MSizeTargAttack_Rand_' Group2 '_' num2str(100*idx)],'Grand2mszT','-v7.3');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating mean and CI of the mean size of the random netwroks 
%Pvalue=0.05;
if isequal(Tail,1)
    
    Pvalue = Alpha;

elseif isequal(Tail,2)
    
    Pvalue = Alpha./2;

end


p_mszT=CL_Pval(Grand2mszT-Grand1mszT,G2mszT-G1mszT,'MeanSizeTrg');

save(['PvalTargAttack_' Group2 'vs' Group1 '_' num2str(100*idx)],'p_mszT');


%#########################################
% plotting relative network size as a function of removed nodes
%#########################################

%targetted attack
figure
plot([1:Sz1]/Sz1,G1mszT,'blackx');hold
plot([1:Sz1]/Sz1,G2mszT,'redo');
G2vsG1T=find(p_mszT<=Pvalue & p_mszT>=0);t1T=ones(1,size(G2vsG1T,1));t1T(:)=1-Pvalue;
plot(G2vsG1T/Sz1,t1T,'red*');
G1vsG2T=find(p_mszT<0 & p_mszT>=-Pvalue);t2T=ones(1,size(G1vsG2T,1));t2T(:)=1-Pvalue;
plot(G1vsG2T/Sz1,t2T,'black*');
xlabel('fraction of removed nodes');ylabel('relative size of the remaining giant component');
legend(Group1,Group2);
hgsave(['MeanSizeDiff_TargAtt_' Group1 '_' Group2  '_' num2str(100*idx) '.fig'])


