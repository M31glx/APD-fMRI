% This function tests the significance of the difference in area under the
% curve (auc) between two groups net measures. the auc is calculated based
% on several density criteria: 1. auc between MinDens and MaxDens; 2. auce
% between MinDens of full connectivity to Max Dens; 3. askes for the
% density interval that you are interested in

%%  

RandorTarg=input('type 1 for "targeted attack", 2 for "random failure", and 3 for both: ');  
data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Net Measures)');%Data1 = spm_select(1,'.mat',['Select mat4SVM .mat (output from concatenate)']);
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


if isequal(RandorTarg,1) || isequal(RandorTarg,3)
    MSizeTarg1 = spm_select(1,['MSizeTargAttack_' Group1],'Select the "MSizeTargAttack_Group1.mat" (output from "Targeted Attack Analysis")');
    MSizeTarg2 = spm_select(1,['MSizeTargAttack_' Group2],'Select the "MSizeTargAttack_Group2.mat" (output from "Targeted Attack Analysis")');
    MSizeTarg1_rand = spm_select(1,['MSizeTargAttack_Rand_' Group1],'Select the "MSizeTargAttack_Rand_Group1.mat" (output from "Targeted Attack Analysis")');
    MSizeTarg2_rand = spm_select(1,['MSizeTargAttack_Rand_' Group2],'Select the "MSizeTargAttack_Rand_Group2.mat" (output from "Targeted Attack Analysis")');
end

if isequal(RandorTarg,2) || isequal(RandorTarg,3)
    MSizeRand1 = spm_select(1,['MSizeRandAttack_' Group1],'Select the "MSizeRandAttack_Group1.mat" (output from "Random Attack Analysis")');
    MSizeRand2 = spm_select(1,['MSizeRandAttack_' Group2],'Select the "MSizeRandAttack_Group2.mat" (output from "Random Attack Analysis")');
    MSizeRand1_rand = spm_select(1,['MSizeRandAttack_Rand_' Group1],'Select the "MSizeRandAttack_Rand_Group1.mat" (output from "Random Attack Analysis")');
    MSizeRand2_rand = spm_select(1,['MSizeRandAttack_Rand_' Group2],'Select the "MSizeRandAttack_Rand_Group2.mat" (output from "Random Attack Analysis")');
end


%% 

if isequal(RandorTarg,1) || isequal(RandorTarg,3)
   
    f1 = load(MSizeTarg1,'G1mszT');G1mszT = f1.G1mszT;
   
    f2 = load(MSizeTarg2,'G2mszT');G2mszT = f2.G2mszT;
   
    f1_rand = load(MSizeTarg1_rand,'Grand1mszT');Grand1mszT = f1_rand.Grand1mszT;
   
    f2_rand = load(MSizeTarg2_rand,'Grand2mszT');Grand2mszT = f2_rand.Grand2mszT;
    xx = size(G1mszT,1);iXax =[1/xx:1/xx:1];
    
end

if isequal(RandorTarg,2) || isequal(RandorTarg,3)

    e1 = load(MSizeRand1,'G1msz');G1msz = e1.G1msz;

    e2 = load(MSizeRand2,'G2msz');G2msz = e2.G2msz;

    e1_rand = load(MSizeRand1_rand,'Grand1msz');Grand1msz = e1_rand.Grand1msz;

    e2_rand = load(MSizeRand2_rand,'Grand2msz');Grand2msz = e2_rand.Grand2msz;
    xx = size(G1msz,1);iXax =[1/xx:1/xx:1];
    
end

%% 

if isequal(RandorTarg,1) || isequal(RandorTarg,3)

    auc_TA_G1 = trapz(iXax,G1mszT);

    auc_TA_G2 = trapz(iXax,G2mszT);

    auc_TA_G1_rand = trapz(iXax,Grand1mszT);

    auc_TA_G2_rand = trapz(iXax,Grand2mszT);
    
end

%% 
if isequal(RandorTarg,2) || isequal(RandorTarg,3)
   
    auc_RA_G1 = trapz(iXax,G1msz);
   
    auc_RA_G2 = trapz(iXax,G2msz);
   
    auc_RA_G1_rand = trapz(iXax,Grand1msz);
   
    auc_RA_G2_rand = trapz(iXax,Grand2msz);

end
%% 

if isequal(RandorTarg,1) || isequal(RandorTarg,3)
	
    p_auc_TargetAttack = CL_Pval(auc_TA_G2_rand - auc_TA_G1_rand, auc_TA_G2 - auc_TA_G1,'AUC_TargetAttack',Tail)
    
end

if isequal(RandorTarg,2) || isequal(RandorTarg,3)

    p_auc_RandomAttack = CL_Pval(auc_RA_G2_rand - auc_RA_G1_rand, auc_RA_G2 - auc_RA_G1,'AUC_RandomAttack',Tail)

end

