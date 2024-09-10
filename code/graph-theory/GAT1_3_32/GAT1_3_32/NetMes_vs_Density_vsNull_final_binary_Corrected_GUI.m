function [NetMes_B1,NetMes_B2,NetMes_B1_rand,NetMes_B2_rand,Thresh_Density1_rand,Thresh_Density2_rand,p_Sigma,p_Gamma,p_Lambda]=NetMes_vs_Density_vsNull_final_binary_Corrected_GUI(Input1,Input2,rand_Input1,rand_Input2,Group1,Group2,MinMes,MaxMes,MesStep,FigOn) 

%'Input1/Input2': input the "unthresholded" correlation matrix for each group (output R1/R2 or P1/P2 from 'Rand_Net_Bootstrap.m')
%'rand_Input1,rand_Input2': Input the randomly generated untheresholded
%                           correlation matrices (R_G1/R_G2 or P_G1/P_G2 output from 'Rand_Net_Bootstrap.m')  each of which is
%                           a (1byN) cell array containing the N randomly generated correlation matrices
%'BiggerIsBetter': enter 1 for r value (correlation) and 0 for p-value;
%'Measure': the measure of interest (e.g. density)
%'MinMes': The minimum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.05 for density)
%'MaxMes': The maximum value of the 'Measure' for which you want to plot
%           the NetMeasure values (default 0.45 for density)


%%
%-- Hadi Hosseini Mar 22, 2011
%   Revised for GUI 
%-- Hadi Hosseini Apr.20,2011

%%
Precision=10;
Step=0.0001;
Measure='density';

%% Input dialog  
GUIinput=0;
if nargin<1
    Data = spm_select(1,'RandNetworksBootstrap*','Select the RandNetworksBootstrap_Results.mat file (output from Generate Random Graphs)');
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;
    Group2 = ff.mat4GAT.g2;
    
    if isfield(ff.mat4GAT,'null')
        NullType = ff.mat4GAT.null;
        RandIter = ff.mat4GAT.nullIter;       
    else
        NullType = 2;
        RandIter = 20;
    end
    
    if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha') 
        
        Alpha = ff.mat4GAT.alpha;
        Tail = ff.mat4GAT.tail;
        
        
    else
        
        Alpha = .05;
        Tail = 2;
        
    
    end
    
    
    
    MinMes=input('please type the lowest density of interest (e.g. 0.1):  ');
    MaxMes=input('please type the highest density of interest (e.g. 0.5):  ');
    MesStep=input('please type the density interval (setps) of interest (e.g. 0.02):  ');
    
    f1 = load(Data,'R1');Input1 = f1.R1;clear f1;Inp1 = Input1;
    f2 = load(Data,'R2');Input2 = f2.R2;clear f2;Inp2 = Input2;
    e1 = load(Data,'R_G1');rand_Input1 = e1.R_G1;clear e1;rInp1 = rand_Input1;
    e2 = load(Data,'R_G2');rand_Input2 = e2.R_G2;clear e2;rInp2 = rand_Input2;
    
    
    if isequal(NullType,3) || isequal(NullType,4)
            
        w1 = load(Data,'data1');d1 = w1.data1;clear w1
        w2 = load(Data,'data2');d2 = w2.data2;clear w2
        
        q1 = load(Data,'dataG1');dr1 = q1.dataG1;clear q1
        q2 = load(Data,'dataG2');dr2 = q2.dataG2;clear q2

    end
    
    GUIinput = 1;
    FigOn = 1;

end


%%

MinMesPlot=MinMes;
MaxMesPlot=MaxMes;
MesStepPlot=MesStep;

N = size(Input1,1);
MinEdge=MinMes*((N^2-N)/2);MinMes=floor(MinEdge);
MaxEdge=MaxMes*((N^2-N)/2);MaxMes=ceil(MaxEdge);
EdgeStep=MesStep*((N^2-N)/2);MesStep=round(EdgeStep);

%% checking input

if isequal(GUIinput,0)%command-based run
    
    if nargin<6
        
        'error...not enough number of input arguments'
        
    elseif isequal(nargin,6)
        
        'default min and max of the measures was applied'
        MinMes=0.05
        MaxMes=0.5
        MesStep=0.01
        NullType = 2;
        RandIter = 20;
        Alpha = .05;
        Tail = 2;
        FigOn = 1;
        rInp1 = rand_Input1;
        rInp2 = rand_Input2;
        
    elseif nargin>10
        
        'error...too many input arguments'
    
    end
    
end


%% Making different threshold

switch Measure
    
    case 'density'
        
        fprintf('%-4s\n',' calculating threshold level. it might take several hours depending on the number of random networks....');
        
        Sz1=size(Input1,1);Sz2=size(Input2,1);
        Input1(1:Sz1+1:end)=0;Input2(1:Sz2+1:end)=0;
        rand_Input1=cellfun(@zero_diag,rand_Input1,'UniformOutput', false);
        rand_Input2=cellfun(@zero_diag,rand_Input2,'UniformOutput', false);             
    
        
        Input1(Input1<0)=0;Input2(Input2<0)=0;
        rand_Input1=cellfun(@neg_zero,rand_Input1,'UniformOutput', false);
        rand_Input2=cellfun(@neg_zero,rand_Input2,'UniformOutput', false);          
        
        R_1=Input1;R_2=Input2;
        Rrand_1=rand_Input1;Rrand_2=rand_Input2;
        Thresh_Density1=[];
        Thresh_Density2=[];
        Thresh_Density1_rand=cell(size(rand_Input1));
        Thresh_Density2_rand=cell(size(rand_Input1));
        
        R1max=max(max(R_1));
        R2max=max(max(R_2));
        RG1max_vector=cellfun(@max,Rrand_1,'UniformOutput', false);RG1max=cellfun(@max,RG1max_vector,'UniformOutput', false);
        RG2max_vector=cellfun(@max,Rrand_2,'UniformOutput', false);RG2max=cellfun(@max,RG2max_vector,'UniformOutput', false);
        
        fprintf('%-4s\n',' calculating threshold level: group 1 network....');
        
        for i=Step:Step:R1max
            
            R_1(R_1<=i)=0;
            Thresh_Density1=[Thresh_Density1;i nnz(triu(R_1))];
            
        end
        
        fprintf('%-4s\n',' calculating threshold level: group 2 network....');
        
        for j=Step:Step:R2max
        
            R_2(R_2<=j)=0;
            Thresh_Density2=[Thresh_Density2;j nnz(triu(R_2))];
            
        end
        
        fprintf('%-4s\n',' calculating threshold level: set 1 random networks....');
        
        for k=1:size(Rrand_1,2)
        
            Rtemp=cell2mat(Rrand_1(k));
            Thresh_Density_temp=[];
            
            for L=Step:Step:cell2mat(RG1max(k))
                
                Rtemp(Rtemp<=L)=0;
                Thresh_Density_temp=[Thresh_Density_temp;L nnz(triu(Rtemp))];
                
            end
            
            Thresh_Density1_rand{k}=Thresh_Density_temp;
            
        end
        
        
        fprintf('%-4s\n',' calculating threshold level: set 2 random networks....');
        
        for p=1:size(Rrand_2,2)
        
            Rtemp=cell2mat(Rrand_2(p));
            Thresh_Density_temp=[];
            
            for q=Step:Step:cell2mat(RG2max(p))
            
                Rtemp(Rtemp<=q)=0;
                Thresh_Density_temp=[Thresh_Density_temp;q nnz(triu(Rtemp))];
            
            end
            
            Thresh_Density2_rand{p}=Thresh_Density_temp;
        
        end

end

save(['ThreshDens_' Group1] ,'Thresh_Density1');
save(['ThreshDens_' Group2],'Thresh_Density2');
save(['ThreshDens_rand_' Group1],'Thresh_Density1_rand','-v7.3');
save(['ThreshDens_rand_' Group2],'Thresh_Density2_rand','-v7.3');

clear R_1 R_2 Rrand_1 Rrand_2

Thresh_Density1 = single(Thresh_Density1);
Thresh_Density2 = single(Thresh_Density2);
Thresh_Density1_rand = cellfun(@single,Thresh_Density1_rand,'UniformOutput', false);
Thresh_Density2_rand = cellfun(@single,Thresh_Density2_rand,'UniformOutput', false);


NetMes_B1=cell(1,length(MinMes:MesStep:MaxMes));
NetMes_B2=cell(1,length(MinMes:MesStep:MaxMes));
NetMes_B1_rand=cell(1,size(rand_Input1,2));
NetMes_B2_rand=cell(1,size(rand_Input1,2));

save NetMes_B1_rand NetMes_B1_rand
save NetMes_B2_rand NetMes_B2_rand

%%

count=0;

for u=MinMes:MesStep:MaxMes
    
    count=count+1;
    
    fprintf('%-4s\n',[' calculating network measures at each density: ' num2str(count) ' out of ' num2str(ceil((MaxMes-MinMes)./MesStep)) '...']);
    
    %%
    Ind=find(Thresh_Density1(:,2) >= u-Precision & Thresh_Density1(:,2) <= u+Precision );
    Dind=Thresh_Density1(Ind,2)-u;
    Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
    Ind=Ind(min(Ind_Ind));
    
    
    if isempty(Ind)
        
        fprintf('%-4s\n',' No same density within the default precision: decreasing the precision for original net 1....');
        Ind=find(Thresh_Density1(:,2) >= u-100*Precision & Thresh_Density1(:,2) <= u+100*Precision );
        Dind=Thresh_Density1(Ind,2)-u;%distance from the u
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
    end
    
    Thresh=Thresh_Density1(Ind,1);
    R_1=Input1;R_1(R_1<=Thresh)=0;
        
    
    if isequal(NullType,3) || isequal(NullType,4)
        
        
        NetMes_Bin = NetMeasures_Binary(R_1,NullType,RandIter,d1,Inp1);
                
    else
        
        NetMes_Bin = NetMeasures_Binary(R_1,NullType,RandIter);
        
    end 

    NetMes_Bin{23,1}=Thresh;
    
    NetMes_B1{count}=NetMes_Bin;clear NetMes_Bin
    
    CorrMat1(:,:,count) = single(R_1);
    
    
    %%
    
    Ind=find(Thresh_Density2(:,2) >= u-Precision & Thresh_Density2(:,2) <= u+Precision );
    Dind=Thresh_Density2(Ind,2)-u;
    Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
    Ind=Ind(min(Ind_Ind));
    
    if isempty(Ind)
        
        fprintf('%-4s\n',' No same density within the default precision: decreasing the precision for original net 2....');
        Ind=find(Thresh_Density2(:,2) >= u-100*Precision & Thresh_Density2(:,2) <= u+100*Precision );
        Dind=Thresh_Density2(Ind,2)-u;
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
    end
    
    Thresh=Thresh_Density2(Ind,1);
    R_2=Input2;R_2(R_2<=Thresh)=0;
    
    if isequal(NullType,3) || isequal(NullType,4)
        
        NetMes_Bin = NetMeasures_Binary(R_2,NullType,RandIter,d2,Inp2);
        
    else
        
        NetMes_Bin=NetMeasures_Binary(R_2,NullType,RandIter);
    
    end
    
    NetMes_Bin{23,1}=Thresh;
    NetMes_B2{count}=NetMes_Bin;clear NetMes_Bin
    CorrMat2(:,:,count) = single(R_2);
    
    %%
    
    for z=1:size(rand_Input1,2)
        
            
        Thresh_Density_temp=cell2mat(Thresh_Density1_rand(z));
        Rtemp=cell2mat(rand_Input1(z));
        
        Ind=find(Thresh_Density_temp(:,2) >= u-Precision & Thresh_Density_temp(:,2) <= u+Precision );
        Dind=Thresh_Density_temp(Ind,2)-u;%distance from the u
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
        if isempty(Ind)
        
            fprintf(['%-4s\n',' No same density within the default precision: decreasing the precision for random nets 1, net ' num2str(z)]);
            Ind=find(Thresh_Density_temp(:,2) >= u-100*Precision & Thresh_Density_temp(:,2) <= u+100*Precision );
            Dind=Thresh_Density_temp(Ind,2)-u;%distance from the u
            Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
            Ind=Ind(min(Ind_Ind));
            
        end
        
        Thresh=Thresh_Density_temp(Ind,1);
        Rtemp(Rtemp<=Thresh)=0;
        
        if isequal(NullType,3) || isequal(NullType,4)
            
            NetMes_Bin = NetMeasures_Binary(Rtemp,NullType,RandIter,dr1{z},rInp1{z});
            
        else
            
            NetMes_Bin = NetMeasures_Binary(Rtemp,NullType,RandIter);
        
        end
        
        NetMes_Bin{23,1}=Thresh;
               
        if z == 1
            load NetMes_B1_rand
        end
        
        NetMes_B1_rand{z}{count} = NetMes_Bin;clear NetMes_Bin
        
        if z == size(rand_Input1,2)
            save('NetMes_B1_rand','NetMes_B1_rand','-v7.3');
            clear NetMes_B1_rand
        
        end
        
        CorrMat1_rand{z}(:,:,count) = single(Rtemp);
        
    end
    
    %%
    
    for z=1:size(rand_Input2,2)
        
        Thresh_Density_temp=cell2mat(Thresh_Density2_rand(z));
        Rtemp=cell2mat(rand_Input2(z));
        
        Ind=find(Thresh_Density_temp(:,2) >= u-Precision & Thresh_Density_temp(:,2) <= u+Precision );
        Dind=Thresh_Density_temp(Ind,2)-u;
        Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
        Ind=Ind(min(Ind_Ind));
        
        if isempty(Ind)
            
            fprintf(['%-4s\n',' No same density within the default precision: decreasing the precision for random nets 2, net ' num2str(z)]);
            Ind=find(Thresh_Density_temp(:,2) >= u-100*Precision & Thresh_Density_temp(:,2) <= u+100*Precision );
            Dind=Thresh_Density_temp(Ind,2)-u;%distance from the u
            Ind_Ind=find(Dind==min(abs(Dind)) | Dind==-min(abs(Dind)));
            Ind=Ind(min(Ind_Ind));
            
        end
        
        Thresh=Thresh_Density_temp(Ind,1);
        Rtemp(Rtemp<=Thresh)=0;
        
        if isequal(NullType,3) || isequal(NullType,4)
        
            NetMes_Bin = NetMeasures_Binary(Rtemp,NullType,RandIter,dr2{z},rInp2{z});
            
        else
            
            NetMes_Bin=NetMeasures_Binary(Rtemp,NullType,RandIter);
        
        end
        
        NetMes_Bin{23,1}=Thresh;
        
        if z == 1
            load NetMes_B2_rand
        end
        
        NetMes_B2_rand{z}{count}=NetMes_Bin;clear NetMes_Bin        
        
        if z == size(rand_Input2,2)
            save('NetMes_B2_rand','NetMes_B2_rand','-v7.3');
            clear NetMes_B2_rand
        end
        
        CorrMat2_rand{z}(:,:,count) = single(Rtemp);

    end
      
end

save(['NetMesBin_' Group1],'NetMes_B1');
save(['NetMesBin_' Group2],'NetMes_B2');

load NetMes_B1_rand;save(['NetMesBin_rand_' Group1],'NetMes_B1_rand','-v7.3');clear NetMes_B1_rand 
load NetMes_B2_rand;save(['NetMesBin_rand_' Group2],'NetMes_B2_rand','-v7.3');clear NetMes_B2_rand 

save(['CorrMat_' Group1],'CorrMat1');
save(['CorrMat_' Group2],'CorrMat2');

save(['CorrMat_rand_' Group1],'CorrMat1_rand','-v7.3');
save(['CorrMat_rand_' Group2],'CorrMat2_rand','-v7.3');


load NetMes_B1_rand 
load NetMes_B2_rand 

%%

fprintf('%-4s\n',' calculating confidence intervals for all network measures...');

Dens_1_rand=[];
MClust_1_rand=[];MDeg_1_rand=[];MStr_1_rand=[];Trans_1_rand=[];Assort_1_rand=[];
GEff_1_rand=[];MLocEff_1_rand=[];Mod_1_rand=[];ModL_1_rand=[];
PathL_1_rand=[];MNodeBetw_1_rand=[];MEdgeBetw_1_rand=[];
Lambda_1_rand=[];Gamma_1_rand=[];Sigma_1_rand=[];

Dens_2_rand=[];
MClust_2_rand=[];MDeg_2_rand=[];MStr_2_rand=[];Trans_2_rand=[];Assort_2_rand=[];
GEff_2_rand=[];MLocEff_2_rand=[];Mod_2_rand=[];ModL_2_rand=[];
PathL_2_rand=[];MNodeBetw_2_rand=[];MEdgeBetw_2_rand=[];
Lambda_2_rand=[];Gamma_2_rand=[];Sigma_2_rand=[];

Dens_1=[];
MClust_1=[];MDeg_1=[];MStr_1=[];Trans_1=[];Assort_1=[];
GEff_1=[];MLocEff_1=[];Mod_1=[];ModL_1=[];
PathL_1=[];MNodeBetw_1=[];MEdgeBetw_1=[];
Lambda_1=[];Gamma_1=[];Sigma_1=[];

Dens_2=[];
MClust_2=[];MDeg_2=[];MStr_2=[];Trans_2=[];Assort_2=[];
GEff_2=[];MLocEff_2=[];Mod_2=[];ModL_2=[];
PathL_2=[];MNodeBetw_2=[];MEdgeBetw_2=[];
Lambda_2=[];Gamma_2=[];Sigma_2=[];

for i=1:size(rand_Input1,2)
    
    Dens_Temp1=[];
    MClust_Temp1=[];MDeg_Temp1=[];MStr_Temp1=[];Trans_Temp1=[];Assort_Temp1=[];
    GEff_Temp1=[];MLocEff_Temp1=[];Mod_Temp1=[];ModL_Temp1=[];
    PathL_Temp1=[];MNodeBetw_Temp1=[];MEdgeBetw_Temp1=[];
    Lambda_Temp1=[];Gamma_Temp1=[];Sigma_Temp1=[];
    
    Dens_Temp2=[];
    MClust_Temp2=[];MDeg_Temp2=[];MStr_Temp2=[];Trans_Temp2=[];Assort_Temp2=[];
    GEff_Temp2=[];MLocEff_Temp2=[];Mod_Temp2=[];ModL_Temp2=[];
    PathL_Temp2=[];MNodeBetw_Temp2=[];MEdgeBetw_Temp2=[];
    Lambda_Temp2=[];Gamma_Temp2=[];Sigma_Temp2=[];
    
    for j=1:count
        
        Dens_Temp1=[Dens_Temp1;NetMes_B1_rand{i}{j}{6}];
        MClust_Temp1=[MClust_Temp1;NetMes_B1_rand{i}{j}{8}];
        Trans_Temp1=[Trans_Temp1;NetMes_B1_rand{i}{j}{9}];
        Assort_Temp1=[Assort_Temp1;NetMes_B1_rand{i}{j}{5}];
        GEff_Temp1=[GEff_Temp1;NetMes_B1_rand{i}{j}{10}];
        MLocEff_Temp1=[MLocEff_Temp1;NetMes_B1_rand{i}{j}{12}];
        Mod_Temp1=[Mod_Temp1;NetMes_B1_rand{i}{j}{13}];
        ModL_Temp1=[ModL_Temp1;NetMes_B1_rand{i}{j}{14}];
        PathL_Temp1=[PathL_Temp1;NetMes_B1_rand{i}{j}{15}];
        MNodeBetw_Temp1=[MNodeBetw_Temp1;NetMes_B1_rand{i}{j}{17}];
        MEdgeBetw_Temp1=[MEdgeBetw_Temp1;NetMes_B1_rand{i}{j}{19}];
        Lambda_Temp1=[Lambda_Temp1;NetMes_B1_rand{i}{j}{20}];
        Gamma_Temp1=[Gamma_Temp1;NetMes_B1_rand{i}{j}{21}];
        Sigma_Temp1=[Sigma_Temp1;NetMes_B1_rand{i}{j}{22}];
        
        Dens_Temp2=[Dens_Temp2;NetMes_B2_rand{i}{j}{6}];
        MClust_Temp2=[MClust_Temp2;NetMes_B2_rand{i}{j}{8}];
        Trans_Temp2=[Trans_Temp2;NetMes_B2_rand{i}{j}{9}];
        Assort_Temp2=[Assort_Temp2;NetMes_B2_rand{i}{j}{5}];
        GEff_Temp2=[GEff_Temp2;NetMes_B2_rand{i}{j}{10}];
        MLocEff_Temp2=[MLocEff_Temp2;NetMes_B2_rand{i}{j}{12}];
        Mod_Temp2=[Mod_Temp2;NetMes_B2_rand{i}{j}{13}];
        ModL_Temp2=[ModL_Temp2;NetMes_B2_rand{i}{j}{14}];
        PathL_Temp2=[PathL_Temp2;NetMes_B2_rand{i}{j}{15}];
        MNodeBetw_Temp2=[MNodeBetw_Temp2;NetMes_B2_rand{i}{j}{17}];
        MEdgeBetw_Temp2=[MEdgeBetw_Temp2;NetMes_B2_rand{i}{j}{19}];
        Lambda_Temp2=[Lambda_Temp2;NetMes_B2_rand{i}{j}{20}];
        Gamma_Temp2=[Gamma_Temp2;NetMes_B2_rand{i}{j}{21}];
        Sigma_Temp2=[Sigma_Temp2;NetMes_B2_rand{i}{j}{22}];
        
        if i<=1
            
            Dens_1=[Dens_1;NetMes_B1{j}{6}];
            MClust_1=[MClust_1;NetMes_B1{j}{8}];
            Trans_1=[Trans_1;NetMes_B1{j}{9}];
            Assort_1=[Assort_1;NetMes_B1{j}{5}];
            GEff_1=[GEff_1;NetMes_B1{j}{10}];
            MLocEff_1=[MLocEff_1;NetMes_B1{j}{12}];
            Mod_1=[Mod_1;NetMes_B1{j}{13}];
            ModL_1=[ModL_1;NetMes_B1{j}{14}];
            PathL_1=[PathL_1;NetMes_B1{j}{15}];
            MNodeBetw_1=[MNodeBetw_1;NetMes_B1{j}{17}];
            MEdgeBetw_1=[MEdgeBetw_1;NetMes_B1{j}{19}];
            Lambda_1=[Lambda_1;NetMes_B1{j}{20}];
            Gamma_1=[Gamma_1;NetMes_B1{j}{21}];
            Sigma_1=[Sigma_1;NetMes_B1{j}{22}];
            
            Dens_2=[Dens_2;NetMes_B2{j}{6}];
            MClust_2=[MClust_2;NetMes_B2{j}{8}];
            Trans_2=[Trans_2;NetMes_B2{j}{9}];
            Assort_2=[Assort_2;NetMes_B2{j}{5}];
            GEff_2=[GEff_2;NetMes_B2{j}{10}];
            MLocEff_2=[MLocEff_2;NetMes_B2{j}{12}];
            Mod_2=[Mod_2;NetMes_B2{j}{13}];
            ModL_2=[ModL_2;NetMes_B2{j}{14}];
            PathL_2=[PathL_2;NetMes_B2{j}{15}];
            MNodeBetw_2=[MNodeBetw_2;NetMes_B2{j}{17}];
            MEdgeBetw_2=[MEdgeBetw_2;NetMes_B2{j}{19}];
            Lambda_2=[Lambda_2;NetMes_B2{j}{20}];
            Gamma_2=[Gamma_2;NetMes_B2{j}{21}];
            Sigma_2=[Sigma_2;NetMes_B2{j}{22}];
        
        end
        
    end
    
    Dens_1_rand(:,i)=Dens_Temp1;
    MClust_1_rand(:,i)=MClust_Temp1;
    Trans_1_rand(:,i)=Trans_Temp1;
    Assort_1_rand(:,i)=Assort_Temp1;
    GEff_1_rand(:,i)=GEff_Temp1;
    MLocEff_1_rand(:,i)=MLocEff_Temp1;
    Mod_1_rand(:,i)=Mod_Temp1;
    ModL_1_rand(:,i)=ModL_Temp1;
    PathL_1_rand(:,i)=PathL_Temp1;
    MNodeBetw_1_rand(:,i)=MNodeBetw_Temp1;
    MEdgeBetw_1_rand(:,i)=MEdgeBetw_Temp1;
    Lambda_1_rand(:,i)=Lambda_Temp1;
    Gamma_1_rand(:,i)=Gamma_Temp1;
    Sigma_1_rand(:,i)=Sigma_Temp1;       
    
    Dens_2_rand(:,i)=Dens_Temp2;
    MClust_2_rand(:,i)=MClust_Temp2;
    Trans_2_rand(:,i)=Trans_Temp2;
    Assort_2_rand(:,i)=Assort_Temp2;
    GEff_2_rand(:,i)=GEff_Temp2;
    MLocEff_2_rand(:,i)=MLocEff_Temp2;
    Mod_2_rand(:,i)=Mod_Temp2;
    ModL_2_rand(:,i)=ModL_Temp2;
    PathL_2_rand(:,i)=PathL_Temp2;
    MNodeBetw_2_rand(:,i)=MNodeBetw_Temp2;
    MEdgeBetw_2_rand(:,i)=MEdgeBetw_Temp2;
    Lambda_2_rand(:,i)=Lambda_Temp2;
    Gamma_2_rand(:,i)=Gamma_Temp2;
    Sigma_2_rand(:,i)=Sigma_Temp2;
    
end

save(['NetMesPerDens_'  Group1],'Dens_1_rand','MClust_1_rand','MDeg_1_rand','MStr_1_rand',...
     'Trans_1_rand','Assort_1_rand','GEff_1_rand','MLocEff_1_rand','Mod_1_rand',...
     'ModL_1_rand','PathL_1_rand','MNodeBetw_1_rand','MEdgeBetw_1_rand',...
     'Lambda_1_rand','Gamma_1_rand','Sigma_1_rand');

save(['NetMesPerDens_'  Group2],'Dens_2_rand','MClust_2_rand','MDeg_2_rand','MStr_2_rand',...
     'Trans_2_rand','Assort_2_rand','GEff_2_rand','MLocEff_2_rand','Mod_2_rand',...
     'ModL_2_rand','PathL_2_rand','MNodeBetw_2_rand','MEdgeBetw_2_rand',...
     'Lambda_2_rand','Gamma_2_rand','Sigma_2_rand');

save(['NetMesPerDens_OrigNet_'  Group1],'Dens_1','MClust_1','MDeg_1','MStr_1',...
     'Trans_1','Assort_1','GEff_1','MLocEff_1','Mod_1',...
     'ModL_1','PathL_1','MNodeBetw_1','MEdgeBetw_1',...
     'Lambda_1','Gamma_1','Sigma_1');

save(['NetMesPerDens_OrigNet_'  Group2],'Dens_2','MClust_2','MDeg_2','MStr_2',...
     'Trans_2','Assort_2','GEff_2','MLocEff_2','Mod_2',...
     'ModL_2','PathL_2','MNodeBetw_2','MEdgeBetw_2',...
     'Lambda_2','Gamma_2','Sigma_2');

%%
 
[mu_MClust1_rand,sigma_MClust1_rand,muCi_MClust1_rand,sigmaCi_MClust1_rand]=normfit(MClust_1_rand',0.05);
[mu_Trans1_rand,sigma_Trans1_rand,muCi_Trans1_rand,sigmaCi_Trans1_rand]=normfit(Trans_1_rand',0.05);
[mu_Assort1_rand,sigma_Assort1_rand,muCi_Assort1_rand,sigmaCi_Assort1_rand]=normfit(Assort_1_rand',0.05);
[mu_GEff1_rand,sigma_GEff1_rand,muCi_GEff1_rand,sigmaCi_GEff1_rand]=normfit(GEff_1_rand',0.05);
[mu_MLocEff1_rand,sigma_MLocEff1_rand,muCi_MLocEff1_rand,sigmaCi_MLocEff1_rand]=normfit(MLocEff_1_rand',0.05);
[mu_Mod1_rand,sigma_Mod1_rand,muCi_Mod1_rand,sigmaCi_Mod1_rand]=normfit(Mod_1_rand',0.05);
[mu_ModL1_rand,sigma_ModL1_rand,muCi_ModL1_rand,sigmaCi_ModL1_rand]=normfit(ModL_1_rand',0.05);
[mu_PathL1_rand,sigma_PathL1_rand,muCi_PathL1_rand,sigmaCi_PathL1_rand]=normfit(PathL_1_rand',0.05);
[mu_MNodeBetw1_rand,sigma_MNodeBetw1_rand,muCi_MNodeBetw1_rand,sigmaCi_MNodeBetw1_rand]=normfit(MNodeBetw_1_rand',0.05);
[mu_MEdgeBetw1_rand,sigma_MEdgeBetw1_rand,muCi_MEdgeBetw1_rand,sigmaCi_MEdgeBetw1_rand]=normfit(MEdgeBetw_1_rand',0.05);
[mu_Lambda1_rand,sigma_Lambda1_rand,muCi_Lambda1_rand,sigmaCi_Lambda1_rand]=normfit(Lambda_1_rand',0.05);
[mu_Gamma1_rand,sigma_Gamma1_rand,muCi_Gamma1_rand,sigmaCi_Gamma1_rand]=normfit(Gamma_1_rand',0.05);
[mu_Sigma1_rand,sigma_Sigma1_rand,muCi_Sigma1_rand,sigmaCi_Sigma1_rand]=normfit(Sigma_1_rand',0.05);

[mu_MClust2_rand,sigma_MClust2_rand,muCi_MClust2_rand,sigmaCi_MClust2_rand]=normfit(MClust_2_rand',0.05);
[mu_Trans2_rand,sigma_Trans2_rand,muCi_Trans2_rand,sigmaCi_Trans2_rand]=normfit(Trans_2_rand',0.05);
[mu_Assort2_rand,sigma_Assort2_rand,muCi_Assort2_rand,sigmaCi_Assort2_rand]=normfit(Assort_2_rand',0.05);
[mu_GEff2_rand,sigma_GEff2_rand,muCi_GEff2_rand,sigmaCi_GEff2_rand]=normfit(GEff_2_rand',0.05);
[mu_MLocEff2_rand,sigma_MLocEff2_rand,muCi_MLocEff2_rand,sigmaCi_MLocEff2_rand]=normfit(MLocEff_2_rand',0.05);
[mu_Mod2_rand,sigma_Mod2_rand,muCi_Mod2_rand,sigmaCi_Mod2_rand]=normfit(Mod_2_rand',0.05);
[mu_ModL2_rand,sigma_ModL2_rand,muCi_ModL2_rand,sigmaCi_ModL2_rand]=normfit(ModL_2_rand',0.05);
[mu_PathL2_rand,sigma_PathL2_rand,muCi_PathL2_rand,sigmaCi_PathL2_rand]=normfit(PathL_2_rand',0.05);
[mu_MNodeBetw2_rand,sigma_MNodeBetw2_rand,muCi_MNodeBetw2_rand,sigmaCi_MNodeBetw2_rand]=normfit(MNodeBetw_2_rand',0.05);
[mu_MEdgeBetw2_rand,sigma_MEdgeBetw2_rand,muCi_MEdgeBetw2_rand,sigmaCi_MEdgeBetw2_rand]=normfit(MEdgeBetw_2_rand',0.05);
[mu_Lambda2_rand,sigma_Lambda2_rand,muCi_Lambda2_rand,sigmaCi_Lambda2_rand]=normfit(Lambda_2_rand',0.05);
[mu_Gamma2_rand,sigma_Gamma2_rand,muCi_Gamma2_rand,sigmaCi_Gamma2_rand]=normfit(Gamma_2_rand',0.05);
[mu_Sigma2_rand,sigma_Sigma2_rand,muCi_Sigma2_rand,sigmaCi_Sigma2_rand]=normfit(Sigma_2_rand',0.05);


N_rand = size(rand_Input1,2);
Pvalue = Alpha;


Ci_MClust=CL_per(MClust_2_rand'-MClust_1_rand',Pvalue);
Ci_Trans=CL_per(Trans_2_rand'-Trans_1_rand',Pvalue);
Ci_Assort=CL_per(Assort_2_rand'-Assort_1_rand',Pvalue);
Ci_GEff=CL_per(GEff_2_rand'-GEff_1_rand',Pvalue);
Ci_MLocEff=CL_per(MLocEff_2_rand'-MLocEff_1_rand',Pvalue);
Ci_Mod=CL_per(Mod_2_rand'-Mod_1_rand',Pvalue);
Ci_ModL=CL_per(ModL_2_rand'-ModL_1_rand',Pvalue);
Ci_PathL=CL_per(PathL_2_rand'-PathL_1_rand',Pvalue);
Ci_MNodeBetw=CL_per(MNodeBetw_2_rand'-MNodeBetw_1_rand',Pvalue);
Ci_MEdgeBetw=CL_per(MEdgeBetw_2_rand'-MEdgeBetw_1_rand',Pvalue);
Ci_Lambda=CL_per(Lambda_2_rand'-Lambda_1_rand',Pvalue);
Ci_Gamma=CL_per(Gamma_2_rand'-Gamma_1_rand',Pvalue);
Ci_Sigma=CL_per(Sigma_2_rand'-Sigma_1_rand',Pvalue);


% p_assort=CL_Pval(Assort_2_rand-Assort_1_rand,Assort_2-Assort_1,'Assort');
% p_GEff=CL_Pval(GEff_2_rand-GEff_1_rand,GEff_2-GEff_1,'GEff');
% p_Gamma=CL_Pval(Gamma_2_rand-Gamma_1_rand,Gamma_2-Gamma_1,'Gamma');
% p_Lambda=CL_Pval(Lambda_2_rand-Lambda_1_rand,Lambda_2-Lambda_1,'Lambda');
% p_MClust=CL_Pval(MClust_2_rand-MClust_1_rand,MClust_2-MClust_1,'MClust');
% p_MEdgeBetw=CL_Pval(MEdgeBetw_2_rand-MEdgeBetw_1_rand,MEdgeBetw_2-MEdgeBetw_1,'MEdgeBetw');
% p_MLocEff=CL_Pval(MLocEff_2_rand-MLocEff_1_rand,MLocEff_2-MLocEff_1,'MLocEff');
% p_MNodeBetw=CL_Pval(MNodeBetw_2_rand-MNodeBetw_1_rand,MNodeBetw_2-MNodeBetw_1,'MNodeBetw');
% p_ModL=CL_Pval(ModL_2_rand-ModL_1_rand,ModL_2-ModL_1,'ModL');
% p_Mod=CL_Pval(Mod_2_rand-Mod_1_rand,Mod_2-Mod_1,'Mod');
% p_PathL=CL_Pval(PathL_2_rand-PathL_1_rand,PathL_2-PathL_1,'PathL');
% p_Sigma=CL_Pval(Sigma_2_rand-Sigma_1_rand,Sigma_2-Sigma_1,'Sigma');
% p_Trans=CL_Pval(Trans_2_rand-Trans_1_rand,Trans_2-Trans_1,'Trans');
% p_Dens=CL_Pval(Dens_2_rand-Dens_1_rand,Dens_2-Dens_1,'Dens');
% p_MDeg=CL_Pval(MDeg_2_rand-MDeg_1_rand,MDeg_2-MDeg_1,'MDeg');
% p_MStr=CL_Pval(MStr_2_rand-MStr_1_rand,MStr_2-MStr_1,'MStr');

p_assort=CL_Pval(Assort_2_rand-Assort_1_rand,Assort_2-Assort_1,'Assort',Tail);
p_GEff=CL_Pval(GEff_2_rand-GEff_1_rand,GEff_2-GEff_1,'GEff',Tail);
p_Gamma=CL_Pval(Gamma_2_rand-Gamma_1_rand,Gamma_2-Gamma_1,'Gamma',Tail);
p_Lambda=CL_Pval(Lambda_2_rand-Lambda_1_rand,Lambda_2-Lambda_1,'Lambda',Tail);
p_MClust=CL_Pval(MClust_2_rand-MClust_1_rand,MClust_2-MClust_1,'MClust',Tail);
p_MEdgeBetw=CL_Pval(MEdgeBetw_2_rand-MEdgeBetw_1_rand,MEdgeBetw_2-MEdgeBetw_1,'MEdgeBetw',Tail);
p_MLocEff=CL_Pval(MLocEff_2_rand-MLocEff_1_rand,MLocEff_2-MLocEff_1,'MLocEff',Tail);
p_MNodeBetw=CL_Pval(MNodeBetw_2_rand-MNodeBetw_1_rand,MNodeBetw_2-MNodeBetw_1,'MNodeBetw',Tail);
p_ModL=CL_Pval(ModL_2_rand-ModL_1_rand,ModL_2-ModL_1,'ModL',Tail);
p_Mod=CL_Pval(Mod_2_rand-Mod_1_rand,Mod_2-Mod_1,'Mod',Tail);
p_PathL=CL_Pval(PathL_2_rand-PathL_1_rand,PathL_2-PathL_1,'PathL',Tail);
p_Sigma=CL_Pval(Sigma_2_rand-Sigma_1_rand,Sigma_2-Sigma_1,'Sigma',Tail);
p_Trans=CL_Pval(Trans_2_rand-Trans_1_rand,Trans_2-Trans_1,'Trans',Tail);
p_Dens=CL_Pval(Dens_2_rand-Dens_1_rand,Dens_2-Dens_1,'Dens',Tail);
p_MDeg=CL_Pval(MDeg_2_rand-MDeg_1_rand,MDeg_2-MDeg_1,'MDeg',Tail);
p_MStr=CL_Pval(MStr_2_rand-MStr_1_rand,MStr_2-MStr_1,'MStr',Tail);


%% plotting difference in net measures

Xax=[MinMesPlot:MesStepPlot:MaxMesPlot];
if ~isequal(size(Xax,2),size(MClust_1,1))
    
    Xax = Xax(1:size(MClust_1,1));
    
end
save DensityInterval Xax


if isequal(FigOn,1)

    %% Plotting difference in net measures between groups
    %clustering
    figure
    plot(Xax,MClust_2'-MClust_1','r*');hold 
    plot(Xax,mu_MClust2_rand-mu_MClust1_rand,'bx');%mean random difference
    plot(Xax,Ci_MClust(1,:)','b--');%mean random difference
    plot(Xax,Ci_MClust(2,:)','b--');%mean random difference
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in clustering coefficient','fontsize',12,'fontweight','b')
    title('Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('Clustering_vs_Null_Nperm.fig')


    %Transitivity
    figure
    plot(Xax,Trans_2'-Trans_1','r*');hold 
    plot(Xax,mu_Trans2_rand-mu_Trans1_rand,'bx');%mean random difference
    plot(Xax,Ci_Trans(1,:)','b--');%mean random difference
    plot(Xax,Ci_Trans(2,:)','b--');%mean random difference
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in transitivity coefficient','fontsize',12,'fontweight','b')
    title('transitivity','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('Transitivity_vs_Null_Nperm.fig')

    %Path length
    figure
    plot(Xax,PathL_2'-PathL_1','r*');hold 
    plot(Xax,mu_PathL2_rand-mu_PathL1_rand,'bx');%mean random difference
    plot(Xax,Ci_PathL(1,:)','b--');%mean random difference
    plot(Xax,Ci_PathL(2,:)','b--');%mean random difference
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in path length','fontsize',12,'fontweight','b')
    title('Char Path length','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('CharPath_vs_Null_Nperm.fig')

    %Assortativity
    figure
    plot(Xax,Assort_2'-Assort_1','r*');hold 
    plot(Xax,mu_Assort2_rand-mu_Assort1_rand,'bx');%mean random difference
    plot(Xax,Ci_Assort(1,:)','b--');%mean random difference
    plot(Xax,Ci_Assort(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in assortativity','fontsize',12,'fontweight','b')
    title('Mean assortativity','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('MeanAssortativity_vs_Null_Nperm.fig')



    %Global Efficciency
    figure
    plot(Xax,GEff_2'-GEff_1','r*');hold 
    plot(Xax,mu_GEff2_rand-mu_GEff1_rand,'bx');%mean random difference
    plot(Xax,Ci_GEff(1,:)','b--');%mean random difference
    plot(Xax,Ci_GEff(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in global efficiency','fontsize',12,'fontweight','b')
    title('Global efficiency','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('GlobalEfficiecny_vs_Null_Nperm.fig')


    %Mean Local Efficciency
    figure
    plot(Xax,MLocEff_2'-MLocEff_1','r*');hold 
    plot(Xax,mu_MLocEff2_rand-mu_MLocEff1_rand,'bx');%mean random difference
    plot(Xax,Ci_MLocEff(1,:)','b--');%mean random difference
    plot(Xax,Ci_MLocEff(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in mean local efficiency','fontsize',12,'fontweight','b')
    title('Mean local efficiency','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('MeanLocalEfficiecny_vs_Null_Nperm.fig')



    %Modularity
    figure
    plot(Xax,Mod_2'-Mod_1','r*');hold 
    plot(Xax,mu_Mod2_rand-mu_Mod1_rand,'bx');%mean random difference
    plot(Xax,Ci_Mod(1,:)','b--');%mean random difference
    plot(Xax,Ci_Mod(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in modularity','fontsize',12,'fontweight','b')
    title('Modularity','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('Modularity_vs_Null_Nperm.fig')


    %Modularity_L
    figure
    plot(Xax,ModL_2'-ModL_1','r*');hold 
    plot(Xax,mu_ModL2_rand-mu_ModL1_rand,'bx');%mean random difference
    plot(Xax,Ci_ModL(1,:)','b--');%mean random difference
    plot(Xax,Ci_ModL(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in modularity_L','fontsize',12,'fontweight','b')
    title('Modularity_L','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('ModularityL_vs_Null_Nperm.fig')



    %Mean Node Betweenness
    figure
    plot(Xax,MNodeBetw_2'-MNodeBetw_1','r*');hold 
    plot(Xax,mu_MNodeBetw2_rand-mu_MNodeBetw1_rand,'bx');%mean random difference
    plot(Xax,Ci_MNodeBetw(1,:)','b--');%mean random difference
    plot(Xax,Ci_MNodeBetw(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in Mean Node Betweenness','fontsize',12,'fontweight','b')
    title('Mean Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('MeanNodeBetweenness_vs_Null_Nperm.fig')



    %Mean Edge Betweenness
    figure
    plot(Xax,MEdgeBetw_2'-MEdgeBetw_1','r*');hold 
    plot(Xax,mu_MEdgeBetw2_rand-mu_MEdgeBetw1_rand,'bx');%mean random difference
    plot(Xax,Ci_MEdgeBetw(1,:)','b--');%mean random difference
    plot(Xax,Ci_MEdgeBetw(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in Mean Edge Betweenness','fontsize',12,'fontweight','b')
    title('Mean Edge Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('MeanEdgeBetweenness_vs_Null_Nperm.fig')


    %Lambda
    figure
    plot(Xax,Lambda_2'-Lambda_1','r*');hold 
    plot(Xax,mu_Lambda2_rand-mu_Lambda1_rand,'bx');%mean random difference
    plot(Xax,Ci_Lambda(1,:)','b--');%mean random difference
    plot(Xax,Ci_Lambda(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in Lambda','fontsize',12,'fontweight','b')
    title('Lambda','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('Lambda_vs_Null_Nperm.fig')

    %Gamma
    figure
    plot(Xax,Gamma_2'-Gamma_1','r*');hold 
    plot(Xax,mu_Gamma2_rand-mu_Gamma1_rand,'bx');%mean random difference
    plot(Xax,Ci_Gamma(1,:)','b--');%mean random difference
    plot(Xax,Ci_Gamma(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in Gamma','fontsize',12,'fontweight','b')
    title('Gamma','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('Gamma_vs_Null_Nperm.fig')


    %Sigma
    figure
    plot(Xax,Sigma_2'-Sigma_1','r*');hold 
    plot(Xax,mu_Sigma2_rand-mu_Sigma1_rand,'bx');%mean random difference
    plot(Xax,Ci_Sigma(1,:)','b--');%mean random difference
    plot(Xax,Ci_Sigma(2,:)','b--');%mean random difference%plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Difference in Sigma','fontsize',12,'fontweight','b')
    title('Sigma','fontsize',14,'fontweight','b','fontangle','italic')
    legend([Group2 ' vs. ' Group1],'Null (mean)','Null (upper bound 95% CI)','Null (lower bound 95% CI)')
    grid on
    hgsave('Sigma_vs_Null_Nperm.fig')


    %% Plotting net measures for each group
    
    %clustering
    figure
    plot(Xax,MClust_1','b*');hold 
    plot(Xax,MClust_2','r*');
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Mean clustering coefficient','fontsize',12,'fontweight','b')
    title('Clustering coefficient','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('Clustering_Value_Nperm.fig')

    %Transitiviy
    figure
    plot(Xax,Trans_1','b*');hold 
    plot(Xax,Trans_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Transitivity','fontsize',12,'fontweight','b')
    title('Transitivity','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('Transitivity_Value_Nperm.fig')


    %Assortativiy
    figure
    plot(Xax,Assort_1','b*');hold 
    plot(Xax,Assort_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Assortativity','fontsize',12,'fontweight','b')
    title('Assortativity','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('Assortativity_Value_Nperm.fig')


    %Global Efficiency
    figure
    plot(Xax,GEff_1','b*');hold 
    plot(Xax,GEff_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Global Efficiency','fontsize',12,'fontweight','b')
    title('Global Efficicency','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('GEfficiency_Value_Nperm.fig')


    %Local Efficciency
    figure
    plot(Xax,MLocEff_1','b*');hold 
    plot(Xax,MLocEff_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Local Efficiency','fontsize',12,'fontweight','b')
    title('Local Efficicency','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('LocalEfficiency_Value_Nperm.fig')


    %Modularity
    figure
    plot(Xax,Mod_1','b*');hold 
    plot(Xax,Mod_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Modularity','fontsize',12,'fontweight','b')
    title('Modularity','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('Modularity_Value_Nperm.fig')


    %Modularity_L
    figure
    plot(Xax,ModL_1','b*');hold 
    plot(Xax,ModL_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Modularity_L','fontsize',12,'fontweight','b')
    title('Modularity_L','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('ModularityL_Value_Nperm.fig')

    %Char path length
    figure
    plot(Xax,PathL_1','b*');hold 
    plot(Xax,PathL_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Characteristic path length','fontsize',12,'fontweight','b')
    title('Characteristic path length','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('CharPathLength_Value_Nperm.fig')


    %Mean Node bEtweenness Centrality
    figure
    plot(Xax,MNodeBetw_1','b*');hold 
    plot(Xax,MNodeBetw_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Mean Node Betweenness','fontsize',12,'fontweight','b')
    title('Mean Node Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('MNodeBetweeness_Value_Nperm.fig')

    %Mean Edge bEtweenness Centrality
    figure
    plot(Xax,MEdgeBetw_1','b*');hold 
    plot(Xax,MEdgeBetw_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Mean Edge Betweenness','fontsize',12,'fontweight','b')
    title('Mean Edge Betweenness','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('MEdgeBetweeness_Value_Nperm.fig')

    %Lambda
    figure
    plot(Xax,Lambda_1','b*');hold 
    plot(Xax,Lambda_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Lambda','fontsize',12,'fontweight','b')
    title('Lambda','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('Lambda_Value_Nperm.fig')



    %Gamma
    figure
    plot(Xax,Gamma_1','b*');hold 
    plot(Xax,Gamma_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Gamma','fontsize',12,'fontweight','b')
    title('Gamma','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('Gamma_Value_Nperm.fig')



    %Sigma
    figure
    plot(Xax,Sigma_1','b*');hold 
    plot(Xax,Sigma_2','r*');
    %plot properties
    xlabel('Density','fontsize',12,'fontweight','b')
    ylabel('Sigma','fontsize',12,'fontweight','b')
    title('Sigma','fontsize',14,'fontweight','b','fontangle','italic')
    legend(Group1,Group2)
    grid on
    hgsave('Sigma_Value_Nperm.fig')

end


fprintf('%-4s\n',' ..........Done ');





