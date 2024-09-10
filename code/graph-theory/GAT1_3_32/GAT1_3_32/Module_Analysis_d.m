function [Com_g1,Coml_g1,Com_g2,Coml_g2] = Module_Analysis_d(Group1,Group2,MinMesPlot,MaxMesPlot,MesStepPlot,NetMes1,NetMes2,ROI,nIter)

%%

GUIinput=0;

if nargin<1

    data_log = spm_select(1,'mat4GATd*','Select mat4GATd .mat (output from Load Data)');
    ff = load(data_log,'mat4GATd');
    Group1 = ff.mat4GATd.g1;Group2 = ff.mat4GATd.g2;ROI = ff.mat4GATd.roi1;
    MinMesPlot=ff.mat4GATd.MinThr;
    MaxMesPlot=ff.mat4GATd.MaxThr;
    MesStepPlot=ff.mat4GATd.MesStep;
    
    Cov = ff.mat4GATd.flagCov;

    if isequal(Cov,0)

        D1 = spm_select(1,['NetMesBin_D_' Group1 '_FinalThrRange'],['Select ' Group1 ' Unadjusted NetMesBin_D_*.mat file that you want to examine']);
        D2 = spm_select(1,['NetMesBin_D_' Group2 '_FinalThrRange'],['Select ' Group2 ' Unadjusted NetMesBin_D_*.mat file that you want to examine']);

    else
        
        D1 = spm_select(1,['NetMesBin_D_Adjusted' Group1 '_FinalThrRange'],['Select ' Group1 ' Adjusted NetMesBin_D_*.mat file that you want to examine']);
        D2 = spm_select(1,['NetMesBin_D_Adjusted' Group2 '_FinalThrRange'],['Select ' Group2 ' Adjusted NetMesBin_D_*.mat file that you want to examine']);

    end
        
        
    if isempty(D1) || isempty(D2)
    
        error('no input identified')
        
    end
    
    nIter = input('number of iteration (e.g. 100): ');
    
    f1=load(D1,'NetMes_Bin');NetMes1=f1.NetMes_Bin;
    f2=load(D2,'NetMes_Bin');NetMes2=f2.NetMes_Bin;
    
       
    GUIinput=1;


end


MidThr=input('target density:');
Xax = [MinMesPlot:MesStepPlot:MaxMesPlot];
MidIdx = find(single(Xax) == single(MidThr));

if isempty(MidIdx)
    
    error('no such target density');
    
end


%% 
j = MidIdx;
count =0;Q1=[];
for i = 1:size(NetMes1,1)
    
        count = count +1; 
        Q1(:,:,count) = NetMes1{i,j}{24,1};

end
O1=mean(Q1,3);


count =0;Q2=[];
for i = 1:size(NetMes2,1)
    
        count = count +1; 
        Q2(:,:,count) = NetMes2{i,j}{24,1};

end
O2=mean(Q2,3);


%% 

for i = 1:nIter
    
    [Cl_g1(i,:),Ql_g1(i)] = modularity_louvain_und(O1);
    Nmodl_g1(i) = max(Cl_g1(i,:));
    
    
    [Cl_g2(i,:),Ql_g2(i)] = modularity_louvain_und(O2);
    Nmodl_g2(i) = max(Cl_g2(i,:));

end

imaxl_g1 = find(Ql_g1==max(Ql_g1));
[All_Nmodl_g1, Idxl_g1] = max(Nmodl_g1(imaxl_g1));
Cl_g1 = Cl_g1(imaxl_g1(Idxl_g1),:);

imaxl_g2 = find(Ql_g2==max(Ql_g2));
[CommonModl_g2, Idxl_g2] = max(Nmodl_g2(imaxl_g2));
Cl_g2 = Cl_g2(imaxl_g2(Idxl_g2),:);

for i=1:max(Cl_g1)
    
    Coml_g1{i}=ROI(Cl_g1==i);
    
end

for i=1:max(Cl_g2)
    
    Coml_g2{i}=ROI(Cl_g2==i);
    
end

%% 

for i = 1:nIter
    
    
    [C_g1(i,:),Q_g1(i)] = modularity_und(O1);
    Nmod_g1(i) = max(C_g1(i,:));
    
    
    [C_g2(i,:),Q_g2(i)] = modularity_und(O2);
    Nmod_g2(i) = max(C_g2(i,:));

end

imax_g1 = find(Q_g1==max(Q_g1));
[All_Nmod_g1, Idx_g1] = max(Nmod_g1(imax_g1));
C_g1 = C_g1(imax_g1(Idx_g1),:);

imax_g2 = find(Q_g2==max(Q_g2));
[CommonMod_g2, Idx_g2] = max(Nmod_g2(imax_g2));
C_g2 = C_g2(imax_g2(Idx_g2),:);

for i=1:max(C_g1)
    
    Com_g1{i}=ROI(C_g1==i);
    
end

for i=1:max(C_g2)
    
    Com_g2{i}=ROI(C_g2==i);
    
end

save(['ModularCommunity_' Group1],'Com_g1','Coml_g1');
save(['ModularCommunity_' Group2],'Com_g2','Coml_g2');

fprintf('%-4s\n','.......DONE.......');
