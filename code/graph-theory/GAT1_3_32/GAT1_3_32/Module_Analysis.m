function [Com_g1,Coml_g1,Com_g2,Coml_g2] = Module_Analysis(Group1,Group2,O1,O2,ROI,nIter)

%%  
GUIinput=0;
if nargin<1
    
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    Data = spm_select(1,'KAN_FixedDens*','Select the KAN_FixedDens_Results.mat (output from "Graph Measures at Minimum Density")');
    
    nIter = input('number of iteration (e.g. 100): ');
    
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2;ROI = ff.mat4GAT.roi1;
    e1 = load(Data,'Output1_Binary');O1 = e1.Output1_Binary;
    e2 = load(Data,'Output2_Binary');O2 = e2.Output2_Binary;

    GUIinput=1;
    
end


if isequal(GUIinput,0)
    
    if nargin<6
    
        error('some input variables are missing')
    
    end
    
end


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


end