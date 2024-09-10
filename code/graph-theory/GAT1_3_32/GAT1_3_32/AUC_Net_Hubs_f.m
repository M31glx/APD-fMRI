function AUC_Net_Hubs_f(Data1,Data2,Group1,Group2,ROIs,SD)

%%
% This function accept the 
% NetMes: the output (NetMes_B1 or NetMes_B2) network measures from NetMEs_vs_Density...
% Dens: The NetMes_B1/B2 index corresponding to minimum density for fully connected graph (output from KAN_FixedDensity) 
% ROIs: name of the ROIs associated with each node (output regionName from rex)
% Param: The criterion used for defining hubs. available options are:
%       -'deg' for Degree
%       -'betw' for Node Betweenness
%SD: the threshold for accepting a node as hub. e.g SD=2 means that nodes
%    whose degrees (or betweenness) are at least 2SD greater than mean will be considered as hub

%Output:
%HubNames: Name of the hubs
%HubInd: index of the hubs (node number)

%%
%-- Hadi Hosseini: Created April 12, 2011
%-- Hadi Hosseini: Updated for GUI April 21,2011
%-- Hadi Hosseini: Updated for functionals Sept. 13,2011

%%

GUIinput=0;

if nargin<1
    
    SD=input('please input the criterion(number) for calculating hubs: n for n*Standard Deviation greater than normal (2 is a good choice)  ');
    data_log = spm_select(1,'mat4GATf*','Select mat4GATf .mat (output from Network Measures)');
    ff = load(data_log,'mat4GATf');
    Group1 = ff.mat4GATf.g1;
    Group2 = ff.mat4GATf.g2;
    
    Data1 = spm_select(1,['AUC_NetMesReg_' Group1 '.mat'],['Select the "AUC_NetMesReg' Group1 '.mat" file ("output from AUC for Regional Data")']);
    Data2 = spm_select(1,['AUC_NetMesReg_' Group2 '.mat'],['Select the "AUC_NetMesReg' Group2 '.mat" file ("output from AUC for Regional Data")']);

    load(Data1);load(Data2);
    
    fprintf('%-4s\n',['loading inputs...']);
    ROIs = ff.mat4GATf.roi1;
    
    GUIinput=1;
    
end

%% 

if isequal(GUIinput,0)

    if nargin<6
    
        error('some input variables are missing')
        
    end

end


dd=pwd;     
mkdir('Hubs');
cd([dd '/Hubs']);

%% 

Mdeg1=mean(AUC_MDeg_1norm);
SDdeg1=std(AUC_MDeg_1norm);
iDeg1=find(AUC_MDeg_1norm >= (Mdeg1+SDdeg1*SD));
HubInd1=iDeg1;   

count=0;
HubNames=[];
for i=HubInd1
    count=count+1;
    HubNames{count}=ROIs{i};
end
if isempty(HubNames)
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Deg_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1');


Mbetw1=mean(AUC_MNodeBetw_1norm);
SDbetw1=std(AUC_MNodeBetw_1norm);
iBetw1=find(AUC_MNodeBetw_1norm >= (Mbetw1+SDbetw1*SD));
HubInd1=iBetw1;   

count=0;
HubNames=[];
for i=HubInd1
    
    count=count+1;
    HubNames{count}=ROIs{i};
    
end

if isempty(HubNames)
    
    HubNames{1}='No Hubs';
    
end

save(['Net_Hubs_Betw_' Group1 '_' num2str(SD) 'SD'],'HubNames','HubInd1');


%%

Mdeg2=mean(AUC_MDeg_2norm);
SDdeg2=std(AUC_MDeg_2norm);
iDeg2=find(AUC_MDeg_2norm >= (Mdeg2+SDdeg2*SD));
HubInd2=iDeg2;   

count=0;
HubNames=[];

for i=HubInd2
    
    count=count+1;
    HubNames{count}=ROIs{i};
    
end

if isempty(HubNames)

    HubNames{1}='No Hubs';
    
end

save(['Net_Hubs_Deg_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2');


%%
Mbetw2=mean(AUC_MNodeBetw_2norm);
SDbetw2=std(AUC_MNodeBetw_2norm);
iBetw2=find(AUC_MNodeBetw_2norm >= (Mbetw2+SDbetw2*SD));
HubInd2=iBetw2;   

count=0;
HubNames=[];

for i=HubInd2
  
    count=count+1;
    HubNames{count}=ROIs{i};
    
end

if isempty(HubNames)
    
    HubNames{1}='No Hubs';
    
end

save(['Net_Hubs_Betw_' Group2 '_' num2str(SD) 'SD'],'HubNames','HubInd2');


fprintf('%-4s\n','.... done ....');


end
