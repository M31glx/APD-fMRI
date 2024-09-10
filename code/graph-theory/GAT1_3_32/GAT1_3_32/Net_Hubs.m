function [ HubNames,HubInd ] = Net_Hubs(NetMes,Group, Dens,ROIs,SD)

%%
% Inputs:
% NetMes: the output (NetMes_B1 or NetMes_B2) network measures from NetMEs_vs_Density...
% Dens: The NetMes_B1/B2 index corresponding to minimum density for fully connected graph (output from KAN_FixedDensity) 
% ROIs: name of the ROIs associated with each node (output regionName from rex)
% Param: The criterion used for defining hubs. available options are:
%       -'deg' for Degree
%       -'betw' for Node Betweenness
% SD: the threshold for accepting a node as hub. e.g SD=2 means that nodes
%    whose degrees (or betweenness) are at least 2SD greater than mean will be considered as hub

% Output:
% HubNames: Name of the hubs
% HubInd: index of the hubs (node number)

%%
%-- Hadi Hosseini: Created April 12, 2011
%-- Hadi Hosseini: Updated for GUI April 21,2011



%% 

GUIinput=0;

if nargin<1

    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;
    Group2 = ff.mat4GAT.g2;
    
    Data = spm_select(1,'KAN_FixedDens*','Select the KAN_FixedDens_Results.mat (output from "Graph Measures at Minimum Density")');
    
    e1 = load(Data,'Output1_Binary');O1 = e1.Output1_Binary;
    e2 = load(Data,'Output2_Binary');O2 = e2.Output2_Binary;
    ROIs=ff.mat4GAT.roi1;
    
    SD=input('please input the criterion(number) for calculating hubs: n for n*Standard Deviation greater than normal (2 is a good choice)  ');
    
    GUIinput=1;
    
end

%% 
if isequal(GUIinput,0)%command-based run
    if nargin<5
        error('some input variables are missing')
    end
end

%% 

g1 = load(Data,'NetMes_Bin1');Allreg1 = g1.NetMes_Bin1;Deg1=Allreg1{1,1};
g2 = load(Data,'NetMes_Bin2');Allreg2 = g2.NetMes_Bin2;Deg2=Allreg2{1,1};

Mdeg1=mean(Deg1);
Mdeg2=mean(Deg2);

SDdeg1=std(Deg1);
SDdeg2=std(Deg2);

iDeg1=find(Deg1>=(Mdeg1+SDdeg1*SD));
iDeg2=find(Deg2>=(Mdeg2+SDdeg2*SD));

HubInd_Deg1=iDeg1;
HubInd_Deg2=iDeg2;

%%

Betw1=Allreg1{16,1};
Betw2=Allreg2{16,1};

Mbetw1 = mean(Betw1);
Mbetw2 = mean(Betw2);

SDbetw1=std(Betw1);
SDbetw2=std(Betw2);

iBetw1=find(Betw1>=(Mbetw1+SDbetw1*SD));
iBetw2=find(Betw2>=(Mbetw2+SDbetw2*SD));

HubInd_Betw1=iBetw1;
HubInd_Betw2=iBetw2;

%% 

count=0;
HubNames_Deg1=[];

for i=HubInd_Deg1

    count=count+1;
    HubNames_Deg1{count}=strrep(ROIs{i},['.rex.data.txt'],'');
    
end

if isempty(HubNames_Deg1)
    
    HubNames_Deg1{1}='No Hubs';
    
end

HubNames_Deg1=HubNames_Deg1'

save(['Net_Hubs_Deg_' Group1 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Deg1','HubInd_Deg1');

%% 

count=0;
HubNames_Betw1=[];

for i=HubInd_Betw1
    
    count=count+1;
    HubNames_Betw1{count}=strrep(ROIs{i},['.rex.data.txt'],'');
    
end

if isempty(HubNames_Betw1)
    
    HubNames_Betw1{1}='No Hubs';
    
end

HubNames_Betw1=HubNames_Betw1'

save(['Net_Hubs_Betw_' Group1 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Betw1','HubInd_Betw1');


%% 
count=0;
HubNames_Deg2=[];

for i=HubInd_Deg2
    
    count=count+1;
    HubNames_Deg2{count}=strrep(ROIs{i},['.rex.data.txt'],'');
    
end

if isempty(HubNames_Deg2)
    
    HubNames_Deg2{1}='No Hubs';
    
end

HubNames_Deg2=HubNames_Deg2'

save(['Net_Hubs_Deg_' Group2 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Deg2','HubInd_Deg2');

%% 

count=0;
HubNames_Betw2=[];

for i=HubInd_Betw2
    
    count=count+1;
    HubNames_Betw2{count}=strrep(ROIs{i},'.rex.data.txt','');
    
end

if isempty(HubNames_Betw2)
    
    HubNames_Betw2{1}='No Hubs';
    
end

HubNames_Betw2=HubNames_Betw2'

save(['Net_Hubs_Betw_' Group2 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Betw2','HubInd_Betw2');

end
