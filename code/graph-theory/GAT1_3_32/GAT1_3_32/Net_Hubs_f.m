function Net_Hubs_f(NetMes,Group, Dens,ROIs,SD)

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

%% Input dialog  
GUIinput=0;
if nargin<1
    WhichGroup=input('network hubs for group1(type 1) or group2 (type 2)? '); 
    SD=input('please input the criterion(number) for calculating hubs: n for n*Standard Deviation greater than normal (2 is a good choice)  ');
    data_log = spm_select(1,'mat4GATf*','Select mat4GATf .mat (output from Network Measures)');
    ff = load(data_log,'mat4GATf');
    switch WhichGroup
        case 1
            Group = ff.mat4GATf.g1;%Group1=input('please type name of the 1st group in single quotation mark ');
        case 2
            Group = ff.mat4GATf.g2;
    end
    cov=ff.mat4GATf.flagCov;%1 for yes, 0 for no
    if isequal(cov,1)
        Data1 = spm_select(1,['NetMesBin_f_Adjusted_' Group '_FinalThrRange.mat'],'Select 1st group residuals (NetMesBin_f_adjusted.mat .mat file)');
    elseif isequal(cov,0)
        Data1 = spm_select(1,['NetMesBin_f_' Group '_FinalThrRange.mat'],'Select 1st group NetMesBin_f.mat file ');
    end
    %%Input conversion
    fprintf('%-4s\n',['loading inputs...']);
    f1 = load(Data1,'NetMes_Bin');NetMes = f1.NetMes_Bin;
    ROIs = ff.mat4GATf.roi1;
    MinMesPlot=ff.mat4GATf.MinThr;
    MaxMesPlot=ff.mat4GATf.MaxThr;
    MesStepPlot=ff.mat4GATf.MesStep;
    %%%%
    % select final density range
    MinThr=input('select your final decision on minimum threshold: ');
    MaxThr=input('select your final decision on maximum threshold: ');
    %%%
    GUIinput=1;%GUI based
end

%% 
if isequal(GUIinput,0)%command-based run
    if nargin<6
        error('some input variables are missing')
    end
    MinThr=input('select your final decision on minimum threshold:');
    MaxThr=input('select your final decision on maximum threshold:');

end

%% select final density range
MinIdx=floor((MinThr-MinMesPlot)/MesStepPlot+1);
MaxIdx=floor((MaxThr-MinMesPlot)/MesStepPlot+1);

%saving the final figs in a regional
dd=pwd;     
mkdir('Hubs');
cd([dd '/Hubs']);

%% unifying each group net by averaging net measures across final density range

%selecting the network measures across specified threshold
NetMes=NetMes(:,MinIdx:MaxIdx);
%%%
MClust=[];MDeg=[];MNodeBetw=[];MLocEff=[];
temp_clust=[];temp_deg=[];temp_nodeb=[];temp_leff=[];

for i=1:size(NetMes,1)%number of subjects in the group
    fprintf('%-4s\n',[' subject ' num2str(i) ' ....']);
    
    for j=1:size(NetMes,2)%thresholds (densities)
    
        temp_deg=[temp_deg;NetMes{i,j}{1,1}];%NThr*90
        temp_nodeb=[temp_nodeb;NetMes{i,j}{16,1}'];%NThr*90

    end

    MDeg=[MDeg;mean(temp_deg)];%Nsbj*90
    MNodeBetw=[MNodeBetw;mean(temp_nodeb)];%Nsbj*90

end

%original
MClust=mean(MClust);%1*90
MDeg=mean(MDeg);%1*90
MNodeBetw=mean(MNodeBetw);%1*90
MLocEff=mean(MLocEff);%1*90

Mdeg=mean(MDeg);%mean degree
SDdeg=std(MDeg);%SD of degree
iDeg=find(MDeg>=(Mdeg+SDdeg*SD));
HubInd=iDeg;   

%saving hubs
%original
count=0;
HubNames=[];
for i=HubInd
    count=count+1;
    HubNames{count}=ROIs{i};
end
if isempty(HubNames)
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Deg_' Group '_' num2str(SD) 'SD'],'HubNames','HubInd');

%% based on betw
%original
Mbetw=mean(MNodeBetw);%mean node betweenness
SDbetw=std(MNodeBetw);%SD of node betweenness
iBetw=find(MNodeBetw>=(Mbetw+SDbetw*SD));
HubInd=iBetw;


%original
count=0;
HubNames=[];
for i=HubInd
    count=count+1;
    HubNames{count}=ROIs{i};
end
if isempty(HubNames)
    HubNames{1}='No Hubs';
end
save(['Net_Hubs_Betw_' Group '_' num2str(SD) 'SD'],'HubNames','HubInd');


fprintf('%-4s\n','.... done ....');

    
end
