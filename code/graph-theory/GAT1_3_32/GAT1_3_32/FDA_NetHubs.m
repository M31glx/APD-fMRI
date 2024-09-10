function FDA_NetHubs(Group1,Group2,fda_MDeg1n,fda_MNodeBetw1n,fda_MDeg2n,fda_MNodeBetw2n,ROI,SD)


%% 
GUIinput=0;

if nargin<1

    SD=input('please input the criterion(number) for calculating hubs: n for n*Standard Deviation greater than normal (2 is a good choice)  ');
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;
    Group2 = ff.mat4GAT.g2;
    ROI=ff.mat4GAT.roi1;
    Data1 = spm_select(1,['fda_NetMesReg_Norm' Group1 '.mat'],['Select the "fda_NetMesReg_Norm' Group1 '.mat" file ("output from fda for Regional Data")']);
    Data2 = spm_select(1,['fda_NetMesReg_Norm' Group2 '.mat'],['Select the "fda_NetMesReg_Norm' Group2 '.mat" file ("output from fda for Regional Data")']);
    
    load(Data1);load(Data2);
    
    GUIinput=1;%GUI based

end

%% 

if isequal(GUIinput,0)
    
    if nargin<5
        error('some input variables are missing')
    end
    
end

%%

Deg1=fda_MDeg1n;
Mdeg1=mean(Deg1);
SDdeg1=std(Deg1);
iDeg1=find(Deg1>=(Mdeg1+SDdeg1*SD));
HubInd_Deg1=iDeg1;

Deg2=fda_MDeg2n;
Mdeg2=mean(Deg2);
SDdeg2=std(Deg2);
iDeg2=find(Deg2>=(Mdeg2+SDdeg2*SD));
HubInd_Deg2=iDeg2;


Betw1=fda_MNodeBetw1n;
Mbetw1=mean(Betw1);
SDbetw1=std(Betw1);
iBetw1=find(Betw1>=(Mbetw1+SDbetw1*SD));
HubInd_Betw1=iBetw1;

Betw2=fda_MNodeBetw2n;
Mbetw2=mean(Betw2);
SDbetw2=std(Betw2);
iBetw2=find(Betw2>=(Mbetw2+SDbetw2*SD));
HubInd_Betw2=iBetw2;


%%

count=0;
HubNames_Deg1=[];

for i=HubInd_Deg1

    count=count+1;
    HubNames_Deg1{count}=strrep(ROI{i},['.rex.data.txt'],[]);
    
end

if isempty(HubNames_Deg1)

    HubNames_Deg1{1}='No Hubs';
    
end

HubNames_Deg1=HubNames_Deg1'
save(['fda_Net_Hubs_Deg_' Group1 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Deg1','HubInd_Deg1');


count=0;
HubNames_Deg2=[];

for i=HubInd_Deg2

    count=count+1;
    HubNames_Deg2{count}=strrep(ROI{i},['.rex.data.txt'],[]);
    
end

if isempty(HubNames_Deg2)

    HubNames_Deg2{1}='No Hubs';
    
end

HubNames_Deg2=HubNames_Deg2'
save(['fda_Net_Hubs_Deg_' Group2 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Deg2','HubInd_Deg2');



%% 

count=0;
HubNames_Betw1=[];

for i=HubInd_Betw1
    
    count=count+1;
    HubNames_Betw1{count}=strrep(ROI{i},['.rex.data.txt'],[]);
    
end

if isempty(HubNames_Betw1)

    HubNames_Betw1{1}='No Hubs';

end

HubNames_Betw1=HubNames_Betw1'
save(['fda_Net_Hubs_Betw_' Group1 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Betw1','HubInd_Betw1');

count=0;
HubNames_Betw2=[];

for i=HubInd_Betw2

    count=count+1;
    HubNames_Betw2{count}=strrep(ROI{i},['.rex.data.txt'],[]);
    
end

if isempty(HubNames_Betw2)

    HubNames_Betw2{1}='No Hubs';
    
end

HubNames_Betw2 = HubNames_Betw2'
save(['fda_Net_Hubs_Betw_' Group2 '_' strrep(num2str(SD),'.','p') 'SD'],'HubNames_Betw2','HubInd_Betw2');


fprintf('%-4s\n','.... done ....');

