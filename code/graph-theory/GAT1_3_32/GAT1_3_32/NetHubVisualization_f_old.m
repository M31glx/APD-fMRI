function NetHubVisualization_f

%% AAL ROI dimension in MNI space (90 Cortical and subcortical regions excluding cerebbelar regions)
aal90 = [-24,-2,-18;27,-1,-19;-44,-62,34;45,-61,37;-7,-80,5;16,-74,8;-12,10,8;14,11,8;-4,34,13;8,36,14;-6,-16,40;8,-10,38;-5,-44,23;7,-43,20;-6,-81,26;13,-81,27;-49,11,18;50,14,20;-36,29,-13;41,31,-13;-46,29,13;50,29,13;-5,53,-9;8,50,-9;-34,31,34;-31,49,-11;33,51,-12;37,32,33;-19,33,41;-5,48,30;9,50,29;-17,46,-15;18,47,-15;22,30,43;-31,-41,-22;34,-40,-22;46,-18,9;-42,-20,9;-25,-22,-11;29,-21,-12;-35,5,2;39,5,1;-15,-69,-6;16,-68,-5;-36,-80,-9;38,-83,-9;-33,-82,15;37,-81,18;-17,-86,27;24,-82,29;-8,14,-13;10,15,-13;-18,-1,-1;21,-1,-1;-8,-27,69;7,-33,67;-21,-17,-22;25,-16,-22;-43,-47,45;46,-48,48;-24,-61,58;26,-60,61;-43,-24,47;41,-27,51;-39,-7,50;41,-10,51;-8,-57,47;10,-57,42;-24,3,1;27,4,1;-5,36,-20;8,34,-19;-47,-10,13;52,-8,13;-6,4,60;8,-1,61;-56,-35,29;57,-33,33;-50,-29,-25;53,-32,-24;-56,-35,-4;57,-39,-3;-37,13,-35;44,13,-34;-40,14,-21;48,13,-18;-53,-22,6;58,-23,5;-11,-19,7;13,-19,7;];%this is the original aal order which some of the rois violate the LR order e.g. MFG, SFG, ...
aal90_LR = [-24,-2,-18;27,-1,-19;-44,-62,34;45,-61,37;-7,-80,5;16,-74,8;-12,10,8;14,11,8;-4,34,13;8,36,14;-6,-16,40;8,-10,38;-5,-44,23;7,-43,20;-6,-81,26;13,-81,27;-49,11,18;50,14,20;-36,29,-13;41,31,-13;-46,29,13;50,29,13;-5,53,-9;8,50,-9;-34,31,34;37,32,33;-31,49,-11;33,51,-12;-19,33,41;22,30,43;-5,48,30;9,50,29;-17,46,-15;18,47,-15;-31,-41,-22;34,-40,-22;46,-18,9;-42,-20,9;-25,-22,-11;29,-21,-12;-35,5,2;39,5,1;-15,-69,-6;16,-68,-5;-36,-80,-9;38,-83,-9;-33,-82,15;37,-81,18;-17,-86,27;24,-82,29;-8,14,-13;10,15,-13;-18,-1,-1;21,-1,-1;-8,-27,69;7,-33,67;-21,-17,-22;25,-16,-22;-43,-47,45;46,-48,48;-24,-61,58;26,-60,61;-43,-24,47;41,-27,51;-39,-7,50;41,-10,51;-8,-57,47;10,-57,42;-24,3,1;27,4,1;-5,36,-20;8,34,-19;-47,-10,13;52,-8,13;-6,4,60;8,-1,61;-56,-35,29;57,-33,33;-50,-29,-25;53,-32,-24;-56,-35,-4;57,-39,-3;-37,13,-35;44,13,-34;-40,14,-21;48,13,-18;-53,-22,6;58,-23,5;-11,-19,7;13,-19,7;];%all the ROis come in LR order

aal_roi={'AMYG-L','AMYG-R','ANG-L','ANG-R','CALC-L','CALC-R','CN-L','CN-R','ACC-L','ACC-R','MCC-L','MCC-R','PCC-L','PCC-R','CUN-L','CUN-R','IFOp-L','IFOp-R','IFOr-L','IFOr-R','IFTr-L','IFTr-R','MedFOr-L','MedFOr-R','MFG-L','MFOr-L','MFOr-R','MFG-R','SFG-L','MedSF-L','MedSF-R','SFOr-L','SFOr-R','SFG-R','FG-L','FG-R','HSHL-L','HSHL-R','HIPP-L','HIPP-R','INS-L','INS-R','LNG-L','LNG-R','IOG-L','IOG-R','MOG-L','MOG-R','SOG-L','SOG-R','OFB-L','OFB-R','PLD-L','PLD-R','PCL-L','PCL-R','PHIP-L','PHIP-R','IPL-L','IPL-R','SPL-L','SPL-R','PoCG-L','PoCG-R','PrCG-L','PrCG-R','PCUN-L','PCUN-R','PUT-L','PUT-R','REC-L','REC-R','RLN-L','RLN-R','SMA-L','SMA-R','SMG-L','SMG-R','ITG-L','ITG-R','MTG-L','MTG-R','MTP-L','MTP-R','STP-L','STP-R','STG-L','STG-R','THL-L','THL-R';};
aal_roi_LR = {'AMYG-L','AMYG-R','ANG-L','ANG-R','CALC-L','CALC-R','CN-L','CN-R','ACC-L','ACC-R','MCC-L','MCC-R','PCC-L','PCC-R','CUN-L','CUN-R','IFOp-L','IFOp-R','IFOr-L','IFOr-R','IFTr-L','IFTr-R','MedFOr-L','MedFOr-R','MFG-L','MFG-R','MFOr-L','MFOr-R','SFG-L','SFG-R','MedSF-L','MedSF-R','SFOr-L','SFOr-R','FG-L','FG-R','HSHL-L','HSHL-R','HIPP-L','HIPP-R','INS-L','INS-R','LNG-L','LNG-R','IOG-L','IOG-R','MOG-L','MOG-R','SOG-L','SOG-R','OFB-L','OFB-R','PLD-L','PLD-R','PCL-L','PCL-R','PHIP-L','PHIP-R','IPL-L','IPL-R','SPL-L','SPL-R','PoCG-L','PoCG-R','PrCG-L','PrCG-R','PCUN-L','PCUN-R','PUT-L','PUT-R','REC-L','REC-R','RLN-L','RLN-R','SMA-L','SMA-R','SMG-L','SMG-R','ITG-L','ITG-R','MTG-L','MTG-R','MTP-L','MTP-R','STP-L','STP-R','STG-L','STG-R','THL-L','THL-R';};

%% freesurfer ROI dimension in MNI space (cortical, subcortical, and cerebellar regions)
fs86 = [-22,-64,-38;24,-64,-37;-15,9,7;-26,0,-2;-21,-5,-2;-12,-20,5;-26,-24,-15;-24,-7,-21;-10,11,-9;-11,-17,-10;16,9,8;27,1,-2;22,-5,-2;12,-19,6;28,-23,-15;25,-6,-21;10,10,-9;12,-16,-10;-54,-47,7;-5,19,28;-36,11,46;-7,-80,19;-23,-8,-34;-35,-43,-22;-40,-70,30;-49,-31,-27;-7,-47,19;-29,-89,-2;-24,30,-17;-14,-68,-5;-6,38,-16;-58,-28,-14;-23,-32,-19;-8,-29,57;-46,13,12;-42,40,-12;-45,30,4;-12,-80,6;-45,-23,42;-6,-19,38;-40,-10,43;-9,-58,35;-5,37,4;-31,45,14;-11,28,42;-22,-63,48;-53,-14,-5;-53,-38,31;-8,63,-10;-32,11,-35;-44,-24,7;-35,5,2;53,-41,6;6,21,28;36,12,46;8,-80,21;25,-9,-33;35,-41,-22;43,-64,30;51,-29,-27;7,-46,18;32,-87,-1;24,29,-18;15,-66,-5;7,37,-16;58,-25,-15;25,-31,-18;8,-27,56;47,13,12;42,40,-11;47,29,5;13,-79,6;44,-22,43;6,-18,38;40,-9,43;10,-57,35;6,36,5;32,46,14;12,30,41;22,-63,49;54,-11,-6;53,-33,31;10,62,-12;33,12,-35;45,-21,7;39,5,1;];
fs_roi={'L.Cerebellum';'R.Cerebellum';'L.Caudate';'L.Putamen';'L.Pallidum';'L.Thalamus';'L.Hippocampus';'L.Amygdala';'L.Accumbens';'L.VentralDC';'R.Thalamus';'R.Caudate';'R.Putamen';'R.Pallidum';'R.Hippocampus';'R.Amygdala';'R.Accumbens';'R.VentralDC';'L.bankssts';'L.caudalanteriorcingulate';'L.caudalmiddlefrontal';'L.cuneus';'L.entorhinal';'L.fusiform';'L.inferiorparietal';'L.inferiortemporal';'L.isthmuscingulate';'L.lateraloccipital';'L.lateralorbitofrontal';'L.lingual';'L.medialorbitofrontal';'L.middletemporal';'L.parahippocampal';'L.paracentral';'L.parsopercularis';'L.parsorbitalis';'L.parstriangularis';'L.pericalcarine';'L.postcentral';'L.posteriorcingulate';'L.precentral';'L.precuneus';'L.rostralanteriorcingulate';'L.rostralmiddlefrontal';'L.superiorfrontal';'L.superiorparietal';'L.superiortemporal';'L.supramarginal';'L.frontalpole';'L.temporalpole';'L.transversetemporal';'L.insula';'R.bankssts';'R.caudalanteriorcingulate';'R.caudalmiddlefrontal';'R.cuneus';'R.entorhinal';'R.fusiform';'R.inferiorparietal';'R.inferiortemporal';'R.isthmuscingulate';'R.lateraloccipital';'R.lateralorbitofrontal';'R.lingual';'R.medialorbitofrontal';'R.middletemporal';'R.parahippocampal';'R.paracentral';'R.parsopercularis';'R.parsorbitalis';'R.parstriangularis';'R.pericalcarine';'R.postcentral';'R.posteriorcingulate';'R.precentral';'R.precuneus';'R.rostralanteriorcingulate';'R.rostralmiddlefrontal';'R.superiorfrontal';'R.superiorparietal';'R.superiortemporal';'R.supramarginal';'R.frontalpole';'R.temporalpole';'R.transversetemporal';'R.insula';};

%% functional ROIs (Shirer 2011)
%func90=
%func_roi=

%%

roi_sch = input('which ROI scheme? type 1 for AAL-LR, 2 for AAL-Orig, 3 for freesurfer,  4 for manual input: ');
hubs = input('size of the nodes corresponds to betweenness(1), degree(2), clustering (3), none(4): '); 
data_log = spm_select(1,'mat4GATf*','Select mat4GATf .mat (output from Load Data)');
ff = load(data_log,'mat4GATf');
Group1 = ff.mat4GATf.g1;Group2 = ff.mat4GATf.g2;R=ff.mat4GATf.roi1;

MinMesPlot=ff.mat4GATf.MinThr;
MaxMesPlot=ff.mat4GATf.MaxThr;
MesStepPlot=ff.mat4GATf.MesStep;

D1 = spm_select(1,'NetMesBin_f',['Select ' Group1 ' Unadjusted NetMesBin_f_*.mat file that you want to examine']);
D2 = spm_select(1,'NetMesBin_f',['Select ' Group2 ' Unadjusted NetMesBin_f_*.mat file that you want to examine']);

if isempty(D1) || isempty(D2)
    error('no input identified')
end
f1=load(D1,'NetMes_Bin');NetMes1=f1.NetMes_Bin;
f2=load(D2,'NetMes_Bin');NetMes2=f2.NetMes_Bin;

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinThr=input('minimum threshold:');MinIdx=find(single(xxx)==single(MinThr));
MaxThr=input('maximum threshold:');MaxIdx=find(single(xxx)==single(MaxThr));
Xax = [MinThr:MesStepPlot:MaxThr];

if MaxIdx > size(NetMes1,2)
    
    MaxIdx = size(NetMes1,2);

end

MidIdx = floor(median(1:1:MaxIdx-MinIdx+1));

NetMes1=NetMes1(:,MinIdx:MaxIdx);
NetMes2=NetMes2(:,MinIdx:MaxIdx);


%%

count =0;Q1=[];
for i = 1:size(NetMes1,1)
    
    for j = 1:size(NetMes1,2)
        
        count = count +1; 
        Q1(:,:,count) = NetMes1{i,j}{24,1};
        
    end

end
O1=mean(Q1,3);


count =0;Q2=[];
for i = 1:size(NetMes2,1)
    
    for j = 1:size(NetMes2,2)
        
        count = count +1; 
        Q2(:,:,count) = NetMes2{i,j}{24,1};
        
    end

end
O2=mean(Q2,3);

%%

DENS1 = density_und(NetMes1{1,MidIdx}{24,1});
DENS2 = density_und(NetMes2{1,MidIdx}{24,1});
Precision = 0.02; 

Dens_OO1 = [];Dens_OO2 = [];

Range = .01:.01:.8;
for OL = Range
    
    OO1 = O1; OO1(OO1 < OL) = 0; Dens_OO1 = [Dens_OO1; density_und(OO1)]; 
    OO2 = O2; OO2(OO2 < OL) = 0; Dens_OO2 = [Dens_OO2; density_und(OO2)]; 
    
end   

%%

Ind1 = find(Dens_OO1 >= DENS1 - Precision & Dens_OO1 <= DENS1 + Precision );
Dind1=Dens_OO1(Ind1)- DENS1;
Ind_Ind1=find(Dind1==min(abs(Dind1)) | Dind1==-min(abs(Dind1)));
Ind1=Ind1(min(Ind_Ind1));

if isempty(Ind1)

    fprintf('%-4s\n',['decreasing the precision for Group ' Group1]);

    Ind1 = find(Dens_OO1 >= DENS1 - 10*Precision & Dens_OO1 <= DENS1 + 10*Precision );
    Dind1=Dens_OO1(Ind1)- DENS1;
    Ind_Ind1=find(Dind1==min(abs(Dind1)) | Dind1==-min(abs(Dind1)));
    Ind1=Ind1(min(Ind_Ind1));

end

O1(O1<Range(Ind1)) = 0; O1(O1>0) = 1;
fprintf('%-4s\n',['Percent Overlap for Group ' Group1 ' is ' num2str(Range(Ind1))]);


%%

Ind2 = find(Dens_OO2 >= DENS2 - Precision & Dens_OO2 <= DENS2 + Precision );
Dind2=Dens_OO2(Ind2)- DENS2;
Ind_Ind2=find(Dind2==min(abs(Dind2)) | Dind2==-min(abs(Dind2)));
Ind2=Ind2(min(Ind_Ind2));

if isempty(Ind2)

    fprintf('%-4s\n',['decreasing the precision for Group ' Group2]);

    Ind2 = find(Dens_OO2 >= DENS2 - 10*Precision & Dens_OO2 <= DENS2 + 10*Precision );
    Dind2=Dens_OO2(Ind2)- DENS2;
    Ind_Ind2=find(Dind2==min(abs(Dind2)) | Dind2==-min(abs(Dind2)));
    Ind2=Ind2(min(Ind_Ind2));

end

O2(O2<Range(Ind2)) = 0; O2(O2>0) = 1;
fprintf('%-4s\n',['Percent Overlap for Group ' Group2 ' is ' num2str(Range(Ind2))]);


%%

reg1=[];reg2=[];

for i=1:size(NetMes1,1)
    
    temp1=[];
    
    for j=1:size(NetMes1,2)
        
        
        switch hubs
        
            case 1
                
                temp1=[temp1;NetMes1{i,j}{16,3}];
            
            case 2
                
                temp1=[temp1;NetMes1{i,j}{1,3}];
                
            case 3
                
                temp1=[temp1;NetMes1{i,j}{7,3}'];
                
        end
        
    end
    
    reg1=[reg1;trapz(Xax,temp1)];

end

reg1=mean(reg1);


for i=1:size(NetMes2,1)
    
    temp2=[];
    
    for j=1:size(NetMes2,2)
        
        
        switch hubs
            
            case 1
                
                temp2=[temp2;NetMes2{i,j}{16,3}];
                
            case 2
                
                temp2=[temp2;NetMes2{i,j}{1,3}];
                
            case 3
                
                temp2=[temp2;NetMes2{i,j}{7,3}'];
                
        end
        
    end
    
    reg2=[reg2;trapz(Xax,temp2)];

end

reg2=mean(reg2);


%%

switch roi_sch

    case 1 
        
        nROI = size(aal90_LR,1);names=aal_roi_LR;
        [rois,ok] = listdlg('PromptString','Select ROIs:',...
                            'SelectionMode','multiple',...
                            'ListString',char(names),...
                            'initialvalue',1:nROI);
        if ~ok,return;end
        names = {names{rois}};xyz = aal90_LR(rois,:);nROI = length(rois);
    
    case 2
        
        nROI = size(aal90,1);names=aal_roi;
        [rois,ok] = listdlg('PromptString','Select ROIs:',...
                            'SelectionMode','multiple',...
                            'ListString',char(names),...
                            'initialvalue',1:nROI);
        if ~ok,return;end
        names = {names{rois}};xyz = aal90(rois,:);nROI = length(rois);
    
    case 3 
    
        nROI = size(fs86,1);names = fs_roi;
        [rois,ok]=listdlg('PromptString','Select ROIs:',...
                            'SelectionMode','multiple',...
                            'ListString',char(names),...
                            'initialvalue',1:nROI);
        if ~ok,return;end
        names = {names{rois}};xyz = fs86(rois,:);nROI = length(rois);
        
%     case 4 
%
%         nROI = size(func90,1);names=func_roi;
%         [rois,ok] = listdlg('PromptString','Select ROIs:',...
%                             'SelectionMode','multiple',...
%                             'ListString',char(names),...
%                             'initialvalue',1:nROI);
%         if ~ok,return;end
%         names = {names{rois}};xyz = func90(rois,:);nROI = length(rois);

    case 4
        
        nROI = size(R,2);names = R;
        
        [DimName,DimPathName] = uigetfile('*.xls','Select the XLS file that contains ROI dimensions');
        dimensions = xlsread([DimPathName DimName]);
        xyz = dimensions;
        
        
        
end


%% 

if ~isequal(hubs,4)

    reg1 = 0.1 + 9 * ((reg1 - min(reg1)) / (max(reg1) - min(reg1)));
    reg2 = 0.1 + 9 * ((reg2 - min(reg2)) / (max(reg2) - min(reg2)));
    Group1_1SD_thresh = mean(reg1) + std(reg1)
    Group1_2SD_thresh = mean(reg1) + 2*std(reg1)
    Group2_1SD_thresh = mean(reg2) + std(reg2)
    Group2_2SD_thresh = mean(reg2) + 2*std(reg2)
    save RegMes1 reg1
    save RegMes2 reg2

else
    
    reg1=ones(1,nROI);
    reg2=ones(1,nROI);
    
end



%% 

ModularYN = input('modular strcture? 1(yes), 2(no): ');

if isequal(ModularYN,1)

    Modules1 = spm_select(1,['ModularCommunity_' Group1 '*.mat'],['Select ' Group1 'ModularCommunity .mat (output from Load Data)']);
    Modules2 = spm_select(1,['ModularCommunity_' Group2 '*.mat'],['Select ' Group2 'ModularCommunity .mat (output from Load Data)']);
    ComType = input('type 1 for ComL, 2 for Com :');
    
    switch ComType
    
        case 1
        
            Mod1 = load(Modules1,'Coml_g1');Mod1=Mod1.Coml_g1;
            Mod2 = load(Modules2,'Coml_g2');Mod2=Mod2.Coml_g2;
        
        case 2
        
            Mod1 = load(Modules1,'Com_g1');Mod1=Mod1.Com_g1;
            Mod2 = load(Modules2,'Com_g2');Mod2=Mod2.Com_g2;
    
    end
    
    Mod_Index1 = [];
    Mod_Index2 = [];
    n1=size(Mod1,2);
    n2=size(Mod2,2);
    
    for i=1:nROI
    
        for j1=1:n1
        
            for k1=1:size(Mod1{j1},2)
            
                if isequal(Mod1{j1}{k1},R{i})
                
                    Mod_Index1=[Mod_Index1;j1];
                
                end
                
            end
            
        end
        
        for j2=1:n2
        
            for k2=1:size(Mod2{j2},2)
            
                if isequal(Mod2{j2}{k2},R{i})
                
                    Mod_Index2=[Mod_Index2;j2];
                    
                end
                
            end
            
        end
        
    end
end

%% 


dlmwrite(['Edge_' Group1 '.edge'],O1,'\t');dlmwrite(['Edge_' Group2 '.edge'],O2,'\t');

nodes1 = cell(nROI,6);nodes2 = nodes1;

for i=1:nROI

    for j=1:3
    
        nodes1{i,j} = xyz(i,j);nodes2{i,j} = xyz(i,j);
    
    end
    
    if isequal(ModularYN,1)
    
        if max(Mod_Index1)<7 || max(Mod_Index1)>11
        
            nodes1{i,4} = Mod_Index1(i);
        
        else
            
            NewColor=[12,15,16,20,13];
            nMod1 = max(Mod_Index1);
            
            for ii=7:nMod1
            
                Mod_Index1(Mod_Index1 == ii) = NewColor(ii-6);
            
            end
            
            nodes1{i,4} = Mod_Index1(i);
        
        end
        
        if max(Mod_Index2)<7 || max(Mod_Index2)>11
        
            nodes2{i,4} = Mod_Index2(i);
        
        else
            
            NewColor=[12,15,16,20,13];
            nMod2 = max(Mod_Index2);
            
            for ii=7:nMod2
            
                Mod_Index2(Mod_Index2 == ii) = NewColor(ii-6);
            
            end
            
            nodes2{i,4} = Mod_Index2(i);
        
        end
        
    else
        
        nodes1{i,4} = 1;nodes2{i,4} = 1;
    
    end
    
    nodes1{i,5} = reg1(i);nodes2{i,5} = reg2(i);
    nodes1{i,6} = names{i};nodes2{i,6} = names{i};

end

filename1 = ['Nodes_' Group1 '.node'];filename2 = ['Nodes_' Group2 '.node'];
fid1 = fopen(filename1, 'w');fid2 = fopen(filename2, 'w');

for row=1:nROI

    fprintf(fid1, '%d %d %d %d %d %s\n', nodes1{row,:});
    fprintf(fid2, '%d %d %d %d %d %s\n', nodes2{row,:});

end

fclose(fid1);fclose(fid2);
    
fprintf('%-4s\n','.......DONE.......');
    
    