%function [p_Sigma, p_Gamma, p_Lambda] = PvalueAtExactMinThreshold(MinDens,NetMesB1,NetMesB2,Net1_rand,Net2_rand,ThreshDens1_rand,ThreshDens2_rand,NullType)
function PvalueAtExactMinThreshold

%% 

%-- Hadi Hosseini, April 2011


%%

GUIinput=0;

if nargin<1 

    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    
    if isfield(ff.mat4GAT,'null')
        NullType = ff.mat4GAT.null;
        Null_Iter = ff.mat4GAT.nullIter;       
    else
        NullType = 2;
        Null_Iter = 20;
    end
    
    if isfield(ff.mat4GAT,'tail') && isfield(ff.mat4GAT,'alpha') 
        
        Alpha = ff.mat4GAT.alpha;
        Tail = ff.mat4GAT.tail;
        
    else
        
        Alpha = .05;
        Tail = 2;
    
    end
    
    
    KAN = spm_select(1,'KAN_FixedDens*','Select KAN_FixedDens_Results (output from Graph Measures at Min Density)');
    RandNets = spm_select(1,'RandNetworksBootstrap*','Select Output from RandNet_Bootstrap');
    ThreshDensity1=spm_select(1,'ThreshDens_rand*','Select Group1 random ThreshDens_rand mat file (output "ThreshDens_rand_GroupName" from compare graphs at a range of densities)');
    ThreshDensity2=spm_select(1,'ThreshDens_rand*','Select Group2 random ThreshDens_rand mat file (output "ThreshDens_rand_GroupName" from compare graphs at a range of densities)');
    
    fprintf('%-4s\n',' loading variables.... ');
    f=load(KAN,'Density1');MinDens=f.Density1;
    f1=load(KAN,'Output1');Net1=f1.Output1;
    f2=load(KAN,'Output2');Net2=f2.Output2;
    c1=load(KAN,'NetMes_Bin1');NetMesB1=c1.NetMes_Bin1;
    c2=load(KAN,'NetMes_Bin2');NetMesB2=c2.NetMes_Bin2;
    
    
    ee1 = load(RandNets,'R1');Input1 = ee1.R1;clear f1;Inp1 = Input1;
    ee2 = load(RandNets,'R2');Input2 = ee2.R2;clear f2;Inp2 = Input2;
    e1=load(RandNets,'R_G1');Net1_rand=e1.R_G1;
    e2=load(RandNets,'R_G2');Net2_rand=e2.R_G2;
    
    d1=load(ThreshDensity1,'Thresh_Density1_rand');ThreshDens1_rand=d1.Thresh_Density1_rand;
    d2=load(ThreshDensity2,'Thresh_Density2_rand');ThreshDens2_rand=d2.Thresh_Density2_rand;
    
    
    if isequal(NullType,3) || isequal(NullType,4)
            
        q1 = load(RandNets,'dataG1');dr1 = q1.dataG1;clear q1
        q2 = load(RandNets,'dataG2');dr2 = q2.dataG2;clear q2
        
        %cxr1 = cellfun(@cov,dr1,'UniformOutput', false);clear dr1
        %cxr2 = cellfun(@cov,dr2,'UniformOutput', false);clear dr2
        
        %RandIter = 10;%computational limit

    end
    
    
    GUIinput=1;
end

if GUIinput == 0
    
    Null_Iter = 20;%defualt
    NullType = 2;
    
end


%%

Intv=10;

Net1_rand=cellfun(@zero_diag,Net1_rand,'UniformOutput', false);
Net2_rand=cellfun(@zero_diag,Net2_rand,'UniformOutput', false);

Net1_rand=cellfun(@neg_zero,Net1_rand,'UniformOutput', false);
Net2_rand=cellfun(@neg_zero,Net2_rand,'UniformOutput', false);          

SigmaR1=[];GammaR1=[];LambdaR1=[];ClustR1=[];CharPathR1=[];TransR1=[];
EdgeBetwR1=[];NodeBetwR1=[];ModLR1=[];LocEffR1=[];GEffR1=[];AssortR1=[];

SigmaR2=[];GammaR2=[];LambdaR2=[];ClustR2=[];CharPathR2=[];TransR2=[];
EdgeBetwR2=[];NodeBetwR2=[];ModLR2=[];LocEffR2=[];GEffR2=[];AssortR2=[];

fprintf('%-4s\n',' Thresholding random networks and calculating p-values at min density of full connectivity (may take several minutes)....  ');


for i=1:size(Net1_rand,2)

   TD1=ThreshDens1_rand{i};TD2=ThreshDens2_rand{i}; 
   Ind1=find(TD1(:,2)>=MinDens-Intv & TD1(:,2)<=MinDens+Intv);
   Ind2=find(TD2(:,2)>=MinDens-Intv & TD2(:,2)<=MinDens+Intv);
   
   Err1=abs(TD1(Ind1,2)-MinDens);Err2=abs(TD2(Ind2,2)-MinDens);
   [MIN1,Imin1]=min(Err1);[MIN2,Imin2]=min(Err2);
   T1=TD1(Ind1(Imin1),1);T2=TD2(Ind2(Imin2),1);
   
   N1=Net1_rand{i};N1(N1<=T1)=0;Net1_rand{i}=N1;
   N2=Net2_rand{i};N2(N2<=T2)=0;Net2_rand{i}=N2;
   
   if isequal(NullType,3) || isequal(NullType,4)
       
       NetMes_Bin1 = NetMeasures_Binary(N1,NullType,Null_Iter,dr1{i},Inp1);
       NetMes_Bin2 = NetMeasures_Binary(N2,NullType,Null_Iter,dr2{i},Inp2);
       
       NetMes_Bin1{23}=T1;
       NetMes_Bin2{23}=T2;
       
       NetMesRand_MinDens1{i}=NetMes_Bin1;
       NetMesRand_MinDens2{i}=NetMes_Bin2;  


   else

       NetMes_Bin1=NetMeasures_Binary(N1,NullType,Null_Iter);
       NetMes_Bin2=NetMeasures_Binary(N2,NullType,Null_Iter);
       
       NetMes_Bin1{23}=T1;
       NetMes_Bin2{23}=T2;
       
       NetMesRand_MinDens1{i}=NetMes_Bin1;
       NetMesRand_MinDens2{i}=NetMes_Bin2;  

    end 

   %%%
   SigmaR1=[SigmaR1;NetMes_Bin1{22}];GammaR1=[GammaR1;NetMes_Bin1{21}];
   LambdaR1=[LambdaR1;NetMes_Bin1{20}];ClustR1=[ClustR1;NetMes_Bin1{8}];
   CharPathR1=[CharPathR1;NetMes_Bin1{15}];TransR1=[TransR1;NetMes_Bin1{9}];
   EdgeBetwR1=[EdgeBetwR1;NetMes_Bin1{19}];NodeBetwR1=[NodeBetwR1;NetMes_Bin1{17}];
   ModLR1=[ModLR1;NetMes_Bin1{14}];LocEffR1=[LocEffR1;NetMes_Bin1{12}];
   GEffR1=[GEffR1;NetMes_Bin1{10}];AssortR1=[AssortR1;NetMes_Bin1{5}];
   
   SigmaR2=[SigmaR2;NetMes_Bin2{22}];GammaR2=[GammaR2;NetMes_Bin2{21}];
   LambdaR2=[LambdaR2;NetMes_Bin2{20}];ClustR2=[ClustR2;NetMes_Bin2{8}];
   CharPathR2=[CharPathR2;NetMes_Bin2{15}];TransR2=[TransR2;NetMes_Bin2{9}];
   EdgeBetwR2=[EdgeBetwR2;NetMes_Bin2{19}];NodeBetwR2=[NodeBetwR2;NetMes_Bin2{17}];
   ModLR2=[ModLR2;NetMes_Bin2{14}];LocEffR2=[LocEffR2;NetMes_Bin2{12}];
   GEffR2=[GEffR2;NetMes_Bin2{10}];AssortR2=[AssortR2;NetMes_Bin2{5}];
   
end


save('NetMesRand_MinDens1','NetMesRand_MinDens1');
save('NetMesRand_MinDens2','NetMesRand_MinDens2');
save('MinDensResults','SigmaR1','SigmaR2','GammaR1','GammaR2',...
    'LambdaR1','LambdaR2','ClustR1','ClustR2','CharPathR1',...
    'CharPathR2','TransR1','TransR2','EdgeBetwR1','EdgeBetwR2',...
    'NodeBetwR1','NodeBetwR2','ModLR1','ModLR2','LocEffR1','LocEffR2',...
    'GEffR1','GEffR2','AssortR1','AssortR2');
save('RandNets_ThreshAtMin','Net1_rand','Net2_rand');

%%

Diff_Sigma_R21=SigmaR2-SigmaR1;
Diff_Gamma_R21=GammaR2-GammaR1;
Diff_Lambda_R21=LambdaR2-LambdaR1;
Diff_Clust_R21=ClustR2-ClustR1;
Diff_CharPath_R21=CharPathR2-CharPathR1;
Diff_Trans_R21=TransR2-TransR1;
Diff_EdgeBetw_R21=EdgeBetwR2-EdgeBetwR1;
Diff_NodeBetw_R21=NodeBetwR2-NodeBetwR1;
Diff_ModL_R21=ModLR2-ModLR1;
Diff_LocEff_R21=LocEffR2-LocEffR1;
Diff_GEff_R21=GEffR2-GEffR1;
Diff_Assort_R21=AssortR2-AssortR1;

Sigma_21=NetMesB2{22}-NetMesB1{22};
Gamma_21=NetMesB2{21}-NetMesB1{21};
Lambda_21=NetMesB2{20}-NetMesB1{20};
Clust_21=NetMesB2{8}-NetMesB1{8};
CharPath_21=NetMesB2{15}-NetMesB1{15};
Trans_21=NetMesB2{9}-NetMesB1{9};
EdgeBetw_21=NetMesB2{19}-NetMesB1{19};
NodeBetw_21=NetMesB2{17}-NetMesB1{17};
ModL_21=NetMesB2{14}-NetMesB1{14};
LocEff_21=NetMesB2{12}-NetMesB1{12};
GEff_21=NetMesB2{10}-NetMesB1{10};
Assort_21=NetMesB2{5}-NetMesB1{5};

p_Sigma = CL_Pval(Diff_Sigma_R21',Sigma_21,'AtMinDens_Sigma',Tail);
p_Gamma = CL_Pval(Diff_Gamma_R21',Gamma_21,'AtMinDens_Gamma',Tail);
p_Lambda = CL_Pval(Diff_Lambda_R21',Lambda_21,'AtMinDens_Lambda',Tail);
p_Clust = CL_Pval(Diff_Clust_R21',Clust_21,'AtMinDens_Clust',Tail);
p_CharPath = CL_Pval(Diff_CharPath_R21',CharPath_21,'AtMinDens_CharPath',Tail);
p_Trans = CL_Pval(Diff_Trans_R21',Trans_21,'AtMinDens_Trans',Tail);

p_EdgeBetw =  CL_Pval(Diff_EdgeBetw_R21',EdgeBetw_21,'AtMinDens_EdgeBetw',Tail);
p_NodeBetw =  CL_Pval(Diff_NodeBetw_R21',NodeBetw_21,'AtMinDens_NodeBetw',Tail);
p_ModL =  CL_Pval(Diff_ModL_R21',ModL_21,'AtMinDens_ModL',Tail);
p_LocEff =  CL_Pval(Diff_LocEff_R21',LocEff_21,'AtMinDens_LocEff',Tail);
p_GEff =  CL_Pval(Diff_GEff_R21',GEff_21,'AtMinDens_GEff',Tail);
p_Assort =  CL_Pval(Diff_Assort_R21',Assort_21,'AtMinDens_Assort',Tail);

%% 
PvalCells = {'Sigma' p_Sigma;'Gamma' p_Gamma;'Lambda' p_Lambda;'Clust' p_Clust;...
             'CharPath' p_CharPath;'Trans' p_Trans;'EdgeBetw' p_EdgeBetw;'NodeBetw' p_NodeBetw;...
             'ModL' p_ModL;'LocEff' p_LocEff;'GEff' p_GEff;'Assort' p_Assort};

[nrows,ncols]= size(PvalCells);

filename = 'PvalsAtMinDens.txt';
fid = fopen(filename, 'w');

for row=1:nrows
    fprintf(fid, '%s,%d \n', PvalCells{row,:});
end

fclose(fid);


fprintf('%-4s\n',' ..........Done ');




