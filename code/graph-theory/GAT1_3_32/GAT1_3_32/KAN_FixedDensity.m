%function [Thresh1,Thresh2,Density1,Density2,NetMes_Bin1,NetMes_Bin2] = KAN_FixedDensity(Input1,Input2,ROIs,func,Group1,Group2,OutputFName,RandIter,varargin)
function KAN_FixedDensity(Input1,Input2,ROIs,func,Group1,Group2,OutputFName,RandIter,varargin)

%  finds the minimum network density at which both networks are connected
%  (no fragmentation in networks) and compute the network measures theresholded at that
%   density

% Input:
%   Input1: first group regional gray matter data (output from "concatenate matrices")
%   Input2: second group regional gray matter data (output from "concatenate matrices")
%   ROIs: a cell array including ROI names (cell(1,nROI))
%   func: 1 for Pearson, 2 for partial correlation analysis
%   Group1: name of group1 (e.g. 'Group1')
%   Group2: name of group2 (e.g. 'Group2')
%   OutputFName: a string containing the output filename
%   RandIter: number of null networks (e.g. 20) 
%   varagin:
%       -optional input arguments are "covariates" for which you need to correct the data:
%       1.'covar1.mat': a .mat file ontaining 'covar1' array with the same number of rows as data1 whose columns contain different covariates(the order of the rows should be consistent with the order of subjects' data in Data1) 
%       2.'covar2.mat': a .mat file ontaining 'covar2' array with the same number of rows as data2 whose columns contain different covariates(the order of the rows should be consistent with the order of subjects' data in Data2)
%

% Output:
%   Thresh1/2: Correlation Threshold at minimum density for network 1/2 
%   Density1/2: minimum density for network 1/2 
%   NetMesBin1/2: network measures at minimm density 

%%
%-- Hadi Hosseini April 19, 2011
%-- updated on Apr 19,2011 for GUI



%%
%clear 
N_repeat=100000;
Precision=20;


%% Input dialog  
GUIinput=0;
Covar1=[];Covar2=[];

if nargin<1 
    
    %FFFFFFFFFFF
    func = input('please type 1 (or 2) for pearson (or partial) correlation analysis  ');
    %FFFFFFFFFFF
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    Group1 = ff.mat4GAT.g1;
    Group2 = ff.mat4GAT.g2;
    NullType = ff.mat4GAT.null;
    RandIter = ff.mat4GAT.nullIter;
    
    OutputFName='KAN_FixedDens_Results.mat';
    
    %Cent = input('centering data? 1 (Yes) 2(No): ');
    Cent = 2;
    
    %%Input conversion
    Input1 = ff.mat4GAT.data1;
    Input2 = ff.mat4GAT.data2;
    ROIs = ff.mat4GAT.roi1;
    
    
    flag_covar=0;
    
    if isequal(ff.mat4GAT.flagCov,1)
        
        
        CovModel = ff.mat4GAT.CovModel;
        
        covar1 = ff.mat4GAT.cov1;
        covar2 = ff.mat4GAT.cov2;
        flag_covar = 1;
        
        Cent = 1;
        
    end
    
    GUIinput=1;
    
end


%% covariates
if isequal(GUIinput,0)

    flag_covar=0;
    
    if ~isempty(varargin)
    
        flag_covar=1;
        covar1=varargin{1};covar2=varargin{2};
        NullType = 1;
        
    end
    
end

%% making original unthresholded correlation matrices
if isequal(flag_covar,0)
    
    if isequal(func,1)
    
        [R_1,P1]=corrcoef(Input1);
        [R_2,P2]=corrcoef(Input2);
    
    elseif isequal(func,2)
    
        [R_1,P1]=PartialOut(Input1);
        [R_2,P2]=PartialOut(Input2);
    
    else
        
        error('no appropriate value for structural/functional analysis');
    
    end
    
elseif isequal(flag_covar,1)

    N1=size(Input1,1);N2=size(Input2,1);
    Input=[Input1;Input2];covar=[covar1;covar2];
    
    if Cent == 1
        
        for ee = 1:size(ROIs,2)
            
            Input(:,ee) = Input(:,ee) - mean(Input(:,ee));
            
        end
    
        for ee = 1:size(covar,2)
            
            covar(:,ee) = covar(:,ee) - mean(covar(:,ee));
            
        end

    end
    
    
    for i=1:size(Input,2)
    
        Stats=regstats(Input(:,i),covar,CovModel);
        res(:,i)=Stats.r;
    
    end
    
    res1=res(1:N1,:);res2=res(N1+1:N1+N2,:);
    
    if isequal(func,1)
    
        [R_1,P1]=corrcoef(res1);
        [R_2,P2]=corrcoef(res2);
    
    elseif isequal(func,2)

        [R_1,P1]=PartialOut(res1);
        [R_2,P2]=PartialOut(res2);
        
    else
        
        error('no appropriate value for structural/functional analysis');
    
    end
    
    Input1 = res1;Input2=res2;
    flag_covar = 0;
    
end

%% Saving the correlation matrices for comparing groups at the end

Cor1 = R_1;Cor2 = R_2;
Pval1 = P1;Pval2 = P2;

save Unthresholded_Correlations R_1 R_2 P1 P2

%% making diagonal zero
Sz1 = size(R_1,1);
Sz2=size(R_2,1);
R_1(1:Sz1+1:end) = 0;
R_2(1:Sz2+1:end) = 0;


%% change negative corrrelations to zero
R_1(R_1<0)=0;R_2(R_2<0)=0;

%% Threshold the arrays

% Here, Density should be accounted as number of edges while Dens is accounted for actual density!

R1 = R_1;R2 = R_2;

Thresh_Density1 = [];
Thresh_Density2 = [];

R1max=nanmax(max(R1));R1min=nanmin(min(R1));StepR1=(R1max-R1min)/N_repeat;
R2max=nanmax(max(R2));R2min=nanmin(min(R2));StepR2=(R2max-R2min)/N_repeat;


fprintf('%-4s\n',' calculating threshold level....');

%for R1
for i=R1min:StepR1:R1max

    R1(R1<=i)=0;

    Thresh_Density1=[Thresh_Density1;i nnz(triu(R1))];

    Sum1C=sum(abs(R1));Stop1C=any(Sum1C==0);
    Sum1R=sum(abs(R1'));Stop1R=any(Sum1R==0);

    if Stop1C==1 || Stop1R==1

        Thresh1=i-StepR1;
        R1=R_1;R1(R1<=Thresh1)=0;
        break

    end

end

%for R2
for j=R2min:StepR2:R2max

    R2(R2<=j)=0;
    
    Thresh_Density2=[Thresh_Density2;j nnz(triu(R2))];
    Sum2C=sum(abs(R2));Stop2C=any(Sum2C==0);
    Sum2R=sum(abs(R2'));Stop2R=any(Sum2R==0);

    if Stop2C==1 || Stop2R==1

        Thresh2=j-StepR2;
        R2=R_2;R2(R2<=Thresh2)=0;
        break

    end

end
  

[Dens1,Nodes1,Edges1] = density_und(R1);
[Dens2,Nodes2,Edges2] = density_und(R2);

Density1=Edges1;
Density2=Edges2;

%% searching for minimum threshold in which both graphs are fully connected

if Density1>Density2
    
    Output1=R1;
    
    Net1_Density=Dens1;
    Density1=Edges1; 
    Thresh1=Thresh1;
    
    Ind=find(Thresh_Density2(:,2)==Density1);
    
    if isempty(Ind)
    
        fprintf('%-4s\n',' no similar density: increasing the interval...');
        Ind=find( Thresh_Density2(:,2)>=Density1-Precision & Thresh_Density2(:,2)<=Density1+Precision);
        
        if size(Ind)>1 
        
            Err=abs(Thresh_Density2(Ind,2)-Density1);
            [MIN,Imin]=min(Err);
            Density2=Thresh_Density2(Ind(Imin),2);
            Thresh2=Thresh_Density2(Ind(Imin),1);
        
            R2=R_2;R2(R2<=Thresh2)=0;Output2=R2;
                     
        elseif isempty(Ind)
            
            fprintf('%-4s\n',' no close density within 1*precision...');
            Ind=find( Thresh_Density2(:,2)>=Density1-2*Precision & Thresh_Density2(:,2)<=Density1+2*Precision);
            
            if size(Ind)>=1 %if more than one density found
            
                Err=abs(Thresh_Density2(Ind,2)-Density1);
                [MIN,Imin]=min(Err);
                Density2=Thresh_Density2(Ind(Imin),2);
                Thresh2=Thresh_Density2(Ind(Imin),1);
                
                R2=R_2;R2(R2<=Thresh2)=0;Output2=R2;
                
                
            elseif isempty(Ind)
            
                fprintf('%-4s\n',' no close density within 2*precision...');
            
            end
            
        end
        
    else
        
        Density2=Thresh_Density2(Ind(1),2);
        Thresh2=Thresh_Density2(Ind(1),1);
        
        R2=R_2;R2(R2<=Thresh2)=0;Output2=R2;
        
    end
    
elseif Density2>Density1

    Output2=R2;
    
    Net2_Density=Dens2;
    Density2=Edges2; 
    Thresh2=Thresh2;
    
    Ind=find(Thresh_Density1(:,2)==Density2);
    
    if isempty(Ind)
    
        fprintf('%-4s\n',' no similar density: increasing the interval...');
        Ind=find( Thresh_Density1(:,2)>=Density2-Precision & Thresh_Density1(:,2)<=Density2+Precision);
        
        if size(Ind)>=1 
        
            Err=abs(Thresh_Density1(Ind,2)-Density2);
            [MIN,Imin]=min(Err);
            Density1=Thresh_Density1(Ind(Imin),2);
            Thresh1=Thresh_Density1(Ind(Imin),1);
            
            R1=R_1;R1(R1<=Thresh1)=0;Output1=R1;
            
            
        elseif isempty(Ind)
        
            fprintf('%-4s\n',' no close density within 1*precision...');
            Ind=find( Thresh_Density1(:,2)>=Density2-2*Precision & Thresh_Density1(:,2)<=Density2+2*Precision);
            
            if size(Ind)>=1 
            
                Err=abs(Thresh_Density1(Ind,2)-Density2);
                [MIN,Imin]=min(Err);
                Density1=Thresh_Density1(Ind(Imin),2);
                Thresh1=Thresh_Density1(Ind(Imin),1);
                
            
                R1=R_1;R1(R1<=Thresh1)=0;Output1=R1;
                               
            elseif isempty(Ind)
            
                fprintf('%-4s\n',' no close density within 2*precision: !!! increase the precision and start over !!!');
            
            end
            
        end
        
    else
        
        Density1=Thresh_Density1(Ind(1),2);
        Thresh1=Thresh_Density1(Ind(1),1);
        
        
        R1=R_1;R1(R1<=Thresh1)=0;Output1=R1;
                
    end
    
end


if isequal(NullType,3) || isequal(NullType,4)
        
    NetMes_Bin1 = NetMeasures_Binary(Output1,NullType,RandIter,Input1,Cor1);
    NetMes_Bin2 = NetMeasures_Binary(Output2,NullType,RandIter,Input2,Cor2);

else
    
    NetMes_Bin1=NetMeasures_Binary(Output1,NullType,RandIter);
    NetMes_Bin2=NetMeasures_Binary(Output2,NullType,RandIter);

end 


Output1_Binary=Output1;Output1_Binary(Output1_Binary>0)=1;
Output2_Binary=Output2;Output2_Binary(Output2_Binary>0)=1;

save(OutputFName,'Output1','Output2','Output1_Binary','Output2_Binary','Thresh1','Thresh2','Density1','Density2','NetMes_Bin1','NetMes_Bin2');

%% plotting the thresholded correlation matrices
figure;imagesc(Output1);title([Group1 'correlation matrix']);hgsave([Group1 '_correlation_matrix.fig']);
figure;imagesc(Output2);title([Group2 'correlation matrix']);hgsave([Group2 '_correlation_matrix.fig']);
figure;imagesc(Output1_Binary);title([Group1 'Binary Graph']);hgsave([Group1 '_Binary_Graph.fig']);
figure;imagesc(Output2_Binary);title([Group2 'Binary Graph']);hgsave([Group2 '_Binary_Graph.fig']);

fprintf('%-4s\n','............................');
MinDens = Density1./(Sz1.*(Sz1-1)./2);
fprintf('%-4s\n', [' The minimum density of full connectivity is : ' num2str(MinDens) ] )
fprintf('%-4s\n','............................');

fprintf('%-4s\n','........... done ............');

    









