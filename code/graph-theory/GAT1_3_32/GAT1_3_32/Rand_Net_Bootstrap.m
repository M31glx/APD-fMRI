function [R1,P1,R2,P2,R_G1,P_G1,R_G2,P_G2]=Rand_Net_Bootstrap(data1,data2,Nperm,OutputFName,func,varargin)

%% 
%   This function makes "Nperm" random connectivity matrices by randomly 
%   assigning each subject's data to a group. It uses sampling with 
%   replacement and the size of the generated samples for each group is
%   the same as the original group size 


%Input:
% data1/2 : group 1/2 data (e.g. regional GMV data for a group with each 
%           row representing one subject's rGMV data)
% Nperm: number of permutation (sampling) (default 100)

% varargin: optional input arguments are those (covariates)for which you 
% need to correct the data:
% Standard format:
% 'var1': an array with the same number of rows as data1 whose columns
%       contin different covariates 
% 'var2': an array with the same number of rows as data1 whose columns
%       contin different covariates

%output
% 'R1'/'P1': Original unthresholded correlation/p-value matrix for group1
% 'R2'/'P2': Original unthresholded correlation/p-value matrix for group2
% 'R_G1/P_G1': Cell array containing Randomly generated unthresholded correlation/p-value matrix for group1
% 'R_G2/P_G2': Cell array containing Randomly generated unthresholded correlation/p-value matrix for group2

%%
%-- Hadi Hosseini, Mar 21, 2011
%-- updated on Apr 20,2011 for GUI

%Ref (Bernhardt et al.,Cereb Cortex 2011)

%% Input dialog  

GUIinput=0;

if nargin<1
    Nperm=input('please type the number of random networks to be generated (default is 1000) and press enter  ');
    if isempty(Nperm)
        Nperm=100;
    end
    OutputFName='RandNetworksBootstrap_Results.mat';

    func=input('please type 1 (or 2) for pearson (or partial) correlation analysis  ');
    
    data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');
    ff = load(data_log,'mat4GAT');
    
    data1 = ff.mat4GAT.data1;
    data2 = ff.mat4GAT.data2;
    
    flag_covar=0;
    
    if isequal(ff.mat4GAT.flagCov,1)
        
        CovModel = ff.mat4GAT.CovModel;
        
        covar1 = ff.mat4GAT.cov1;
        covar2 = ff.mat4GAT.cov2;
        
        flag_covar = 1;
    end
    
    GUIinput=1;
    
end


%% covariates

if isequal(GUIinput,0)
    
    if nargin <3
        NPerm=1000;
        OutputFName='RandNetworksBootstrap_Results.mat';
        func=1;
    
    elseif nargin <2
        
        'error: please provide one data set for each group)' 
        
    end
    
    flag_covar=0;
    
    if ~isempty(varargin)
        
        flag_covar=1;
        covar1=varargin{1};covar2=varargin{2};
        
    end
    
end

%% making original unthresholded correlation matrices

if isequal(flag_covar,0)
    
    if isequal(func,1)
         
        [R1,P1]=corrcoef(data1);
        [R2,P2]=corrcoef(data2);
         
    elseif isequal(func,2)
         
        [R1,P1]=PartialOut(data1);
        [R2,P2]=PartialOut(data2);
    
    else
        
        error('no appropriate value for structural/functional analysis');
    
    end
    
elseif isequal(flag_covar,1)
    
    n1=size(data1,1);n2=size(data2,1);
    Input=[data1;data2];covar=[covar1;covar2];
    
    for i=1:size(Input,2)
    
        Stats=regstats(Input(:,i),covar,CovModel);
        res(:,i)=Stats.r;
        
    end
    
    res1=res(1:n1,:);res2=res(n1+1:n1+n2,:);
    
    if isequal(func,1)
 
        [R1,P1]=corrcoef(res1);
        [R2,P2]=corrcoef(res2);
    
    elseif isequal(func,2)
         
        [R1,P1]=PartialOut(res1);
         [R2,P2]=PartialOut(res2);
    
    else
        
        error('no appropriate value for structural/functional analysis');
    
    end
    
    data1 = res1;data2 = res2;
    flag_covar = 0;
    
end
    
%% mixing the data for random sampling

Sz1=size(data1,1);Sz2=size(data2,1);
data=[data1;data2];
RandIndex=randperm(Sz1+Sz2);
RandData(1:Sz1+Sz2,:)=data(RandIndex(1:Sz1+Sz2),:);

%% Random sampling and generating Nperm unthresholded correlation matrices,  

dataG1=cell(1,Nperm);dataG2=cell(1,Nperm);
R_G1=cell(1,Nperm);R_G2=cell(1,Nperm);
P_G1=cell(1,Nperm);P_G2=cell(1,Nperm);

for i=1:Nperm
    
    fprintf('%-4s\n',['generating random network #' num2str(i) '...']);
    Samp1=randsample(Sz1+Sz2,Sz1,'true');
    Samp2=randsample(Sz1+Sz2,Sz2,'true');
    dataG1{i}=RandData(Samp1,:);
    dataG2{i}=RandData(Samp2,:);
    
    if isequal(func,1)
        
        [R_G1{i},P_G1{i}]=corrcoef(dataG1{i});
        [R_G2{i},P_G2{i}]=corrcoef(dataG2{i});
        
    elseif isequal(func,2)
        
        [R_G1{i},P_G1{i}]=PartialOut(dataG1{i});
        [R_G2{i},P_G2{i}]=PartialOut(dataG2{i}); 
        
    else
        
        error('no appropriate value for structural/functional analysis');
        
    end
    
end

save(OutputFName,'R1','P1','R2','P2','R_G1','P_G1','R_G2','P_G2','dataG1','dataG2','data1','data2');


end

