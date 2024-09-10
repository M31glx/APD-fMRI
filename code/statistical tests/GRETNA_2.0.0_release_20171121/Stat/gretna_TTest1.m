function [TTest1_T, TTest1_P] = gretna_TTest1(DependentFiles, Covariates, Base)
% TTest1_T = gretna_TTest1(DependentDirs, OutputName, OtherCovariates, Base)
% Perform one sample t test.
% Input:
%   DependentFiles - the metric file of dependent variable. 1 by 1 cell
%   OutputName - the output name.
%   Covariates - The other covariates. 1 by 1 cell 
%   Base - the base of one sample T Test. 0: default.
% Output:
%   TTest1_T - the T value, also write image file out indicated by OutputName
%___________________________________________________________________________
% Written by Wang Xindi 140615.
% State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University, China, 100875
%___________________________________________________________________________
% Modified from DPABI's y_TTest1_Image
% Written by YAN Chao-Gan 100411.
% State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University, China, 100875
% Core function re-written by YAN Chao-Gan 14011.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com

if ~exist('Base', 'var')
    Base=0;
end

DependentMatrix=[];
CovariatesMatrix=[];

for i=1:1
    if ischar(DependentFiles{i})
        Matrix=load(DependentFiles{i});
        fprintf('\n\tGroup %.1d: Load File %s:\n', i, DependentFiles{i});        
    else
        Matrix=DependentFiles{i};
    end
    DependentMatrix=cat(1, DependentMatrix, Matrix);
    
    if exist('Covariates','var') && ~isempty(Covariates)
        CovariatesMatrix=cat(1, CovariatesMatrix, Covariates{i});
    end
end

[NumOfSample, NumOfMetric]=size(DependentMatrix);
DependentMatrix=DependentMatrix-Base;

Regressors=ones(NumOfSample, 1);
Regressors=[Regressors, CovariatesMatrix];

Contrast=zeros(1, size(Regressors, 2));
Contrast(1)=1;

[b_OLS_metric, t_OLS_metric, TTest1_T, r_OLS_metric] = gretna_GroupAnalysis(DependentMatrix, Regressors, Contrast, 'T');

DOF=NumOfSample-size(Regressors, 2);
TTest1_P=2*(1-tcdf(abs(TTest1_T), DOF));

fprintf('\n\tOne Sample T Test Calculation finished.\n');