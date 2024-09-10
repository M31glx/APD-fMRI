function Load_data

%%
Group1 = input('please type name of the 1st group in single quotation mark: ');
mat4GAT.g1 = Group1;
Group2 = input('please type name of the 2nd group in single quotation mark: ');
mat4GAT.g2 = Group2;


AnyCOV = input( 'any nuisance covariates? type 1 for YES, 2 for NO :');

if AnyCOV == 1

    CovM = input('modeling covariates with(type 1) or wo(type2) interactions: ');
        
    if isequal(CovM,1)

        mat4GAT.CovModel = 'interaction';

    else

        mat4GAT.CovModel = 'linear';

    end
    
    COV = input( 'make covariate mat files from excel sheets? type 1 for YES, 2 for NO: ');
   
   if isequal(COV,1)
       [CovName,CovPathName] = uigetfile('*.xls',['Select the Covariates XLS file for' Group1]);
       covar = xlsread([CovPathName CovName]);
       mat4GAT.cov1=covar;
       [CovName,CovPathName] = uigetfile('*.xls',['Select the Covariates XLS file for' Group2]);
       covar = xlsread([CovPathName CovName]);
       mat4GAT.cov2 = covar;
       mat4GAT.flagCov = 1;

   elseif isequal(COV,2)
       covar1 = spm_select(1,'.mat',['Select ' Group1 ' covar.mat file']);
       covar2 = spm_select(1,'.mat',['Select ' Group2 ' covar.mat file']);
       c1 = load(covar1,'covar');covar1 = c1.covar;
       mat4GAT.cov1 =covar1;
       c2 = load(covar2,'covar');covar2 = c2.covar;
       mat4GAT.cov2 = covar2;
       mat4GAT.flagCov = 1;
   
   end

else
    
    mat4GAT.flagCov = 0;
    
end
   

RexOrDirect = questdlg('Output from REX?','','Yes','No (input data directly)','Yes');
if strcmp(RexOrDirect,'Yes') 
    
    Data1 = spm_select(1,'.mat',['Select ' Group1 ' REX.mat (output from REX)']);
    Data2 = spm_select(1,'.mat',['Select ' Group2 ' REX.mat (output from REX)']);
    f1 = load(Data1,'params');mat4SVM = f1.params.ROIdata;regionName = f1.params.ROInames;
    mat4GAT.data1=mat4SVM;mat4GAT.roi1=regionName;
    f2 = load(Data2,'params');mat4SVM = f2.params.ROIdata;regionName = f2.params.ROInames;
    mat4GAT.data2=mat4SVM;mat4GAT.roi2=regionName;

else 

    MatOrXls = input( 'Type 1(or 2) for .mat (or XLS) input: ');
    
    switch MatOrXls
    
        case 1
            
             Data1 = spm_select(1,'.mat',['Select ' Group1 ' mat4SVM.mat file']);
             Data2 = spm_select(1,'.mat',['Select ' Group2 ' mat4SVM.mat file']);
             f1 = load(Data1,'mat4SVM');mat4GAT.data1 = f1.mat4SVM;
             f2 = load(Data2,'mat4SVM');mat4GAT.data2 = f2.mat4SVM;
             h1 = load(Data1,'regionName');mat4GAT.roi1 = h1.regionName;
             h2 = load(Data2,'regionName');mat4GAT.roi2 = h2.regionName;
        
        case 2
            
            [FileName,PathName] = uigetfile('*.xls',['Select the Excel data file for' Group1]);
            [mat4SVM,regionName] = xlsread([PathName FileName]);
            mat4GAT.data1=mat4SVM;mat4GAT.roi1=regionName;
            [FileName,PathName] = uigetfile('*.xls',['Select the Excel data file for' Group2]);
            [mat4SVM,regionName] = xlsread([PathName FileName]);
            mat4GAT.data2 = mat4SVM;mat4GAT.roi2 = regionName;
    
    end
    
end


null_flag = input('2 (same deg dist), 3 (H-Q-S), 4 (RandCorr): '); 
mat4GAT.null = null_flag;

Null_Iter = input('number of null networks (e.g. 20):  ');
mat4GAT.nullIter = Null_Iter;

Tail = input('one (type 1) or two (type 2) -tailed test:  ');
mat4GAT.tail = Tail;


Alpha = input('significance threshold level (e.g. 0.05):  ');
mat4GAT.alpha = Alpha;


save(['mat4GAT_' Group1 '_' Group2],'mat4GAT');
