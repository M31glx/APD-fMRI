function RegionalMapping

%%

%-- Hadi Hosseini: June, 2011
%-- updated April, 2014

% Reference (Hosseini et al.,Plos One 2012)

%%
Met = input('type 1(Betw), 2(Deg), or 3(Clustering) :');
FDR = input('type 1(uncorrected), 2(fdr corrected): '); 


data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');

ff = load(data_log,'mat4GAT');

Group1 = ff.mat4GAT.g1;
Group2 = ff.mat4GAT.g2;
ROI = ff.mat4GAT.roi1;
pthr = ff.mat4GAT.alpha;
Tail = ff.mat4GAT.tail;

switch FDR
    
    case 1 %uncorrected maps


        switch Met

            case 1 

                Data = spm_select(1,'RegNodeBetwNorm.*_pval.mat','Select the RegNodeBetwNorm*_pval.mat (output from "regional measures")');
                ee = load(Data,'p');
                pval = ee.p;
                Name = [Group2 '_vs_' Group1 '_Betw_RegionalMap.img'];

            case 2

                Data = spm_select(1,'RegDegNorm.*_pval.mat','Select the RegDegNorm*_pval.mat (output from "regional measures")');
                ee = load(Data,'p');
                pval = ee.p;
                Name = [Group2 '_vs_' Group1 '_Deg_RegionalMap.img'];

            case 3 

                Data = spm_select(1,'RegClustNorm.*_pval.mat','Select the RegClustNorm*_pval.mat (output from "regional measures")');
                ee = load(Data,'p');
                pval = ee.p;
                Name = [Group2 '_vs_' Group1 '_Clust_RegionalMap.img'];

        end
        
    case 2%fdr corrected maps
        
        if isequal(Tail,1)
        
            switch Met

                case 1 

                    Data = spm_select(1,'fdr_pval_Group1gdanGroup2.*NodeBetw','Select the pval mat file: fdr_pval_Group1gdanGroup2*_RegNodeBetwNorm.mat');
                    ee = load(Data,'pp1');
                    pval_g1g2 = ee.pp1;

                    Data = spm_select(1,'fdr_pval_Group2gdanGroup1.*NodeBetw','Select the pval mat file: fdr_pval_Group2gdanGroup1*_RegNodeBetwNorm.mat');
                    ee = load(Data,'pp2');
                    pval_g2g1 = ee.pp2;

                    Name = [Group2 '_vs_' Group1 '_Betw_RegionalMap_fdr.img'];


                case 2

                    Data = spm_select(1,'fdr_pval_Group1gdanGroup2.*Deg','Select the pval mat file: fdr_pval_Group1gdanGroup2*_RegDegNorm.mat');
                    ee = load(Data,'pp1');
                    pval_g1g2 = ee.pp1;

                    Data = spm_select(1,'fdr_pval_Group2gdanGroup1.*Deg','Select the pval mat file: fdr_pval_Group2gdanGroup1*_RegDegNorm.mat');
                    ee = load(Data,'pp2');
                    pval_g2g1 = ee.pp2;

                    Name = [Group2 '_vs_' Group1 '_Deg_RegionalMap_fdr.img'];

                case 3

                    Data = spm_select(1,'fdr_pval_Group1gdanGroup2.*Clust','Select the pval mat file: fdr_pval_Group1gdanGroup2*_RegClustNorm.mat');
                    ee = load(Data,'pp1');
                    pval_g1g2 = ee.pp1;

                    Data = spm_select(1,'fdr_pval_Group2gdanGroup1.*Clust','Select the pval mat file: fdr_pval_Group2gdanGroup1*_RegClustNorm.mat');
                    ee = load(Data,'pp2');
                    pval_g2g1 = ee.pp2;

                    Name = [Group2 '_vs_' Group1 '_Clust_RegionalMap_fdr.img'];

            end
            
            
        else
            
            switch Met

                case 1 

                    Data = spm_select(1,'fdr_pval.*NodeBetw','Select the pval mat file: fdr_pval*_RegNodeBetwNorm.mat');
                    ee = load(Data,'pp');
                    pval = ee.pp;

                    Name = [Group2 '_' Group1 '_Diff_Betw_RegionalMap_fdr.img'];


                case 2

                    Data = spm_select(1,'fdr_pval.*Deg','Select the pval mat file: fdr_pval*_RegDegNorm.mat');
                    ee = load(Data,'pp');
                    pval = ee.pp;

                    Name = [Group2 '_' Group1 '_Diff_Deg_RegionalMap_fdr.img'];

                case 3

                    Data = spm_select(1,'fdr_pval.*Clust','Select the pval mat file: fdr_pval*_RegClustNorm.mat');
                    ee = load(Data,'pp');
                    pval = ee.pp;

                    Name = [Group2 '_' Group1 '_Diff_Clust_RegionalMap_fdr.img'];

            end
            
            
            
        end
        
        
end
        
        
     

%% 

switch FDR
    
    
    case 1%uncorrected
        
        if isequal(Tail,1)
        
        
            [i2,j2] = find(pval <= pthr & pval > 0);
            roi_g2g1 = ROI(i2);
            scale_g2g1 = log ( 1 ./ abs(pval(i2)) );


            [i1,j1] = find(pval >= -pthr & pval < 0);
            roi_g1g2 = ROI(i1);
            scale_g1g2 = -log ( 1 ./ abs(pval(i1)) );

            roi_all = cat(2,roi_g2g1, roi_g1g2);
            scale_all = cat(1,scale_g2g1, scale_g1g2);
            
        else
            
            [i2,j2] = find(pval <= pthr & pval > 0);
            roi_g2g1 = ROI(i2);
            scale_g2g1 = log ( 1 ./ abs(pval(i2)) );

            roi_all = roi_g2g1;
            scale_all = scale_g2g1;
            
            
        end
        
    case 2%fdr corrected
        
        if isequal(Tail,1)
        
            [i2,j2] = find(pval_g2g1 <= pthr & pval_g2g1 > 0);
            roi_g2g1 = ROI(i2);
            scale_g2g1 = log ( 1 ./ abs(pval_g2g1(i2)) );


            [i1,j1] = find(pval_g1g2 <= pthr & pval_g1g2 > 0);
            roi_g1g2 = ROI(i1);
            scale_g1g2 = -log ( 1 ./ abs(pval_g1g2(i1)) );

            roi_all = cat(2,roi_g2g1, roi_g1g2);
            scale_all = cat(1,scale_g2g1, scale_g1g2);
            
        else
            
            [ii,jj] = find(pval <= pthr);
            roi = ROI(ii);
            scale = log ( 1 ./ abs(pval(ii)) );


            roi_all = roi;
            scale_all = scale;
            
            
           
            
        end
        
end
      


%% 


if ~isempty(scale_all)

    for i = 1:size(scale_all,1)
        
        if i==1
            
            Eval = [num2str(scale_all(i)) '*i' num2str(i)];
        
        else
            
            if scale_all(i) > 0
                
                str_scale = ['+' num2str(scale_all(i))];
                
            else
                
                str_scale = num2str(scale_all(i));
            end
            
                
            Eval = [Eval , str_scale '*i' num2str(i)];
        
        end
        
    end
    
else
    
    error('No regions survived')
    return
    
end


%%

fprintf('%-4s\n',' Select the IMAGES corresponding to following ROIs "in the same order" printed here ....');
    
rrr = roi_all'
Name2 = strrep(Name,'.img','');
cell2csv(Name2,rrr);
open(Name2)

Q = spm_imcalc_ui({},Name,Eval);

fprintf('%-4s\n',' ..........Done ');





