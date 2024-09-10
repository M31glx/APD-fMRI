function R_null = rand_null_hh(d,nr,r)

% d:   original data (or residuals), nSbj*nROI
% r:   original correlation matrix(unthresholded)
% nr:  number of random networks to be generated, scalar

nroi = size(d,2);

for jj = 1:nr
    
    for kk=1:size(d,1)
        
        rand_ind = randperm(nroi);
        dr(kk,:) = d(kk,rand_ind);
        
    end
 
    r_n = corrcoef(dr);
    r_n(1:nroi+1:end)=0;
    r(1:nroi+1:end)=0;
    
    
    
    if mad(r_n(:),1) <= mad(r(:),1)
        
        Dt = [mad(r_n(:),1)./mad(r(:),1)].*[median(r(:)) - median(r_n(:))];
        
        r_null = r_n + Dt;
        
        r_null(r_null == Dt) = 0;
        
        if max(r_null(:)) >= max(r(:))
            
            r_null(r_null >= max(r(:))) = max(r(:));
        
        end
            
        
    elseif mad(r_n(:),1) > mad(r(:),1)
        
        Dt = [mad(r(:),1)./mad(r_n(:),1)].*[median(r(:)) - median(r_n(:))];
        
        r_null = r_n + Dt;
        
        r_null(r_null == Dt) = 0;

        
        if max(r_null(:)) >= max(r(:))
            
            r_null(r_null >= max(r(:))) = max(r(:));
        
        end
            
        
    end
    
    r_null(r_null<0) = 0;
    
    R_null(:,:,jj) = r_null;

end


