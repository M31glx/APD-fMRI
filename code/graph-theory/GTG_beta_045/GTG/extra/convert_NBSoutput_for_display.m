function convert_NBSoutput_for_display(nbs,outname,disp_type,node_size_type,node_color_type,edge_size_type,node_colors,cmap,circos_disp_type)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 01.15.15
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2015. Software can be modified and redistributed, but modifed, 
% redistributed versions must have the same rights

if ~exist('outname','var') || isempty(outname)
    outname = 'temp';
end
if ~exist('disp_type','var') || isempty(disp_type)
    disp_type = 1;
end
if ~exist('node_size_type','var') || isempty(node_size_type)
    node_size_type = 3;
end
if ~exist('node_color_type','var') || isempty(node_color_type)
    node_color_type = 4;
end
if ~exist('edge_size_type','var') || isempty(edge_size_type)
    edge_size_type = 3;
end

for net = 1:nbs.NBS.n
    labels = nbs.NBS.node_label(logical(sum(full(nbs.NBS.con_mat{1,net})+full(nbs.NBS.con_mat{1,net})') > 0));
    labels = labels(:)';
    mask   = (full(nbs.NBS.con_mat{1,net})+full(nbs.NBS.con_mat{1,net})');
    
    if edge_size_type==1 % Edge strength reflects effect size
        pvals = cdf('t',nbs.GLM.contrast(logical(nbs.GLM.contrast~=0))*nbs.NBS.test_stat,size(nbs.GLM.X,1)-size(nbs.GLM.X,2));
        zmat  = icdf('norm',pvals,0,1).*mask;
        zmat  = zmat(sum(mask,2)>0,sum(mask,2)>0);
        if maxn(abs(zmat))>1
            zmat = (zmat./maxn(abs(zmat)));
        end
        zmat = weight_conversion(zmat,'autofix');
    elseif edge_size_type==2 % Edge strength reflects mean strength across participants
        tempcoords              = repmat(logical(triu(ones(length(nbs.NBS.node_label))-diag(ones(size(nbs.NBS.node_label))),1)),[1,1,size(nbs.GLM.y,1)]);
        full_conmat             = zeros(length(nbs.NBS.node_label),length(nbs.NBS.node_label),size(nbs.GLM.y,1));
        full_conmat(tempcoords) = shiftdim(nbs.GLM.y,1);
        full_conmat             = full_conmat+permute(full_conmat,[2 1 3]);
        fish_full_conmat        = atanh(full_conmat);
        full_mean_mat           = atan(mean(fish_full_conmat,3)).*mask;
        zmat                    = full_mean_mat(sum(mask,2)>0,sum(mask,2)>0);
        if maxn(abs(zmat))>1
            zmat = zmat./maxn(abs(zmat));
        end
        zmat = weight_conversion(zmat,'autofix');
    elseif edge_size_type==3 % All significant edges equal to 1
        zmat = mask(sum(mask,2)>0,sum(mask,2)>0);
    end
    
    if edge_size_type~=3 % Re-scale range to 1:10
        red_mask = zmat~=0;
        zmat(red_mask) = zmat(red_mask)-minn(zmat(red_mask));
        zmat = (zmat./maxn(abs(zmat))).*9;
        zmat(red_mask) = zmat(red_mask)+1;
    end
    
    if node_size_type==1 % Node size reflects mean Node Strength across participants
        if edge_size_type~=2
            tempcoords              = repmat(logical(triu(ones(length(nbs.NBS.node_label))-diag(ones(size(nbs.NBS.node_label))),1)),[1,1,size(nbs.GLM.y,1)]);
            full_conmat             = zeros(length(nbs.NBS.node_label),length(nbs.NBS.node_label),size(nbs.GLM.y,1));
            full_conmat(tempcoords) = shiftdim(nbs.GLM.y,1);
            full_conmat             = full_conmat + permute(full_conmat,[2 1 3]);
            fish_full_conmat        = atanh(full_conmat);
            full_mean_mat           = atan(mean(fish_full_conmat,3));
        end
        full_mean_mat = weight_conversion(full_mean_mat,'autofix');
        % full_mean_mat = abs(full_mean_mat);     % Uncomment if you want to use the absolute value of weights instead of just positive weights
        node_sizes    = strengths_und_sign(full_mean_mat);
        node_sizes    = node_sizes(sum(mask,2)>0);
    elseif node_size_type==2 % Node size reflects the magnitude of NBS effects for all differential links for that node
        if edge_size_type==1
            node_sizes = strengths_und_sign(abs(zmat));
        else
            pvals = cdf('t',nbs.GLM.contrast(logical(nbs.GLM.contrast~=0))*nbs.NBS.test_stat,size(nbs.GLM.X,1)-size(nbs.GLM.X,2));
            pzmat = icdf('norm',pvals,0,1).*mask;
            pzmat = pzmat(sum(mask,2)>0,sum(mask,2)>0);
            if maxn(abs(pzmat))>1
                pzmat = pzmat./maxn(abs(pzmat));
            end
            node_sizes = strengths_und_sign(abs(pzmat));
        end
    elseif node_size_type==3 % All nodes have equal size
        node_sizes = ones(length(labels),1);
    end
    if node_size_type~=3 % Re-scale range to 1:10
        node_sizes = node_sizes-min(node_sizes);
        node_sizes = (node_sizes./max(node_sizes)).*9;
        node_sizes = node_sizes'+1;
    end
    
    if node_color_type==1 % Node color reflects mean Node Strength across participants
        if edge_size_type~=2 && node_size_type~=1
            tempcoords              = repmat(logical(triu(ones(length(nbs.NBS.node_label))-diag(ones(size(nbs.NBS.node_label))),1)),[1,1,size(nbs.GLM.y,1)]);
            full_conmat             = zeros(length(nbs.NBS.node_label),length(nbs.NBS.node_label),size(nbs.GLM.y,1));
            full_conmat(tempcoords) = shiftdim(nbs.GLM.y,1);
            full_conmat             = full_conmat + permute(full_conmat,[2 1 3]);
            fish_full_conmat        = atanh(full_conmat);
            full_mean_mat           = atan(mean(fish_full_conmat,3));
        end
        full_mean_mat = weight_conversion(full_mean_mat,'autofix');
        % full_mean_mat = abs(full_mean_mat);     % Uncomment if you want to use the absolute value of weights instead of just positive weights
        node_colors   = strengths_und_sign(full_mean_mat);
        node_colors   = node_colors(sum(mask,2)>0);
    elseif node_color_type==2 % Node value reflects the magnitude of NBS effects for all differential links for that node
        if node_size_type==2
            node_colors = node_sizes;
        elseif edge_size_type==1
            node_colors = strengths_und_sign(abs(zmat));
        else
            pvals = cdf('t',nbs.GLM.contrast(logical(nbs.GLM.contrast~=0))*nbs.NBS.test_stat,size(nbs.GLM.X,1)-size(nbs.GLM.X,2));
            pzmat = icdf('norm',pvals,0,1).*mask;
            pzmat = pzmat(sum(mask,2)>0,sum(mask,2)>0);
            if maxn(abs(pzmat))>1
                pzmat = pzmat./maxn(abs(pzmat));
            end
            node_colors = strengths_und_sign(abs(pzmat));
        end
    elseif node_color_type==3 % Node color values have been entered
        if length(node_colors)>length(node_sizes)
            node_colors = node_colors(sum(mask,2)>0);
        end
        [~,~,node_colors] = unique(node_colors);
    elseif node_color_type==4 % All nodes have equal value
        node_colors = ones(length(labels),1);
    end
    
    if node_color_type~=3 && node_color_type~=4% Re-scale range to 1:10
        node_colors = node_colors-min(node_colors);
        node_colors = (node_colors./max(node_colors)).*9;
        node_colors = node_colors'+1;
    end
    
    if exist('cmap','var') && ~isempty(cmap)
        if ischar(cmap)
            eval(['cmap = ',cmap,'(round(1.3*size(zmat,1)));']);
            cmap(size(zmat,1)+1:end,:) = [];
        end
    else
        cmap = parula(size(zmat,1));
    end
    
    if disp_type==1 % Create output for BrainNet Viewer
        coords = nbs.NBS.node_coor(logical(sum(full(nbs.NBS.con_mat{1,net})+full(nbs.NBS.con_mat{1,net})') > 0),:);
        labels   = strrep(labels,'-','');
        labels   = strrep(labels,'_','');
        labels   = strrep(labels,' ','');
        labels   = strrep(labels,'.','');
        nodeinfo = [num2cell([coords,node_colors,node_sizes]');labels];
        if nbs.NBS.n > 1
            dlmwrite([outname,'_',num2str(net),'.edge'],zmat,'\t');
            fid = fopen([outname,'_',num2str(net),'.node'],'w');
        else
            dlmwrite([outname,'.edge'],zmat,'\t');
            fid = fopen([outname,'.node'],'w');
        end
        fprintf(fid,'%i\t%i\t%i\t%i\t%i\t%s\n',nodeinfo{:});
        fclose(fid);
    elseif disp_type==2 % Create output for SONIA
        [node_colors, ~, node_color_inds] = unique(node_colors);
        cc_mod                          = hsv(numel(node_colors)+1);
        nodeinfo                        = [num2cell(1:length(labels));labels;num2cell(10*ones(size(node_sizes')));repmat({'ellipse'},1,length(labels));num2cell(cc_mod(node_color_inds,1))';num2cell(cc_mod(node_color_inds,2))';num2cell(cc_mod(node_color_inds,3))';repmat({'black'},1,length(labels));repmat({'10'},1,length(labels));repmat({'black'},1,length(labels));repmat({'0'},1,length(labels));repmat({'1'},1,length(labels))];
        [r,c]                           = find(triu(zmat~=0));
        arc_colors                      = (cc_mod(node_color_inds(r),:)+cc_mod(node_color_inds(c),:))/2;
        arcinfo                         = [num2cell([r,c,(100*abs(zmat(triu(zmat)~=0)))]');num2cell(zeros(size(arc_colors(:,1))))';num2cell(zeros(size(arc_colors(:,2))))';num2cell(zeros(size(arc_colors(:,3))))';repmat({'0'},1,length(r));repmat({'1'},1,length(r))];
        if nbs.NBS.n > 1
            fid = fopen([outname,'_',num2str(net),'.son'],'w');
        else
            fid = fopen([outname,'.son'],'w');
        end
        fprintf(fid,'NodeId\tLabel\tNodeSize\tNodeShape\tRedRGB\tGreenRGB\tBlueRGB\tLabelColor\tLabelSize\tBorderColor\tStartTime\tEndTime\n');
        fprintf(fid,'%i\t%s\t%i\t%s\t%f\t%f\t%f\t%s\t%s\t%s\t%s\t%s\n',nodeinfo{:});
        fprintf(fid,'FromId\tToId\tArcWeight\tRedRGB\tGreenRGB\tBlueRGB\tStartTime\tEndTime\n');
        fprintf(fid,'%i\t%i\t%f\t%f\t%f\t%f\t%s\t%s\n',arcinfo{:});
        fclose(fid);
    elseif disp_type==3 % Create simulated anealing circle graph
        labels = strrep(labels,'_','');
        if nbs.NBS.n>1
            outname_temp = [outname,'_mod_aneal_',num2str(net)];
        else
            outname_temp = [outname,'_mod_aneal'];
        end
        [zmat_reord_simanneal,simanneal_ind] = reorder_matrix(abs(zmat),'circ',0);
        zmat_reord_simanneal                 = zmat_reord_simanneal.*sign(zmat(simanneal_ind,simanneal_ind));
        simanneal_labels                     = labels(simanneal_ind);
        circ_graph_jms(zmat_reord_simanneal,simanneal_labels,outname_temp);
    elseif disp_type==4 % Create diagonal organization circle graph
        labels = strrep(labels,'_','');
        if nbs.NBS.n>1
            outname_temp = [outname,'_mod_diag_',num2str(net)];
        else
            outname_temp = [outname,'_mod_diag'];
        end
        [zmat_reord_mat,mat_ind] = reorderMAT(abs(zmat),100,'circ');
        zmat_reord_mat           = zmat_reord_mat.*sign(zmat(mat_ind,mat_ind));
        mat_labels               = labels(mat_ind);
        circ_graph_jms(zmat_reord_mat,mat_labels,outname_temp);
    elseif disp_type==5 % Create modular circle graph
        labels = strrep(labels,'_','');
        if nbs.NBS.n>1
            outname_temp = [outname,'_mod_circ_',num2str(net)];
        else
            outname_temp = [outname,'_mod_circ'];
        end
        [mod_ind,zmat_reord_mod] = reorder_mod(abs(zmat),node_colors);
        zmat_reord_mod           = abs(zmat_reord_mod).*sign(zmat(mod_ind,mod_ind));
        mod_labels               = labels(mod_ind);
        node_colors              = node_colors(mod_ind);
        figure;
        circ_graph_jms(zmat_reord_mod,mod_labels,outname_temp,node_colors);
    elseif disp_type==6 % gephi
        labels = strrep(labels,'-','');
        labels = strrep(labels,'_','');
        labels = strrep(labels,' ','');
        labels = strrep(labels,'.','');
        
        [r,c] = find(triu(zmat~=0));
        arcinfo = [{'source','target','weight'};labels(r)',labels(c)',num2cell(abs(zmat(triu(zmat)~=0)))];
        %arcinfo = [labels(r)',labels(c)'];
        %czmat   = num2cell(double(zmat));
        %czmat   = [labels;czmat];
        %czmat   = [[cell(1,1);labels'],czmat];
        %czmat   = [[' ';labels'],czmat];
        if nbs.NBS.n > 1
            outname2 = [outname,'_',num2str(net),'_gephi.csv'];
        else
            outname2 = [outname,'_gephi.csv'];
        end
        cell2csv(outname2,arcinfo,';')
        %cell2csv(outname2,czmat,';',2001)
    elseif disp_type==7 % cytoscape sif file
        labels  = strrep(labels,'__','/');
        labels  = strrep(labels,'-',' ');
        labels  = strrep(labels,'_',' ');
        labels  = strrep(labels,'  ',' ');
        labels  = strrep(labels,'  ',' ');
        labels  = strrep(labels,'.',' ');
        newlabs = cell(length(labels),1);
        for q = 1:length(labels)
            newlabs{q} = strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(labels{q},'_10',''),'_9',''),'_8',''),'_7',''),'_6',''),'_5',''),'_4',''),'_3',''),'_2',''),'_1','');
        end
        for q = 1:length(labels)
            for w = 1:length(labels)
                if q ~=w
                    if strcmp(newlabs{q},newlabs{w})
                        newlabs{q} = labels{q};
                    end
                end
            end
        end
        labels  = newlabs{:}';
        [r,c]   = find(triu(zmat~=0));
        arcinfo = [labels(r);labels(c)];
        if nbs.NBS.n > 1
            fid = fopen([outname,'_',num2str(net),'_cyto.sif'],'w');
        else
            fid = fopen([outname,'_cyto.sif'],'w');
        end
        fprintf(fid,'%s fc %s\n',arcinfo{:});
        fclose(fid);
    elseif disp_type==8 % cytoscape file with weight info
        labels  = strrep(labels,'__','/');
        labels  = strrep(labels,'-',' ');
        labels  = strrep(labels,'_',' ');
        labels  = strrep(labels,'  ',' ');
        labels  = strrep(labels,'  ',' ');
        labels  = strrep(labels,'.',' ');
        newlabs = cell(length(labels),1);
        for q = 1:length(labels)
            newlabs{q} = strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(labels{q},'_10',''),'_9',''),'_8',''),'_7',''),'_6',''),'_5',''),'_4',''),'_3',''),'_2',''),'_1','');
        end
        for q = 1:length(labels)
            for w = 1:length(labels)
                if q ~=w
                    if strcmp(newlabs{q},newlabs{w})
                        newlabs{q} = labels{q};
                    end
                end
            end
        end
        labels  = newlabs(:)';
        [r,c]   = find(triu(zmat~=0));
        arcinfo = [labels(r);labels(c);num2cell(abs(zmat(triu(zmat)~=0)))'];
        if nbs.NBS.n > 1
            fid = fopen([outname,'_',num2str(net),'_cyto.csv'],'w');
        else
            fid = fopen([outname,'_cyto.csv'],'w');
        end
        fprintf(fid,'source,target,weight\n');
        fprintf(fid,'%s,%s,%f\n',arcinfo{:});
        %fprintf(fid,'%s (fc) %s = %f\n',arcinfo{:});
        fclose(fid);
    elseif disp_type==9 % circos-like circular graph
        newlabs = cell(length(labels),1);
        for q = 1:length(labels)
            newlabs{q} = strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(labels{q},'_10',''),'_9',''),'_8',''),'_7',''),'_6',''),'_5',''),'_4',''),'_3',''),'_2',''),'_1','');
        end
        for q = 1:length(labels)
            for w = 1:length(labels)
                if q ~=w
                    if strcmp(newlabs{q},newlabs{w})
                        newlabs{q} = labels{q};
                    end
                end
            end
        end
        labels = newlabs;
        labels = strrep(labels,'__','/');
        labels = strrep(labels,'-',' ');
        labels = strrep(labels,'_',' ');
        labels = strrep(labels,'  ',' ');
        labels = strrep(labels,'  ',' ');
        labels = strrep(labels,'.',' ');
        
        if ~exist('circos_disp_type','var') || isempty(circos_disp_type)
            circos_disp_type = 'deg_sort';
        end
        switch circos_disp_type
            case 'deg_sort'
                [~,ind] = sort(degrees_und(zmat),'descend');
                labels  = labels(ind);
                zmat    = zmat(ind,ind);
            case 'sim_anneal'
                [zmat_reord,ind] = reorder_matrix(abs(zmat),'circ',0);
                zmat             = zmat_reord.*sign(zmat(ind,ind));
                labels           = labels(ind);
%                 num_ann_reps = 10;
%                 zmat_reord = zeros(size(zmat,1),size(zmat,2),num_ann_reps);
%                 ind        = zeros(size(zmat,1),num_ann_reps);
%                 cost       = zeros(num_ann_reps,1);
%                 for q = 1:num_ann_reps
%                     [zmat_reord(:,:,q),ind(:,q),cost(q)] = reorder_matrix(abs(zmat),'circ',0);
%                 end
%                 best_fit = find(cost==min(cost),1);
%                 disp(['Lowest cost = ',num2str(min(cost))])
%                 zmat     = zmat_reord(:,:,best_fit).*sign(zmat(ind(:,best_fit),ind(:,best_fit)));
%                 labels   = labels(ind(:,best_fit));
            case 'diag'
                [zmat_reord,ind] = reorderMAT(abs(zmat),100,'circ');
                zmat             = zmat_reord.*sign(zmat(ind,ind));
                labels           = labels(ind);
        end
        circularGraph(zmat,'Label',labels,'Colormap',cmap);
    elseif disp_type==10 % binarized circos-like circular graph
        newlabs = cell(1,length(labels));
        for q = 1:length(labels)
            newlabs{q} = strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(labels{q},'_10',''),'_9',''),'_8',''),'_7',''),'_6',''),'_5',''),'_4',''),'_3',''),'_2',''),'_1','');
        end
        for q = 1:length(labels)
            for w = 1:length(labels)
                if q ~=w
                    if strcmp(newlabs{q},newlabs{w})
                        newlabs{q} = labels{q};
                    end
                end
            end
        end
        labels = newlabs;
        labels = strrep(labels,'__','/');
        labels = strrep(labels,'-',' ');
        labels = strrep(labels,'_',' ');
        labels = strrep(labels,'  ',' ');
        labels = strrep(labels,'  ',' ');
        labels = strrep(labels,'.',' ');
        
        if ~exist('circos_disp_type','var')
            circos_disp_type = 'deg_sort';
        end
        switch circos_disp_type
            case 'deg_sort'
                [~,ind] = sort(degrees_und(zmat),'descend');
                labels  = labels(ind);
                zmat    = zmat(ind,ind);
            case 'sim_anneal'
                [zmat_reord,ind] = reorder_matrix(abs(zmat),'circ',0);
                zmat             = zmat_reord.*sign(zmat(ind,ind));
                labels           = labels(ind);
            case 'diag'
                [zmat_reord,ind] = reorderMAT(abs(zmat),100,'circ');
                zmat             = zmat_reord.*sign(zmat(ind,ind));
                labels           = labels(ind);
        end
        circularGraph(double(zmat>0),'Label',labels,'Colormap',cmap);
    end
end
