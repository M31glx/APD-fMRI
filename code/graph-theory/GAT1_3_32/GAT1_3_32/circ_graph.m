function circ_graph(graph,roi,G)

% input: 
%       graph: nroi*nroi binary matrix including network connections
%       roi: 1*nroi cell array including roi names
%       G: group name (string)

nroi = size(roi,2);

theta = linspace(0,2*pi,nroi+1);

theta = theta(1:end-1);

[x,y] = pol2cart(theta,1);

[ind1,ind2] = find(graph);

figure;
plot(x,y,'.k','markersize',20);hold on
% plot([x(ind1); x(ind2)],[y(ind1); y(ind2)],'b'); 
xlim([min(x)-.5 max(x)+.5]) 
ylim([min(y)-.5 max(y)+.5]) 

text(x(1:round(nroi/4)), y(1:round(nroi/4)), roi(1,1:round(nroi/4)), 'VerticalAlignment','bottom', 'HorizontalAlignment','left','fontsize',11); 
text(x(round(nroi/4)+1:round(nroi/2)), y(round(nroi/4)+1:round(nroi/2)), roi(1,round(nroi/4)+1:round(nroi/2)), 'VerticalAlignment','bottom', 'HorizontalAlignment','right','fontsize',11); 
text(x(round(nroi/2)+1:round(3*nroi/4)), y(round(nroi/2)+1:round(3*nroi/4)), roi(1,round(nroi/2)+1:round(3*nroi/4)), 'VerticalAlignment','top', 'HorizontalAlignment','right','fontsize',11); 
text(x(round(3*nroi/4)+1:nroi), y(round(3*nroi/4)+1:nroi), roi(1,round(3*nroi/4)+1:nroi), 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',11); 


graphNorm = abs(graph)/max(max(abs(graph)));
line_widths = graphNorm.*4;

% Move the range of values so they start at 0.001 (if we set them to start at 0, the lowest correlation value would be seen as no connection)
graphNorm = (graph + abs(min(min(graph))) + 0.001);
% Normalise the correlation values, but keep the sign
graphNorm2 = graphNorm/max(max(abs(graphNorm)));
% Set the 0 connections back to 0
graphNorm2(find(graph==0)) = 0;

% There cannot be 0 values for line width, so set 0 to 0.0000000000001
line_widths(find(line_widths==0)) = 0.0000000000001;


% Draw lines depicting each ROI-ROI connection
[unique_colors, tmp, unique_color_inds] = unique(graph);
unique_color_inds=reshape(unique_color_inds, nroi, nroi);
cc=jet(numel(unique_colors)+1);
colormap(jet);
for i=1:nroi
    for j=i+1:nroi
        
        if (graph(i,j) ~= 0)
            
            h2=plot([x(1,i) x(1,j)],[y(1,i) y(1,j)],'color', cc(unique_color_inds(i,j),:), 'LineWidth', line_widths(i,j));
            
        end
            
    end
end

axis equal off
set(gca,'YTick',[])
set(gca,'YColor','w')


FIG = ['Circular_Graph_' G ];
hgsave(FIG)
print('-depsc',FIG)
print('-djpeg','-r300',FIG)
print('-dtiff','-r600',FIG)
print('-dpng','-r300',FIG)


