function ADJ2NetGrh(CIJ, CIJ_Coordinate, fname, arcs)
% function ADJ2NetGrh(CIJ, CIJ_Coordinate, fname, arcs, threshold)
% writes a Pajek .net files from a MATLAB matrix
% E.g. ADJ2NetGrh(CIJ, CIJ_Coordinate, fname, arcs, threshold)
% Inputs:
%           CIJ = adjacency matrix
%           CIJ_Coordinate = coordinate matrix of nodes
%           fname = filename minus .net extension
%           arcs = 1 produces a directed net, arcs = 0 an undirected net
%           threshold = threshold of the extraction
%           
% Edited by Lijie Huang, 2010/4/13
%
%==========================================================================
N = length(CIJ);
CIJ_temp = triu(CIJ - eye(N));

% commented by Ashkan, no needed for threshold.
% for i = 1:N
%     for j = i:N
%         if CIJ_temp(i,j) < threshold
%             CIJ_temp(i,j) = 0;
%         end
%     end
% end

N = size(CIJ_temp,1);
fid = fopen(cat(2,fname,'.net'), 'w');

%%%VERTICES
fprintf(fid, '*vertices %6i \r', N);
for i = 1:N
    fprintf(fid, '%6i "%6i" ', [i i]);
    fprintf(fid, '%6f %6f %6f \r',[CIJ_Coordinate(i,2) CIJ_Coordinate(i,1) CIJ_Coordinate(i,3)]);
end

%%%ARCS/EDGES
if arcs
    fprintf(fid, '*arcs \r');
else
    fprintf(fid, '*edges \r');
end

for i = 1:N
    for j = 1:N
        if CIJ_temp(i,j) ~= 0
            fprintf(fid, '%6i %6i %6f \r', [i j CIJ_temp(i,j)]);
        end
    end
end

fclose(fid)

currdir = pwd;
netname = strcat(currdir,'\',cat(2,fname,'.net'));
%{
flog = fopen('C:\pajek\Pajek\Pajek.log','w');
fprintf(flog, 'NETBEGIN 1 \r');
fprintf(flog, 'CLUBEGIN 1 \r');
fprintf(flog, 'PERBEGIN 1 \r');
fprintf(flog, 'CLSBEGIN 1 \r');
fprintf(flog, 'HIEBEGIN 1 \r');
fprintf(flog, 'VECBEGIN 1 \r\r');
fprintf(flog, 'Msg Reading Network   ---    %s \r', netname);
fprintf(flog, 'N 1 RDN "%s" (111) \r', netname);
fprintf(flog, 'E 1 DRAW 0 0 0 0 0 \r');
fprintf(flog, 'E 1 BITMAP 0 0 0 0 0 "%s" \r', strcat(currdir,'\','X.bmp'));
fprintf(flog, 'EXIT \r');
fclose(flog)
dos('runPajek.bat');
%}
end

