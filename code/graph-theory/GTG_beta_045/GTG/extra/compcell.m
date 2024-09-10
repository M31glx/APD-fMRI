function Acomp = compcell(A1,A2)

% Compare contents of two cell arrays. 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 03.07.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if any(size(A1)~=size(A2))
    error('Array sizes must be equal')
elseif length(size(A1))>4
    error('Arrays with greater than 4 dimensions are not supported')
elseif iscell(A1{1,1,1,1})
    error('Cell arrays within the array are not supported')
end

Acomp = false(size(A1));
for dim1 = 1:size(A1,1)
    for dim2 = 1:size(A1,2)
        for dim3 = 1:size(A1,3)
            for dim4 = 1:size(A1,4)
                Acomp(dim1,dim2,dim3,dim4) = ~any(A1{dim1,dim2,dim3,dim4}~=A2{dim1,dim2,dim3,dim4});
            end
        end
    end
end
