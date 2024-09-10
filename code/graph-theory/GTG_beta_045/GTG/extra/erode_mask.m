function [mask_out] = erode_mask(mask_in,rad)

% Erode binary mask_in with a sphere of radius = rad
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 10.06.14
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

%diam = (rad*2)-1;
diam = rad*2;
mask_out = imerode(mask_in,strel('ball',diam,diam,diam));
mask_out(mask_out<0.5)=0;
mask_out(mask_out>=0.5)=1;
