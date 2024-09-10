function [DVARS,DVARS_norm,DVARS_precnorm] = clct_DVARS(data,mask,varargin)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 09.02.15
% 
% Calculates normalized DVARS, based on Tom Nichols' DVARS.sh script. 
% 
% Input:
% data     = 4d (time is 4th dimension) timeseries matrix
% mask     = 3d binary mask indicating what is brain
% varargin = enter 1 to use a robust estimator for the standard deviation,
%            enter 2 to use the traditional estimator
% 
% Output:
% DVARS          = unstandardized DVARS
% DVARS_norm     = standardized DVARS (the default from the Nichols script)
% DVARS_precnorm = precision-normalized DVARS 
%                  From the Nichols script: "Before taking the SD, temporal
%                                           difference images are 
%                                           standardized voxel-wise giving 
%                                           a more precisely normalized 
%                                           DVARS measure. A side effect, 
%                                           however, is that high-variance 
%                                           parts of the image are 
%                                           down-weighted relative to 
%                                           low-variance areas."
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if nargin==3
    std_type = varargin{1};
elseif nargin>3
    error('Too many input arguments')
else
    std_type = 1;
end

data = double(data);
mask = logical(mask);

if std_type==1
    sigm = iqr(data,4)./1.349;
elseif std_type==2
    sigm = var(data,0,4);
end

tdim = size(data,4);
d1mc = data(:,:,:,1:(tdim-1))-repmat(mean(data(:,:,:,1:(tdim-1)),4),[1,1,1,tdim-1]);
d2mc = data(:,:,:,2:tdim)-repmat(mean(data(:,:,:,2:tdim),4),[1,1,1,tdim-1]);
d    = (d1mc-d2mc);

ar1             = sum(d1mc.*d2mc,4)./sqrt(sum(d1mc.^2,4).*sum(d2mc.^2,4));
ar1(isnan(ar1)) = 0;
vdiff           = sqrt(2*(1-ar1)).*sigm;
svdiff          = mean(vdiff(mask));
prec_d          = (d./repmat(vdiff,[1,1,1,tdim-1])).^2;

DVARS          = nan(tdim-1,1);
DVARS_norm     = nan(tdim-1,1);
DVARS_precnorm = nan(tdim-1,1);
for t = 1:tdim-1
    ctp               = d(:,:,:,t).^2;
    prec_ctp          = prec_d(:,:,:,t);
    DVARS(t)          = sqrt(mean(ctp(mask)));
    DVARS_norm(t)     = DVARS(t)/svdiff;
    DVARS_precnorm(t) = sqrt(mean(prec_ctp(mask)));
end
