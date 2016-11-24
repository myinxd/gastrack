function [DensityMat] = getMosaic(cords,data,step,direction)
% [DensityMat] = getMosaic_new(cords,data,peakflag,wlim,step,direction)
% This code translate the density vector into mosaic matrix.
%
% Input
% data
% cords
% wlim: limitation of width, e.g. 2 [Unit: Mpc]
% step: moscaic width, e.g. 0.1 Mpc
% direction: direction of projection, default as 'z'
%
% Output:
% DensityMat: sparse matrix of the field
% AxisLim: limitation of axes
%
% Version: 1.0
% Author: Zhixian MA <zxma_sjtu@qq.com>
% Date: 2016/11/15

if nargin < 3
    step = 0.02;
    direction = 'z';
elseif nargin < 4
   direction = 'z';
end

% Init
cord_x = cords(1,:);
cord_y = cords(2,:);
cord_z = cords(3,:);

% Scarlarize
wlim = 2 / 2;
id_x = (cord_x > -wlim) .* (cord_x <= wlim);
id_y = (cord_y > -wlim) .* (cord_y <= wlim);
id_z = (cord_z > -wlim) .* (cord_z <= wlim);
idx = id_x .* id_y .* id_z;

idx = find(idx==1);
cord_x = cord_x(idx);
cord_y = cord_y(idx);
cord_z = cord_z(idx);
data = data(idx);

if direction == 'z'
    % Binned
    cord_col = cord_x; 
    cord_row = cord_y;
elseif direction == 'x'
    % Binned
    cord_col = cord_y; 
    cord_row = cord_z;
elseif direction == 'y'
    % Binned
    cord_col = cord_x; 
    cord_row = cord_z;
else
    printf('Error in direction mode.\n')
    return
end

% Mosaic and projection
bin = -wlim:step:wlim;
scale = length(bin)-1;
DensityMat = zeros(scale,scale);
for i = 1 : scale
    bin_idx = find(cord_col >= bin(i) & cord_col < bin(i+1));
    temp_row = cord_row(bin_idx);
    temp_data = data(bin_idx);
    for j = 1 : scale
        bin_idy = find(temp_row >= bin(j) & temp_row < bin(j+1));
        DensityMat(j,i) = sum(temp_data(bin_idy))/(step^2);
    end
end