function [partlist,maj_idx,min_idx,maj_cls,min_cls] = getParticles(hdf5file,params,numhalo,numgas,step)
% [partlist,maj_idx,min_idx,maj_cls,min_cls] = getParticles(hdf5file,params,numhalo,numgas,step)
% A tool to get particles based on the provided region parameters.
% Input
% hdf5file: filepath of the hdf5.
% params: sector parameters, [radius_low,radius_high,angle_low,angle_high],
%         angles should be larger than 0, and with a unit of degree.
% numhalo: number of halo particles in the major cluster
% numgas: number of gas particles in the major cluster
% step: box length of the mosaic 
%
% Output
% maj_idx: indices of major particles in the region
% min_idx: indices of minor particles in the region
% maj_cls: the major cluster
% min_cls: the minor cluster
% 
% Version: 1.0
% Author: Zhixian MA <zxma_sjtu@qq.com>
% Date: 2016/11/24 Happy thanksgiving!

if nargin < 5
    step = 0.01;
end

%% Init
unit_len = 3.08568e21;
cm2Mpc = 3.240779289469756e-25;
radius_low = params(1); radius_high = params(2);
angle_low = params(3)/180*pi; angle_high = params(4)/180*pi;

if angle_high < angle_low
    disp('Wrong angles.')
    return
end

if angle_low >= 2 * pi
    angle_high = angle_high - (2 * pi)*floor(angle_low/(2*pi));
    angle_low = mod(angle_low,2 * pi);
end

if angle_high > 2 * pi
    angle_low = [angle_low,0];
    angle_high = [2*pi, angle_high-2*pi];
end

part_maj = 0;
part_min = 0;
part_all = 0;

% gas
gas_cords_path = getPath({'PartType0','Coordinates'});
gas_den_path = getPath({'PartType0','Density'});
gas_id_path = getPath({'PartType0','ParticleIDs'});
% halo
halo_cords_path = getPath({'PartType1','Coordinates'}) ;
halo_mass_path = getPath({'PartType1','Masses'});
halo_id_path = getPath({'PartType1','ParticleIDs'});

%% Clusters
gas_den = h5read(hdf5file,gas_den_path);
gas_cords = h5read(hdf5file,gas_cords_path);
gas_id = h5read(hdf5file,gas_id_path);
% halo
halo_mass = double(h5read(hdf5file,halo_mass_path)*10e8)/10e8;
halo_cords = h5read(hdf5file,halo_cords_path);
halo_id = h5read(hdf5file,halo_id_path);
% Transform
gas_cords = gas_cords * unit_len * cm2Mpc;
halo_cords = halo_cords * unit_len * cm2Mpc;

% cluster_major
gas_idx = (gas_id <= numgas-1);
maj_cls.gas_id = gas_id(gas_idx);
maj_cls.gas_den = gas_den(gas_idx);
maj_cls.gas_cords = gas_cords(:,gas_idx);
halo_idx = (halo_id <= numhalo+length(gas_id)-1);
maj_cls.halo_id = halo_id(halo_idx);
maj_cls.halo_mass = halo_mass(halo_idx);
maj_cls.halo_cords = halo_cords(:,halo_idx);
% cluster_minor
gas_idx = (gas_id > numgas);
min_cls.gas_id = gas_id(gas_idx);
min_cls.gas_den = gas_den(gas_idx);
min_cls.gas_cords = gas_cords(:,gas_idx);
halo_idx = (halo_id > numhalo+length(gas_id)-1);
min_cls.halo_idx = halo_id(halo_idx);
min_cls.halo_mass = halo_mass(halo_idx);
min_cls.halo_cords = halo_cords(:,halo_idx);

% Get mosaic and find peak in three directions
halo_z = getMosaic(maj_cls.halo_cords,maj_cls.halo_mass,step,'z');
halo_y = getMosaic(maj_cls.halo_cords,maj_cls.halo_mass,step,'y');
halo_x = getMosaic(maj_cls.halo_cords,maj_cls.halo_mass,step,'x');
% Find peak
[y_p1,x_p1] = find(halo_z==max(halo_z(:)));
y_p1 = (y_p1-1) * step - 1;
x_p1 = (x_p1-1) * step - 1;
[z_p1,x_p2] = find(halo_y==max(halo_y(:)));
z_p1 = (z_p1-1) * step - 1;
x_p2 = (x_p2-1) * step - 1;
[z_p2,y_p2] = find(halo_x==max(halo_x(:)));
y_p2 = (y_p2-1) * step - 1;
z_p2 = (z_p2-1) * step - 1;
x_p = (x_p1+x_p2)/2;
y_p = (y_p1+y_p2)/2;
z_p = (z_p1+z_p2)/2;

%% Statistic particles
% distances
maj_idx = (zeros(1,length(maj_cls.gas_id)));
min_idx = (zeros(1,length(min_cls.gas_id)));
% major
maj_x = maj_cls.gas_cords(1,:) - x_p;
maj_y = maj_cls.gas_cords(2,:) - y_p;
maj_z = maj_cls.gas_cords(3,:) - z_p;
maj_dist = sqrt(maj_x.^2 + maj_y.^2 + maj_z.^2);
% minor
min_x = min_cls.gas_cords(1,:) - x_p;
min_y = min_cls.gas_cords(2,:) - y_p;
min_z = min_cls.gas_cords(3,:) - z_p;
min_dist = sqrt(min_x.^2 + min_y.^2 + min_z.^2);

% angles
maj_ang = asin(abs(maj_y)./maj_dist);
min_ang = asin(abs(min_y)./min_dist);

for i = 1 :length(angle_low)
    % major
    maj_idx_dist = (maj_dist <= radius_high) .* (maj_dist>=radius_low);
    % Quarant1
    idx_q1 = (maj_x >= 0).*(maj_y>=0);
    idx_ang1 = (maj_ang <= angle_high(i)) .* (maj_ang >= angle_low(i));
    % Quarant2
    idx_q2 = (maj_x < 0).*(maj_y>=0);
    idx_ang2 = ((pi-maj_ang) <= angle_high(i)) .* ((pi-maj_ang) >= angle_low(i));
    % Quarant3
    idx_q3 = (maj_x < 0).*(maj_y< 0);
    idx_ang3 = ((pi+maj_ang) <= angle_high(i)) .* ((pi+maj_ang) >= angle_low(i));
    % Quarant4
    idx_q4 = (maj_x >= 0).*(maj_y<=0);
    idx_ang4 = ((2*pi-maj_ang) <= angle_high(i)) .* ((2*pi-maj_ang) >= angle_low(i));
    % Combine idx
    maj_idx_t = (idx_ang1.*idx_q1) + (idx_ang2.*idx_q2) + (idx_ang3.*idx_q3) + (idx_ang4.*idx_q4);
    maj_idx_t = logical(maj_idx_dist) .* logical(maj_idx_t);
    maj_idx = logical(maj_idx_t) + logical(maj_idx);
    maj_idx = logical(maj_idx);
    part_maj = part_maj + sum(maj_idx_t);
    
    % minor
    min_idx_dist = (min_dist <= radius_high) .* (min_dist>=radius_low);
    % Quarant1
    idx_q1 = (min_x >= 0).*(min_y>=0);
    idx_ang1 = (min_ang <= angle_high(i)) .* (min_ang >= angle_low(i));
    % Quarant2
    idx_q2 = (min_x < 0).*(min_y>=0);
    idx_ang2 = ((pi-min_ang) <= angle_high(i)) .* ((pi-min_ang) >= angle_low(i));
    % Quarant3
    idx_q3 = (min_x < 0).*(min_y< 0);
    idx_ang3 = ((pi+min_ang) <= angle_high(i)) .* ((pi+min_ang) >= angle_low(i));
    % Quarant4
    idx_q4 = (min_x >= 0).*(min_y<=0);
    idx_ang4 = ((2*pi-min_ang) <= angle_high(i)) .* ((2*pi-min_ang) >= angle_low(i));
    % Combine idx
    min_idx_t = (idx_ang1.*idx_q1) + (idx_ang2.*idx_q2) + (idx_ang3.*idx_q3) + (idx_ang4.*idx_q4);
    min_idx_t = logical(min_idx_dist) .* logical(min_idx_t);
    min_idx = logical(min_idx_t) + logical(min_idx);
    min_idx = logical(min_idx);
    part_min = part_min + sum(min_idx_t);
end

part_all = part_maj + part_min;
partlist = [part_all,part_maj,part_min];
