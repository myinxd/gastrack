function [peaklist,maj_gas,min_gas] = findpeaks(hdf5file,NumGas,NumHalo,step,wlim,direction,SaveFlag)
% [peaklist,cls_maj,cls_min] = findpeaks(hdf5file,NumGas,NumHalo,SaveFlag)
% Find peaks of the two clusters
%
% Input
% hdf5file: file name of the hdf5
% NumGas: Number of gas particles in the major cluster
% NumHalo: Number of halo particles in the minor cluster
% SaveFlag: if equals to one, then save the information of clusters
%
% Output
% peaklist: list of peaks
% cls_maj: particles of major
% cls_min: particles of minor
%
% Version: 1.0
% Author: Zhixian MA <zxma_sjtu@qq.com>
% Date: 2016/11/15

if nargin < 7
    SaveFlag = 0;
elseif nargin < 6
    SaveFlag = 0;
    direction = 'z';
end

%% Init
% get id
Id = regexp(hdf5file,'[0-9][0-9][0-9]'); 
FrameId = hdf5file(Id:Id+2);
unit_mass = 1.989E43;
unit_len = 3.08568e21;
cm2Mpc = 3.240779289469756e-25;
% gas
gas_cords_path = getPath({'PartType0','Coordinates'});
gas_den_path = getPath({'PartType0','Density'});
gas_id_path = getPath({'PartType0','ParticleIDs'});
% halo
halo_cords_path = getPath({'PartType1','Coordinates'}) ;
halo_mass_path = getPath({'PartType1','Masses'});
halo_id_path = getPath({'PartType1','ParticleIDs'});

%% Load
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

%% Split
% cluster_major
cls_maj.gas_id = find(gas_id <= NumGas-1);
cls_maj.gas_den = gas_den(cls_maj.gas_id);
cls_maj.gas_cords = gas_cords(:,cls_maj.gas_id);
cls_maj.halo_id = find(halo_id <= NumHalo+length(gas_id)-1);
cls_maj.halo_mass = halo_mass(cls_maj.halo_id);
cls_maj.halo_cords = halo_cords(:,cls_maj.halo_id);
% cluster_minor
cls_min.gas_id = find(gas_id > NumGas);
cls_min.gas_den = gas_den(cls_min.gas_id);
cls_min.gas_cords = gas_cords(:,cls_min.gas_id);
cls_min.halo_id = find(halo_id > NumHalo+length(gas_id)-1);
cls_min.halo_mass = halo_mass(cls_min.halo_id);
cls_min.halo_cords = halo_cords(:,cls_min.halo_id);
% save
if SaveFlag == 1
    fname = ['Snap_',FrameId,'_cls'];
    save(fname,'cls_maj','cls_min');
end

%% Get mosaic
% cluster major
% step = 0.02;
maj_gas = getMosaic(cls_maj.gas_cords,cls_maj.gas_den,step,direction);
maj_halo = getMosaic(cls_maj.halo_cords,cls_maj.halo_mass,step,direction);
% cluster minor
min_gas = getMosaic(cls_min.gas_cords,cls_min.gas_den,step,direction);
min_halo = getMosaic(cls_min.halo_cords,cls_min.halo_mass,step,direction);

%% Find peaks
[peak_maj] = max(maj_halo(:));
[maj_row,maj_col] = find(maj_halo==peak_maj);
maj_row = maj_row * step - 1;
maj_col = maj_col * step - 1;

[peak_min] = max(min_halo(:));
[min_row,min_col] = find(min_halo==peak_min);
min_row = min_row * step - 1;
min_col = min_col * step - 1;

p1 = [maj_col,maj_row];
p2 = [min_col,min_row];
peaklist = [p1;p2];

%% Cut
if wlim < 2
    idx_diff = fix((1 - wlim/2)/step);
    maj_gas = maj_gas(idx_diff+1:end-idx_diff,idx_diff+1:end-idx_diff);
    min_gas = min_gas(idx_diff+1:end-idx_diff,idx_diff+1:end-idx_diff);
end























