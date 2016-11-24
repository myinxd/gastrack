function [maj_same,min_same] = cmpRegions(reffile,cmpfile,params_ref,params_cmp,params_sim)
% [maj_same,min_same] = cmpRegions(reffile,cmpfile,parms_ref,parms_cmp,params_sim)
% The combined function to track particles between different times
%
% Input
% reffile: path of the referred time
% cmpfile: path of the compared time
% params_ref: sector params of the referred one
% params_cmp: sector params of the compared one
% params_sim: simulation parameters, [numhalo,numgas,step]
%
% output
% maj_same: number of same particles in the major cluster
% min_same: number of same particles in the minor cluster
%
% Version: 1.0
% Date: 2016/11/24
% Author: Zhixian MA <zxma_sjtu@qq.com>

% Init
numhalo = params_sim(1);
numgas = params_sim(2);
if length(params_sim) == 3
    step = params_sim(3);
else
    step = 0.01;
end

% Times
idx_ref = regexp(reffile,'[0-9][0-9][0-9]');
snap_ref = str2num(reffile(idx_ref:idx_ref+2));
idx_cmp = regexp(cmpfile,'[0-9][0-9][0-9]');
snap_cmp = str2num(cmpfile(idx_cmp:idx_cmp+2));

% Search particles
[part_ref,maj_idx_ref,min_idx_ref,maj_cls_ref,min_cls_ref] = getParticles(reffile,params_ref,numhalo,numgas,step);
[part_cmp,maj_idx_cmp,min_idx_cmp,maj_cls_cmp,min_cls_cmp] = getParticles(cmpfile,params_cmp,numhalo,numgas,step);

% Compare
maj_same = cmp_id(maj_cls_ref,maj_cls_cmp,maj_idx_ref,maj_idx_cmp);
min_same = cmp_id(min_cls_ref,min_cls_cmp,min_idx_ref,min_idx_cmp);

% Disp
fprintf('Total particles at %.2f Gyr: %d\n',snap_ref * 0.02, part_ref(1));
fprintf('Major particles at %.2f Gyr: %d\n',snap_ref * 0.02, part_ref(2));
fprintf('Minor particles at %.2f Gyr: %d\n',snap_ref * 0.02, part_ref(3));
fprintf('\n');

fprintf('Total particles at %.2f Gyr: %d\n',snap_cmp * 0.02, part_cmp(1));
fprintf('Major particles at %.2f Gyr: %d\n',snap_cmp * 0.02, part_cmp(2));
fprintf('Minor particles at %.2f Gyr: %d\n',snap_cmp * 0.02, part_cmp(3));
fprintf('Major particles from %.2f Gyr: %d\n',snap_ref * 0.02, maj_same);
fprintf('Minor particles from %.2f Gyr: %d\n',snap_ref * 0.02, min_same);
end

function numsame = cmp_id(cls1,cls2,idx1,idx2)
% function cmp_id(cls1,cls2,idx1,idx2)

partId1 = cls1.gas_id(idx1);
partId2 = cls2.gas_id(idx2);
sameId = intersect(partId1,partId2);

numsame = length(sameId);

end