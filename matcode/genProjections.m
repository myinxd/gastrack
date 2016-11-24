function genProjections(snap_path,save_path,params,direction)
% genProjections(snap_path,save_path,params)
% This code generate all of the frames of the simulation
% 
% Input
% snap_path: path that snaps are saved
% save_path: path to save the images
% params: parameters of the simulation, [numhalo,numgas,step,wlim]
% direction: projection direction, 'x','y' or 'z'
% 
% Version 1.0
% Author: Zhixian MA <zxma_sjtu@qq.com>
% Date: 2016/11/15

% Init
warning off
files = dir(snap_path);

numhalo = params(1);
numgas = params(2);
step = params(3);
wlim = params(4);

NumSamples = length(files);
% Circulation
i = 3;
while i <= NumSamples
    % get id
    snap = files(i).name;
    if ~strcmp(snap(end-3:end),'hdf5')
        i = i + 1;
        continue;
    end
    disp(snap);
    temp_id = regexp(snap,'[0-9][0-9][0-9]');
    frameId = snap(temp_id:temp_id+2);
    
    % find peaks and the sub gas maps
    [peaklist,maj_gas,min_gas] = findpeaks([snap_path,snap],numgas,numhalo,step,wlim,direction);
    
    % combine
    scale = size(maj_gas,1);
    % major cluster
    gca1 = figure(1);
    set (gca1,'Position',[0,0,scale,scale]);
    imshow(flipud(log(maj_gas)),[],'border','tight');
    axis normal;
    colormap('cool')
    f1 = getframe(gca1);
    im1 = f1.cdata;
    % minor cluster
    gca2 = figure(2);
    set (gca2,'Position',[0,0,scale,scale]);
    imshow(flipud(log(min_gas)),[],'border','tight');
    axis normal;
    colormap('gray')
    f2 = getframe(gca2);
    im2 = f2.cdata;

    % resize
    try
        im = im1+im2;
    catch
        disp('Error')
        i  = i;
        continue
    end
    im = imresize(im,[800,800],'bilinear');

    figure(3);
    %set (gca,'Position',[0,0,800,800]);
    imshow(im,'XData',[-wlim/2,wlim/2],'YData',[wlim/2,-wlim/2])
    xlabel('x (Mpc)','fontsize',12)
    ylabel('y (Mpc)','fontsize',12)
    axis on
    hold on
    scatter(peaklist(1,1),peaklist(1,2),300,'bx','linewidth',1.5);
    hold on
    scatter(peaklist(2,1),peaklist(2,2),200,'r+','linewidth',1.5);
    set(gca,'YDir','normal')
    % Mark time
    text(-0.35,0.35,[num2str((str2double(frameId) * 0.02)),' Gyr'],'fontsize',12)
    
    % save
    fname = [save_path,'gas_',frameId,'_width_',num2str(wlim),'.png']; 
    saveas(gca,fname);
    
    cla(gca)
    i = i + 1;
end

close all