# Copyright (C) 2017 Zhixian MA <zxma_sjtu@qq.com>

"""
Update
======
[2017-09-08]: Hide labels and ticks of both x and y axes.
[2017-09-08]: Hide label of colorbar
[2017-09-08]: Add scale of the image at the southwest.

Reference
=========
http://matplotlib.org/examples/axes_grid/demo_edge_colorbar.html?highlight=cbar
"""
import os
import shutil
import yt
import numpy as np
import yt.units as units
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


# Custom packages
import gastrack
import partsplit

import sys
sys.setrecursionlimit(1000000)


def get_peak(fn,gt,ps,snap):
    ps.get_single_file(fn)
    # Init save names
    haloname = ("halo_Projection_z_density_%f_%d_0.png" %
                (gt.width[0], gt.buff_size[0]))
    gasname = ("snap_Projection_z_density_%f_%d_0.png" %
                (gt.width[0], gt.buff_size[0]))
    # get peaks
    gt.save = True
    width_old = gt.width
    gt.width = (2, 'Mpc')
    # cluster1
    f_c1 = ('subs/snap_%03d_c1.hdf5' % (snap))
    ds1 = gt.load_hdf5(f_c1)
    halo_c1 = ("halo_Projection_z_density_%f_%d_c1.png" %
                        (gt.width[0], gt.buff_size[0]))
    halo_map_c1 = gt.get_halo_map(ds1, halo_c1)
    peak_maj = gt.locate_peak(halo_map_c1)
    os.remove(os.path.join(gt.input_dir, f_c1))
    # cluster2
    f_c2 = ('subs/snap_%03d_c2.hdf5' % (snap))
    ds2 = gt.load_hdf5(f_c2)
    halo_c2 = ("halo_Projection_z_density_%f_%d_c2.png" %
               (gt.width[0], gt.buff_size[0]))
    halo_map_c2 = gt.get_halo_map(ds2, halo_c2)
    peak_min = gt.locate_peak(halo_map_c2)
    os.remove(os.path.join(gt.input_dir, f_c2))

    return peak_maj[0,:], peak_min[0,:]


def main(argv):

    # Parameters
    folder = ['./snaps','./snaps','./snaps']
    snap = int(argv[0])
    fns = ['snap_%03d_none.hdf5' % (snap),'snap_%03d_sfr.hdf5' % (snap),'snap_%03d_cool.hdf5' % (snap)]

    t = argv[1]
    case_flag = argv[2]
    # Defination of the axesgrid
    fig = plt.figure(figsize=(3,1))
    grid = AxesGrid(fig,rect=(0.1,0.1,0.8,0.8),
                    nrows_ncols = (1, 3),
                    axes_pad = 0.03,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="3%",
                    cbar_pad="0%")


    # case dict
    case_dict = {0:"Case A", 1: "Case B", 2: "Case C"}

    for i, fn in enumerate(fns):
        # p = None
        # Load the data and create a single plot
        print("Processing on %s" % (fn))
        fn = os.path.join(folder[i],fn)
        # fname = yt.load(fn) # load data
        unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
                    'UnitMass_in_g'            :   1.989e+43,
                    'UnitVelocity_in_cm_per_s' :      100000}
        bbox_lim = 1e5   # kpc
        bbox = [[-bbox_lim,bbox_lim],
            [-bbox_lim,bbox_lim],
            [-bbox_lim,bbox_lim]]  # limits
        ds = yt.load(fn,unit_base=unit_base,bounding_box=bbox)
        ds.index
        ad= ds.all_data()
        p = yt.ProjectionPlot(ds, 'z', ('deposit', 'PartType1_density'), width=(0.8, 'Mpc'),fontsize=10)

        # Ensure the colorbar limits match for all plots
        p.set_zlim('PartType1_density', 1e-3, 0.2)
        p.set_buff_size(3200)


        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[('deposit', 'PartType1_density')]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]

        # Finally, this actually redraws the plot.
        p._setup_plots()


    # Detect peaks and mark them
    gt = gastrack.GasTrack(input_dir='./',output_dir='./tmp')
    ps_output = os.path.join(gt.input_dir, 'subs')

    for i, fn in enumerate(fns):
        # Load the data and create a single plot
        gt = gastrack.GasTrack(input_dir=folder[i],output_dir='./tmp')
        if not os.path.exists(gt.output_dir):
            os.mkdir(gt.output_dir)
        ps_output = os.path.join(gt.input_dir,'subs')
        if os.path.exists(ps_output):
            shutil.rmtree(ps_output)
            os.mkdir(ps_output)
        else:
            os.mkdir(ps_output)
        ps = partsplit.PartSplit(gt.input_dir, ps_output)

        print("Processing on %s" % (fn))
        peak_maj,peak_min = get_peak(fn,gt,ps,snap)
        peak_maj.shape

        x = np.linspace(0.235,0.31,20)
        y = -0.33 * np.ones(x.shape)
        
        grid[i].axes.text(x=-0.35,y=0.3, s='%s=%.2f Gyr' %(t, snap*0.02),fontdict={"size":10,"color":"white","weight":'bold'})
        if case_flag == '1':
            grid[i].axes.text(x= 0.2,y=0.3, s=case_dict[i], fontdict={"size":10,"color":"white","weight":'bold'})
        grid[i].axes.text(x= 0.19, y=-0.3, s="100 kpc", fontdict={"size":10,"color":"white","weight":'bold'})
        grid[i].axes.plot(x,y,'w-')
        grid[i].scatter(peak_maj[1],peak_maj[2],s=100, c='b', marker='x')
        grid[i].scatter(peak_min[1],peak_min[2],s=80, c='r', marker='+')
        grid[i].set_xticks([])
        grid[i].set_ylabel("")
        grid[i].set_xlabel("")
        grid[i].set_yticks([])

    for cax in grid.cbar_axes:
        #cax.toggle_label(True)
        cax.axis[cax.orientation].set_label('')    
    # save image
    plt.savefig('./images/snap_%3d_comb.png' % snap,dpi=300)

if __name__ == "__main__":
    main(sys.argv[1:])
