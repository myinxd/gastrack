# gastrack

This repository is the backup, in which a simple tool for tracking gas motion in galaxies collision (merging) is provided. In our work, the [Gadget](http://www.gadgetcode.org/) is utlized to simulate galaxies collision in the cluster namely `HCG 62`, and simulated results are saved as [hdf5](https://support.hdfgroup.org/HDF5/) structure. 

Two *PartTypes* are in this simulation, i.e. the gas and halo particles. We simulated the revolution of their massive density and motions, 
so as to analyze its merging staffs. Thus, to precisedly track the motion lines of the two merging gas, a gas tracker is needed.

## Method
Previously, the strategy of peak detection in two dimensional images is utlized. However, there were three parameters should be tuned, which were the `neighbors`,`peak threshold` and `number of peaks`. After several times of experiments, we can't find any group for all of them to get the optimum result, especially for the time that two galxies were nearly merged to one.

Since that, we change our method. We separate the whole simulated particles into two subsets, in which each contained the particles belonging to its galaxy. Then, the location of paticles with the largest density is detected, and combined into the peak list we required. Finally, the two singly detected peak are marked onto the gas map, and an example image is as follows. 

<center>
<img src="https://github.com/myinxd/gastrack/blob/master/example/gas_example.png" width=400 />
Fig1. Gas density projection example
</center>

## New tool 
In order to compare and track particles of different clusters between different snaps (i.e. times), a new tool is added. By means of this, particles in regions of interested (ROI) are tracked. Thus, the detailed motivations and revolutions of the merging can be analyzed.  

## Matlab code
Scripts written in matlab style is provided, which can be fetched from [matcode](https://github.com/myinxd/gastrack/tree/master/matcode). And usage of these codes is embedded in these scripts.
You can type `help func_name` at matlab's `commind window` for details.

```matlab
   >> help func_name
```

## How to use it?
- To use the code, some python packages are required. 
  - numpy
  - h5py: processing the `hdf5` files
  - yt: processing simulated particles, making projection and plotting results.

- Run the tool
  - For gas tracking
  ```sh
  $ python3 gastrack.py <input_dir> <output_dir>
  ```
  - For particles tracking
  ```sh
  $ python3 countParticles.py <snap_target> <snap_ref> <output_path> <numhalo> <num gas> \\
    <radius_low> <radius_high> <angle_low> <angle_high>
  ```

## References
- Collette, A.,"Python and HDF5", O'reilly, 2013.
- h5py documents, https://docs.h5py.org
- The yt project, https://yt-project.org
- matplotlib, http://matplotlib.org/

## Author
- Zhixian MA <`zxma_sjtu(at)qq.com`>

## License
Unless otherwise declared:

- Codes developed are distributed under the [MIT license](https://opensource.org/licenses/mit-license.php);
- Documentations and products generated are distributed under the [Creative Commons Attribution 3.0 license](https://creativecommons.org/licenses/by/3.0/us/deed.en_US);
- Third-party codes and products used are distributed under their own licenses.
