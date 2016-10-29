# gastrack

This repository is the backup, in which a simple tool for tracking gas motion in galaxies collision (merging) is provided. In our work, the [Gadget](http://www.gadgetcode.org/) is utlized to simulate galaxies collision in the cluster namely `HCG 62`. 

Two *PartTypes* are in this simulation, i.e. the gas and halo particles. We simulated the revolution of their massive density and motions, 
so as to analyze its merging staffs. Thus, to precisedly track the motion lines of the two merging gas, a gas tracker is needed.

## Method
Previously, the strategy of peak detection in two dimensional images is utlized. However, there were three parameters should be tuned, which were the `neighbors`,`peak threshold` and `number of peaks`. After several times of experiments, we can't find any group for all of them to get the optimum result, especially for the time that two galxies were nearly merged to one.

Since that, we change our method. We separate the whole simulated particles into two subsets, in which each contained the particles belonging to its galaxy. Then, the location of paticles with the largest density is detected, and combined into the peak list we required. Finally, the two singly detected peak are marked onto the gas map. 

## How to use it?
- To use the code, some python packages are required. 
  - numpy
  - h5py: processing the `hdf5` files
  - yt: processing simulated particles, making projection and plotting results.

- Run the tool
  ```sh
  $ python3 gastrack.py <input_dir> <output_dir>
  ```

## References
- Collette, A.,"Python and HDF5", O'reilly, 2013.
- h5py documents, https://docs.h5py.org
- The yt project, https://yt-project.org

## Author
- Zhixian MA <`zxma_sjtu(at)qq.com`>

## License
Unless otherwise declared:

- Codes developed are distributed under the [MIT license](https://opensource.org/licenses/mit-license.php);
- Documentations and products generated are distributed under the [Creative Commons Attribution 3.0 license](https://creativecommons.org/licenses/by/3.0/us/deed.en_US);
- Third-party codes and products used are distributed under their own licenses.
