# ISR_vol_img
Developers: Juan Araújo (UmU), Francisco López (UmU). UmU: Umeå University.

## Supporting code for the manuscript:
Efficient Computation and visualization of Ionospheric Volumetric Images for the enhanced interpretation of Incoherent Scatter Radar Data,
Juan Araújo<sup>c</sup>, Francisco López, Stefan Johansson, Assar Westman, Madelen Bodin\
<sup>c</sup> Correspondance: juan.araujo@umu.se 

### Abstract
Incoherent scatter radar (ISR) techniques provide reliable measurements for the analysis of ionospheric plasma, which are derived by employing antennas that transmit and receive radio waves.
Recent developments in ISR technologies are capable of generating high-resolution 3D data. Examples of such technologies employ the so-called phased-array antenna like the AMISR in north
America or the upcoming EISCAT 3D (E3D) in the northern Fennoscandia region. E3D will be capable of generating the highest resolution ISR datasets that have been ever measured. However, running E3D experiments will be costly in terms of energy consumption and staffing. To use these resources in the most efficient manner, we present a novel computational strategy for the generation of high-resolution and smooth volumetric ionospheric images that represent ISR data. Through real-time visualization, our computational framework will enable a fast decision-making during the monitoring process, where the experimental parameters are adapted in real time as the radars monitor specific phenomena. The proposed strategy employs an effective mesh generator along with an efficient interpolator specialized for ISR technologies. An interactive framework is designed for the exploration and visualization of spatio-temporal and volumetric ionospheric images. The strategy is targeted at offering the analyst a wider range of alternatives in order to interpret ISR data. The proposed interpolation strategy is explicit, it takes account of the challenging cases of missing data and multi-scalar interpolation. Furthermore, in our computations, we have structured data allocation to minimize cache misses and to promote efficient vectorized operations. Finally, we efficiently parallelized our implementation for multi-core architectures with OpenMP and computations confirm outstanding speed-ups of our strategy with execution times that allows for real-time usage. We describe our implementation and demonstrate the interactive visual exploration of ionospheric data supplemented by interactive data transformations and filters. The visualization strategy is generic in the sense that can be applied to a large variety of data sets, supports the interactive visual analysis, and allows for effective navigation in the volume image.

## Software requirements (Linux systems)
+ g++ version 9.4.0, https://gcc.gnu.org/
+ HDF5 C++ API 1.12.2, https://docs.hdfgroup.org/archive/support/HDF5/doc1.8/cpplus_RM/index.html
+ OpenGL Mathematics (GLM), https://www.opengl.org/sdk/libs/GLM/

## ISR data
The data used for this application is in .hdf5 format and must be stored under the directory data/\
Download AMISR data from http://cedar.openmadrigal.org/, following the instructions provided there.

Examples of data files used in the manuscript:\
pfa150317.004.hdf5, pfa170128.002.hdf5, pfa200116.002.hdf5, ras161121.002.hdf5, ras190510.004.hdf5, ras200113.004.hdf5

## Compilation  (Linux systems):
+ g++ -std=c++17 -Wall -O3 -mavx -mfma -mavx2 -march=native -c interpolation_multi_parallel.cpp -o interpolation_multi_parallel.o -fopenmp;

+ g++ -std=c++17 -Wall -I/usr/lib/x86_64-linux-gnu/hdf5/serial/lib/include -I/src -I . -I ./src interpolation_multi_parallel.o src/*.cpp main.cpp -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -o main -lhdf5_cpp -lhdf5 -fopenmp;



## Execution:
./main -F [filename] -[options]

### Options:
+ n: number of 1D interpolation points. Total interpolation points is the cube of this parameter.
+ v: verbose on interpolation process
+ t: number of threads in parallel interpolation
+ p: Shepard parameter

### Examples:
+ ./main -F ras161121.002.hdf5
+ ./main -F ras161121.002.hdf5 -n 32
+ ./main -F ras161121.002.hdf5 -n 32 -p 3 
+ ./main -F ras161121.002.hdf5 -n 32 -p 3 -t 8

## Output:
The process of the computation of the volumetric images is output on the terminal prompt.
The result of the computation of the volumetric images is written in .vtk files and stored under the directory output/ \
Output .vtk files can be visualized on the Paraview software (https://www.paraview.org/) using volume rendering.

### Acknowledgment:
The computation of the volumetric images utilizes a specialized version of the so-called Shepard interpolation \
Donald Shepard. A two-dimensional interpolation function for irregularly-spaced data. In Proceedings of the 1968 23rd ACM national conference, pages 517–524, 1968. \
In the project, we extended the Shepard algorithm for including: gradient estimations, filtering of missing data, efficiency and parallelization as well as vectorized operations at processor level. All these are novel contributions which are briefly described in our manuscript.

+ In our implementation, we modified the Shepard implementation from John Burkardt,  
https://people.sc.fsu.edu/~jburkardt/

+ For Convex-Hull computations, we used the implementation by Miguel Vieira, 
https://github.com/MiguelVieira

