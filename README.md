# ISR_vol_img
Developers: Juan Araújo (UmU), Francisco López (UmU).\
UmU: Umeå University.

## Supporting code for the manuscript:
Efficient Computation and visualization of Ionospheric Volumetric Images for the enhanced interpretation of Incoherent Scatter Radar Data,
Juan Araújo<sup>c</sup>, Francisco López, Stefan Johansson, Assar Westman, Madelen Bodin, 2024.\
<sup>c</sup> Correspondance: juan.araujo@umu.se 

### Abstract
Incoherent scatter radar (ISR) techniques provide reliable measurements for the analysis of ionospheric plasma. 
Recent developments in ISR technologies allow the generation of high-resolution 3D data. Examples of such technologies employ the so-called phased-array antenna systems like the AMISR in North America or the upcoming EISCAT_3D in the Northern Fennoscandia region. 
EISCAT3D will be capable of generating the highest resolution ISR datasets that have ever been measured. However, running EISCAT3D experiments will be costly in terms of energy consumption and staffing. To use these resources in the most efficient manner, we present a novel computational strategy for the generation of high-resolution and smooth volumetric ionospheric images that represent ISR data. Through real-time processing, our computational framework will enable a fast decision-making during the monitoring process, where the experimental parameters are adapted in real time as the radars monitor specific phenomena. We describe our
strategy, which implements a flexible mesh generator along with an efficient interpolator specialized for ISR technologies. The proposed strategy is generic in the sense that it can be applied to a large variety of data sets and supports the interactive visual analysis and exploration of ionospheric data supplemented by interactive data transformations and filters.

## Software requirements (Linux systems)
+ g++ version 9.4.0, https://gcc.gnu.org/
+ HDF5 C++ API 1.12.2, https://docs.hdfgroup.org/archive/support/HDF5/doc1.8/cpplus_RM/index.html
+ OpenGL Mathematics (GLM), https://www.opengl.org/sdk/libs/GLM/

## ISR data
The data used for this application is in .hdf5 format and must be stored under the directory data/\
Download AMISR data from http://cedar.openmadrigal.org/, following the instructions provided there.\
In the project, we illustrate our strategy applied to data from the Resolute Bay ISR (RISR-C) and Poker Flat (PFISR).

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
The computation of the volumetric images utilizes a specialized version of the so-called Shepard interpolation:
+ Donald Shepard. _A two-dimensional interpolation function for irregularly-spaced data_. In Proceedings of the 1968 23rd ACM national conference, pages 517–524, 1968.

In the project, we extended the Shepard algorithm for including: gradient estimation plus correction, filtering of missing data, multi-scalar interpolation, efficiency and parallelization as well as vectorized operations at processor level. All these are novel contributions which are briefly described in our manuscript.

+ In our implementation, we modified the Shepard implementation from John Burkardt,\
https://people.sc.fsu.edu/~jburkardt/

+ For Convex-Hull computations, we used the implementation by Miguel Vieira,\
https://github.com/MiguelVieira

