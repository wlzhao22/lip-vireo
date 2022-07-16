# lip-vireo
This is a C++ project including traditional image keypoint detectors and descriptors. 

# Overview
In nowadays, traditional image keypoint detectors and descriptors are no longer in fashion. They have given place to the deep features, which become increasingly popular since 2012. However, the theories upon which the traditional detectors and descriptors are built are still shiny in nowadays. The scale-space theory from Dr. Tony Lindberg is beautiful. I learned a lot from his papers and book. It took me hundreds of weekends to implement all the detectors and descriptors in this project. This project also saw the birth of 'Flip invariant SIFT', which is my proposal. It well addressed the sensitivity issue of SIFT to flip transformation. The project name was given by my Ph.D supervisor as this project was started during my PhD in VIREO research group in City University of Hong Kong. 

The included detectors are Hessian, laplacian-of-Gaussian, Harris and SURF. The included descriptors are SIFT, LJET, Flip-invariant SIFT, RIFT, and SPIN. Please find out more details from the manual.

# Detectors and Descriptors
#### Detectors: Difference-of-Gaussian (DoG), Laplacian-of-Gaussian (LoG), Hessian, Hessian-Laplacian, Harris, Harris-Laplacian, SURF, and Hessian-Affine
#### Descriptors: SIFT [1], F-SIFT [2], SURF [3], FIND, SPIN, LJET, and RIFT

###### [1] David G. Lowe: Distinctive Image Features from Scale-Invariant Keypoints, IJCV'04.
###### [2] Wan-Lei Zhao, Chong-Wah Ngo: Flip-invariant SIFT for Copy and Object Detection, TIP'13. 
###### [3] Herbert Bay, Tinne Tuytelaars, Luc Van Gool: SURF: Speeded up robust features, ECCV'06

# Compile
##### This project can be smoothly compiled under Ubuntu and MacOS, while it is also possible to be compiled successfully (by MinGW) under Windows with only small changes. For Windows, user should install MinGW first, then edit "iotool.cpp", comment out the following lines in the file

```code
#ifndef LINUX
#define LINUX
#endif
```

#### Step 1. Download libjpeg, zlib and libpng packages and compile them respectively, copy libjpeg.a, zlib.a and libpng.a to "libs/" folder
##### for MacOS, please use jpeg_MacOS.zip, for Ubuntu 16.x or later please use jpeg_MacOS.zip

#### Step 2. Run make install for zlib and libpng on your system
#### Step 3. Run "make release" under the folder as follows

```shell
$ cd lip-vireo
$ make release
```

##### The prerequisite packages could be found under the "packages" folder or from the following URLs
###### 1. https://zlib.net/zlib-1.2.12.tar.xz
###### 2. https://udomain.dl.sourceforge.net/project/libpng/libpng16/1.6.37/libpng-1.6.37.tar.xz
###### 3. http://sources.buildroot.net/libjpeg/

# Manual
##### A PDF manual of the overview about the detectors and descriptors, as well as the call of lip-vireo from commandline could be found under "manual" folder.

# Copyright Issues
##### 1. DoG+SIFT has been patented by David G. Lowe in US (Patent No. US 6,711,293B1 [20]). The patent holder is University of British Columbia, Canada
##### 2. Copyright holder for Fast Hessian and SURF descriptor are Tinne Tuytelaars and et al.;
##### 3. Copyright holder for FIND descriptor are Mr. Xiao-Jie Guo and et al. from University of Tianjin, China.


# Author
##### Wan-Lei Zhao
###### He was born in a village in Yunnan province, China. He loves to share his work.
