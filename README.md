# lip-vireo
This is a C++ project including traditional image keypoint detectors and descriptors. 

# Overview
In nowadays, traditional image keypoint detectors and descriptors are no longer in fashion. They have given place to the deep features, which become increasingly popular since 2012. However, the theories upon which the traditional detectors and descriptors are built are still shiny in nowadays. The scale-space theory from Dr. Tony Lindberg is beautiful. I learned a lot from his papers and book. It took me hundreds of weekends to implement all the detectors and descriptors in this project. This project also saw the birth of 'Flip invariant SIFT', which is my proposal. It well addressed the sensitivity issue of SIFT to flip transformation.

# Compile
##### This project can be smoothly compiled under Ubuntu and MacOS, while it is also possible to be compiled successfully (by MingW) under Windows with only small changes.
### Step 1. download libjpeg, zlib and libpng and compile them respectively, copy libjpeg.a, zlib.a and libpng.a to "libs/" folder
### Step 2. run "make release"

# Manual
##### A PDF manual about the overview about the detectors and descriptors, as well as the call of lip-vireo from commandline could be found under "manual" folder.

# Important Notice
##### DoG+SIFT has been patented by David G. Lowe in US (Patent No. US 6,711,293B1 [20]). The patent holder is University of British Columbia, Canada
##### Copyright holder for Fast Hessian and SURF descriptor are Tinne Tuytelaars and et al.;
##### Copyright holder for FIND descriptor are Mr. Xiao-Jie Guo and et al. from University of Tianjin, China.
