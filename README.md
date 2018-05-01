# Diffusional Kurtosis Estimator
The Diffusional Kurtosis Estimator (DKE) is software for processing 
Diffusional Kurtosis Imaging (DKI) data. 

## Downloading DKE
DKE is written primarily in MATLAB. It can be run either from executables 
or from source code. An advantage of running DKE from executables is that 
it can be run under the free MATLAB Compiler Runtime 2012a (MCR), whereas 
MATLAB is needed for running DKE from source code. An advantage of running 
DKE from source code is that you can modify your copy of the source code 
to suit your own purposes.

### Executables
The DKE executables for 32- and 64-bit Windows, 64-bit Linux, and 64-bit macOS 
can be downloaded from the MUSC CBI website. If you want to run DKE this way 
[click here to register and download DKE executables](http://academicdepartments.musc.edu/cbi/dki/dke-swreg.html).

### Source Code
The source code can be downloaded from this site using either the Downloads 
menu or via "git clone" (if you have a git client installed). For example, 
at the Linux or macOS command prompt:
```
$ git clone git@git.musc.edu:dal220/dke.git
```
If you download DKE via "git clone", you can later update to the latest 
version on this site with "git pull":
```
$ git pull
```

## Running DKE from source code in MATLAB
Most of the source code for DKE is in the mfiles subdirectory, which should be 
added to your MATLAB path. DKE makes use of MEX files (MATLAB executable files) 
for speed. Previously-compiled MEX files are in the subdirectories that end in 
"*mex". To use these files, add the appropriate subdirectory to your MATLAB 
path. For example, if you want to run DKE in MATLAB on macOS:
```
>> addpath /path/to/mfiles
>> addpath /path/to/mac64mex
>> dke /path/to/dke_parameters.dat
```

Note that the previously-compiled MEX files might not work on your computer, 
depending on the specific version of the operating system that you are using. 
In that case, you need to recompile the MEX files yourself with MATLAB Coder.

You can create MEX files for augment_constraints, compute_rd, and 
compute_rf by opening the corresponding .prj files in MATLAB. These files are 
in separate directories alongside the mfiles directory. MATLAB Coder will 
convert the .m files into C code and compile the C code into MEX files.

You can also create the MEX file for mlsei using MATLAB Coder. The source 
code is in C (converted from Fortran). Instructions are in mexlsei_setup.m in 
the mlsei/mfiles directory, but some of the steps have been done. Make the 
f2c library, move it to the "c" directory, and then create the MEX file. 
For example, on a 64-bit Linux computer (computer type = glnxa64), in MATLAB:
```
>> pwd

ans =

/home/username/src/dke/mlsei/c/libf2c

>> !make -f makefile.lx64
cc -c f77vers.c
cc -c i77vers.c
cc -c -DSkip_f2c_Undefs -O -DNON_UNIX_STDIO -fPIC main.c
...
ar: creating libf2cx64.a
ranlib libf2cx64.a

>> !mv libf2cx64.a ..

>> cd ../..

>> pwd

ans =

/home/username/src/dke/mlsei

>> mex -L "./c" -lf2cx64 -outdir . -v ./c/mexlsei.c ./c/dlsei.c
...
```

The MEX files can then be moved to the appropriate directory ending in "mex", 
for example linux64mex if you are using 64-bit Linux.


## Useful Links for DKE
[Introduction to DKI](http://academicdepartments.musc.edu/cbi/dki/index.html)

[DKI Protocols](http://academicdepartments.musc.edu/cbi/dki/protocols.html)

[DKI Data Processing](http://academicdepartments.musc.edu/cbi/dki/dke.html)

[Results](http://academicdepartments.musc.edu/cbi/dki/results.html)

[FAQs](http://academicdepartments.musc.edu/cbi/dki/faq.html)

[KIN](http://academicdepartments.musc.edu/cbi/dki/kin.html)

[References](http://academicdepartments.musc.edu/cbi/dki/references.html)

[Forum](https://www.nitrc.org/forum/?group_id=652)
