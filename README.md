# DeepRepeat: Estimation of short tandem repeats by deep learning on Oxford Nanopore sequencing signals data 

In DeepRepeat, we reasonably assume that directly adjacent repeats share similar signal distribution. And then, we convert a repeat and its upstream and downstream repeats into RGB channels of a color image, feed images of repeats and non-repeats into a deep convolutional neural network to learn whether an event of signals are in a repeat or not, and determine repeat counts for individuals by aligning long reads against a reference genome and by summarizing and modeling the repeat prediction for a certain allele from multiple reads with Gaussian mixture distribution. 

# System requirements
## Hardware requirements
There is specific hardware requirements to use DeepRepeat if you can successfully install all dependent packages. If you want to train your own DeepRepeat models, you need a computing node with tens of GB memory.

## Software requirements
Please refer to `environment.yml` for detail. For your quick reference, DeepRepeat needs
```
  - python=3.6
  - hdf5=1.10.1
  - htslib=1.9
  - scikit-learn
  - tensorflow=1.9
  - samtools
  - minimap2
```

# Installation
It is easy to install the dependent packages of DeepRepeat using `annoconda`. Thus, please install `annoconda` first, and then follow the commands below to install DeepRepeat.

```
git clone https://github.com/WGLab/DeepRepeat
cd DeepRepeat
conda env create -f environment.yml
source activate py36deeprepeat   #if you change conda env name, please replace `py36deeprepeat`
cd bin/scripts
export DR_conda_base="./" #"/home/liuq1/anaconda2"  #replace this folder for your own annoconda folder
g++ -O3 -std=c++11 -o IndexF5files ComFunction.c Fast5Index.c IndexF5files.c
h5c++ -O3 -std=c++11 -I $DR_conda_base/include -L$DR_conda_base/lib -lhts -o genomic1FE ComFunction.c ComOption.c BamReader.c Fast5Index.c Fast5Reader.c RepeatFeatExtract.c genomic1FE.c
cd ../../
```

Then, you can run `python DeepRepeat.py`


# General Usage

# Revision history
For release history, please visit [here](https://github.com/WGLab/DeepRepeat/releases). 

# Getting help
Please refer to the [DeepRepeat issue pages](https://github.com/WGLab/DeepRepeat/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing DeepRepeat
***Please cite the publication below if you use our tool***

Qian Liu, Li Fang, Alex Mas Monteys, Beverly L. Davidson, Kai Wang. Direct detection and quantification of short tandem repeats 
on ionic signal data from Nanopore sequencing.


