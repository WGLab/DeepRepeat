# DeepRepeat: Estimation of short tandem repeats by deep learning on Oxford Nanopore sequencing signals data 

In DeepRepeat, we reasonably assume that directly adjacent repeats share similar signal distribution. And then, we convert a repeat and its upstream and downstream repeats into RGB channels of a color image, feed images of repeats and non-repeats into a deep convolutional neural network to learn whether an event of signals are in a repeat or not, and determine repeat counts for individuals by aligning long reads against a reference genome and by summarizing and modeling the repeat prediction for a certain allele from multiple reads with Gaussian mixture distribution. 

# System requirements
## Hardware requirements

## Software requirements

# Installation

# General Usage

# Revision history
For release history, please visit [here](https://github.com/WGLab/DeepRepeat/releases). 

# Getting help
Please refer to the [DeepRepeat issue pages](https://github.com/WGLab/DeepRepeat/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing DeepRepeat
***Please cite the publication below if you use our tool***
Qian Liu, Li Fang, Alex Mas Monteys, Shauna A. Ebanks, Beverly L. Davidson, Kai Wang. Estimation of short tandem repeats by deep learning on Oxford Nanopore sequencing signals data. 
