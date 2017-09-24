## Synopsis

Matrix completion with entries over the set of reals has been popular in recommendation systems. One of the most efficient methods for solving this problem has been 
alternating minimization. However, when the entries are over discrete alphabets, then real valued relaxation may not lead to the correct underlying archetypes. This leads us
to the framework of Discrete Matrix Completion where an EM-like algorithm iterates over the real-valued matrix completion to further refine the results and give superior performance 
than naive rounding techniques. 

## Motivation

One of the main motivations for this problem is application to haplotype phasing in genomics. The input consists of noisy reads which provide insufficient information about the haplotypes.
The haplotypes can be considered to be similar to the varied types of feature vector combinations, representing categories of individual in a recommendation setting. 

## Stage-wise Alternating Minimization

There is a separate folder than implements another algorithm for matrix completion known as the stage-wise alternating minimization. This algorithm is not scalable for higher rank, however it 
has theoretical foundations and provable guarantees.  **Low-rank matrix completion using alternating minimization ** (STOC '13 Jain et al.)  

## DataSet

A sample dataset is provided in the read data incomplete matrix format. Considerable pre-processing is required for obtaining this data matrix from the reads, but is not the focus of the 
present work. In this project, we stress on the algorithmic techniques for haplotype reconstruction.

## Citation

If you use this code, please cite ** Resolving multicopy duplications de novo using polyploid phasing ** (RECOMB 2017, Chaisson et. al.)

## Running the Code

- Download the git folder and cd to `DiscreteMatCom`. Run `main_DMC('sim.0.0.0005')` in MATLAB. 
Or, if you want to run from terminal : ` $ matlab -nodesktop -nosplash -nojvm -r "main_DMC('sim.0.0.0005'); exit"  `

- cd to folder `StageAltMin`. Run `main_StageAlt('sim.0.0.0005')` in MATLAB. 
Or, if you want to run from terminal : ` $ matlab -nodesktop -nosplash -nojvm -r "main_StageAlt('sim.0.0.0005'); exit"  `



