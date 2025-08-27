# User’s Manual of CoInfecProfiler

## Preamble
CoInfecProfiler is a novel tool that reconstructs haplotypes within-host from sequencing data and identifies virus co-infection by combining regularized regression and paired-end linkage. The workflow of CoInfecProfiler is briefed as follows: (A) Mutations are detected by another tool, e.g., GATK. (B) A regression equation is constructed using a reference matrix (which needs to be built by user) and mutations and paired-end linkage information. (C) The frequency of each haplotype is estimated by regularized regression.

## Installation
CoInfecProfiler.jar is a batteries-included JAR executable. All needed external jar packages are included in the downloadable, `CoInfecProfiler.jar`. To download all necessary files, users can use the command `git clone https://github.com/XXX/CoInfecProfiler.git` As we used an R package L0Learn, the users have to install R and L0Learn (https://cran.r-project.org/web/packages/L0Learn/index.html). The versions of R and R package L0Learn that we have used on our platform are: version 2.1.0 for L0Learn and version 4.3.1 for R. Users are also expected to have java (version: 1.8) on their platform.
Several other tools are prerequisites for running CoInfecProfiler. Users can download and install them from the websites:

- bwa: https://github.com/lh3/bwa
- samtools 1.4: http://www.htslib.org/download/
- GATK 4.3: https://software.broadinstitute.org/gatk/download/index

## Function
SARS-CoV-2-Reconstruction-VCF: Using SARS-CoV-2-Reconstruction-VCF to reconstruct the haplotypes and their frequency information in the samples.

## Quick start with included example data
Example data is provided in the "Example" folder. After updating absolute paths of executable (such as bwa etc) and parent folder in the param.txt file, users can run CoInfecProfiler by a simple commandline:

### Usage

``` bash
java -jar CoInfecProfiler.jar SARS-CoV-2-Reconstruction-VCF param.txt
```

CoInfecProfiler will output the final haplotype results in the **"output"** folder: `SARS-CoV-2-Variants.txt`

### CoInfecProfiler params file
```
BWA=/path/to/bwa
Samtools=/path/to/samtools
GATK=/path/to/gatk
Rscript_Path=/path/to/Rscript
SNP_Matrix=/path/to/reference_matrix
Fastq_1=/path/to/test_1.fastq
Fastq_2=/path/to/test_2.fastq
Output_Path=/path/to/out
Proj_Name=test
Reference=/path/to/ref.fasta
Min_Hap_Freq=0.01
SNV_Cutoff=0.01
Reconstruction_Start=1
Reconstruction_End=29903
Regression_Gamma_Min=0.001
Regression_Gamma_Max=1.0
Regression_n_Gamma=10
Regression_Lambda=0.0005
Sum_Weight=20.0
Regularized_Regression=L0L1
MAF_Weight=2.0
LD_Weight=0.5
AF_Zero_Weight=2.0
Min_LD_Reads=10
```
## Citation

## Contact

## Copyright License (MIT Open Source)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
