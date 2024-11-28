# Module xChIPPy - X chromosome ChIPseq analysis in Python
Analysis of ChIPSeq and ATACseq data sets using genome annotation GTF and read alignment data.

Find and set the two path variables:

`DATA_FILE_PATH` needs to be set to a path for loading analysis data from and writing output to.
`ANNOTATION_FILE_PATH` needs to be set to a folder containing files with genomic coordinates of features and regions from the literature.

For performing a ChIPseq / ATACseq analysis download and unzip the compressed wig files from the link below into the corresponding folder.

https://polybox.ethz.ch/index.php/s/SGzXOB8JYQAt6m1

To load the precomputed analysis data, find the following line in the script and set LOAD_ANALYSIS:

`LOAD_CURRENT_ANALYSIS = True
`

To perform ChIPseq analysis set the following variables in the script:

`CHIPSEQ_ANALYSIS = True`

`ATACSEQ_ANALYSIS = False`

A ChIPseq analysis is represented in an `XchipAnalysis` object that contains read coverage data and annotation of the experiments. Core functions are implemented as members of this object. Separate functions are implemented for extensive analysis of ChIPseq and ATACseq data.

The following functions will load an existing analysis:

`analysis = XchipAnalysis()`

`analysis.loadAnalysis(DATA_FILE_PATH)`

To save an analysis use:

`analysis.saveAnalysis(DATA_FILE_PATH)`

This version of xpRNApy is based on Python 2.7 and was used with the following packages for performing the analysis:

    pandas.__version__
    u'0.24.2'
    numpy.__version__
    '1.16.6'
    seaborn.__version__
    '0.9.1'
    matplotlib.__version__
    '2.2.5'
    scipy.__version__
    '1.2.3'
    HTSeq.__version__
    '0.11.1'

