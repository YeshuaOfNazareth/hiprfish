# HiPR-FISH Probe Design and Image Analysis 

## Acknowledgement
This suite of code makes use of open source packages, including `numpy`, `pandas`, `biopython`, `bioformats`, `javabridge`, `scikit-image`, `scikit-learn`, and `scipy`.

## HiPR-FISH Image Analysis
Image analysis pipelines and scripts for HiPR-FISH experiments

### Overview

This pipeline enables automated image analysis for highly multiplexed FISH experiments on microbial communities. In most cases, the main pipeline is a snakemake workflow. There are also standalone scripts used for specific analyses presented in our paper.

## HiPR-FISH Probe Design
Probe design pipeline for HiPR-FISH experiments

### Acknowledgements
We would like to thank Jakob Wirbel for their help with testing the probe design pipeline. 

### Overview

This pipeline enables design of complex oligo probe sets used for highly multiplexed FISH experiments on microbial communities. The main pipeline is a snakemake workflow. There are two versions of the pipeline. The `hiprfish-probe-design-consensus` version uses the consensus approach by designing probes from the taxon consensus sequence for each taxon. The `hiprfish-probe-design-molecule` version designs probes from each individual 16S molecule from PacBio sequencing datasets and pool all unique probes for subsequent evaluation. The probe evaluation and selection is identifical in either version.  

### Required resources

The pipeline requires a local copy of the 16SMicrobial database from NCBI.

### Before running the pipeline
0. This pipeline was run on Ubuntu.

1. Download [Miniconda](https://repo.anaconda.com/miniconda/Miniconda3-py38_4.11.0-Linux-x86_64.sh). Install it via the following commands and type yes if needed:
      `cd downloads`\
      `bash Miniconda3-py38_4.11.0-Linux-x86_64.sh (if you get errors here, verify the file name.)`\
      `echo "PATH=\$PATH:/home/TYPE YOUR USERNAME HERE/.local/bin" >> ~/.bashrc`

2. Install the environment by running the following commands and typing "yes" if you have a choice:
     `conda create -n hiprfish python=3.8`\
     `conda activate hiprfish`\
     `conda install pandas=1.3.1`\
     `conda install -c bioconda primer3 -y`\
     `conda install -c anaconda joblib -y`\
     `conda install -c anaconda biopython -y`\
     `pip install snakemake`\
     `pip install ete3`\
     `pip install SetCoverPy`\
     `pip install tables`\
     `pip install openpyxl`\
     `pip install matplotlib`\
     `python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose`
     
3. Create a folder on the desktop and label it “hiprfish”. Download archive from github and extract contents of hiprfish-master to newly created hiprfish folder.

4. Open terminal and type:
     `cd hiprfish`\
     `chmod +x usearch11.0.667_i86linux32`

5. Install blast via terminal:
     `sudo apt install ncbi-blast+`

6. [Download](https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz) a copy of the NCBI 16S rRNA database. Move database contents (taxdb, etc.) inside the  16S_ribosomal_RNA folder, which is found in hiprfish>RNA. Make sure to delete text.txt file.

7. Install primer3 via terminal:
     `sudo apt-get install -y build-essential g++ cmake git-all`\
     `git clone https://github.com/primer3-org/primer3.git primer3`\
     `cd primer3/src`\
     `make`\
     `make test`

8. Edit the `hiprfish_config.json file` to point the pipeline to the correct directories.
   a.	Go to hiprfish/hiprfish_code/probe-design/hiprfish-probe-design-consensus/hiprfish_config.json and open hiprfish_config via text editor
   b.	If you have been following this guide, then your config parameters should look like this:
{
    "__default__" :
    {
        "SCRIPTS_PATH" : "/home/TYPE YOUR USERNAME/Desktop/hiprfish/hiprfish_code/probe-design/hiprfish-probe-design-consensus",
        "DATA_DIR" : "/home/ TYPE YOUR USERNAME /Desktop/hiprfish/data"
    },
    "blast":
    {
        "16s_db" : "/home/ TYPE YOUR USERNAME /Desktop/hiprfish/RNA/16S_ribosomal_RNA/16S_ribosomal_RNA"
    },
    "primer3":
    {
        "primer3_exec_dir" : /home/TYPE YOUR USERNAME/primer3/src/primer3_core",
        "primer3_config_dir" : "/home/ TYPE YOUR USERNAME /Desktop/hiprfish/primer3_config/primer3_config"
    },
    "usearch":
    {
        "usearch_dir" : "/home/ TYPE YOUR USERNAME /Desktop/hiprfish/usearch11.0.667_i86linux32"
    },
    "simulations" :
    {
        "simulation_table" : "/home/ TYPE YOUR USERNAME /Desktop/hiprfish/simulation_table/simulation_table_example.csv"
    }


### Input
1. Simulation summary file (simulation_table_test.csv)
   - A csv file containing all the designs to be run.
      * `DESIGN_ID`: identifier for each design
      * `SAMPLE`: name of the input FASTA file without file extension
      * `TARGET_RANK`: desired taxonomic rank for the probe design. Availabel options: phylum, class, order, family, genus, and species
      * `SIMILARITY`: similarity cut off for grouping 16S sequences. A low cut off (e.g. 0.1) essentially means 16S sequences will be grouped by their lineage information in the NCBI 16SMicrobial database. A higher cut off can be used to subdivide sequences within a given taxon. Higher cut off values generally leads to longer run time.
      * `MAX_CONTINUOUS_HOMOLOGY`: maximum continuous homology (measured in bp) for a probe-target hit to be considered significant. Lower values leads to more stringent designs. Default is 14 bp.
      * `MIN_TM`: minimum melting temperature threhold
      * `MAX_TM`: maximum melting temperature threhold
      * `GC`: minimum probe GC content threhold
      * `INCLUDE_START`: number of nucleotides to exclude at the beginning of the 16S sequences
      * `INCLUDE_END`: number of nucleotides to exclude at the end of the 16S sequences
      * `PROBE_SELECTION_METHOD`: method for selecting probes. Available options are
         1. `SingleBestProbe`: select the top probe for each taxa, if available
         2. `AllSpecific`: select all probes that are specific and only specific to its target taxon
         3. `AllSpecificPStartGroup`: select all probes that are specific and only specific to its target taxon within each segment of the 16S sequences. By default the 16S sequences are dividied into block resolutions of 100bp regions. If there are less than 15 probes available (average one probe per block), the block resolution is modified in 20bp decrements until there are 15 probes or the block resolution is zero, whichever happens first.
         4. `MinOverlap`: select all probes that are specific and only specific to its target taxon with minimum overlap in their target coverage
         5. `TopN`: select the top *n* probes for each taxa
       * `PRIMERSET`: primer sets to include in the final probes. There are three sets (A, B, and C) availble in the current version. User specific primer sets can also be added if necessary.
       * `OTU`: boolean to indicate whether to group 16S sequences only by their similarity. Generally set to `F` for ease of taxonomic interpretation of the probe designs, but could be useful if very high taxonomic resolution is desired.
       * `TPN`: number of top probes to select for each taxon, if the probe selection method is set to `TopN`
       * `FREQLL`: minimum abundance threshold. Default is zero, and is generally left at zero. Can be increased in situations where the in silico taxonomic coverage is not as good as desired. A higher value means increasing the probe design space for the more abundance sequences at the risk of those probes mishybridizing to the lower abundance taxa in the experiment.
       * `BOT`: minimum blast on target rate threshold. Probes with blast on target values lower than this value is considered *promiscuous*, and is not included in the final probe pool.
       * `BITSCORE_THRESH`: blast bitscore cutoff. Any blast hits (between probe and target sequence) with a score higher than this number will be considered significant and used for evaluation of probe specificity.
       * `BARCODESELECTION`: method for barcode assignment to taxa. Available options are:
         1. MostSimple: assign barcodes by barcode complexity, starting with the simplest ones. Barcodes with more bits are considered more complex.
         2. Random: randomly assign barcodes to taxa
         3. MostComplex: assign barcodes by barcode complexity, starting with the most complex ones. Barcodes with more bits are considered more complex.
       * `BPLC`: minimum blocking probe length threhold. Blocking probes with length lower than this threshold is considered likely to be washed off and do not need to be included in the final probe pool. Default is 15 bp.
       * `HELPER_PROBE_REPEAT`: number of times to repeat helper sequences in the final complex oligo pool. Default is 14 so that the any inidividual helper sequence is roughly at the same concentration as any individual encoding probe.
       * `SOD`: assumed sodium concentration for caluclation of melting temperatures. Default is 390.
       * `DNACONC`: assumed probe concentration for calculation of melting temperatures. Default is 5.
       * `MT_CUTOFF`: off target melting temperature cutoff. Probes would need to have a maximum off target melting temp smaller than this number to be considered specific. Default is 60. Note that this seems high because there seem to be a constant offset between melting temperatures calculated by primer3 and biopython built-in melting temperature calculation. This parameter refers to the calculation from the biopython implementation, which generally are higher than the primer3 calculations.
       * `OT_GC_CUTOFF`: off target maximum GC count. A probe would only be considered specific if any of its off-target binding sites have less than this many bases of G or C. Default is 7.
       * `THEME_COLOR`: overall theme color for axes and labels of the generated plots. Available options are:
          1. black: plots will have black axes and labels - works well against light background slides
          2. white: plots will have white axes and labels - works well against dark background slides
2. FASTA file
   - A FASTA file containing full length 16S sequences of the community to be probed. This file can be curated from public databases, or it can come from your own long read sequencing datasets, such as those from PacBio. The input file should be placed in `DATA_DIR/[SAMPLE]/input/[SAMPLE].fasta`.


### Output

1. Simulation results file
   - A csv file containing all the parameters for all the designs, as well as some summary statistics for each design
2. Probe folder
   - A folder containing selected probe summary files for each taxa, a concatenated file containing all selected probes, a file containing information for all the blocking probes, as well as text files that can be sent as is to array synthesis vendors for complex oligo pool synthesis. View your results inside hiprfish>data>simulation>dsgn0000 folder and find useful statistics in simulation_table folder.

### Before you start
•	Go to hiprfish>data and delete simulation folder
•	Go to hiprfish>data>mock_  community and delete species and utilities folders
•	Go to hiprfish>data>mock_  community>input and delete all files except mock_community.fasta 
•	Input desired sequences inside mock_community.fasta file
•	Set desired variables inside simulation_table_example.csv

### Running the pipeline
1.	Go to hiprbish-probe-design-consesus folder (/home/YOUR USERNAME/Desktop/hiprfish/hiprfish_code/probe-design/hiprfish-probe-design-consensus/) and right click on the background. Select open in terminal.
2.	Type conda activate hiprfish
3.	Type snakemake --configfile hiprfish_config.json -j n (where n is the number of cores you are willing to utilize, I personally use 8)

