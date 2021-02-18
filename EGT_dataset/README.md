# General information
This directory contains the results and scripts of the article "Gene Tree and Species Tree Reconciliation with Endosymbiotic Gene Transfer", Anselmetti et al., submitted.

# Results reproduction
To reproduce the results of this article -> read and execute the lines of the file "data_preprocessing.sh".
Some files have to be manually produced or modified:
- The 11 Plantae species tree file (MitoCOGs/dataset/11_Plantae_dataset/11_Plantae_species_tree.nwk) has been produced manually based on the topology in Kannan et al., 2014
- 2 manual modifications have to be done on the data preprocessing pipeline (see l.123 & l.127) due to a format to the end line of the Kannan et al., 2014 dataset files that does not allow automatic modifications

# Directory structure and content
```
.
├── bin
├── data_preprocessing.sh
├── MitoCOGs
│   ├── 11_Plantae_dataset
│   └── raw_data
├── MitoCOGs_dataset_statistics.ods
├── pipeline_for_endorex.py
└── README.md

```

- "bin/" directory: contains bin of the software/tool used in the pipeline to produce input data for EndoRex on the 11 Plantae dataset (pipeline implemented in script "pipeline_for_endorex.py")
- "data_preprocessing.sh": bash script to get input data for the EndoRex pipeline (pipeline_for_endorex.py) from the Kannan et al., 2014 dataset
- "MitoCOGs/" directory: contains 2 sub-directories
	- "11_Plantae_dataset/": contains all data and results files of the 11 Plantae dataset
	- "raw_data/" subdirectory: contains the raw data from Kannan et al., 2014 obtained at ftp://ftp.ncbi.nih.gov/pub/koonin/MitoCOGs
- "pipeline_for_endorex.py": pyhton script to produce input data for EndoRex software from data preprocessed by the "data_preprocessing.sh" script+manual modifications 
- "MitoCOGs_dataset_statistics.ods": Statistics on the Kannan et al., 2014 and the 11 Plantae dataset
- README.md: current file
