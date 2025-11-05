# Modified Peptide Analysis Instructions

## Downloads

[Anaconda](https://www.anaconda.com/download)
[VSCode](https://code.visualstudio.com/download)

## Installation

1) Open Anaconda Powershell Prompt

2) Navigate to the directory where the scripts are located

3) Run following command: conda env create -f environment.yml

## How to run

1) Export mascot output as mzIdentML file. Select all Optional Protein Hit Information settings (Description, Length in residues and protein coverage, Taxonomy,
Taxonomy ID, Protein sequence.) and Optional Peptide Match Information (Show duplicate peptides, Export site analysis data)

2) Create folder with name of project (e.g. 240924oraimi). All files go in this folder (.mzid, .txt and .toml)

3) Open proteome discoverer for each search results. Click File, Export, To Text.

4) Click PSMs box only and save.

5) Move PSMs.txt file to the same folder as mzIdentML files.

6) Fill in config.toml file (example are here: Z:\proteinchem\IoloSquires\mascot-phosphopeptide\00-current\mascot_ppa_input)

7) Activate the conda environment and run main.py. Run main.py by writing in the command line "python main.py [directory\to\input\folder]. The results will appear in a folder called "output".
