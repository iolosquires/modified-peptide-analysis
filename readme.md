1) Export mascot output as mzIdentML file. Select all Optional Protein Hit Information settings (Description, Length in residues and protein coverage, Taxonomy,
Taxonomy ID, Protein sequence.) and Optional Peptide Match Information (Show duplicate peptides, Export site analysis data)

2) Create folder in input folder with name of project (e.g. 240924oraimi). Put .mzid files into the folder you have created.

3) Export PSMs file from Proteome Discoverer and add to project folder.

4) Modify config.toml file in the project folder. Here's an example config file:
            analysis_name = "240924oraimi"
            mascot_filename = ["F297454.mzid","F297455.mzid"]
            pd_filename = ["ORaimi_20240924_01_PSMs.txt","ORaimi_20240924_02_PSMs.txt"]
            sample_name = ["240924oraimi_1","240924oraimi_2"]
            score_cutoff = 19.0
            search_protein = "2146"
            mrc_db = "PA_MRCdatabase-240910.fasta"

5) Create conda environment from environment.yml file. Activate the conda environment and run main.py. Run main.py by writing in the command line "python main.py [directory\to\input\folder]. The results will appear in a folder called "output".
