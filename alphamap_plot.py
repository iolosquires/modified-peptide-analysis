import tomli as tomllib
import argparse
from pathlib import Path
from alphamap.organisms_data import import_fasta
from alphamap.importing import import_data
from alphamap.preprocessing import format_input_data
from alphamap.sequenceplot import plot_peptide_traces, uniprot_color_dict
from alphamap.uniprot_integration import uniprot_feature_dict
from alphamap.organisms_data import import_uniprot_annotation
import glob
import re

parser = argparse.ArgumentParser(
    prog="mascot-ppa", description="Mascot phosphopeptide analysis (PPA) tool"
)
parser.add_argument("input_directory")
args = parser.parse_args()
input_directory = Path(args.input_directory)

with open(str(input_directory / "config.toml"), "rb") as f:
    config = tomllib.load(f)

# list of output directories
filename_dict = dict(zip(config["mascot_filename"], config["sample_name"]))
analysis_name = config["analysis_name"]

species = config["species"]
uniprot_id = config["uniprot_for_plot"]  # Nek1 isoform 3
mod_search = config["mod_search"]

output_dir = f"{input_directory}\\output\\"  # \\alphamap_{uniprot_id}.html'
output_savename = f"{input_directory}\\output\\alphamap_{uniprot_id}.html"
for_alphamap_path = glob.glob(output_dir + "*_for_alphamap.tsv")
for_alphamap_mascot = [
    re.search(r"F\d{6}", item).group() + ".mzid" for item in for_alphamap_path
]
for_alphamap_samplename = [filename_dict[item] for item in for_alphamap_mascot]
human_uniprot = import_uniprot_annotation(species)
human_fasta = import_fasta(species)

regex_mod_dict = {"Phospho": r"\[Phospho\s*\(STY\)\]", "GlyGly": r"\[GlyGly\s*\(K\)\]"}

df_to_plot = []

for i, file in enumerate(for_alphamap_path):
    df = import_data(file)
    formatted_df = format_input_data(
        df=df, fasta=human_fasta, modification_exp=regex_mod_dict[mod_search]
    )
    df_to_plot.append(formatted_df)

fig = plot_peptide_traces(
    df_to_plot,
    name=for_alphamap_samplename,
    protein=uniprot_id,
    fasta=human_fasta,
    uniprot=human_uniprot,
    selected_features=["DOMAIN", "REGION", "MOTIF", "MOD_RES", "BINDING", "ACT_SITE"],
    uniprot_feature_dict=uniprot_feature_dict,
    uniprot_color_dict=uniprot_color_dict,
)

fig.write_html(output_savename)
