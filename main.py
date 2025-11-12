import tomllib
import pandas as pd
import argparse
from pathlib import Path

import functions.ppa_functions as iolo
import functions.ppa_plots as ioloplot
import functions.ppa_checks as iolotest
import functions.ppa_dataclasses as iolodata

from functions.alphamap_functions import mascot_script_to_alphamap

# parse arguments
parser = argparse.ArgumentParser(
    prog="mascot-ppa", description="Mascot phosphopeptide analysis (PPA) tool"
)
parser.add_argument("input_directory")
args = parser.parse_args()
input_directory = Path(args.input_directory)

with open(str(input_directory / "config.toml"), "rb") as f:
    config_toml = tomllib.load(f)

iolo.save_everything_to_output(input_directory, "output/code")
ed = iolo.open_read_experimental_design(input_directory)

# create dataclasses from config
config = iolodata.create_config_from_toml(config_toml,ed)

paths = iolodata.create_paths_from_toml(
    config_toml,
    input_directory,
    r"z:/proteinchem/IoloSquires/mascot-phosphopeptide/00-current/mascot_ppa_script/alphamap_plot.py",
    r"C:/Users/ISquires001/AppData/Local/anaconda3/envs/alphamap2/python.exe",
)

paths.output_directory.mkdir(parents=False, exist_ok=True)
if config.create_old_plot:
    paths.plot_path.mkdir(parents=False, exist_ok=True)

file_contains_phospho = []
pd_file_contains_phospho = []
merged_data = []
merge_data_pdf = []

###loop over masccot files (mzid files) and pd output (txt files)
for mascot_filename, pd_filename, sample_name, search_record in zip(config.mascot_filenames, 
                                                                    config.pd_filenames, 
                                                                    config.sample_names, 
                                                                    config.search_protein_list):
    
    mascot_savename = str(mascot_filename).replace(".mzid", "")
    POI_record, mrc_db = iolo.check_uniprot_id(config, paths,search_record)
    
    print(f"Processing sample: {sample_name}")

    mascot_df = iolo.import_mascot_file(mascot_filename, input_directory)
    df_wanted = iolo.filter_accession_and_score(mascot_df, config, search_record)

    if len(df_wanted) == 0:
        
        file_contains_phospho.append(False)
        pd_file_contains_phospho.append(False)
        merged_data.append("No peptides found")
        merge_data_pdf.append("No peptides found")
        continue

    df_wanted_phospho = iolo.find_phosphopeptides_in_mascot(df_wanted, config, search_record)

    

    if len(df_wanted_phospho) == 0:
        
        file_contains_phospho.append(False)
        pd_file_contains_phospho.append(False)
        merged_data.append("No phosphopeptides found")
        merge_data_pdf.append("No phosphopeptides found")
        continue

    file_contains_phospho.append(True)
    df_phospho_grouped, phos_pep_lengths, phos_start_pos = (
        iolo.process_mascot_phospho_dataframe(df_wanted_phospho)
    )

    if pd_filename == "none":
        size_pd_df = 0
    else:
        pd_output, size_pd_df, pd_savename = iolo.load_and_process_pd_file(
            pd_filename, input_directory, config, search_record
        )

    if size_pd_df > 0:
        pd_file_contains_phospho.append(True)
        pd_grouped, pd_peptide_site_dict, pd_peptide_is_dict, pd_conf_dict = (
            iolo.process_pd_phospho_dataframe(pd_output, df_phospho_grouped, POI_record)
        )
        df_phospho_grouped = iolo.add_pd_info_to_mascot(
            df_phospho_grouped, pd_peptide_site_dict, pd_peptide_is_dict
        )

        mascot_for_merge = df_phospho_grouped[
            [
                "Mascot:score",
                "pd_is",
                "chargeState",
                "experimentalMassToCharge",
                "pd_peptide",
                "Mascot:PTM site assignment confidence",
                "pd_conf",
                "phos_in_protein",
                "poi_start_stop",
            ]
        ].copy()

        pd_for_merge = pd_grouped[
            [
                "Ions Score",
                "Charge",
                "m/z [Da]",
                "lc_pep",
                "phosphors_conf",
                "phos_in_protein",
                "poi_start_stop",
            ]
        ].copy()

        mascot_for_merge.columns = [
            "Mascot Score",
            "PD Score",
            "Charge State",
            "M/Z",
            "Peptide",
            "PTM Confidence Mascot",
            "PTM Confidence PD",
            "Position in Protein",
            "Start Stop",
        ]
        pd_for_merge.columns = [
            "PD Score",
            "Charge State",
            "M/Z",
            "Peptide",
            "PTM Confidence PD",
            "Position in Protein",
            "Start Stop",
        ]

        # processing merged
        merged = pd.concat([mascot_for_merge, pd_for_merge], axis=0)
        merged = iolo.process_merged_dataframe(mascot_for_merge)

        conf_dict = iolo.create_conf_dict(df_phospho_grouped)
        merge_conf_dict = iolo.combine_pd_mascot_confs(pd_conf_dict, conf_dict)

        # Create phosphopeptide plot
        if config.create_old_plot:
            ioloplot.create_phospho_peptide_plot(
                merged, POI_record, merge_conf_dict, paths, sample_name
            )

        merged_filter = iolo.filter_merged_dataframe_by_localisation_confidence(
            merged, config
        )
        # check start stop positions are correct
        iolotest.check_position_correct_aa(
            merged_filter["Position in Protein"], merged_filter["Peptide"], POI_record
        )
        iolotest.check_peptides_correct_position(
            merged_filter["Peptide"], merged_filter["Start Stop"], POI_record
        )

        # Store things from loop
        merged_data.append(merged_filter)
        data = iolo.create_pdf_data(merged_filter)
        merge_data_pdf.append(data)

        mascot_script_to_alphamap(
            merged,
            config.uniprot_id,
            paths.output_directory / (mascot_savename + "_for_alphamap.tsv"),
        )

    # case where no pd file is given or no phosphopeptides in pd file
    elif size_pd_df == 0:
        pd_file_contains_phospho.append(False)

        df_phospho_grouped["pd_peptide"] = df_phospho_grouped.apply(
            lambda x: iolo.lowercase_modified_residue(
                x["PeptideSequence"], x["phos_positions_list"]
            ),
            axis=1,
        )

        conf_dict = iolo.create_conf_dict(df_phospho_grouped)

        mascot_for_merge = df_phospho_grouped[
            [
                "Mascot:score",
                "chargeState",
                "experimentalMassToCharge",
                "pd_peptide",
                "Mascot:PTM site assignment confidence",
                "phos_in_protein",
                "poi_start_stop",
            ]
        ].copy()

        mascot_for_merge.columns = [
            "Mascot Score",
            "Charge State",
            "M/Z",
            "Peptide",
            "PTM Confidence Mascot",
            "Position in Protein",
            "Start Stop",
        ]
        if config.create_old_plot:
            phos_peptides = list(
                mascot_for_merge["Peptide"].apply(iolo.capitalise_peptides)
            )
            phos_start_pos = mascot_for_merge["Start Stop"].to_list()
            pep_pos_strings = [str(item) for item in phos_start_pos]
            pep_dict = dict(zip(pep_pos_strings, phos_peptides))
            phos_confs = [conf_dict.get(i, 0) for i in range(1, len(POI_record.seq) + 1)]
            list_colors = ["red"] * len(conf_dict)
            mod_color_dict = dict(zip(conf_dict.keys(), list_colors))
            phos_start_pos_unique = list(set(phos_start_pos))

            ioloplot.master_phosphopeptide_plot(
                phos_start_pos_unique,
                pep_dict,
                POI_record,
                mod_color_dict,
                phos_confs,
                100,
                output_path=paths.plot_path,
                sample_name=sample_name,
            )

        mascot_for_merge["Start Stop"] = mascot_for_merge["Start Stop"].apply(
            lambda x: (x[0] - 1, x[1])
        )
        mascot_for_merge = iolo.process_merged_dataframe(mascot_for_merge)

        iolotest.check_position_correct_aa(
            mascot_for_merge["Position in Protein"],
            mascot_for_merge["Peptide"],
            POI_record,
        )
        iolotest.check_peptides_correct_position(
            mascot_for_merge["Peptide"], mascot_for_merge["Start Stop"], POI_record
        )

        merged_data.append(mascot_for_merge)
        data = iolo.create_pdf_data(mascot_for_merge)
        merge_data_pdf.append(data)
        mascot_script_to_alphamap(
            mascot_for_merge,
            config.uniprot_id,
            paths.output_directory / (mascot_savename + "_for_alphamap.tsv"),
        )

iolo.create_excel_report(config, merged_data, file_contains_phospho, paths)
iolo.create_pdf_report(
    config, merge_data_pdf, file_contains_phospho, pd_file_contains_phospho, paths
)

