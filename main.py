
from itertools import compress
import pyteomics.mzid as mzid
import tomllib
import subprocess
import pandas as pd
import argparse
from pathlib import Path
import logging
from fpdf import FPDF
from tqdm import tqdm

import functions.ppa_functions as iolo
import functions.ppa_plots as ioloplot
import functions.ppa_checks as iolotest
import functions.ppa_dataclasses as iolodata

from functions.alphamap_functions import mascot_script_to_alphamap

#parse arguments
parser = argparse.ArgumentParser(prog='mascot-ppa',
                                description='Mascot phosphopeptide analysis (PPA) tool')
parser.add_argument("input_directory")
args = parser.parse_args()
input_directory = Path(args.input_directory)

with open(str(input_directory / "config.toml"), "rb") as f:
    config_toml = tomllib.load(f)

#create dataclasses from config
config = iolodata.create_config_from_toml(config_toml)
assert len(config.mascot_filenames) == len(config.sample_names), "Number of mascot files and sample names do not match"

paths = iolodata.create_paths_from_toml(config_toml, 
                                        input_directory,
                                        r"z:/proteinchem/IoloSquires/mascot-phosphopeptide/00-current/mascot_ppa_script/alphamap_plot.py",
                                        r"C:/Users/ISquires001/AppData/Local/anaconda3/envs/alphamap2/python.exe"
                                        )

paths.output_directory.mkdir(parents=False, exist_ok=True)
paths.plot_path.mkdir(parents=False, exist_ok=True)

POI_record, mrc_db = iolo.check_uniprot_id(config,
                                      paths)


logging.basicConfig(filename=paths.log_file, 
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

### starting checks complete, start the analysis
logging.info("Starting the analysis")

file_contains_phospho = []
pd_file_contains_phospho = []
merged_data = []
merge_data_pdf = []

###loop over masccot files (mzid files) and pd output (txt files)
for mascot_filename,pd_filename, sample_name in tqdm(zip(config.mascot_filenames,
                                                         config.pd_filenames,
                                                         config.sample_names),
                                                     total=len(config.mascot_filenames)):
    
    search_func_dict = {True: iolo.contains_search_term_uniprot, 
                        False: iolo.contains_search_term}

    mascot_savename = str(mascot_filename).replace('.mzid','')
    ## create mascot file dataframe, find phosphorylated peptides, find position in protein, create coverage plot, create phosphosite plot,
    logging.info("Starting analysis for %s", mascot_filename)

    mascot_file = iolo.mascot_file_check(mascot_filename,input_directory)

    logging.info("Opening mzid file: %s", str(mascot_file))

    mascot_df= mzid.DataFrame(str(mascot_file))
    logging.info("Mascot file: %s loaded with %s rows", mascot_filename,len(mascot_df))
    
    df_wanted = mascot_df[(mascot_df['Mascot:score'] >= float(config['score_cutoff'])) & 
          (mascot_df['accession'].apply(lambda x: search_func_dict[config.uniprot_id_check](x, config.search_protein)))].copy()
    
    #print(df_wanted.head())

    if len(df_wanted) == 0:
        logging.info("No peptides found for %s", mascot_filename)
        file_contains_phospho.append(False)
        pd_file_contains_phospho.append(False)
        merged_data.append('No peptides found')
        merge_data_pdf.append('No peptides found')
        continue

    df_wanted['start_end']= df_wanted.apply(lambda x: list(zip(x.start,x.end)), axis = 1)

    if config.uniprot_id_check:
        df_wanted['accession'] = df_wanted['accession'].apply(iolo.extract_between_pipes)

    df_wanted['accid_start_end_dict'] = df_wanted.apply(lambda x: dict(zip(x.accession, x.start_end)), axis = 1)
    df_wanted['poi_start_stop']= df_wanted.apply(lambda x: x.accid_start_end_dict[config.search_protein], axis = 1)
    df_wanted['has_phospho'] = df_wanted['Modification'].apply(iolo.find_phospho_mod)
    
    df_wanted_phospho = df_wanted[df_wanted['has_phospho'] == True].copy()

    logging.info("%s phospho peptides", len(df_wanted_phospho))
    
    if len(df_wanted_phospho) == 0:
        logging.info("No phosphopeptides found for %s", mascot_filename)
        file_contains_phospho.append(False)
        pd_file_contains_phospho.append(False)
        merged_data.append('No phosphopeptides found')
        merge_data_pdf.append('No phosphopeptides found')
        continue
    
    file_contains_phospho.append(True)
    df_wanted_phospho['phos_positions'] = df_wanted_phospho['Modification'].apply(iolo.get_phospho_positions)
    phospho_cols_for_grouping = ['PeptideSequence','chargeState','phos_positions',
                                 'Mascot:score','Mascot:PTM site assignment confidence',
                                 'poi_start_stop','experimentalMassToCharge']
    df_phospho_grouped = df_wanted_phospho[phospho_cols_for_grouping].copy()
    
    
    df_phospho_grouped['single_acceptor']=df_phospho_grouped['PeptideSequence'].apply(iolo.find_phospho_single_acceptor_site)
    df_phospho_grouped['double_acceptor']=df_phospho_grouped['PeptideSequence'].apply(iolo.find_phospho_double_acceptor_site)
    df_phospho_grouped['missing_conf'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].isna()

    condition_single = (df_phospho_grouped['single_acceptor'] == True) & (df_phospho_grouped['missing_conf'] == True)
    condition_double = (df_phospho_grouped['double_acceptor'] == True) & (df_phospho_grouped['missing_conf'] == True)
    condition_missing = (df_phospho_grouped['single_acceptor'] == False) & (df_phospho_grouped['missing_conf'] == True) & (df_phospho_grouped['double_acceptor'] == False)
    

    df_phospho_grouped['Mascot:PTM site assignment confidence'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].fillna(condition_single.map({True: 100}))
    df_phospho_grouped['Mascot:PTM site assignment confidence'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].fillna(condition_double.map({True: 100}))
    df_phospho_grouped['Mascot:PTM site assignment confidence'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].fillna(condition_missing.map({True: 0}))
    
    group_has_conf_site = df_phospho_grouped.groupby(['PeptideSequence','chargeState','phos_positions'])['Mascot:PTM site assignment confidence'].transform(lambda x: x.max() >= 80)
    df_phospho_grouped['top_mascot_and_over80'] = (df_phospho_grouped.groupby(['PeptideSequence','chargeState','phos_positions'])['Mascot:score'].transform(lambda x: x == x.max())) & group_has_conf_site 
    df_phospho_grouped=df_phospho_grouped.round({'Mascot:score': 0})
    
    phos_pep_lengths = df_phospho_grouped['PeptideSequence'].str.len().to_list()
    phos_start_pos = [item[0] for item in df_phospho_grouped['poi_start_stop']]
    df_phospho_grouped['phos_positions_list'] = df_phospho_grouped['phos_positions'].apply(iolo.phospho_position_string_to_list)
    df_phospho_grouped['phos_in_protein'] = df_phospho_grouped.apply(lambda x: iolo.phos_site_in_protein(x['poi_start_stop'],x['phos_positions_list']),axis=1)
    
    if pd_filename == "none":
        size_pd_df = 0
    else:
        pd_savename = str(pd_filename).replace('.txt','')
        logging.info("Starting analysis for %s", pd_filename)
        pd_file = iolo.mascot_file_check(pd_filename,input_directory)
        logging.info("Opening PD txt file: %s", str(pd_file))
        pd_output = pd.read_csv(str(pd_file),dtype={'Protein Accessions': 'str'},sep='\t')
        logging.info("PD file: %s loaded with %s rows", pd_file,len(pd_output))
        if "ptmRS: Best Site Probabilities" in pd_output.columns:
            pd_output['PhosphoRS: Best Site Probabilities'] = pd_output['ptmRS: Best Site Probabilities']
        pd_output = pd_output[pd_output['Protein Accessions'].apply(lambda x: config.search_protein in x)].copy() 
        pd_output = pd_output.dropna(subset=['PhosphoRS: Best Site Probabilities']).copy()
        pd_output = pd_output[pd_output['Ions Score'] >= config['score_cutoff']].copy()
        pd_output['phospho_rs'] = pd_output['PhosphoRS: Best Site Probabilities'].apply(iolo.mods_to_list_pd_output)
        pd_output['phospho_rs'] = pd_output['phospho_rs'].apply(iolo.remove_non_phospho_mods)
        pd_output = pd_output[pd_output['phospho_rs'] != "No Phospho"].copy()
        logging.info("PD file: %s has %s rows after filtering", pd_file,len(pd_output))
        size_pd_df = len(pd_output)
    if size_pd_df > 0:
        pd_file_contains_phospho.append(True)
        logging.info("%s phosphopeptides found for %s", len(pd_output), pd_filename)
        
        pd_output['peptide_trim'] = pd_output['Annotated Sequence'].apply(iolo.get_peptide_from_pd_output)
        pd_output['phosphors_conf'] = pd_output['phospho_rs'].apply(iolo.get_conf_pd)
        pd_output['phos_pos_pep'] = pd_output['phospho_rs'].apply(iolo.get_pos_pd)
        pd_output['peptide_cap'] = pd_output['peptide_trim'].apply(iolo.capitalise_peptides)
        pd_output['poi_start_stop'] = pd_output['peptide_cap'].map(dict(zip(df_phospho_grouped['PeptideSequence'],
                                                                            df_phospho_grouped['poi_start_stop'])))
        
        pd_output['mean_conf'] = pd_output['phosphors_conf'].apply(iolo.mean_conf_pd)
        
        pd_grouped = pd_output.sort_values('mean_conf',ascending=False).groupby(['peptide_trim','Charge','phos_pos_pep'],as_index=False).first()
        pd_grouped['phos_positions_list'] = pd_grouped['phos_pos_pep'].apply(iolo.phos_pos_to_ints)
        
        # If nans in poi_start_stop, align peptide to protein using own function.
        
        if pd_grouped['poi_start_stop'].isnull().sum() > 0:
            #split intwo two dataframes, one with nans and one without
            pd_grouped_nan = pd_grouped[pd_grouped['poi_start_stop'].isnull()].copy()
            pd_grouped_not_nan = pd_grouped.dropna(subset=['poi_start_stop']).copy()
            pd_grouped_nan['poi_start_stop'] = pd_grouped_nan.apply(lambda x: iolo.align_peptide_to_protein(x['peptide_cap'],POI_record),axis=1)
            pd_grouped = pd.concat([pd_grouped_nan,pd_grouped_not_nan],axis=0)

        pd_grouped['phos_in_protein'] = pd_grouped.apply(lambda x: [item + x['poi_start_stop'][0] -2 for item in x['phos_positions_list']],axis=1)
        pd_grouped['lc_pep'] = pd_grouped.apply(lambda x: iolo.lowercase_modified_residue(x['peptide_cap'],x['phos_positions_list']),axis=1)
        pd_grouped['lc_pep_cs'] = pd_grouped['Charge'].astype(str) + '-' + pd_grouped['lc_pep']  
        pd_peptide_site_dict = dict(zip(pd_grouped['lc_pep_cs'],pd_grouped['phosphors_conf']))
        pd_peptide_is_dict = dict(zip(pd_grouped['lc_pep_cs'],pd_grouped['Ions Score']))
        pd_grouped_explode = pd_grouped.explode(['phosphors_conf','phos_positions_list','phos_in_protein'])
        pd_output_ge_grouped = pd_grouped_explode.sort_values('phosphors_conf',ascending=False).groupby(['phos_in_protein'],as_index=False).first()
        pd_conf_dict = dict(zip(pd_output_ge_grouped['phos_in_protein'],pd_output_ge_grouped['phosphors_conf']))
        df_phospho_grouped['pd_peptide'] = df_phospho_grouped.apply(lambda x: iolo.lowercase_modified_residue(x['PeptideSequence'],x['phos_positions_list']),axis=1)
        df_phospho_grouped['pd_peptide_cs'] = df_phospho_grouped['chargeState'].astype(str) + '-'  + df_phospho_grouped['pd_peptide']
        df_phospho_grouped['pd_conf'] = df_phospho_grouped['pd_peptide_cs'].map(pd_peptide_site_dict)                                                                                                                
        df_phospho_grouped['pd_is'] = df_phospho_grouped['pd_peptide_cs'].map(pd_peptide_is_dict)

        # Merge data from PD and Mascot and save to list.

        mascot_for_merge = df_phospho_grouped[['Mascot:score','pd_is','chargeState','experimentalMassToCharge',
                                               'pd_peptide', 'Mascot:PTM site assignment confidence',
                                               'pd_conf','phos_in_protein','poi_start_stop']].copy()
        
        
        pd_for_merge = pd_grouped[['Ions Score','Charge','m/z [Da]',
                                   'lc_pep','phosphors_conf','phos_in_protein',
                                   'poi_start_stop']].copy()
        
        
        mascot_for_merge.columns = ['Mascot Score','PD Score','Charge State',
                                    'M/Z','Peptide','PTM Confidence Mascot',
                                    'PTM Confidence PD','Position in Protein','Start Stop']
        pd_for_merge.columns= ['PD Score','Charge State','M/Z',
                               'Peptide','PTM Confidence PD','Position in Protein',
                               'Start Stop']
        merged = pd.concat([mascot_for_merge,pd_for_merge],axis=0)

        
        merged = merged.round({'M/Z': 1})
        merged['Mascot Group Conf'] = merged.groupby(['Peptide','Charge State','M/Z'])['PTM Confidence Mascot'].transform(lambda x : [x.tolist()]*len(x))
        merged = merged.sort_values(by=['Mascot Score'], ascending=False).drop_duplicates(subset=['Peptide','Charge State','M/Z'],keep='first')
        group_sizes = mascot_for_merge.groupby(['Peptide', 'Charge State']).size()
        
        phos_peptides = list(merged['Peptide'].apply(iolo.capitalise_peptides))
        phos_start_pos = merged['Start Stop'].to_list()
        pep_pos_strings = [str(item) for item in phos_start_pos]
        pep_dict = dict(zip(pep_pos_strings,phos_peptides))
        df_phospho_conf_dict = df_phospho_grouped.sort_values('Mascot:PTM site assignment confidence',ascending=False)
        df_phospho_conf_dict2 = df_phospho_conf_dict.explode('phos_in_protein')
        df_phospho_conf_dict2 = df_phospho_conf_dict2.drop_duplicates(subset='phos_in_protein')
        conf_dict = dict(zip(df_phospho_conf_dict2['phos_in_protein'],df_phospho_conf_dict2['Mascot:PTM site assignment confidence']))
        
        merge_conf_dict = iolo.combine_pd_mascot_confs(pd_conf_dict,conf_dict)
        
        phos_confs = [merge_conf_dict.get(i,0) for i in range(1,len(POI_record.seq)+1)]
        list_colors = ['red']*len(merge_conf_dict)
        mod_color_dict = dict(zip(merge_conf_dict.keys(),list_colors))
        phos_start_pos_unique = list(set(phos_start_pos))
        
        ioloplot.master_phosphopeptide_plot(phos_start_pos_unique,
                                            pep_dict,
                                            POI_record,
                                            mod_color_dict,
                                            phos_confs,
                                            100,
                                            output_path=paths.plot_path,
                                            sample_name=sample_name)

        # Filter merged data on localisation confidence and add to list for pdf creation.

        m = merged.apply(lambda x: iolo.filter_both_confidences(x['PTM Confidence Mascot'],x['PTM Confidence PD'],0),axis=1)
        merged_filter = merged[m].copy()
        merged_filter['Start Stop'] = merged_filter['Start Stop'].apply(lambda x: (x[0]-1,x[1]))
        
        #check start stop positions are correct
        iolotest.check_position_correct_aa(merged_filter['Position in Protein'],merged_filter['Peptide'],POI_record)
        iolotest.check_peptides_correct_position(merged_filter['Peptide'],merged_filter['Start Stop'],POI_record)
        
        #Store things from loop
        
        merged_data.append(merged_filter)
        df_pdf = merged_filter.map(str)
        columns =[['Mascot Score',
                   'PD Score','Charge State',
                   'M/Z','Peptide','Conf Mascot',
                   'Conf PD','Position in Protein',
                   'Start Stop','Mascot Group Conf']]
        rows = df_pdf.values.tolist()
        data = columns + rows
        merge_data_pdf.append(data)

        mascot_script_to_alphamap(merged,
                                 config['uniprot_for_plot'],
                                 paths.output_path / (mascot_savename + '_for_alphamap.tsv'))

    elif size_pd_df == 0:
        pd_file_contains_phospho.append(False)
        df_phospho_grouped['pd_peptide'] = df_phospho_grouped.apply(lambda x: iolo.lowercase_modified_residue(x['PeptideSequence'],x['phos_positions_list']),axis=1)
        mascot_for_merge = df_phospho_grouped[['Mascot:score',
                                               'chargeState',
                                               'experimentalMassToCharge',
                                               'pd_peptide',
                                               'Mascot:PTM site assignment confidence',
                                               'phos_in_protein',
                                               'poi_start_stop']].copy()
        
        mascot_for_merge.columns = ['Mascot Score','Charge State',
                                    'M/Z','Peptide','PTM Confidence Mascot',
                                    'Position in Protein','Start Stop']
        #Create phosphopeptide plot.
        merged = mascot_for_merge
        phos_peptides = list(merged['Peptide'].apply(iolo.capitalise_peptides))
        phos_start_pos = merged['Start Stop'].to_list()
        pep_pos_strings = [str(item) for item in phos_start_pos]
        pep_dict = dict(zip(pep_pos_strings,phos_peptides))
        df_phospho_conf_dict = df_phospho_grouped.sort_values('Mascot:PTM site assignment confidence',ascending=False)
        df_phospho_conf_dict2 = df_phospho_conf_dict.explode('phos_in_protein')
        df_phospho_conf_dict2 = df_phospho_conf_dict2.drop_duplicates(subset='phos_in_protein')
        conf_dict = dict(zip(df_phospho_conf_dict2['phos_in_protein'],df_phospho_conf_dict2['Mascot:PTM site assignment confidence']))
        phos_confs = [conf_dict.get(i,0) for i in range(1,len(POI_record.seq)+1)]
        list_colors = ['red']*len(conf_dict)
        mod_color_dict = dict(zip(conf_dict.keys(),list_colors))
        phos_start_pos_unique = list(set(phos_start_pos))
        
        ioloplot.master_phosphopeptide_plot(phos_start_pos_unique,
                                            pep_dict,
                                            POI_record,
                                            mod_color_dict,
                                            phos_confs,
                                            100,
                                            output_path=paths.plot_path,
                                            sample_name=sample_name)

        mascot_for_merge['Start Stop'] = mascot_for_merge['Start Stop'].apply(lambda x: (x[0]-1,x[1]))
        mascot_for_merge = mascot_for_merge.round({'M/Z': 1})
        mascot_for_merge['Mascot Group Conf'] = mascot_for_merge.groupby(['Peptide','Charge State','M/Z'])['PTM Confidence Mascot'].transform(lambda x : [x.tolist()]*len(x))
        mascot_for_merge = mascot_for_merge.sort_values(by=['Mascot Score'], ascending=False).drop_duplicates(subset=['Peptide','Charge State','M/Z'],keep='first')
        
        iolotest.check_position_correct_aa(mascot_for_merge['Position in Protein'],mascot_for_merge['Peptide'],POI_record)
        iolotest.check_peptides_correct_position(mascot_for_merge['Peptide'],mascot_for_merge['Start Stop'],POI_record)
        
        merged_data.append(mascot_for_merge)
        df_pdf = mascot_for_merge.map(str)
        columns =[['Mascot Score','Charge State','M/Z','Peptide',
                   'PTM Confidence Mascot','Position in Protein','Start Stop',
                   'Mascot Group Conf']]
        rows = df_pdf.values.tolist()
        data = columns + rows
        merge_data_pdf.append(data)
        mascot_script_to_alphamap(mascot_for_merge,
                                 config['uniprot_for_plot'],
                                 paths.output_path / (mascot_savename + '_for_alphamap.tsv'))
        

logging.info("Creating pdf report. Have %s dataframes to add to pdf", len(merged_data))
pdf = FPDF()
excel_dir = paths.output_path / (config['analysis_name'] + '_phosphopeptide_report.xlsx')

# Create pdf report and excel

data_for_excel_filter = list(compress(merged_data, file_contains_phospho))
sample_names_filter = list(compress(config.sample_names, file_contains_phospho))

with pd.ExcelWriter(excel_dir,engine="xlsxwriter")as writer:
    for data,mascot_file in zip(data_for_excel_filter,sample_names_filter):
        logging.info("Adding data to Excel for %s which has phosphopeptides", mascot_file)
        data.to_excel(writer, sheet_name=mascot_file,index=False)

col_width_dict = {True: (5,5,5,3,12,5,5,5,5,5), False: (5,5,3,12,5,5,5,5)}

for data,mascot_file,phospho_indicator,pd_phospho in zip(merge_data_pdf,config.sample_names,file_contains_phospho,pd_file_contains_phospho):
    
    
    


    if phospho_indicator:
        logging.info("Adding data to pdf for %s which has phosphopeptides", mascot_file)
        pdf.add_page()
        pdf.set_font("Times", style="B", size=18)
        pdf.cell(0, 10, mascot_file, align="C")
        pdf.ln(10) 
        pdf.set_font("Times", size=6)
        #add logo image
        pdf.image(str(paths.logo_file), x=20, y=5, w=40, h=15)
        with pdf.table(
            borders_layout="SINGLE_TOP_LINE",
            cell_fill_color=200,  # grey
            cell_fill_mode="ROWS",
            line_height=pdf.font_size * 2.5,
            text_align="CENTER",
            width=160,
            col_widths=col_width_dict[pd_phospho]
        ) as table:
            for data_row in data:
                row = table.row()
                for datum in data_row:
                    row.cell(datum)
        
    else:

        logging.info("Adding data to pdf for %s which does not have phosphopeptides", mascot_file)
        pdf.add_page()
        pdf.set_font("Times", style="B", size=18)
        pdf.cell(0, 10, mascot_file, align="C")
        pdf.ln(10) 
        pdf.set_font("Times", size=7)
        pdf.image(str(paths.logo_file), x=20, y=10, w=40, h=15)

pdf_dir = paths.output_path / (config['analysis_name'] + '_phosphopeptide_report.pdf')
pdf.output(pdf_dir)
logging.info("Analysis finished")

print('Creating Alphamap....')
#call alphamap script

# Build command
cmd = [paths.alphamap_python_path, paths.alphamap_python_script] + [args.input_directory]

# Run the script
result = subprocess.run(cmd,
                        text=True)

print('Analysis finished')