import functions.ppa_functions as iolo
import functions.ppa_plots as ioloplot
import functions.ppa_checks as iolotest
from itertools import compress
import pyteomics.mzid as mzid
import tomllib
import re
import pandas as pd
import argparse
from pathlib import Path
import logging
from fpdf import FPDF
from tqdm import tqdm
import sys
from functions.alphamap_functions import mascot_script_to_alphamap

#parse arguments
parser = argparse.ArgumentParser(prog='mascot-ppa',
                                description='Mascot phosphopeptide analysis (PPA) tool')
parser.add_argument("input_directory")
args = parser.parse_args()
input_directory = Path(args.input_directory)

with open(str(input_directory / "config.toml"), "rb") as f:
    config = tomllib.load(f)

###do config file, directory, file checks, create log file and output directory
assert len(config['mascot_filename']) == len(config['sample_name']), "Number of mascot files and sample names do not match"

search_protein = config['search_protein']
sample_names = config['sample_name']
uniprot_id_check = False

if re.search('[A-Z]', search_protein):
    uniprot_id_check = True

print("uniprot id search:", uniprot_id_check)
logo_file,output_path = iolo.file_check(config['analysis_name'])

log_filename = output_path / "run.log"
plot_path = output_path / "plots"
if not plot_path.exists():
    plot_path.mkdir(parents=False, exist_ok=True)

logging.basicConfig(filename=log_filename, level=logging.INFO,format='%(asctime)s - %(levelname)s - %(message)s')

### starting checks complete, start the analysis
logging.info("Starting the analysis")

#get the accession ids sequence from MRC database
if uniprot_id_check == False:
    mrc_db_path = "mrc_db\\fasta_formated_version\\" + config['mrc_db']
    mrc_db = iolo.mrc_db_to_dict(mrc_db_path)
    POI_record = mrc_db[search_protein]

elif uniprot_id_check == True:
    #make class called POI_record with attributes id and seq
    class POI_record:
        def __init__(self, seq, name):
            self.seq = seq
            self.name = name

    poi_record_instance = POI_record(seq = config['seq'],
                                     name = search_protein)


#mascot_data = []
search_func_dict = {True: iolo.contains_search_term_uniprot, False: iolo.contains_search_term}

file_contains_phospho = []
pd_file_contains_phospho = []
merged_data = []
merge_data_pdf = []

###loop over masccot files (mzid files) and pd output (txt files)
for mascot_filename, sample_name in tqdm(zip(config['mascot_filename'],config['sample_name']),total=len(config['mascot_filename'])):
    
    mascot_savename = str(mascot_filename).replace('.mzid','')
    ## create mascot file dataframe, find phosphorylated peptides, find position in protein, create coverage plot, create phosphosite plot,
    logging.info("Starting analysis for %s", mascot_filename)

    mascot_file = iolo.mascot_file_check(mascot_filename,input_directory)

    logging.info("Opening mzid file: %s", str(mascot_file))

    mascot_df= mzid.DataFrame(str(mascot_file))
    logging.info("Mascot file: %s loaded with %s rows", mascot_filename,len(mascot_df))
    
    
    df_wanted = mascot_df[(mascot_df['Mascot:score'] >= float(config['score_cutoff'])) & 
          (mascot_df['accession'].apply(lambda x: search_func_dict[uniprot_id_check](x, search_protein)))].copy()
    
    if len(df_wanted) == 0:
        logging.info("No peptides found for %s", mascot_filename)
        file_contains_phospho.append(False)
        pd_file_contains_phospho.append(False)
        merged_data.append('No peptides found')
        merge_data_pdf.append('No peptides found')
        continue

    #df_wanted.to_csv("debug_start_stop.csv")
    df_wanted['start_end']= df_wanted.apply(lambda x: list(zip(x.start,x.end)), axis = 1)

    if uniprot_id_check:
        df_wanted['accession'] = df_wanted['accession'].apply(iolo.extract_between_pipes)

    df_wanted['accid_start_end_dict'] = df_wanted.apply(lambda x: dict(zip(x.accession, x.start_end)), axis = 1)
    #df_wanted.to_csv("debug_key_error")
    df_wanted['poi_start_stop']= df_wanted.apply(lambda x: x.accid_start_end_dict[search_protein], axis = 1)
    df_wanted['has_ubi'] = df_wanted['Modification'].apply(iolo.find_ubiquitin_mod)
    
    df_wanted_phospho = df_wanted[df_wanted['has_ubi'] == True].copy()

    logging.info("%s ubi peptides", len(df_wanted_phospho))
    
    if len(df_wanted_phospho) == 0:
        logging.info("No ubipeptides found for %s", mascot_filename)
        file_contains_phospho.append(False)
        pd_file_contains_phospho.append(False)
        merged_data.append('No ubipeptides found')
        merge_data_pdf.append('No ubipeptides found')
        continue
    
    file_contains_phospho.append(True)
    df_wanted_phospho['ubi_positions'] = df_wanted_phospho['Modification'].apply(iolo.get_ubi_positions)
    
    phospho_cols_for_grouping = ['PeptideSequence','chargeState','ubi_positions',
                                 'Mascot:score','Mascot:PTM site assignment confidence',
                                 'poi_start_stop','experimentalMassToCharge']
    df_phospho_grouped = df_wanted_phospho[phospho_cols_for_grouping].copy()
    
    df_phospho_grouped['single_acceptor']=df_phospho_grouped['PeptideSequence'].apply(iolo.find_ubi_single_acceptor_site)
    df_phospho_grouped['double_acceptor']=df_phospho_grouped['PeptideSequence'].apply(iolo.find_ubi_double_acceptor_site)
    df_phospho_grouped['missing_conf'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].isna()

    condition_single = (df_phospho_grouped['single_acceptor'] == True) & (df_phospho_grouped['missing_conf'] == True)
    condition_double = (df_phospho_grouped['double_acceptor'] == True) & (df_phospho_grouped['missing_conf'] == True)
    condition_missing = (df_phospho_grouped['single_acceptor'] == False) & (df_phospho_grouped['missing_conf'] == True) & (df_phospho_grouped['double_acceptor'] == False)
    

    df_phospho_grouped['Mascot:PTM site assignment confidence'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].fillna(condition_single.map({True: 100}))
    df_phospho_grouped['Mascot:PTM site assignment confidence'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].fillna(condition_double.map({True: 100}))
    df_phospho_grouped['Mascot:PTM site assignment confidence'] = df_phospho_grouped['Mascot:PTM site assignment confidence'].fillna(condition_missing.map({True: 0}))
    
    group_has_conf_site = df_phospho_grouped.groupby(['PeptideSequence','chargeState','ubi_positions'])['Mascot:PTM site assignment confidence'].transform(lambda x: x.max() >= 80)
    df_phospho_grouped['top_mascot_and_over80'] = (df_phospho_grouped.groupby(['PeptideSequence','chargeState','ubi_positions'])['Mascot:score'].transform(lambda x: x == x.max())) & group_has_conf_site 
    df_phospho_grouped=df_phospho_grouped.round({'Mascot:score': 0})
    
    phos_pep_lengths = df_phospho_grouped['PeptideSequence'].str.len().to_list()
    phos_start_pos = [item[0] for item in df_phospho_grouped['poi_start_stop']]
    df_phospho_grouped.to_csv("debug_ubi_positions.csv")
    df_phospho_grouped['ubi_positions_list'] = df_phospho_grouped['ubi_positions'].apply(iolo.ubi_position_string_to_list)
    df_phospho_grouped['ubi_in_protein'] = df_phospho_grouped.apply(lambda x: iolo.phos_site_in_protein(x['poi_start_stop'],x['ubi_positions_list']),axis=1)
    
    pd_file_contains_phospho.append(False)
    df_phospho_grouped['pd_peptide'] = df_phospho_grouped.apply(lambda x: iolo.lowercase_modified_residue(x['PeptideSequence'],x['ubi_positions_list']),axis=1)
    mascot_for_merge = df_phospho_grouped[['Mascot:score',
                                            'chargeState',
                                            'experimentalMassToCharge',
                                            'pd_peptide',
                                            'Mascot:PTM site assignment confidence',
                                            'ubi_in_protein',
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
    df_phospho_conf_dict2 = df_phospho_conf_dict.explode('ubi_in_protein')
    df_phospho_conf_dict2 = df_phospho_conf_dict2.drop_duplicates(subset='ubi_in_protein')
    conf_dict = dict(zip(df_phospho_conf_dict2['ubi_in_protein'],df_phospho_conf_dict2['Mascot:PTM site assignment confidence']))
    
    phos_confs = [conf_dict.get(i,0) for i in range(1,len(poi_record_instance.seq)+1)]
    
    list_colors = ['red']*len(conf_dict)
    mod_color_dict = dict(zip(conf_dict.keys(),list_colors))
    phos_start_pos_unique = list(set(phos_start_pos))
    
    ioloplot.master_phosphopeptide_plot(phos_start_pos_unique,
                                        pep_dict,
                                        poi_record_instance,
                                        mod_color_dict,
                                        phos_confs,
                                        100,
                                        output_path=plot_path,
                                        sample_name=sample_name)

    
    mascot_for_merge['Start Stop'] = mascot_for_merge['Start Stop'].apply(lambda x: (x[0]-1,x[1]))
    mascot_for_merge = mascot_for_merge.round({'M/Z': 1})
    mascot_for_merge['Mascot Group Conf'] = mascot_for_merge.groupby(['Peptide','Charge State','M/Z'])['PTM Confidence Mascot'].transform(lambda x : [x.tolist()]*len(x))
    mascot_for_merge = mascot_for_merge.sort_values(by=['Mascot Score'], ascending=False).drop_duplicates(subset=['Peptide','Charge State','M/Z'],keep='first')
    
    
    #check start stop positions are correct
    #mascot_for_merge.to_csv("check_pos_prot")
    #print("length of protein",len(poi_record_instance.seq))
    iolotest.check_position_correct_aa(mascot_for_merge['Position in Protein'],mascot_for_merge['Peptide'],poi_record_instance)
    iolotest.check_peptides_correct_position(mascot_for_merge['Peptide'],mascot_for_merge['Start Stop'],poi_record_instance)
    
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
                              output_path / (mascot_savename + '_for_alphamap.tsv'))
    

logging.info("Creating pdf report. Have %s dataframes to add to pdf", len(merged_data))
pdf = FPDF()
excel_dir = output_path / (config['analysis_name'] + '_ubipeptide_report.xlsx')

# Create pdf report and excel

data_for_excel_filter = list(compress(merged_data, file_contains_phospho))
sample_names_filter = list(compress(sample_names, file_contains_phospho))

with pd.ExcelWriter(excel_dir,engine="xlsxwriter")as writer:
    for data,mascot_file in zip(data_for_excel_filter,sample_names_filter):
        logging.info("Adding data to Excel for %s which has ubipeptides", mascot_file)
        data.to_excel(writer, sheet_name=mascot_file,index=False)

col_width_dict = {True: (5,5,5,3,12,5,5,5,5,5), False: (5,5,3,12,5,5,5,5)}

for data,mascot_file,phospho_indicator,pd_phospho in zip(merge_data_pdf,sample_names,file_contains_phospho,pd_file_contains_phospho):
    
    if phospho_indicator:
        logging.info("Adding data to pdf for %s which has ubipeptides", mascot_file)
        pdf.add_page()
        pdf.set_font("Times", style="B", size=18)
        pdf.cell(0, 10, mascot_file, align="C")
        pdf.ln(10) 
        pdf.set_font("Times", size=6)
        #add logo image
        pdf.image(str(logo_file), x=20, y=5, w=40, h=15)
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

        logging.info("Adding data to pdf for %s which does not have ubipeptides", mascot_file)
        pdf.add_page()
        pdf.set_font("Times", style="B", size=18)
        pdf.cell(0, 10, mascot_file, align="C")
        pdf.ln(10) 
        pdf.set_font("Times", size=7)
        pdf.image(str(logo_file), x=20, y=10, w=40, h=15)

pdf_dir = output_path / (config['analysis_name'] + '_ubipeptide_report.pdf')
pdf.output(pdf_dir)
logging.info("Analysis finished")
print('Analysis finished')