import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from functions.ppa_functions import capitalise_peptides 


def create_sequences_for_plot (phos_confs,POI_record,pep_strings,chunk_size):
    
    sequences_for_chop = [phos_confs] + [list(range(1,len(phos_confs)+1))] + [POI_record.seq] + pep_strings

    assert type(phos_confs) == list, 'phos_confs must be a list'
    assert type(pep_strings) == list, 'pep_strings must be a list'

    chopped_sequences = []
    for seq in sequences_for_chop:
        chopped = chop_strings(seq,chunk_size)
        chopped_sequences.append(chopped)
    return chopped_sequences

def find_chunks_with_phosphosites (seqs, mod_color_dict):
    skip_list = []
    for num_seqs in range(0,len(seqs[0])):
        aa_num = seqs[1][num_seqs]
        skip = True
        #check if there are any phosphosites in the chunk
        for aa_number in aa_num:
            if aa_number in mod_color_dict.keys():
                skip= False
                break
        skip_list.append(skip)
    chunk_for_plot = [item for item,skip in zip(range(0,len(seqs[0])),skip_list) if not skip]
    return chunk_for_plot
def phosphopep_coverage_plot (chunk_for_plot,chopped_sequences,mod_color_dict,POI_record,output_path,sample_name):

    fig,axs = plt.subplots(len(chunk_for_plot),1,figsize=(20,4*len(chunk_for_plot)))

    # Ensure axs is always iterable
    if len(chunk_for_plot) == 1:
        axs = [axs]  

    fig.suptitle(f'Modified peptide plot for {POI_record.name}', fontsize=16)
    
    font = {'family': 'monospace',
                'color':  'black',
                'weight': 'normal',
                'size': 8,
                }

    for ax_num,num_seqs in enumerate(chunk_for_plot):
        
        aa_num = chopped_sequences[1][num_seqs]
        phospho_confidences = np.array(chopped_sequences[0][num_seqs])
        protein = chopped_sequences[2][num_seqs]
        pep_strings = [item[num_seqs] for item in chopped_sequences[3:]]
        sequences = [protein] + pep_strings

        categories = ['Poor (0-80%)','High (80-100%)']
        cutoffs = [(0,80),(80,100)]
        colours = ['grey','red']
    
        for cat, cut, col in zip(categories,cutoffs,colours):
            pc_for_mask = phospho_confidences.copy()
            mask = ma.masked_outside(pc_for_mask, cut[0],cut[1] )
            pc_mask2 = pc_for_mask.copy()
            pc_mask2[mask.mask] = 0
            axs[ax_num].bar(aa_num, pc_mask2, align='edge', linewidth=0, color=col,label=cat)

        if ax_num == 0:
            axs[ax_num].legend()
            
        axs[ax_num].set_ylim([0, 100])
        ypos_xlbl = -25

        for j in range(0,len(sequences)):
            
            for i in range(0, len(protein)):
                    
                    if j == 0:
                        
                        axs[ax_num].text(float(aa_num[i]-1),ypos_xlbl, sequences[j][i],
                        bbox=dict(facecolor=mod_color_dict.get(aa_num[i]-1,'grey'),alpha=0.5),fontdict=font)

                    else:
                        
                        axs[ax_num].text(float(aa_num[i]-1),ypos_xlbl, sequences[j][i],fontdict=font)
                        
            ypos_xlbl -= 10    
                
        
    fig.tight_layout()
    plot_dir = output_path / (sample_name + '_cov_plot.png')
    plt.savefig(plot_dir, dpi=600)
    plt.close()

def create_pep_groups_plot (groups,pep_dict):
    pep_groups = []
    for nonoverlap_group in groups:
        group = []
        for start_stop in nonoverlap_group:
            group.append(pep_dict[str(start_stop)])
        
        pep_groups.append(group)
    return pep_groups

def fill_pep_groups (pep_groups,groups,POI_record):
   
    pep_strings = []
    for group_coords, group_seqs in zip(groups,pep_groups):
        if len(group_coords) == 0:
            continue
        pep_strings.append(create_peptide_string(POI_record.seq, group_seqs, group_coords))
    
    return pep_strings

def master_phosphopeptide_plot (phos_start_pos_unique,pep_dict,POI_record,mod_color_dict,phos_confs,chunk_size,output_path,sample_name):

    groups = all_non_overlap(phos_start_pos_unique)
    pep_groups = create_pep_groups_plot(groups,pep_dict)
    pep_strings = fill_pep_groups(pep_groups,groups,POI_record)
    chopped_sequences = create_sequences_for_plot(phos_confs,POI_record,pep_strings,chunk_size)
    chunk_for_plot = find_chunks_with_phosphosites(chopped_sequences,mod_color_dict)
    phosphopep_coverage_plot(chunk_for_plot,chopped_sequences,mod_color_dict,POI_record,output_path,sample_name)

def chop_strings (protein,chunk_size):
    string_parts = []
    continue_chopping = True
    i = 0
    while continue_chopping:
        shifter = i*chunk_size
        i += 1
        chopper = (0+shifter,chunk_size+shifter)
        sliced_protein = protein[chopper[0]:chopper[1]]
        if len(sliced_protein) == 0:
            continue_chopping = False
            break
        string_parts.append(sliced_protein)
    return string_parts

def all_non_overlap (test):
    
    all_groups = []
    groups, not_in_group = find_non_overlap(test)
    all_groups.append(groups)
    while len(not_in_group) > 0:
        groups, not_in_group = find_non_overlap(not_in_group)
        all_groups.append(groups)
    return all_groups

def find_non_overlap (test):
    test.sort()
    groups = []
    not_in_group = []
    group_end = test[0][1]
    for start, end in test[1:]:
        if start >= group_end:
            group_end = end
            groups.append((start, end))
        else:
            not_in_group.append((start, end))
    groups.append(test[0]) 
    groups.sort()

    return groups, not_in_group

def create_peptide_string (protein, peptides, pep_pos):
    pep_string = '_'*len(protein)
    first_loop = True
    for pep, pos in zip(peptides, pep_pos):
        if first_loop:
            pep2 = pep_string[:pos[0]-1] + pep + pep_string[pos[1]:]
            first_loop = False
        else:
            pep2 = pep2[:pos[0]-1] + pep + pep2[pos[1]:]

    #assert len(pep2) == len(pep_string), 'Lengths do not match between protein and peptides string'
    return pep2  


def create_phospho_peptide_plot(merged,POI_record,merge_conf_dict,paths,sample_name):

    phos_peptides = list(merged['Peptide'].apply(capitalise_peptides))
    phos_start_pos = merged['Start Stop'].to_list()
    pep_pos_strings = [str(item) for item in phos_start_pos]
    pep_dict = dict(zip(pep_pos_strings,phos_peptides))
    phos_confs = [merge_conf_dict.get(i,0) for i in range(1,len(POI_record.seq)+1)]
    list_colors = ['red']*len(merge_conf_dict)
    mod_color_dict = dict(zip(merge_conf_dict.keys(),list_colors))
    phos_start_pos_unique = list(set(phos_start_pos))
    
    master_phosphopeptide_plot(phos_start_pos_unique,
                                        pep_dict,
                                        POI_record,
                                        mod_color_dict,
                                        phos_confs,
                                        100,
                                        output_path=paths.plot_path,
                                        sample_name=sample_name)
