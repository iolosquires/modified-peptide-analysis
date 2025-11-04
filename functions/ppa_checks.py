import functions.ppa_functions as iolo

def check_position_correct_aa (position_column,peptide_column,protein_seq):
    for df_pos,peptide in zip(position_column, peptide_column):
        pos = [(char,index+1) for index,char in enumerate(peptide) if char.islower()]
        for position,output_pos in zip(pos,df_pos):
            residue = protein_seq.seq[output_pos]
            peptide_residue = position[0].upper()
            assert residue == peptide_residue, f"{residue} in protein is not equal to {peptide_residue,position[1]} in peptide {peptide} Is the mrc db uptodate?"

def check_peptides_correct_position (peptide_column,position_column,protein_seq):
    for pep,df_ss in zip(peptide_column,position_column):
        
        pred_pep = protein_seq.seq[df_ss[0]:df_ss[1]]
        pep = iolo.capitalise_peptides(pep)
        assert pred_pep == pep, f"{pred_pep} is not equal to {pep}"