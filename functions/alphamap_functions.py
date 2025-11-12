from pathlib import Path
import pandas as pd

def create_assigned_mods(test):
    mod_dict = {'S': '(79.9663)',
                'T': '(79.9663)',
                'Y': '(79.9663)',
                'K': '(114.0429)'
                }
    
    test_split = list(test)
    #find positions of lowercase letters
    lowercase_positions = [i for i, letter in enumerate(test_split) if letter.islower()]
    residue = [test_split[i].upper() for i in lowercase_positions]

    assigned_mods = []

    for i, res in zip(lowercase_positions, residue):
        mod_string = str(i+1) + res + mod_dict[res]
        assigned_mods.append(mod_string)
    
    return ",".join(assigned_mods)

def open_read_experimental_design (input_directory):
    
    experimental_design_file = Path(input_directory) / "ExperimentalDesign.txt"
    assert experimental_design_file.exists(), "experimental design file does not exist"

    ed = pd.read_csv(experimental_design_file, sep="\t")
    return ed


def mascot_script_to_alphamap(df, protein_id, output_dir):


    df['Protein ID'] = protein_id
    df['Assigned Modifications'] = df['Peptide'].apply(create_assigned_mods)
    df['Peptide'] = df['Peptide'].str.upper()

    df = df.groupby(["Peptide","Assigned Modifications"], as_index = False).agg({
        "Protein ID" : "first",
        "Charge State": (lambda x: ",".join(map(str, x))),
        # for all other columns: turn into lists
        **{col: "first" for col in df.columns if col not in ["Protein ID", "Peptide","Assigned Modifications","Charge State"]},
    })
    df.to_csv(output_dir, index=False, sep = "\t")
