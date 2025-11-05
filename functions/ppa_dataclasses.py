import dataclasses
import functions.ppa_functions as iolo
import re


@dataclasses.dataclass
class proteinRecord:
    name: str
    seq: str

@dataclasses.dataclass
class configInfo:
    analysis_name: str
    mascot_filenames: list
    pd_filenames: list
    sample_names: list
    mrc_db: str
    search_protein: str
    seq : str
    score_cutoff: float
    species: str
    uniprot_id: str
    mod_search: str
    uniprot_id_check: bool = dataclasses.field(default=False)

    def __post_init__(self):
        if re.search('[A-Z]', self.search_protein):
            self.uniprot_id_check = True



@dataclasses.dataclass
class configPathInfo:
    input_directory: str
    mrc_db: str
    alphamap_python_script: str
    alphamap_python_path: str
    mrc_db_path: str = dataclasses.field(init=False)
    output_directory: str = dataclasses.field(init=False)
    logo_file: str = dataclasses.field(init=False)
    log_file : str = dataclasses.field(init=False)
    plot_path: str = dataclasses.field(init=False)

    def __post_init__(self):
        self.output_directory = self.input_directory + "/output"
        self.logo_file = iolo.get_logo_file()
        self.log_file = self.output_directory + "/run.log"
        self.plot_path = self.output_directory + "/plots"
        self.mrc_db_path = "Z:/proteinchem/CURRENT MRC DATABASE/" + self.mrc_db
    

def create_config_from_toml(config: dict) -> configInfo:
    
    return configInfo(
        analysis_name = config['analysis_name'],
        mascot_filenames = config['mascot_filename'],
        pd_filenames = config['pd_filename'],
        sample_names = config['sample_name'],
        search_protein = config['search_protein'],
        seq = config['seq'],
        score_cutoff = config['score_cutoff'],
        species = config['species'],
        uniprot_id = config['uniprot_for_plot'],
        mod_search = config['mod_search']
    )

def create_paths_from_toml(config: dict, 
                           input_directory: str,
                           alphamap_python_script: str,
                           alphamap_python_path: str) -> configPathInfo:
    
    return configPathInfo(
        input_directory = input_directory,
        mrc_db = config['mrc_db'],
        alphamap_python_script = alphamap_python_script,
        alphamap_python_path = alphamap_python_path
    )