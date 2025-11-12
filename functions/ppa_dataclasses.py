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
    search_protein_list: list
    seq: str
    score_cutoff: float
    species: str
    uniprot_id: str
    mod_search: str
    uniprot_id_check: bool
    site_localisation_cutoff: float = dataclasses.field(default=0.0)
    create_old_plot: bool = dataclasses.field(default=False)

    def __post_init__(self):
        # turn search_protein_list into list of strings
        self.search_protein_list = [str(item) for item in self.search_protein_list]
     


@dataclasses.dataclass
class configPathInfo:
    input_directory: str
    mrc_db: str
    alphamap_python_script: str
    alphamap_python_path: str
    mrc_db_path: str = dataclasses.field(init=False)
    output_directory: str = dataclasses.field(init=False)
    logo_file: str = dataclasses.field(init=False)
    log_file: str = dataclasses.field(init=False)
    plot_path: str = dataclasses.field(init=False)

    def __post_init__(self):
        # self.input_directory = str(self.input_directory)
        self.output_directory = self.input_directory / "output"
        self.logo_file = iolo.get_logo_file()
        self.log_file = self.output_directory / "run.log"
        self.plot_path = self.output_directory / "plots"
        self.mrc_db_path = "Z:/proteinchem/CURRENT MRC DATABASE/" + self.mrc_db


def create_config_from_toml(config: dict, ExperimentalDesign) -> configInfo:
    return configInfo(
        
        mascot_filenames=ExperimentalDesign["mascot_filename"],
        pd_filenames=ExperimentalDesign["pd_filename"],
        sample_names=ExperimentalDesign["sample_name"],
        search_protein_list=ExperimentalDesign["search_protein"],

        analysis_name=config["analysis_name"],
        mrc_db=config["mrc_db"],
        seq=config["seq"],
        score_cutoff=config["score_cutoff"],
        species=config["species"],
        uniprot_id=config["uniprot_for_plot"],
        mod_search=config["mod_search"],
        create_old_plot = config["create_old_plot"],
        uniprot_id_check=config["uniprot_id_check"]

    )


def create_paths_from_toml(
    config: dict,
    input_directory: str,
    alphamap_python_script: str,
    alphamap_python_path: str,
) -> configPathInfo:
    return configPathInfo(
        input_directory=input_directory,
        mrc_db=config["mrc_db"],
        alphamap_python_script=alphamap_python_script,
        alphamap_python_path=alphamap_python_path,
    )
