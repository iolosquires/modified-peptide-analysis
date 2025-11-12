from itertools import compress
import re
from Bio import SeqIO
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
from typing import Union
import numpy as np
import matplotlib.patches as patches
from pathlib import Path
import os
import matplotlib.gridspec as gridspec
import statistics
import pandas as pd
import pyteomics.mzid as mzid
from fpdf import FPDF
import shutil

from functions.ppa_dataclasses import proteinRecord

file_check_dictionary = {True: "three_file_check", False: "two_file_check"}

PHOSPHO_ST_PATTERN = "3"
PHOSPHO_Y_PATTERN = "4"



def save_everything_to_output(source_dir, save_dir):
    """
    Copy everything (files and folders) from source_dir to an 'output' subdirectory.

    Parameters
    ----------
    source_dir : str
        Path to the directory whose contents should be copied.
    """
    # Ensure the source directory exists
    if not os.path.isdir(source_dir):
        raise FileNotFoundError(f"Source directory not found: {source_dir}")

    # Create the output directory path
    output_dir = os.path.join(source_dir, save_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Copy all contents
    for item in os.listdir(source_dir):
        src_path = os.path.join(source_dir, item)
        dst_path = os.path.join(output_dir, item)

        # Skip copying the output folder into itself
        if src_path == output_dir:
            continue

        # Copy files or directories
        if os.path.isdir(src_path):
            shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
        else:
            shutil.copy2(src_path, dst_path)

    print(f"✅ All contents from '{source_dir}' copied to '{output_dir}'")



def load_and_process_pd_file(pd_filename, input_directory, config, search_record):
    pd_savename = str(pd_filename).replace(".txt", "")
    pd_file = mascot_file_check(pd_filename, input_directory)
    pd_output = pd.read_csv(str(pd_file), dtype={"Protein Accessions": "str"}, sep="\t")
    if "ptmRS: Best Site Probabilities" in pd_output.columns:
        pd_output["PhosphoRS: Best Site Probabilities"] = pd_output[
            "ptmRS: Best Site Probabilities"
        ]
    pd_output = pd_output[
        pd_output["Protein Accessions"].apply(lambda x: search_record in x)
    ].copy()
    pd_output = pd_output.dropna(subset=["PhosphoRS: Best Site Probabilities"]).copy()
    pd_output = pd_output[pd_output["Ions Score"] >= config.score_cutoff].copy()
    pd_output["phospho_rs"] = pd_output["PhosphoRS: Best Site Probabilities"].apply(
        mods_to_list_pd_output
    )
    pd_output["phospho_rs"] = pd_output["phospho_rs"].apply(remove_non_phospho_mods)
    pd_output = pd_output[pd_output["phospho_rs"] != "No Phospho"].copy()
    size_pd_df = len(pd_output)
    return pd_output, size_pd_df, pd_savename


def check_uniprot_id(config, paths, search_record):
    if not config.uniprot_id_check:
        mrc_db = mrc_db_to_dict(paths.mrc_db_path)
        POI_record = mrc_db[str(search_record)]
        return POI_record, mrc_db
    else:
        POI_record = proteinRecord(name=search_record, seq=config.seq)
        return POI_record, None


def mrc_db_to_dict(path_to_mrc):
    mrc_db = SeqIO.to_dict(SeqIO.parse(path_to_mrc, "fasta"))
    mrc_dict = {k.split("|")[0]: v for k, v in mrc_db.items()}
    return mrc_dict


def coverage_plot_v1(protein, peptide_position, output_path, title):
    sample_name = protein.name

    """ generates plot and saves to file """
    fig, ax = plt.subplots(figsize=(6, 1.75))  ## 6, 1.75## 10, 2.5
    patches = []
    ax.set_xlim(0, len(protein.seq))
    ax.set_ylim(0.5, 1.0)  ## 0.5, 1.0

    ax.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=True,  # ticks along the bottom edge are off
        top=False,  # True,         # ticks along the top edge are off
        labelbottom=True,
    )  # labels along the bottom edge are off
    ax.tick_params(  # plt.tick_params(
        axis="y",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        right=False,  # ticks along the bottom edge are off
        left=False,  # True,         # ticks along the top edge are off
        labelright=False,  # labels along the bottom edge are off
        labelleft=False,
    )  ##True)# labels along the bottom edge are off

    plot_list = []
    for n, item in enumerate(peptide_position):
        pep_len = item[1] - item[0]
        plot_list.append([item[0], 0.5, pep_len, 0.5])

    coordinates = plot_list

    for i in coordinates:
        rect = plt.Rectangle((i[0], i[1]), i[2], i[3], color="grey", alpha=0.5)  # 0.4
        patches.append(rect)

    p = PatchCollection(patches, cmap=mpl.cm.jet, color="grey", alpha=0.5)  # 0.3
    ax.add_collection(p)
    ax.set_title("peptide coverage   [sample: " + title + "]", fontsize=10)  # , pad=5
    ax.set_xlabel("AA position   [protein: " + str(sample_name) + "]", fontsize=10)
    ax.set_ylabel("-P", fontsize=10, rotation=0, ha="left", labelpad=20)  #
    ax.tick_params(axis="x", which="major", labelsize=10)

    if len(protein.seq) <= 300:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))  ## 200  ## 250
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))  ##  20  ##  50
    elif len(protein.seq) <= 500:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(100))  ## 200  ## 250
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))  ##  20  ##  50
    elif len(protein.seq) <= 1000:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(100))  ## 200  ## 250
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))  ##  20  ##  50
    elif len(protein.seq) >= 1000:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(250))  ## 200  ## 250
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))  ##  20  ##  50
    else:
        ax.xaxis.set_major_locator(ticker.AutoLocator())
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    plt.subplots_adjust(left=0.1, right=0.9, top=0.8, bottom=0.3)
    plot_dir = output_path / (title + "cov_plot.png")
    plt.savefig(plot_dir, dpi=300)
    plt.close()
    return fig


def get_phospho_position(peptide):
    position = []
    for match in re.finditer(PHOSPHO_ST_PATTERN + "|" + PHOSPHO_Y_PATTERN, peptide):
        position.append((match.start() - 1))

    return position


def _vis_psites_init(domain_position, domain_color, length):
    # from autoprot https://github.com/ag-warscheid/autoprot/blob/dev/autoprot/visualization/annotation.py
    if domain_position is None:
        domain_position = []
    # check if domain_color is a cmap name
    try:
        cm = plt.get_cmap(domain_color)
        color = cm(np.linspace(0, 1, len(domain_position)))
    except ValueError as e:
        if isinstance(domain_color, str):
            color = [
                domain_color,
            ] * len(domain_position)
        elif isinstance(domain_color, list):
            if len(domain_color) != len(domain_position):
                raise TypeError("Please provide one domain colour per domain") from e
            else:
                color = domain_color
        else:
            raise TypeError(
                "You must provide a colormap name, a colour name or a list of colour names"
            ) from e

    lims = (1, length)
    height = lims[1] / 25

    return lims, height, color, domain_position


def vis_psites(
    output_path,
    title,
    length: int,
    domain_position: Union[list[tuple[int]], None] = None,
    ps: Union[list[int], None] = None,
    pl: Union[list[str], None] = None,
    plc: Union[list[str], None] = None,
    pls: int = 4,
    ax: plt.Axes = None,
    domain_color: str = "tab10",
    ret_fig: bool = False,
):
    # noinspection PyUnresolvedReferences
    # noinspection PyShadowingNames
    # from autoprot https://github.com/ag-warscheid/autoprot/blob/dev/autoprot/visualization/annotation.py
    """
    Visualize domains and phosphosites on a protein of interest.

    Parameters
    ----------
    name : str
        Name of the protein.
        Used for plot title.
    length : int
        Length of the protein.
    domain_position : list of tuples of int
        Each element is a tuple of domain start and end postiions.
    ps : list of int
        position of phosphosites.
    pl : list of str
        label for ps (has to be in same order as ps).
    plc : list of colours
        optionally one can provide a list of colors for the phosphosite labels.
    pls : int, optional
        Fontsize for the phosphosite labels. The default is 4.
    ax: matplotlib axis, optional
        To draw on an existing axis
    domain_color: str
        Either a matplotlib colormap (see https://predictablynoisy.com/matplotlib/gallery/color/colormap_reference.html)
        or a single color
    ret_fig : bool
        Return fig as element.
        Default set to False.

    Returns
    -------
    matplotlib.figure
        The figure object.

    Examples
    --------
    Draw an overview on the phosphorylation of AKT1S1.

    >>> name = "AKT1S1"
    >>> length = 256
    >>> domain_position = [(35,43),
    ...                    (77,96)]
    >>> ps = [88, 92, 116, 183, 202, 203, 211, 212, 246]
    >>> pl = ["pS88", "pS92", "pS116", "pS183", "pS202", "pS203", "pS211", "pS212", "pS246"]

    colors (A,B,C,D (gray -> purple), Ad, Bd, Cd, Dd (gray -> teal) can be used to indicate regulation)

    >>> plc = ['C', 'A', 'A', 'C', 'Cd', 'D', 'D', 'B', 'D']
    >>> autoprot.visualization.vis_psites(name, length, domain_position, ps, pl, plc, pls=12)

    .. plot::
        :context: close-figs

        name = "AKT1S1"
        length = 256
        domain_position = [(35,43),
                           (77,96)]
        ps = [88, 92, 116, 183, 202, 203, 211, 212, 246]
        pl = ["pS88", "pS92", "pS116", "pS183", "pS202", "pS203", "pS211", "pS212", "pS246"]
        plc = ['C', 'A', 'A', 'C', 'Cd', 'D', 'D', 'B', 'D']
        vis.vis_psites(name, length, domain_position, ps, pl, plc, pls=12)
        plt.show()

    """

    lims, height, color, domain_position = _vis_psites_init(
        domain_position, domain_color, length
    )

    if ax is None:
        fig1 = plt.figure(figsize=(15, 2))
        ax1 = fig1.add_subplot(111, aspect="equal")
    else:
        ax1 = ax

    # background of the whole protein in grey
    ax1.add_patch(patches.Rectangle((0, 0), length, height, color="lightgrey"))

    for idx, (start, end) in enumerate(domain_position):
        width = end - start
        ax1.add_patch(patches.Rectangle((start, 0), width, height, color=color[idx]))

    # only plot phospho site if there are any
    if ps is not None:
        text_color = {
            "A": "gray",
            "Ad": "gray",
            "B": "#dc86fa",
            "Bd": "#6AC9BE",
            "C": "#aa00d7",
            "Cd": "#239895",
            "D": "#770087",
            "Dd": "#008080",
        }

        for idx, site in enumerate(ps):
            plt.axvline(site, 0, 1, color="red")
            plt.text(
                site - 1,
                height - (height + height * 0.15),
                pl[idx] if pl is not None else "",
                fontsize=pls,
                rotation=90,
                color=text_color[plc[idx]] if plc is not None else "black",
            )

    plt.subplots_adjust(left=0.25)
    plt.ylim(height)
    plt.xlim(lims)

    ax1.axes.get_yaxis().set_visible(False)

    plt.tight_layout()

    plot_dir = output_path / (title + "_phos_vis_plot.png")
    plt.savefig(plot_dir, dpi=300)
    plt.close()
    if ret_fig:
        return ax1


def mascot_file_check(mascot_filename, input_dir):
    mascot_file = str(input_dir) + "//" + mascot_filename

    assert Path(mascot_file).exists(), "mascot file %s does not exist" % mascot_filename

    return mascot_file


def coverage_plot_v3(
    protein,
    sample_name,
    start_position_list_1,
    pep_len_list_1,
    start_position_list_2,
    pep_len_list_2,
    output_path,
):
    """generates plot and saves to file"""

    mpl.rcParams["figure.figsize"] = (6, 2)  # 6, 1.75 ## 6, 2.5
    mpl.rcParams["savefig.dpi"]
    fig = plt.figure()
    # ---------------------------------------
    gs = gridspec.GridSpec(
        2,
        1,  ##3,2
        # width_ratios=[7, 1],
        height_ratios=[1, 1],
    )
    # ---------------------------------------
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    # ax3 = fig.add_subplot(gs[4])
    # ax4 = fig.add_subplot(gs[3])
    plt.subplots_adjust(
        left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.0
    )
    # ---------------------------------------
    ax1.grid(False)  # True, alpha=0.25, zorder=1)
    ax2.grid(False)  # True, alpha=0.25, zorder=1)
    # ax3.grid(True, alpha=0.25, zorder=1)
    # ax4.grid(True, alpha=0.25, zorder=1)
    # ---------------------------------------
    ax1.set_xlim(0, len(protein.seq))
    ax1.set_ylim(0.5, 1.0)  ## 0.5, 1.0
    # ax1.set_title("["+str(extra_data)+"]", fontsize=10) ## -> subtitle #12
    # ax1.axvline(x=-5, linewidth =1, linestyle='dashed' ,color='r', alpha=0.75, zorder=2)
    # ax1.axvline(x=5, linewidth =1, linestyle='dashed', color='r', alpha=0.75, zorder=2)
    ax1.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # True,      # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False,
    )  # labels along the bottom edge are off
    ax1.tick_params(  # plt.tick_params(
        axis="y",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        right=False,  # ticks along the bottom edge are off
        left=False,  # True,         # ticks along the top edge are off
        labelright=False,  # labels along the bottom edge are off
        labelleft=False,
    )  ##True)# labels along the bottom edge are off
    # ---------------------------------------
    ax2.set_xlim(0, len(protein.seq))
    ax2.set_ylim(0.5, 1.0)  ## 0.5, 1.0
    # ax1.set_title("["+str(extra_data)+"]", fontsize=10) ## -> subtitle #12
    # ax1.axvline(x=-5, linewidth =1, linestyle='dashed' ,color='r', alpha=0.75, zorder=2)
    # ax1.axvline(x=5, linewidth =1, linestyle='dashed', color='r', alpha=0.75, zorder=2)
    ax2.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=True,  # ticks along the bottom edge are off
        top=False,  # True,         # ticks along the top edge are off
        labelbottom=True,
    )  # labels along the bottom edge are off
    ax2.tick_params(  # plt.tick_params(
        axis="y",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        right=False,  # ticks along the bottom edge are off
        left=False,  # True,         # ticks along the top edge are off
        labelright=False,  # labels along the bottom edge are off
        labelleft=False,
    )  ##True)# labels along the bottom edge are off
    #
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontsize(10)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label2.set_fontsize(10)
    # ax1.xaxis.set_minor_locator(plt.MaxNLocator(50))
    # ax2.xaxis.set_minor_locator(plt.MaxNLocator(50))
    # ax1.xaxis.set_major_locator(ticker.MultipleLocator(200))  ## 200  ## 250
    # ax1.xaxis.set_minor_locator(ticker.MultipleLocator(50))   ##  20  ##  50
    # ax2.xaxis.set_major_locator(ticker.MultipleLocator(200))  ## 200  ## 250
    # ax2.xaxis.set_minor_locator(ticker.MultipleLocator(50))   ##  20  ##  50
    if len(protein.seq) <= 300:
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))  ## 200  ## 250
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(10))  ##  20  ##  50
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(50))  ## 200  ## 250
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(10))  ##  20  ##  50
    elif len(protein.seq) <= 500:
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))  ## 200  ## 250
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(10))  ##  20  ##  50
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(100))  ## 200  ## 250
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(10))  ##  20  ##  50
    elif len(protein.seq) <= 1000:
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))  ## 200  ## 250
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(50))  ##  20  ##  50
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(100))  ## 200  ## 250
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(50))  ##  20  ##  50
    elif len(protein.seq) >= 1000:
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(250))  ## 200  ## 250
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(50))  ##  20  ##  50
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(250))  ## 200  ## 250
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(50))  ##  20  ##  50
    else:
        ax1.xaxis.set_major_locator(ticker.AutoLocator())
        ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax2.xaxis.set_major_locator(ticker.AutoLocator())
        ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    #
    ax1.set_title(
        "peptide coverage   [sample: " + sample_name + "]", fontsize=10
    )  # , pad=5
    ##ax1.set_xlabel('[AA position]', fontsize=12)
    ##ax2.set_title('peptide coverage', fontsize=16)#, pad=5
    ax2.set_xlabel("AA position" + protein.name, fontsize=10)
    ax1.set_ylabel("-P", fontsize=10, rotation=0, ha="left", labelpad=20)  #
    ax2.set_ylabel("+P", fontsize=10, rotation=0, ha="left", labelpad=20)  #
    # ---------------------------------------
    patches_1 = []
    patches_2 = []
    plot_list_1 = []
    plot_list_2 = []

    for n in range(len(start_position_list_1)):
        start = start_position_list_1[n]
        pep_len = pep_len_list_1[n]
        plot_list_1.append([start, 0.5, pep_len, 0.5])
    #
    for n in range(len(start_position_list_2)):
        start = start_position_list_2[n]
        pep_len = pep_len_list_2[n]
        plot_list_2.append([start, 0.5, pep_len, 0.5])
    #
    coordinates_1 = plot_list_1
    coordinates_2 = plot_list_2
    #
    for i in coordinates_1:
        rect = plt.Rectangle((i[0], i[1]), i[2], i[3], color="grey", alpha=0.5)  # 0.4
        patches_1.append(rect)
    #
    for i in coordinates_2:
        rect = plt.Rectangle((i[0], i[1]), i[2], i[3], color="red", alpha=0.5)  # 0.4
        patches_2.append(rect)
    #
    p1 = PatchCollection(
        patches_1, cmap=mpl.cm.jet, alpha=0.5, color="grey"
    )  # 0.3 ## matplotlib
    p2 = PatchCollection(
        patches_2, cmap=mpl.cm.jet, alpha=0.5, color="red"
    )  # 0.3 ## matplotlib

    ax1.add_collection(p1)
    ax2.add_collection(p2)

    plt.subplots_adjust(left=0.1, right=0.9, top=0.85, bottom=0.25)
    plot_dir = output_path / (sample_name + "_cov_plot.png")
    plt.savefig(plot_dir, dpi=600)
    plt.close()


def contains_search_term(accession_list, search_term):
    return any(search_term == str(item) for item in accession_list)


def contains_search_term_uniprot(accession_list, search_term):
    return any(search_term in str(item) for item in accession_list)


def get_phospho_positions(mod_dict_list):
    phos_position = ""
    for mod_dict in mod_dict_list:
        if mod_dict["name"] == "Phospho":
            phos_position += str(mod_dict["location"]) + "p"
    return phos_position


def get_logo_file():
    logo_file = Path(os.getcwd()) / "logo" / "ppu_logo.png"

    assert logo_file.exists(), "logo file does not exist"

    return logo_file


def find_phospho_mod(mod_dict_list):
    try:
        for mod_dict in mod_dict_list:
            if mod_dict["name"] == "Phospho":
                found_phospho = True
                break
            else:
                found_phospho = False
    except:  # noqa: E722
        found_phospho = False

    return found_phospho


def find_phospho_single_acceptor_site(peptide):
    if peptide.count("S") + peptide.count("T") + peptide.count("Y") == 1:
        return True
    else:
        return False


def find_phospho_double_acceptor_site(peptide):
    if peptide.count("S") + peptide.count("T") + peptide.count("Y") == 2:
        return True
    else:
        return False


def mods_to_list_pd_output(string):
    test2 = "".join(string.split()).split(";")
    test4 = [item.split(":") for item in test2]
    test5 = [item[0].split("(") for item in test4]

    new_list = []
    for item1, item2 in zip(test4, test5):
        merge = [item1[1]] + item2

        new_list.append(merge)
    return new_list


def remove_non_phospho_mods(mod_list):
    new_list = []
    for item in mod_list:
        if item[2] == "Phospho)":
            new_list.append(item)

    if new_list == []:
        new_list = "No Phospho"
    return new_list


def find_ubiquitin_mod(mod_dict_list):
    try:
        for mod_dict in mod_dict_list:
            if mod_dict["name"] == "GG":
                found_phospho = True
                break
            else:
                found_phospho = False
    except:  # noqa: E722
        found_phospho = False

    return found_phospho


def find_ubi_single_acceptor_site(peptide):
    if peptide.count("K") == 1:
        return True
    else:
        return False


def find_ubi_double_acceptor_site(peptide):
    if peptide.count("K") == 2:
        return True
    else:
        return False


def get_ubi_positions(mod_dict_list):
    phos_position = ""
    for mod_dict in mod_dict_list:
        if mod_dict["name"] == "GG":
            phos_position += str(mod_dict["location"]) + "u"
    return phos_position


def get_peptide_from_pd_output(string):
    pep_search = re.search(r"\.\w+\.", string)

    pep_search = pep_search.group(0)[1:-1]
    return pep_search


def capitalise_peptides(pep):
    replace_dict = {"s": "S", "t": "T", "y": "Y", "m": "M", "c": "C", "k": "K"}

    pep2 = "".join([replace_dict.get(item, item) for item in pep])

    return pep2


def get_position(test):
    return int(test[-1])


def combine_pd_mascot_confs(pd_conf_dict, conf_dict):
    merge_conf_dict = {}

    for key in pd_conf_dict.keys():
        merge_conf_dict[key] = max(conf_dict.get(key, 0), pd_conf_dict[key])

    return conf_dict | merge_conf_dict


def get_conf_pd(conf_list):
    return [float(item[0]) for item in conf_list]


def get_pos_pd(position_list):
    pos = [int(item[1][1:]) for item in position_list]
    pos_string = str(pos)
    return pos_string


def mean_conf_pd(conf_list):
    return statistics.mean(conf_list)


def phos_pos_to_ints(string):
    test = "".join(string.split())
    # turn test back into a list
    test2 = test[1:-1].split(",")
    test3 = [int(item) for item in test2]
    return test3


def phos_site_in_protein(start_stop, phos_in_peptide):
    phos_positions = []
    for site in phos_in_peptide:
        phos_positions.append(site + start_stop[0] - 2)
    return phos_positions


def lowercase_modified_residue(peptide, pos):
    replace_dict = {"S": "s", "T": "t", "Y": "y", "m": "M", "C": "c", "K": "k"}
    new = []
    for index, letter in enumerate(peptide):
        if index + 1 in pos:
            new.append(replace_dict.get(letter, letter))
        else:
            new.append(letter)

    return "".join(new)


def phospho_position_string_to_list(string):
    test_split = string.split("p")
    phospho_positions = []
    for item in test_split:
        if item == "":
            continue
        else:
            phospho_positions.append(int(item))
    return phospho_positions


def ubi_position_string_to_list(string):
    test_split = string.split("u")
    phospho_positions = []
    for item in test_split:
        if item == "":
            continue
        else:
            phospho_positions.append(int(item))
    return phospho_positions


def filter_both_confidences(mascot, phosphors, cutoff):
    mascot_check = False
    phosphors_check = False

    # if tm == True:
    # return True

    if mascot >= cutoff:
        mascot_check = True
    if isinstance(phosphors, list):
        if any(item >= cutoff for item in phosphors):
            phosphors_check = True
    if mascot_check | phosphors_check:
        return True
    else:
        return False


def align_peptide_to_protein(peptide, protein):
    searching = True
    end = 0
    while searching:
        start = protein.seq[end:].find(peptide)
        if start == -1:
            searching = False

        else:
            end = start + len(peptide)

            position = (start + 1, end)

    return position


def extract_between_pipes(lst):
    return [re.search(r"\|(.*?)\|", item).group(1) for item in lst]


def filter_accession_and_score(df, config, search_record):
    return df[
        (df["Mascot:score"] >= float(config.score_cutoff))
        & (
            df["accession"].apply(
                lambda x: search_func_dict[config.uniprot_id_check](
                    x, search_record
                )
            )
        )
    ]


search_func_dict = {True: contains_search_term_uniprot, False: contains_search_term}


def import_mascot_file(mascot_filename, input_directory):
    mascot_file = mascot_file_check(mascot_filename, input_directory)
    return mzid.DataFrame(str(mascot_file))


def find_phosphopeptides_in_mascot(df_wanted, config,search_record):
    df_wanted["start_end"] = df_wanted.apply(
        lambda x: list(zip(x.start, x.end)), axis=1
    )
    if config.uniprot_id_check:
        df_wanted["accession"] = df_wanted["accession"].apply(extract_between_pipes)
    df_wanted["accid_start_end_dict"] = df_wanted.apply(
        lambda x: dict(zip(x.accession, x.start_end)), axis=1
    )
    df_wanted["poi_start_stop"] = df_wanted.apply(
        lambda x: x.accid_start_end_dict[search_record], axis=1
    )
    df_wanted["has_phospho"] = df_wanted["Modification"].apply(find_phospho_mod)
    df_wanted_phospho = df_wanted[df_wanted["has_phospho"]].copy()

    return df_wanted_phospho


def process_mascot_phospho_dataframe(df_wanted_phospho):
    df_wanted_phospho["phos_positions"] = df_wanted_phospho["Modification"].apply(
        get_phospho_positions
    )
    phospho_cols_for_grouping = [
        "PeptideSequence",
        "chargeState",
        "phos_positions",
        "Mascot:score",
        "Mascot:PTM site assignment confidence",
        "poi_start_stop",
        "experimentalMassToCharge",
    ]
    df_phospho_grouped = df_wanted_phospho[phospho_cols_for_grouping].copy()

    df_phospho_grouped["single_acceptor"] = df_phospho_grouped["PeptideSequence"].apply(
        find_phospho_single_acceptor_site
    )
    df_phospho_grouped["double_acceptor"] = df_phospho_grouped["PeptideSequence"].apply(
        find_phospho_double_acceptor_site
    )
    df_phospho_grouped["missing_conf"] = df_phospho_grouped[
        "Mascot:PTM site assignment confidence"
    ].isna()

    condition_single = (df_phospho_grouped["single_acceptor"]) & (
        df_phospho_grouped["missing_conf"]
    )
    condition_double = (df_phospho_grouped["double_acceptor"]) & (
        df_phospho_grouped["missing_conf"]
    )
    condition_missing = (
        (~df_phospho_grouped["single_acceptor"])
        & (df_phospho_grouped["missing_conf"])
        & (~df_phospho_grouped["double_acceptor"])
    )

    df_phospho_grouped["Mascot:PTM site assignment confidence"] = df_phospho_grouped[
        "Mascot:PTM site assignment confidence"
    ].fillna(condition_single.map({True: 100}))
    df_phospho_grouped["Mascot:PTM site assignment confidence"] = df_phospho_grouped[
        "Mascot:PTM site assignment confidence"
    ].fillna(condition_double.map({True: 100}))
    df_phospho_grouped["Mascot:PTM site assignment confidence"] = df_phospho_grouped[
        "Mascot:PTM site assignment confidence"
    ].fillna(condition_missing.map({True: 0}))

    group_has_conf_site = df_phospho_grouped.groupby(
        ["PeptideSequence", "chargeState", "phos_positions"]
    )["Mascot:PTM site assignment confidence"].transform(lambda x: x.max() >= 80)
    df_phospho_grouped["top_mascot_and_over80"] = (
        df_phospho_grouped.groupby(
            ["PeptideSequence", "chargeState", "phos_positions"]
        )["Mascot:score"].transform(lambda x: x == x.max())
    ) & group_has_conf_site
    df_phospho_grouped = df_phospho_grouped.round({"Mascot:score": 0})

    phos_pep_lengths = df_phospho_grouped["PeptideSequence"].str.len().to_list()
    phos_start_pos = [item[0] for item in df_phospho_grouped["poi_start_stop"]]
    df_phospho_grouped["phos_positions_list"] = df_phospho_grouped[
        "phos_positions"
    ].apply(phospho_position_string_to_list)
    df_phospho_grouped["phos_in_protein"] = df_phospho_grouped.apply(
        lambda x: phos_site_in_protein(x["poi_start_stop"], x["phos_positions_list"]),
        axis=1,
    )

    return df_phospho_grouped, phos_pep_lengths, phos_start_pos


def process_pd_phospho_dataframe(pd_output, df_phospho_grouped, POI_record):
    pd_output["peptide_trim"] = pd_output["Annotated Sequence"].apply(
        get_peptide_from_pd_output
    )
    pd_output["phosphors_conf"] = pd_output["phospho_rs"].apply(get_conf_pd)
    pd_output["phos_pos_pep"] = pd_output["phospho_rs"].apply(get_pos_pd)
    pd_output["peptide_cap"] = pd_output["peptide_trim"].apply(capitalise_peptides)
    pd_output["poi_start_stop"] = pd_output["peptide_cap"].map(
        dict(
            zip(
                df_phospho_grouped["PeptideSequence"],
                df_phospho_grouped["poi_start_stop"],
            )
        )
    )

    pd_output["mean_conf"] = pd_output["phosphors_conf"].apply(mean_conf_pd)

    pd_grouped = (
        pd_output.sort_values("mean_conf", ascending=False)
        .groupby(["peptide_trim", "Charge", "phos_pos_pep"], as_index=False)
        .first()
    )
    pd_grouped["phos_positions_list"] = pd_grouped["phos_pos_pep"].apply(
        phos_pos_to_ints
    )

    # If nans in poi_start_stop, align peptide to protein using own function.

    if pd_grouped["poi_start_stop"].isnull().sum() > 0:
        # split intwo two dataframes, one with nans and one without
        pd_grouped_nan = pd_grouped[pd_grouped["poi_start_stop"].isnull()].copy()
        pd_grouped_not_nan = pd_grouped.dropna(subset=["poi_start_stop"]).copy()
        pd_grouped_nan["poi_start_stop"] = pd_grouped_nan.apply(
            lambda x: align_peptide_to_protein(x["peptide_cap"], POI_record), axis=1
        )
        pd_grouped = pd.concat([pd_grouped_nan, pd_grouped_not_nan], axis=0)

    pd_grouped["phos_in_protein"] = pd_grouped.apply(
        lambda x: [
            item + x["poi_start_stop"][0] - 2 for item in x["phos_positions_list"]
        ],
        axis=1,
    )
    pd_grouped["lc_pep"] = pd_grouped.apply(
        lambda x: lowercase_modified_residue(
            x["peptide_cap"], x["phos_positions_list"]
        ),
        axis=1,
    )
    pd_grouped["lc_pep_cs"] = (
        pd_grouped["Charge"].astype(str) + "-" + pd_grouped["lc_pep"]
    )
    pd_peptide_site_dict = dict(
        zip(pd_grouped["lc_pep_cs"], pd_grouped["phosphors_conf"])
    )
    pd_peptide_is_dict = dict(zip(pd_grouped["lc_pep_cs"], pd_grouped["Ions Score"]))
    pd_grouped_explode = pd_grouped.explode(
        ["phosphors_conf", "phos_positions_list", "phos_in_protein"]
    )
    pd_output_ge_grouped = (
        pd_grouped_explode.sort_values("phosphors_conf", ascending=False)
        .groupby(["phos_in_protein"], as_index=False)
        .first()
    )
    pd_conf_dict = dict(
        zip(
            pd_output_ge_grouped["phos_in_protein"],
            pd_output_ge_grouped["phosphors_conf"],
        )
    )
    return pd_grouped, pd_peptide_site_dict, pd_peptide_is_dict, pd_conf_dict


def add_pd_info_to_mascot(df_phospho_grouped, pd_peptide_site_dict, pd_peptide_is_dict):
    df_phospho_grouped["pd_peptide"] = df_phospho_grouped.apply(
        lambda x: lowercase_modified_residue(
            x["PeptideSequence"], x["phos_positions_list"]
        ),
        axis=1,
    )
    df_phospho_grouped["pd_peptide_cs"] = (
        df_phospho_grouped["chargeState"].astype(str)
        + "-"
        + df_phospho_grouped["pd_peptide"]
    )
    df_phospho_grouped["pd_conf"] = df_phospho_grouped["pd_peptide_cs"].map(
        pd_peptide_site_dict
    )
    df_phospho_grouped["pd_is"] = df_phospho_grouped["pd_peptide_cs"].map(
        pd_peptide_is_dict
    )
    return df_phospho_grouped


def process_merged_dataframe(merged):
    merged = merged.round({"M/Z": 1})
    merged["Mascot Group Conf"] = merged.groupby(["Peptide", "Charge State", "M/Z"])[
        "PTM Confidence Mascot"
    ].transform(lambda x: [x.tolist()] * len(x))
    merged = merged.sort_values(by=["Mascot Score"], ascending=False).drop_duplicates(
        subset=["Peptide", "Charge State", "M/Z"], keep="first"
    )
    return merged


def create_conf_dict(df_phospho_grouped):
    df_phospho_conf_dict = df_phospho_grouped.sort_values(
        "Mascot:PTM site assignment confidence", ascending=False
    )
    df_phospho_conf_dict2 = df_phospho_conf_dict.explode("phos_in_protein")
    df_phospho_conf_dict2 = df_phospho_conf_dict2.drop_duplicates(
        subset="phos_in_protein"
    )
    conf_dict = dict(
        zip(
            df_phospho_conf_dict2["phos_in_protein"],
            df_phospho_conf_dict2["Mascot:PTM site assignment confidence"],
        )
    )
    return conf_dict


def filter_merged_dataframe_by_localisation_confidence(merged, config):
    m = merged.apply(
        lambda x: filter_both_confidences(
            x["PTM Confidence Mascot"],
            x["PTM Confidence PD"],
            config.site_localisation_cutoff,
        ),
        axis=1,
    )
    merged_filter = merged[m].copy()
    merged_filter["Start Stop"] = merged_filter["Start Stop"].apply(
        lambda x: (x[0] - 1, x[1])
    )

    return merged_filter


def create_pdf_data(merged_filter):
    df_pdf = merged_filter.map(str)
    columns = [
        [
            "Mascot Score",
            "PD Score",
            "Charge State",
            "M/Z",
            "Peptide",
            "Conf Mascot",
            "Conf PD",
            "Position in Protein",
            "Start Stop",
            "Mascot Group Conf",
        ]
    ]
    rows = df_pdf.values.tolist()
    data = columns + rows
    return data


def create_excel_report(config, merged_data, file_contains_phospho, paths):
    excel_dir = paths.output_directory / (
        config.analysis_name + "_phosphopeptide_report.xlsx"
    )
    data_for_excel_filter = list(compress(merged_data, file_contains_phospho))
    sample_names_filter = list(compress(config.sample_names, file_contains_phospho))

    with pd.ExcelWriter(excel_dir, engine="xlsxwriter") as writer:
        for data, mascot_file in zip(data_for_excel_filter, sample_names_filter):
            data.to_excel(writer, sheet_name=mascot_file, index=False)


def create_pdf_report(
    config, merge_data_pdf, file_contains_phospho, pd_file_contains_phospho, paths
):
    pdf = FPDF()
    col_width_dict = {
        True: (5, 5, 5, 3, 12, 5, 5, 5, 5, 5),
        False: (5, 5, 3, 12, 5, 5, 5, 5),
    }

    for data, mascot_file, phospho_indicator, pd_phospho in zip(
        merge_data_pdf,
        config.sample_names,
        file_contains_phospho,
        pd_file_contains_phospho,
    ):
        if phospho_indicator:
            pdf.add_page()
            pdf.set_font("Times", style="B", size=18)
            pdf.cell(0, 10, mascot_file, align="C")
            pdf.ln(10)
            pdf.set_font("Times", size=6)
            # add logo image
            pdf.image(str(paths.logo_file), x=20, y=5, w=40, h=15)
            with pdf.table(
                borders_layout="SINGLE_TOP_LINE",
                cell_fill_color=200,  # grey
                cell_fill_mode="ROWS",
                line_height=pdf.font_size * 2.5,
                text_align="CENTER",
                width=160,
                col_widths=col_width_dict[pd_phospho],
            ) as table:
                for data_row in data:
                    row = table.row()
                    for datum in data_row:
                        row.cell(datum)

        else:
            pdf.add_page()
            pdf.set_font("Times", style="B", size=18)
            pdf.cell(0, 10, mascot_file, align="C")
            pdf.ln(10)
            pdf.set_font("Times", size=7)
            pdf.image(str(paths.logo_file), x=20, y=10, w=40, h=15)

    pdf_dir = paths.output_directory / (
        config.analysis_name + "_phosphopeptide_report.pdf"
    )
    pdf.output(pdf_dir)

def open_read_experimental_design (input_directory):
    
    experimental_design_file = Path(input_directory) / "ExperimentalDesign.txt"
    assert experimental_design_file.exists(), "experimental design file does not exist"

    ed = pd.read_csv(experimental_design_file, sep="\t")
    return ed

