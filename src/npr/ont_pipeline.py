import glob
import json
import os
import re
import sys
from pathlib import Path
import pandas as pd
import requests

from npr.communication import send_email
from npr.snakehelper import get_seqdir, glob2reports


def analysis_done(flowcell, config):
    """
    Determine whether flowcell has already been analysed
    Notice that such an analysis could have been done at different locations
    """
    loc1 = os.path.join(
        config["paths"]["outputDir"], os.path.basename(flowcell), "analysis.done"
    )
    if os.path.exists(loc1):
        return True

    return False


def find_new_flowcell(config):
    """
    look for new flowcells inf offloadDir
    baseDir flowcells are flagged with 'transfer.final'
    """
    #################################### offload ####################################
    offload_path = os.path.join(config["paths"]["offloadDir"])
    dirs = []
    if offload_path:
        # Assume that report*.json mark 'ready' flow cells
        # depending on the offload_path
        pattern = os.path.join(offload_path, "**", "report*.json")
        jsons = glob.glob(pattern, recursive=True)
        for j in jsons:
            dirs.append(os.path.dirname(j))

    # filter flowcells based on config['target_flowcell']
    if config["target_flowcell"]:
        dirs = [d for d in dirs if config["target_flowcell"] in d]
        if not dirs and config["force_processing"]:
            pattern = os.path.join(offload_path, "*")
            # I guess we can skip filtering here ?
            dirs = [d for d in glob.glob(pattern) if config["target_flowcell"] in d]

    # Iterate over all flowcell in dir
    for flowcell in dirs:
        if config["options"]["verbosity"]:
            print(f"Working with {flowcell}")

        # test and continue if 'analysis.done' exists for this flowcell
        if analysis_done(flowcell, config):
            if config["options"]["verbosity"]:
                print(f" {flowcell} found")
            continue

        # needed here to communicate flowcell with send_email
        config["info_dict"]["base_path"] = flowcell

        ss = os.path.join(flowcell, "SampleSheet.csv")

        flowcell_path = os.path.join(
            config["paths"]["outputDir"],
            os.path.basename(flowcell),
        )

        ss2 = os.path.join(
            flowcell_path,
            "reports",
            "SampleSheet.csv",
        )
        Path(flowcell_path, 'reports').mkdir(parents=True, exist_ok=True)

        if not os.path.isfile(ss):
            _ss_msg = f"No SampleSheet.csv file found in {flowcell_path}\n"
            if not get_samplesheet_from_parkour(flowcell.split("_")[-2], config, ss2):
                _ss_msg += f"Failed to retrieve sampleSheet via parkour for {flowcell.split("_")[-2]}\n"
                if not os.path.isfile(ss2):
                    _ss_msg += f"No SampleSheet found at {ss2}\n"
                    send_email("Error for flowcell:", _ss_msg, config, failure=True)
                    sys.exit("no sampleSheet.")

        # return flowcell to ont()
        msg = "SampleSheet.csv file found.\n"
        config["input"] = {"name": os.path.basename(flowcell)}

        return (os.path.basename(flowcell), msg, flowcell)

    return (None, None, None)


def get_samplesheet_from_parkour(flowcell, config, output_csv_path):
    """
    Saves the samplesheet from parkour, directly to the outputDir
    """
    postfixes = ["-5", "-4", "-3", "-2", "-1", ""]

    for pf in postfixes:
        flowcell_id = flowcell + pf
        params = {"flowcell_id": flowcell_id}

        try:
            res = requests.get(
                config["parkour"]["url"] + "/api/flowcells/retrieve_samplesheet/",
                auth=(config["parkour"]["user"], config["parkour"]["password"]),
                params=params,
                verify=config["parkour"]["pem"],
            )
            if res.status_code == 200:
                with open(output_csv_path, "wb") as file:
                    file.write(res.content)
                print(
                    f"Samplesheet for flowcell {flowcell_id} available at {output_csv_path}"
                )

                return True

        except requests.exceptions.RequestException as e:
            print(f"Request failed for flowcell {flowcell_id}: {e}")
            return False

        except OSError as e:
            print(f"File operation failed for flowcell {flowcell_id}: {e}")
            return False


def read_flowcell_info(config, info_dict, base_path):
    """
    - set paths
    - copy json, summaries & samplesheet
    - parse HTML and json file for flowcell, kit & barcoding
    """

    flowcell = config["input"]["name"]
    info_dict["base_path"] = base_path
    flowcell_path = os.path.join(config["paths"]["outputDir"], flowcell)
    if not os.path.exists(flowcell_path):
        os.mkdir(flowcell_path)
    info_dict["flowcell_path"] = flowcell_path
    # Copy over reports.
    for globStr in ["*json", "*html", "*txt", "SampleSheet.csv"]:
        glob2reports(globStr, base_path, flowcell_path)
    json_file = glob.glob(os.path.join(flowcell_path, "reports", "*json"))
    html_file = glob.glob(os.path.join(flowcell_path, "reports", "*html"))
    if html_file:
        match_html_key_val = re.compile(r'"title": "(.+?)", "value": "(.+?)"')

        print("[green]Reading info from html[/green]")
        html_cont = ""
        with open(html_file[0]) as f:
            html_cont = f.read()

        info_dict["software"] = {}
        matches = re.finditer(match_html_key_val, html_cont)
        for match in matches:
            name = match.group(1)
            value = match.group(2)
            if name == "Flow cell type":
                info_dict["flowcell"] = value
            elif name == "Kit type":
                info_dict["kit"] = value
            elif name == "MinKNOW":
                info_dict["software"]["MinKNOW"] = value
            elif name == "Bream":
                info_dict["software"]["Bream"] = value
            elif name == "MinKNOW Core":
                info_dict["software"]["MinKNOW Core"] = value
            elif name == "Configuration":
                info_dict["software"]["Configuration"] = value
            elif name == "Dorado":
                info_dict["software"]["Dorado"] = value

    if json_file:
        print("[green]Reading info from json[/green]")
        with open(json_file[0]) as f:  # assume only 1 file
            jsondata = json.load(f)

        if "flowcell" not in info_dict:
            info_dict["flowcell"] = jsondata["protocol_run_info"]["meta_info"]["tags"][
                "flow cell"
            ]["string_value"]
        if "kit" not in info_dict:
            info_dict["kit"] = jsondata["protocol_run_info"]["meta_info"]["tags"][
                "kit"
            ]["string_value"]
        # HTML file is not reporting barcoding, we need it from the json
        info_dict["barcoding"] = bool(
            jsondata["protocol_run_info"]["meta_info"]["tags"]["barcoding"][
                "bool_value"
            ]
        )

        for rg in jsondata["protocol_run_info"]["args"]:
            if rg == "--barcoding" and not info_dict["barcoding"]:
                print("Json bool for barcoding wrong ! Override !")
                info_dict["barcoding"] = True
            if rg.startswith("barcoding_kits"):
                start = rg.find("[") + 1
                stop = rg.find("]")
                info_dict["barcode_kit"] = rg[start:stop].replace('"', "")
        if "barcode_kit" not in info_dict:
            if info_dict["barcoding"]:
                print(
                    "[red] Not barcoding kit found in json. Default to flowcell kit.[/red]"
                )
                info_dict["barcode_kit"] = info_dict["kit"]
            else:
                info_dict["barcode_kit"] = "no_bc"

    print(f'flowcell_info_parsed: flowcell = {info_dict["flowcell"]}')
    print(f'flowcell_info_parsed: kit = {info_dict["kit"]}')
    print(f'flowcell_info_parsed: barcode_kit = {info_dict["barcode_kit"]}')
    return info_dict


def read_samplesheet(config):
    """
    read samplesheet
    """
    sample_sheet = pd.read_csv(
        os.path.join(
            config["info_dict"]["flowcell_path"], "reports", "SampleSheet.csv"
        ),
        sep=",",
        skiprows=1,  # first line is: [Data],,,,,,,,,,
    )

    # convert: No_index? --> no_bc
    if "I7_Index_ID" not in sample_sheet.columns:
        print("[red] Column I7_Index_ID not defined in sample_sheet [/red]")
        exit(1)
    # legacy fix: replace missing barcode (NaN or 'No_index.*') by 'no_bc'
    sample_sheet["I7_Index_ID"] = sample_sheet["I7_Index_ID"].fillna("no_bc")
    sample_sheet["I7_Index_ID"] = sample_sheet["I7_Index_ID"].replace(
        "No_index.*", "no_bc", regex=True
    )

    data = dict()
    data["projects"] = []
    data["samples"] = []
    for index, row in sample_sheet.iterrows():
        assert row["Sample_ID"] not in data.keys()

        # Remove any empty space from the sample project when the PI has 2 last names
        sample_project = row["Sample_Project"].replace(" ", "")

        if row["Sample_Project"] not in data["projects"]:
            # data["projects"].append(row["Sample_Project"])
            data["projects"].append(sample_project)
        if row["Sample_ID"] not in data["samples"]:
            data["samples"].append(row["Sample_ID"])
        data[row["Sample_ID"]] = dict(
            {
                "Sample_Name": row["Sample_Name"],
                # "Sample_Project": row["Sample_Project"],
                "Sample_Project": sample_project,
                "barcode_kits": config["info_dict"]["barcode_kit"],
                "index_id": row["I7_Index_ID"]
                .replace("BP", "barcode")
                .replace("NB", "barcode"),
                "Sample_ID": row["Sample_ID"],
            }
        )
    if len(data["samples"]) > 1 and not config["info_dict"]["barcoding"]:
        print(
            "[red] Danger, barcoding inferred as false, but more then 1 sample in samplesheet. [/red]"
        )
        print("[red] Overwrite barcoding boolean to true. [/red]")
        config["info_dict"]["barcoding"] = True
        print("[red] Defaulting to flowcell kit as barcode kit ! [/red]")
        config["info_dict"]["barcode_kit"] = config["info_dict"]["kit"]
    print(
        "[green] Barcode kit determined as: {}".format(
            config["info_dict"]["barcode_kit"]
        )
    )
    return config["info_dict"]["barcode_kit"], data


def get_periphery(config):
    """
    Return the full pathname to the flowcell in the periphery
    as it should be there after transfer
    Warning1: assumes a _single_ project per flowcell
    Warning2: assumes regular group directory - not applicable to 'external' projects
    """
    group = config["data"]["projects"][0].split("_")[-1].lower()
    group = remove_spaces(group)
    print("Group ", group)
    groupdir = os.path.join(config["paths"]["groupDir"], group)
    # Warning: get_seqdir _creates_ directories as side effect
    groupONT = get_seqdir(groupdir, "sequencing_data")
    fc_base = os.path.basename(config["info_dict"]["flowcell_path"])
    periphery = os.path.join(groupONT, fc_base)
    return periphery


def get_dest_path(config, dir):
    """
    for a given directory, regex: Project_\\d+)_(\\w+)_(\\w+)
    get target destination based on 'pi_name' in 'dir'
    in contrast to get_periphery, this function can handle "unknown" or
    undeterminable PIs and will invoke an external_groupDir in such cases
    """

    # default pi_name if it cannot be determined from basename(dir)
    pi_name = "unknown"
    # try to get pi_name from basename(dir)
    match = re.match(r"(Project_\d+)_([A-Za-z-]+)_([A-Za-z]+)", os.path.basename(dir))

    if match:
        p_id, username, pi_name = match.groups()
        pi_name = pi_name.lower()

    # assume PI has directories all set up
    groupdir = os.path.join(config["paths"]["groupDir"], pi_name)

    if not os.path.exists(groupdir):
        # group pi_name has no volumne so transfer somewhere else (external runs)
        dest_path = os.path.join(
            config["paths"]["external_groupDir"],
            pi_name,
            "sequencing_data/OxfordNanopore",
        )
        print(f"No group directory {groupdir} for PI {pi_name}")
        print(f"Use {dest_path}")
    else:
        # get the most recent destination path: e.g. "{groupdir}/sequencing_data6/OxfordNanopore"
        dest_path = get_seqdir(groupdir, "sequencing_data")

    dest_path = os.path.join(
        dest_path, os.path.basename(config["info_dict"]["flowcell_path"])
    )
    return dest_path


def remove_spaces(input_string):
    # Check if there are empty spaces
    if " " in input_string:
        # Remove empty spaces
        return input_string.replace(" ", "")
    return input_string
