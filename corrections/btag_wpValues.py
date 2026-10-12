from .general import pog_folder_names,period_names,periods
import json
import os
import correctionlib

JSON_PATH = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/{}/latest/btagging.json.gz"
ALTERNATIVE_JSON_PATH = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/{}/add_2025_b_and_c_WPs/btagging.json.gz"

def getBTagWPValues(config):
    era = config.get("era")
    period_unc = period_names[era]
    tagger= config.get("bTagAlgo", "particleNet")


    pog_folder = pog_folder_names["BTV"][period_unc]
    json_path = JSON_PATH.format(pog_folder)
    if not os.path.exists(json_path):
        json_path = ALTERNATIVE_JSON_PATH.format(pog_folder)
    cset = correctionlib.CorrectionSet.from_file(json_path)

    wp_values_name = f"{tagger}_wp_values"
    wp_correction = cset[wp_values_name]

    wp_dict = {}

    target_wps = ["L", "M", "T"]

    for wp in target_wps:
        try:
            # Evaluate the correction by passing the WP name string
            threshold = wp_correction.evaluate(wp)
            wp_dict[wp] = threshold
        except Exception as e:
            print(f"Could not evaluate WP '{wp}' for {tagger_name}: {e}")

    # print("Resulting WP Dictionary:")
    # print(wp_dict)
    return wp_dict


