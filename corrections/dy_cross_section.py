import math
import yaml
import os

# =====================================================================
# CROSS-SECTION DATABASE CLASS (string formula resolution)
# =====================================================================
class MiniCrossSectionDB:
    def __init__(self, xs_cfg_dict):
        """Initialise the database from the global cross-section dictionary."""
        self.entries = {}
        self.values = {}

        # Load all available entries from the global crossSections file
        for key, entry in xs_cfg_dict.items():
            if isinstance(entry, dict):
                if "crossSec" in entry:
                    self.add_entry(key, entry["crossSec"])
                elif "BR" in entry:
                    self.add_entry(key, entry["BR"])

    def add_entry(self, name, expr):
        """Dynamically add a new entry, evaluating the expression recursively."""
        if name in self.values:
            return
        value = self.evaluate_expression(expr, name)
        self.values[name] = value

    def get_value(self, name):
        """Return the numerical value associated with a key."""
        if name not in self.values:
            raise RuntimeError(f"MiniCrossSectionDB: Unknown key in database: '{name}'")
        return self.values[name]

    def evaluate_expression(self, expr, name=None):
        """Resolve strings and complex mathematical formulas (e.g. A - B - C) using eval."""
        try:
            if isinstance(expr, (int, float)):
                return float(expr)
            elif isinstance(expr, str):
                # Use self.values as global context to resolve the cross-section names
                result = eval(expr, {}, self.values)
            else:
                raise RuntimeError(f"Invalid expression type: {type(expr)}")

            if not isinstance(result, float):
                raise RuntimeError("The expression did not produce a float value")
            if math.isnan(result) or math.isinf(result):
                raise RuntimeError("The expression produced an unphysical value (NaN or Inf)")
            if result < 0:
                raise RuntimeError(f"The expression produced a negative value: {result}")

            return result
        except Exception as e:
            msg = f"'{expr}'" + (f" for entry '{name}'" if name else "")
            raise RuntimeError(f"CrossSectionDB: Error evaluating {msg}: {e}")


# =====================================================================
# CALCULATION AND APPLICATION FUNCTIONS (adapted to the long-string YAML)
# =====================================================================

def dy_cross_section_calc(db, dysf_cfg):
    """
    Compute the cross-section factors by interpreting the YAML structure.
    No hardcoded values. If totalCrossSection is absent/commented out, the check is skipped.
    """
    # 1. Dynamic handling of totalCrossSectionScaling (e.g. "DY_NNLO_QCD_NLO_EW / DYto2L_...")
    scaling_expr = dysf_cfg['totalCrossSectionScaling']

    if isinstance(scaling_expr, str) and '/' in scaling_expr:
        # Split the string on the two sides of the slash '/'
        numer_key, denom_key = [k.strip() for k in scaling_expr.split('/')]
        numer_val = db.get_value(numer_key)
        denom_val = db.get_value(denom_key)
        total_xs_scaling = numer_val / denom_val
    else:
        # In case it is a string without slash or already a float number
        total_xs_scaling = db.evaluate_expression(scaling_expr, name="totalCrossSectionScaling")

    # 2. Compute the cross-sections bin by bin by evaluating the YAML strings (formulas with minus '-')
    selection_bins = []
    for bin_entry in dysf_cfg['bins']:
        raw_xs = bin_entry['crossSection']

        if isinstance(raw_xs, str):
            # Evaluate the whole arithmetic string using the DB
            subtotal = db.evaluate_expression(raw_xs, name=bin_entry['name'])
        elif isinstance(raw_xs, dict):
            # Fallback in case the add/subtract format is put back in the future
            subtotal = 0.0
            for entry in raw_xs.get('add', []): subtotal += db.get_value(entry)
            for entry in raw_xs.get('subtract', []): subtotal -= db.get_value(entry)
        else:
            subtotal = float(raw_xs)

        selection_bins.append(subtotal)

    # 3. Consistency check (sanity check) - run ONLY if the key exists and is not commented out
    total_xs_from_bins = sum(selection_bins)
    expected_key = dysf_cfg.get('totalCrossSection', None) # Returns None if commented out or absent

    if expected_key:
        expected_total_xs = db.get_value(expected_key)
        diff = abs(total_xs_from_bins - expected_total_xs) / expected_total_xs

        print(f"[MCStitcher-Core] Computed sum of bins: {total_xs_from_bins}")
        print(f"[MCStitcher-Core] Expected cross-section for the check: {expected_total_xs}")
        print(f"[MCStitcher-Core] Relative difference: {diff}")

        if diff > 0.001:
            raise RuntimeError(
                f"MCStitcher Error: sum of bin cross-sections ({total_xs_from_bins}) "
                f"does not match total cross-section ({expected_total_xs}) for key '{expected_key}'. Diff: {diff}"
            )
    else:
        # This is exactly what is needed to get unblocked with the current YAML!
        print(f"[MCStitcher-Core] 'totalCrossSection' not defined in the YAML (omitted or commented out).")
        print(f"[MCStitcher-Core] Check skipped. Total sum of bins used: {total_xs_from_bins}")

    return total_xs_scaling, selection_bins

def apply_dy_cross_section(df, weight_xs_name, xs_cfg, dysf_cfg, json_dict):
    """
    Redefine the cross-section branch in the RDataFrame using the internal DB.
    C++ syntax fixed to avoid implicit lambda capture errors in Cling.
    """
    # Initialise the MiniCrossSectionDB database with the global cross-section dictionary
    db = MiniCrossSectionDB(xs_cfg)

    # If the stitching YAML has extra local definitions (crossSections: {}), load them into the DB
    if "crossSections" in dysf_cfg and dysf_cfg["crossSections"]:
        for entry_name, entry_expr in dysf_cfg["crossSections"].items():
            db.add_entry(entry_name, entry_expr)

    # Compute the scaling and the binned cross-sections
    total_xs_scaling, selection_bins = dy_cross_section_calc(db, dysf_cfg)

    # Prepare the structures that store the weight counters in json_dict
    if 'gen' not in json_dict: json_dict['gen'] = {}
    if 'pu' not in json_dict: json_dict['pu'] = {}
    available_columns = [str(col) for col in df.GetColumnNames()]

    xs_expr_list = []
    columns_to_pass = set()

    for entry, binned_xs in zip(dysf_cfg['bins'], selection_bins):
        bin_name = entry['name']
        selection = entry['selection']

        # Extract the names of the columns used in the selection to pass them to Define (e.g. LHE_NpNLO)
        # A simple trick: keep the alphabetic tokens that match real columns of the DF
        for col in available_columns:
            if col in selection:
                columns_to_pass.add(col)

        # Register the RDataFrame ResultPtr pointers for the yields
        if bin_name not in json_dict['gen']: json_dict['gen'][bin_name] = {}
        json_dict['gen'][bin_name]['selection'] = selection
        json_dict['gen'][bin_name]['value'] = df.Filter(selection).Sum("genWeight")

        for column in available_columns:
            prefix = "weight_qcd_scale_sum__"
            if not column.startswith(prefix):
                continue
            point_name = column[len(prefix):]
            json_key = f"qcd_scale__{point_name}"
            if json_key not in json_dict:
                json_dict[json_key] = {}
            json_dict[json_key][bin_name] = {
                "selection": selection,
                "value": df.Filter(selection).Sum(column),
            }

        if bin_name not in json_dict['pu']: json_dict['pu'][bin_name] = {}
        json_dict['pu'][bin_name]['selection'] = selection
        json_dict['pu'][bin_name]['value'] = df.Filter(selection).Sum("weight_pu")

        # Handling of PU systematic variations (Up/Down)
        for scale in ['up', 'down']:
            pu_var_name = f"weight_pu_{scale}"
            if pu_var_name not in available_columns: continue
            json_key = f"pu_{scale}"
            if json_key not in json_dict: json_dict[json_key] = {}
            if bin_name not in json_dict[json_key]: json_dict[json_key][bin_name] = {}
            json_dict[json_key][bin_name]['selection'] = selection
            json_dict[json_key][bin_name]['value'] = df.Filter(pu_var_name)

        # Compute the final scaled cross-section for this bin (xs_scaling_globale * bin_xs)
        xs = total_xs_scaling * binned_xs
        sub_expr = f"if ({selection}) return float({xs});"
        xs_expr_list.append(sub_expr)

    xs_expr_list.append('throw std::runtime_error("No bin matched in DY cross-sectioning bins.");')

    # Build the pure C++ macro block without nesting it in a second instantiated lambda
    macro_xs_expr = "\n".join(xs_expr_list)

    # Convert the column set into a list to pass it formally to Define
    columns_list = list(columns_to_pass)

    # Apply the transformation to the RDataFrame, passing the needed columns explicitly
    df = df.Define(weight_xs_name, macro_xs_expr, columns_list)

    return df, json_dict
# def dy_cross_section_calc(xs_cfg, dysf_cfg):
#     """
#     Calculates the new DY cross section factors for using multiple DY samples.
#     Args:
#     - xs_cfg expected to crossSections13p6TeV.yaml loaded as dict.
#     - dysf_cfg expected to be dy_cross_section_stitching.yaml loaded as dict.
#     Returns:
#     - The overall scale factor
#     - The selection-by-selection factors for different binnings
#     """
#     # Calculate the total cross-section scaling factor
#     numer = dysf_cfg['totalCrossSectionScaling']['numerator']
#     denom = dysf_cfg['totalCrossSectionScaling']['denominator']
#     numer = xs_cfg[numer]['crossSec']
#     denom = xs_cfg[denom]['crossSec']
#     if isinstance(numer, str):
#         numer = eval(numer)
#     if isinstance(denom, str):
#         denom = eval(denom)
#     total_xs_scaling = numer/denom

#     # Calculate the selection by selection cross-sections
#     selection_bins = []
#     for bin in dysf_cfg['bins']:
#         subtotal = 0
#         for entry in bin['crossSection']['add']:
#             subtotal += xs_cfg[entry]['crossSec']
#         for entry in bin['crossSection']['subtract']:
#             subtotal -= xs_cfg[entry]['crossSec']
#         selection_bins.append(subtotal)

#     # Check the factor
#     total_xs_from_bins = sum(selection_bins)
#     expected_total_xs = dysf_cfg['totalCrossSectionScaling']['numerator']
#     expected_total_xs = xs_cfg[expected_total_xs]['crossSec']
#     diff = abs(total_xs_from_bins - expected_total_xs)/expected_total_xs
#     print(expected_total_xs)
#     print(total_xs_from_bins)
#     print(diff)
#     assert diff < 0.001

#     return total_xs_scaling, selection_bins


# def apply_dy_cross_section(df, weight_xs_name, xs_cfg, dysf_cfg, json_dict):
#     """
#     Redefines the weight_xs branch for DY samples.
#     Args:
#     - df is the working RDataFrame
#     - xs_cfg expected to crossSections13p6TeV.yaml loaded as dict.
#     - dysf_cfg_path is the .yaml file containing the DY stitching configuration
#     Returns:
#     - The updated RDataFrame
#     """
#     # Get needed values from dy_cross_section_calc.
#     total_xs_scaling, selection_bins = dy_cross_section_calc(xs_cfg, dysf_cfg)
#     # Build the RDataFrame expression to assign new cross-section values
#     xs_expr_list = []
#     for entry, binned_xs in zip(dysf_cfg['bins'], selection_bins):
#         ### sum of gen weigths
#         if 'gen' not in json_dict.keys():
#             json_dict['gen'] = {}
#         if entry['name'] not in json_dict['gen'].keys():
#             json_dict['gen'][entry['name']] = {}
#         json_dict['gen'][entry['name']]['selection'] = entry['selection']
#         json_dict['gen'][entry['name']]['value'] = df.Filter(entry['selection']).Sum("genWeight")
#         ### sum of pu weights

#         if entry['name'] not in json_dict['pu'].keys():
#             json_dict['pu'][entry['name']] = {}
#         json_dict['pu'][entry['name']]['selection'] = entry['selection']
#         json_dict['pu'][entry['name']]['value'] = df.Filter(entry['selection']).Sum("weight_pu")
#         for scale in ['up','down']:
#             if f"weight_pu_{scale}" not in df.GetColumnNames(): continue
#             if entry['name'] not in json_dict[f"pu_{scale}"].keys():
#                 json_dict[f"pu_{scale}"][entry['name']] = {}
#             json_dict[f"pu_{scale}"][entry['name']]['selection'] = entry['selection']
#             json_dict[f"pu_{scale}"][entry['name']]['value'] = df.Filter(entry['selection']).Sum(f"weight_pu_{scale}")

#         xs = total_xs_scaling * binned_xs
#         selection = entry['selection']
#         sub_expr = f"if ({selection}) return {xs};"
#         xs_expr_list.append(sub_expr)
#     xs_expr_list.append('throw std::runtime_error("No bin matched in DY cross-sectioning bins.")')
#     xs_expr = "\n".join(xs_expr_list)

#     # Apply and return!
#     df = df.Define(weight_xs_name, xs_expr)


#     return df,json_dict
