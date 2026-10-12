def jet_variation_suffixes(df, want_variations, syst_cfg):
    """Include every defined JES p4 variation in the skim's derived columns."""
    suffixes = [""]
    if not want_variations:
        return suffixes
    scales = syst_cfg.get("scales", ["up", "down"])
    for source in ("JER", "JES_Total"):
        template = syst_cfg["systematics"][source]["jet_suffix"]
        suffixes.extend(template.format(scale=scale) for scale in scales)
    suffixes.extend(
        name[len("Jet_p4"):]
        for name in sorted(str(column) for column in df.GetColumnNames())
        if name.startswith("Jet_p4_JES") and any(name.endswith(scale) for scale in scales)
    )
    return list(dict.fromkeys(suffixes))
