"""Resolve jet horn mitigation from the effective skim selection."""


def horn_mitigation_enabled(era, selections):
    if era not in ("Run3_2024", "Run3_2025", "Run3_2026", "2024", "2025", "2026"):
        return True
    expression = "".join(str(selections.get("jet_horn_veto_expr", "false")).split())
    # Both forms are used to disable the veto in configurations/CLI overrides.
    return expression not in ("false", "(false)", "0", "(0)",
                              "(abs(v_ops::eta(Jet_p4))<0)", "")


def configure_horn_veto(selections, era, mode):
    if mode == 'without':
        selections['jet_horn_veto_expr'] = '(abs(v_ops::eta(Jet_p4)) < 0)'
    elif mode == 'with':
        upper = ' && abs(v_ops::eta(Jet_p4)) < 3.0' if str(era).removeprefix('Run3_') in ('2024', '2025', '2026') else ''
        selections['jet_horn_veto_expr'] = (
            '(abs(v_ops::eta(Jet_p4)) >= 2.5' + upper +
            ' && v_ops::pt(Jet_p4) <= 50.0)')
