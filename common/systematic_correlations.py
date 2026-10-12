"""One nuisance naming policy for production, aggregation and datacards."""
import re


def short_era(era):
    return str(era).removeprefix('Run3_')


def correlation_tag(info, era):
    """None preserves legacy naming; an empty tag means fully correlated."""
    if 'isCorrelated' not in info:
        if 'Eras_correlated' in info:
            raise ValueError('Eras_correlated requires isCorrelated')
        return None
    enabled = info['isCorrelated']
    if not isinstance(enabled, bool):
        raise ValueError('isCorrelated must be a YAML boolean')
    era = short_era(era)
    if not enabled:
        return era
    groups = info.get('Eras_correlated', 'all')
    if groups == 'all':
        return ''
    if not isinstance(groups, list) or not groups:
        raise ValueError('Eras_correlated must be all, a nonempty era list, or a list of groups')
    groups = groups if all(isinstance(g, list) for g in groups) else [groups]
    seen = set()
    selected = era
    for group in groups:
        if not group or any(isinstance(e, (list, dict, bool)) for e in group):
            raise ValueError('Invalid correlation group')
        names = [short_era(e) for e in group]
        if any(not re.fullmatch(r'20\d\d(?:EE|BPix)?', e) for e in names):
            raise ValueError('Correlation groups must contain physical Run3 eras')
        if len(set(names)) != len(names) or seen.intersection(names):
            raise ValueError('Correlation groups must be disjoint, without duplicate eras')
        seen.update(names)
        if era in names:
            selected = '_'.join(sorted(names))
    return selected


def nuisance_name(info, era, **labels):
    template = info['name']
    tag = correlation_tag(info, era)
    if tag is None:
        return template.format(era=short_era(era), **labels)
    if '{era}' in template:
        if tag:
            template = template.replace('{era}', tag)
        else:
            template = template.replace('_{era}', '').replace('{era}_', '').replace('{era}', '')
    elif tag:
        template += '_' + tag
    return template.format(**labels)


def configured_nuisances(cfg, era):
    """Yield (configuration key, resolved name), including JER and process labels."""
    for section in ('systematics', 'weights', 'derived_systematics'):
        for key, info in cfg.get(section, {}).items():
            if key == 'Central' or not info.get('name') or info.get('derived_envelope'):
                continue
            entries = [dict(info, name=f'{component}{{era}}') for component in info.get('components', [])] or [info]
            for entry in entries:
                processes = cfg.get('pdf', {}).get('process_labels', {}).values() if '{pdf_process}' in entry['name'] else []
                processes = set(processes) or {''}
                for process in processes:
                    if '{process}' in entry['name']:
                        for label in set(cfg.get('qcd_scale', {}).get('process_labels', {}).values()):
                            yield key, nuisance_name(entry, era, process=label, scale='', pdf_process=process)
                    else:
                        yield key, nuisance_name(entry, era, scale='', pdf_process=process)
    qcd = cfg.get('qcd_scale', {})
    if qcd.get('enabled'):
        for variation in qcd.get('variations', []):
            for label in set(qcd.get('process_labels', {}).values()):
                yield 'QCDScale', nuisance_name(dict(qcd, **variation), era, process=label)
