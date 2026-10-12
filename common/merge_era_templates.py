"""Complete era-specific shapes with the nominal from unaffected eras."""

def histogram_keys(directory, prefix=''):
    import ROOT
    result = set()
    for key in directory.GetListOfKeys():
        cls = ROOT.TClass.GetClass(key.GetClassName())
        name = prefix + key.GetName()
        if cls and cls.InheritsFrom('TDirectory'):
            result.update(histogram_keys(directory.Get(key.GetName()), name + '/'))
        elif cls and cls.InheritsFrom('TH1'):
            result.add(name)
    return result


def nominal_for(key, nominal_keys):
    """Longest nominal key k in the same directory with key == k + '_' + suffix.

    The same answer as scanning every nominal key, found by cutting the name at
    each '_' from the right: files with ~10^5 histograms make the scan O(N^2)."""
    parent, _, name = key.rpartition('/')
    head = parent + '/' if parent else ''
    cut = len(name)
    while True:
        cut = name.rfind('_', 0, cut)
        if cut <= 0:
            return None
        candidate = head + name[:cut]
        if candidate in nominal_keys:
            return candidate


def complete_era_shapes(output, inputs, nominal_inputs=None):
    """After hadd, pad each missing shifted key with that era's nominal.

    Inputs contain one process per era; absent processes contribute zero.
    nominal_inputs supplies Central files for variations-only productions.
    """
    import ROOT
    nominal_inputs = nominal_inputs or inputs
    if len(inputs) != len(nominal_inputs):
        raise ValueError('One nominal file is required per era input')
    opened = []
    def open_root(path, mode='READ'):
        handle = ROOT.TFile.Open(str(path), mode)
        if not handle or handle.IsZombie():
            raise RuntimeError(f'Cannot open ROOT file: {path}')
        opened.append(handle)
        return handle
    try:
        sources = [open_root(p) for p in inputs]
        source_keys = [histogram_keys(f) for f in sources]
        varied = set().union(*source_keys)
        varied = {k for k in varied if k.endswith(('Up', 'Down'))}
        if not varied:
            return
        nominals = [open_root(p) for p in nominal_inputs]
        nominal_keys = set().union(*(histogram_keys(f) for f in nominals))
        nominal_keys = {k for k in nominal_keys if not k.endswith(('Up', 'Down'))}
        target = open_root(output, 'UPDATE')
        for key in sorted(varied):
            missing = [i for i, keys in enumerate(source_keys) if key not in keys]
            if not missing:
                continue
            nominal_key = nominal_for(key, nominal_keys)
            if nominal_key is None:
                raise RuntimeError(f'No nominal histogram for {key}; supply --nominal-dir')
            hist = target.Get(key).Clone()
            hist.SetDirectory(0)
            for i in missing:
                nominal = nominals[i].Get(nominal_key)
                if not nominal:
                    # Region/category/variable selections may differ by era.
                    # An observable absent also from Central contributes zero.
                    continue
                if hist.GetNcells() != nominal.GetNcells() or hist.GetDimension() != nominal.GetDimension():
                    raise RuntimeError(f'Incompatible binning for {key}')
                for axis in ('GetXaxis', 'GetYaxis', 'GetZaxis')[:hist.GetDimension()]:
                    left, right = getattr(hist, axis)(), getattr(nominal, axis)()
                    if left.GetNbins() != right.GetNbins() or any(
                        left.GetBinLowEdge(b) != right.GetBinLowEdge(b)
                        for b in range(1, left.GetNbins() + 2)
                    ):
                        raise RuntimeError(f'Incompatible bin edges for {key}')
                if not hist.Add(nominal):
                    raise RuntimeError(f'Cannot add nominal for {key}')
            parent, _, name = key.rpartition('/')
            directory = target.GetDirectory(parent) if parent else target
            directory.WriteTObject(hist, name, 'Overwrite')
    finally:
        for handle in reversed(opened):
            handle.Close()
