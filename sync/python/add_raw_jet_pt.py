#!/usr/bin/env python3
"""Add the raw (uncorrected for JEC and JER) jet pT to a sync tuple.

skim_v4 does not store Jet_rawFactor, so it is read back from the NanoAOD files
of the datasets (config/Run3_2024/samples_withfiles.yaml) for the tuple events.
Events are matched on run:luminosityBlock:event and, since MC extensions reuse
the same numbers, on the jet collection itself: the NanoAOD Jet_pt must equal
the tuple's Jet_pt_nocorr jet by jet. Any tuple event without exactly one
matching NanoAOD event stops the script.

Added branches (MISSING = -1000 where the jet does not exist):

  Jet_rawFactor, Jet_pt_raw      all jets, NanoAOD order (= skim order)
  SelectedJet_pt_raw             selected jets, in the SelectedJet_* order
  leadingjet_pt_raw, subleadingjet_pt_raw, thirdleadingjet_pt_raw
  vbfjet1_pt_raw, vbfjet2_pt_raw

with Jet_pt_raw = Jet_pt_nocorr * (1 - Jet_rawFactor), the raw pT from which
corrections/jets.h reapplies the JEC. The NanoAOD reads are cached per file in
--cache-dir, so an interrupted run resumes.

  python3 sync/python/add_raw_jet_pt.py TUPLE.root --dataset Muon0_Run2024H --dataset Muon1_Run2024H
"""
import argparse
import hashlib
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import yaml

REPOSITORY = Path(__file__).resolve().parents[2]
ERA = "Run3_2024"
MISSING = -1000.0
KEYS = ["run", "luminosityBlock", "event"]

# Events of the tuple, set before the worker processes fork.
_TUPLE_EVENTS = None
_TUPLE_KEYS = None


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tuple", type=Path, help="sync tuple from prepare_tuple_forSync.py (updated in place)")
    parser.add_argument("--dataset", action="append", required=True,
                        help="dataset of samples_withfiles.yaml whose NanoAOD holds the events (repeatable)")
    parser.add_argument("--cache-dir", type=Path, help="per-NanoAOD-file cache (default: <tuple stem>_nano_cache)")
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--tree", default="Events")
    return parser.parse_args()


def nano_files(datasets):
    catalogue = yaml.safe_load((REPOSITORY / "config" / ERA / "samples_withfiles.yaml").read_text())
    files = []
    for dataset in datasets:
        if dataset not in catalogue:
            raise KeyError(f"{dataset} is not in config/{ERA}/samples_withfiles.yaml")
        files += catalogue[dataset]["filelist"]
    return list(dict.fromkeys(files))


def scan_nano(url, cache):
    """Jets of the NanoAOD events of one file that are in the tuple."""
    if cache.is_file():
        return cache
    import awkward as ak
    import uproot
    if url.startswith("/eos/"):
        url = "root://eoscms.cern.ch/" + url
    parts = []
    for attempt in range(3):
        try:
            with uproot.open(url, handler=uproot.source.xrootd.XRootDSource, timeout=600) as handle:
                tree = handle["Events"]
                for chunk in tree.iterate(KEYS + ["Jet_pt", "Jet_rawFactor"], step_size="200 MB"):
                    events = ak.to_numpy(chunk["event"]).astype(np.uint64)
                    pos = np.clip(np.searchsorted(_TUPLE_EVENTS, events), 0, len(_TUPLE_EVENTS) - 1)
                    candidate = np.flatnonzero(_TUPLE_EVENTS[pos] == events)
                    if not len(candidate):
                        continue
                    run = ak.to_numpy(chunk["run"])[candidate]
                    lumi = ak.to_numpy(chunk["luminosityBlock"])[candidate]
                    keep = candidate[[(int(r), int(l), int(e)) in _TUPLE_KEYS
                                      for r, l, e in zip(run, lumi, events[candidate])]]
                    if len(keep):
                        parts.append(chunk[keep])
            break
        except OSError as error:  # transient XRootD/EOS errors
            parts = []
            if attempt == 2:
                raise RuntimeError(f"{url}: {error}") from error
    found = ak.concatenate(parts) if parts else None
    tmp = cache.with_suffix(".tmp.npz")
    if found is None or len(found) == 0:
        np.savez(tmp, run=np.zeros(0, np.uint32), luminosityBlock=np.zeros(0, np.uint32),
                 event=np.zeros(0, np.uint64), counts=np.zeros(0, np.int64),
                 Jet_pt=np.zeros(0, np.float32), Jet_rawFactor=np.zeros(0, np.float32))
    else:
        np.savez(tmp, run=ak.to_numpy(found["run"]), luminosityBlock=ak.to_numpy(found["luminosityBlock"]),
                 event=ak.to_numpy(found["event"]).astype(np.uint64),
                 counts=ak.to_numpy(ak.num(found["Jet_pt"])),
                 Jet_pt=ak.to_numpy(ak.flatten(found["Jet_pt"])),
                 Jet_rawFactor=ak.to_numpy(ak.flatten(found["Jet_rawFactor"])))
    os.replace(tmp, cache)
    return cache


def raw_factors(tuple_arrays, caches):
    """Jet_rawFactor of every tuple event, in tuple order (awkward array)."""
    import awkward as ak
    candidates = {}
    for cache in caches:
        with np.load(cache) as handle:
            stored = {name: handle[name] for name in handle.files}
        offsets = np.concatenate([[0], np.cumsum(stored["counts"])])
        for i, key in enumerate(zip(stored["run"].tolist(), stored["luminosityBlock"].tolist(),
                                    stored["event"].tolist())):
            jets = slice(offsets[i], offsets[i + 1])
            candidates.setdefault(key, []).append((stored["Jet_pt"][jets], stored["Jet_rawFactor"][jets]))
    factors, problems = [], {"not_found": 0, "ambiguous": 0}
    for run, lumi, event, pt_nocorr in zip(*(tuple_arrays[k] for k in KEYS), tuple_arrays["Jet_pt_nocorr"]):
        pt_nocorr = np.asarray(pt_nocorr, dtype=np.float32)
        same = [raw for pt, raw in candidates.get((int(run), int(lumi), int(event)), [])
                if len(pt) == len(pt_nocorr) and np.array_equal(pt, pt_nocorr)]
        if len(same) != 1:
            problems["not_found" if not same else "ambiguous"] += 1
            factors.append(np.full(len(pt_nocorr), np.nan, np.float32))
            continue
        factors.append(same[0])
    if problems["not_found"] or problems["ambiguous"]:
        raise RuntimeError(f"{problems['not_found']} tuple events without a matching NanoAOD event, "
                           f"{problems['ambiguous']} with more than one: check the --dataset list")
    return ak.Array(factors)


def write_raw_columns(tuple_path, tree, raw_factor, event):
    """Rewrite the tuple with the raw-pT branches, keeping its metadata."""
    import awkward as ak
    import uproot
    import ROOT
    friend_path = tuple_path.with_name(f"{tuple_path.stem}_rawfactor_friend.root")
    with uproot.recreate(friend_path) as handle:
        handle["RawFactor"] = {"RawJet": ak.zip({"rawFactor": ak.values_astype(raw_factor, np.float32)}),
                               "RawJet_event": np.asarray(event, dtype=np.uint64)}

    ROOT.gROOT.SetBatch(True)
    source = ROOT.TFile.Open(str(tuple_path))
    metadata = json.loads(source.Get("metadata").GetTitle()) if source.Get("metadata") else None
    events = source.Get(tree)
    # The multithreaded Snapshot marks the tuple as reshuffled. The friend was
    # written in the tuple's stored entry order, so the entries are aligned.
    events.ResetBit(ROOT.TTree.kEntriesReshuffled)
    events.AddFriend("RawFactor", str(friend_path))
    df = ROOT.RDataFrame(events)
    misaligned = df.Filter("static_cast<ULong64_t>(event) != RawJet_event").Count().GetValue()
    if misaligned:
        raise RuntimeError(f"{misaligned} entries of the raw-factor friend are not aligned with the tuple")
    m = f"{MISSING}f"
    df = (df.Define("Jet_rawFactor", "ROOT::VecOps::RVec<float>(RawJet_rawFactor)")
            .Define("Jet_pt_raw", "ROOT::VecOps::RVec<float>(Jet_pt_nocorr * (1.f - Jet_rawFactor))")
            .Define("SelectedJet_pt_raw", "ROOT::VecOps::RVec<float>(ROOT::VecOps::Take(Jet_pt_raw, SelectedJet_idx))"))
    for index, name in enumerate(("leadingjet", "subleadingjet", "thirdleadingjet")):
        df = df.Define(f"{name}_pt_raw", f"SelectedJet_pt_raw.size() > {index} ? SelectedJet_pt_raw[{index}] : {m}")
    for leg in (1, 2):
        df = df.Define(f"vbfjet{leg}_pt_raw",
                       f"HasVBF && VBFJetIdx_{leg} >= 0 && static_cast<size_t>(VBFJetIdx_{leg}) < "
                       f"SelectedJet_pt_raw.size() ? SelectedJet_pt_raw[VBFJetIdx_{leg}] : {m}")
    added = ["Jet_rawFactor", "Jet_pt_raw", "SelectedJet_pt_raw", "leadingjet_pt_raw", "subleadingjet_pt_raw",
             "thirdleadingjet_pt_raw", "vbfjet1_pt_raw", "vbfjet2_pt_raw"]
    original = [str(b.GetName()) for b in events.GetListOfBranches() if str(b.GetName()) not in added]
    output = tuple_path.with_name(f"{tuple_path.stem}.rawjet.tmp.root")
    options = ROOT.RDF.RSnapshotOptions()
    options.fMode = "RECREATE"
    df.Snapshot(tree, str(output), original + added, options)
    source.Close()
    if metadata is not None:
        metadata["branches"] = sorted(set(metadata.get("branches", [])) | set(added))
        metadata["conventions"]["jets"]["Jet_pt_raw, SelectedJet_pt_raw, *_pt_raw"] = (
            "raw jet pT (no JEC, no JER) = Jet_pt_nocorr * (1 - Jet_rawFactor), "
            "Jet_rawFactor read from the NanoAOD (sync/python/add_raw_jet_pt.py)")
        text = json.dumps(metadata, indent=2, default=str)
        handle = ROOT.TFile.Open(str(output), "UPDATE")
        ROOT.TNamed("metadata", text).Write()
        handle.Close()
        tuple_path.with_name(f"{tuple_path.stem}_metadata.json").write_text(text)
    os.replace(output, tuple_path)
    friend_path.unlink()
    print(f"[raw-jet] added {', '.join(added)} to {tuple_path}")


def main():
    global _TUPLE_EVENTS, _TUPLE_KEYS
    args = arguments()
    import uproot
    if not args.tuple.is_file():
        raise FileNotFoundError(args.tuple)
    with uproot.open(args.tuple) as handle:
        arrays = handle[args.tree].arrays(KEYS + ["Jet_pt_nocorr"], library="np")
    _TUPLE_EVENTS = np.unique(arrays["event"].astype(np.uint64))
    _TUPLE_KEYS = set(zip(arrays["run"].tolist(), arrays["luminosityBlock"].tolist(), arrays["event"].tolist()))
    files = nano_files(args.dataset)
    cache_dir = args.cache_dir or args.tuple.with_name(f"{args.tuple.stem}_nano_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    caches = {url: cache_dir / f"{hashlib.sha1(url.encode()).hexdigest()[:16]}.npz" for url in files}
    todo = [url for url, cache in caches.items() if not cache.is_file()]
    print(f"[raw-jet] {len(arrays['event'])} tuple events, {len(files)} NanoAOD files "
          f"({len(files) - len(todo)} already cached)")
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(scan_nano, url, caches[url]): url for url in todo}
        for done, future in enumerate(as_completed(futures), 1):
            future.result()
            print(f"[raw-jet] {done}/{len(todo)} NanoAOD files read", flush=True)
    write_raw_columns(args.tuple, args.tree, raw_factors(arrays, list(caches.values())), arrays["event"])


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
