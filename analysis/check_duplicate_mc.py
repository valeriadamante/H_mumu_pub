#!/usr/bin/env python3

import ROOT
import argparse
import sys
import os
from collections import Counter

ROOT.gROOT.SetBatch(True)


def get_lumi_branch(tree):
    for name in ["luminosityBlock", "lumi", "lumiblock"]:
        if tree.GetBranch(name):
            return name
    return None


def open_tree(path, tree_name):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open file: {path}")

    tree = f.Get(tree_name)
    if not tree:
        raise RuntimeError(f"Tree '{tree_name}' not found in file: {path}")

    return f, tree


def build_event_key(entry, branches):
    return tuple(int(getattr(entry, b)) for b in branches)


def main():
    parser = argparse.ArgumentParser(
        description="Check duplicated events inside one ROOT file or between two files."
    )

    parser.add_argument("file1")
    parser.add_argument("file2", nargs="?", default=None)
    parser.add_argument("--tree", default="Events")
    parser.add_argument(
        "--branches",
        nargs="+",
        default=None,
        help="Branches used to define event identity. Default: run lumi event",
    )
    parser.add_argument(
        "--use-full-event-id",
        action="store_true",
        help="Use FullEventId only",
    )
    parser.add_argument(
        "--max-print",
        type=int,
        default=20,
        help="Maximum number of duplicated events to print",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Optional output txt report",
    )

    args = parser.parse_args()

    f1, t1 = open_tree(args.file1, args.tree)
    f2, t2 = open_tree(args.file2, args.tree) if args.file2 else (None, None)

    if args.use_full_event_id:
        branches = ["FullEventId"]
    elif args.branches:
        branches = args.branches
    else:
        lumi = get_lumi_branch(t1)
        if lumi is None:
            raise RuntimeError(
                "Could not find lumi branch. Use --branches or --use-full-event-id."
            )
        branches = ["run", lumi, "event"]

    for b in branches:
        if not t1.GetBranch(b):
            raise RuntimeError(f"Branch '{b}' not found in file1")
        if t2 and not t2.GetBranch(b):
            raise RuntimeError(f"Branch '{b}' not found in file2")

    print("============================================================")
    print("[INFO] Duplicate event checker")
    print(f"[INFO] file1    = {args.file1}")
    if args.file2:
        print(f"[INFO] file2    = {args.file2}")
    print(f"[INFO] tree     = {args.tree}")
    print(f"[INFO] branches = {branches}")
    print(f"[INFO] entries1 = {t1.GetEntries()}")
    if t2:
        print(f"[INFO] entries2 = {t2.GetEntries()}")
    print("============================================================")

    seen_file1 = set()
    internal_dups_file1 = Counter()

    print("[INFO] Reading file1...")

    for i, entry in enumerate(t1):
        key = build_event_key(entry, branches)

        if key in seen_file1:
            internal_dups_file1[key] += 1
        else:
            seen_file1.add(key)

        if i > 0 and i % 500000 == 0:
            print(f"[INFO] file1 processed entries: {i}")

    duplicates_between = []
    internal_dups_file2 = Counter()
    seen_file2 = set()

    if t2:
        print("[INFO] Scanning file2...")
        for i, entry in enumerate(t2):
            key = build_event_key(entry, branches)
            if key in seen_file2:
                internal_dups_file2[key] += 1
            else:
                seen_file2.add(key)
            if key in seen_file1:
                duplicates_between.append(key)
            if i > 0 and i % 500000 == 0:
                print(f"[INFO] file2 processed entries: {i}")

    n_dup_between = len(duplicates_between)
    n_unique_dup_between = len(set(duplicates_between))

    print("\n============================================================")
    print("[RESULT]")
    print(f"Unique events in file1:              {len(seen_file1)}")
    if t2:
        print(f"Unique events in file2:              {len(seen_file2)}")
        print(f"Duplicate entries between files:     {n_dup_between}")
        print(f"Unique duplicated events between:    {n_unique_dup_between}")
    print(f"Internal duplicate events in file1:  {len(internal_dups_file1)}")
    if t2:
        print(f"Internal duplicate events in file2:  {len(internal_dups_file2)}")
    print("============================================================")

    if n_unique_dup_between > 0:
        print(f"\n[INFO] First {min(args.max_print, n_unique_dup_between)} duplicated events:")
        for key in list(dict.fromkeys(duplicates_between))[:args.max_print]:
            print("  " + " ".join(f"{b}={v}" for b, v in zip(branches, key)))

    if args.report:
        with open(args.report, "w") as out:
            out.write("# Duplicate event report\n")
            out.write(f"file1: {args.file1}\n")
            if args.file2:
                out.write(f"file2: {args.file2}\n")
            out.write(f"tree: {args.tree}\n")
            out.write(f"branches: {branches}\n\n")

            out.write(f"unique_file1: {len(seen_file1)}\n")
            if t2:
                out.write(f"unique_file2: {len(seen_file2)}\n")
                out.write(f"duplicate_entries_between: {n_dup_between}\n")
                out.write(f"unique_duplicated_events_between: {n_unique_dup_between}\n")
            out.write(f"internal_duplicate_events_file1: {len(internal_dups_file1)}\n")
            if t2:
                out.write(f"internal_duplicate_events_file2: {len(internal_dups_file2)}\n")
            out.write("\n")

            out.write("# duplicated events between files\n")
            for key in sorted(set(duplicates_between)):
                out.write(" ".join(f"{b}={v}" for b, v in zip(branches, key)) + "\n")

        print(f"\n[INFO] Report written to: {args.report}")

    f1.Close()
    if f2:
        f2.Close()

    if internal_dups_file1 or internal_dups_file2 or n_unique_dup_between > 0:
        sys.exit(2)


if __name__ == "__main__":
    main()
