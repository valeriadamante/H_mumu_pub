#!/usr/bin/env python3
"""Stage 1: validate one dataset's ROOT ntuples and JSON reports."""

import argparse
import json
import multiprocessing
import os
import re
import sys
import time
from pathlib import Path

sys.path.append(os.environ.get("ANALYSIS_PATH", str(Path(__file__).resolve().parents[1])))

import common.utilities as utilities
from common.utilities import (
    write_manifest,
    discover_root_files,
    validate_file,
)


def discover_json_files(path):
    if path.endswith(".json"):
        return [os.path.abspath(path)]
    return sorted(
        os.path.abspath(os.path.join(root, name))
        for root, _, names in os.walk(path)
        for name in names
        if name.endswith(".json")
    )


def validate_json(path, retries=3, retry_delay=2.0):
    last_reason = "unknown validation error"
    for attempt in range(1, retries + 1):
        try:
            with open(path) as handle:
                payload = json.load(handle)
            if not isinstance(payload, dict):
                return path, False, "top-level JSON value is not an object"
            else:
                return path, True, ""
        except json.JSONDecodeError as error:
            # Syntax errors are deterministic; retry only transient I/O errors.
            return path, False, repr(error)
        except Exception as error:
            last_reason = repr(error)
        if attempt < retries and retry_delay:
            time.sleep(retry_delay)
    return path, False, f"{last_reason} (failed after {retries} attempts)"


def pairing_key(path, is_json=False):
    stem = Path(path).stem
    indexed_match = re.fullmatch(
        r"report_(\d+)" if is_json else r"skim_(\d+)",
        stem,
    )
    if indexed_match:
        return indexed_match.group(1)
    if is_json:
        for suffix in ("_skim_report", "_report"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
    elif stem.endswith("_skim"):
        stem = stem[: -len("_skim")]
    return stem


def is_empty_root_result(result):
    return not result[1] and result[2].startswith("empty tree ")


def format_invalid_files(invalid_roots, invalid_jsons):
    """Return actionable diagnostics for the final validation error."""
    details = []
    for label, items in (("ROOT", invalid_roots), ("JSON", invalid_jsons)):
        for item in items:
            details.append(f"{label} {item['path']}: {item['reason']}")
    return details


def pair_results(root_results, json_results):
    roots_by_key = {}
    jsons_by_key = {}
    for result in root_results:
        roots_by_key.setdefault(pairing_key(result[0]), []).append(result)
    for result in json_results:
        jsons_by_key.setdefault(pairing_key(result[0], is_json=True), []).append(result)

    valid_roots = []
    valid_jsons = []
    invalid_roots = []
    invalid_jsons = []
    ignored_empty_roots = []
    for key in sorted(set(roots_by_key) | set(jsons_by_key)):
        root_items = roots_by_key.get(key, [])
        json_items = jsons_by_key.get(key, [])
        good_roots = [item for item in root_items if item[1]]
        empty_roots = [item for item in root_items if is_empty_root_result(item)]
        good_jsons = [item for item in json_items if item[1]]
        pair_is_valid = len(good_roots) == 1 and len(good_jsons) == 1
        empty_pair_is_valid = (
            not good_roots and len(empty_roots) == 1 and len(good_jsons) == 1
        )
        if pair_is_valid:
            root_path = good_roots[0][0]
            json_path = good_jsons[0][0]
            valid_roots.append(root_path)
            valid_jsons.append(json_path)
        elif empty_pair_is_valid:
            empty_path = empty_roots[0][0]
            valid_jsons.append(good_jsons[0][0])
            ignored_empty_roots.append(
                {
                    "path": empty_path,
                    "pairing_key": key,
                    "reason": empty_roots[0][2],
                }
            )
        for path, intrinsically_valid, intrinsic_reason in root_items:
            if (pair_is_valid and intrinsically_valid) or (
                empty_pair_is_valid and is_empty_root_result(
                    (path, intrinsically_valid, intrinsic_reason)
                )
            ):
                continue
            if not intrinsically_valid:
                reason = intrinsic_reason
            elif not json_items:
                reason = "missing JSON counterpart"
            elif not good_jsons:
                reason = "JSON counterpart is invalid"
            else:
                reason = "ambiguous ROOT/JSON pairing"
            invalid_roots.append({"path": path, "pairing_key": key, "reason": reason})

        for path, intrinsically_valid, intrinsic_reason in json_items:
            if (pair_is_valid or empty_pair_is_valid) and intrinsically_valid:
                continue
            if not intrinsically_valid:
                reason = intrinsic_reason
            elif not root_items:
                reason = "missing ROOT counterpart"
            elif not good_roots:
                reason = "ROOT counterpart is invalid"
            else:
                reason = "ambiguous ROOT/JSON pairing"
            invalid_jsons.append({"path": path, "pairing_key": key, "reason": reason})

    return (
        valid_roots,
        valid_jsons,
        invalid_roots,
        invalid_jsons,
        ignored_empty_roots,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--era", required=True)
    parser.add_argument("--dataset-name", required=True)
    parser.add_argument("--root-input", required=True)
    parser.add_argument("--json-input", required=True)
    parser.add_argument("--output-manifest", required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-delay", type=float, default=2.0)
    parser.add_argument("--progress-every", type=int, default=25)
    args = parser.parse_args()
    if args.workers < 1 or args.retries < 1 or args.progress_every < 1:
        parser.error("--workers, --retries and --progress-every must be >= 1")
    if args.retry_delay < 0:
        parser.error("--retry-delay must be >= 0")

    analysis_path = os.environ.get(
        "ANALYSIS_PATH", str(Path(__file__).resolve().parents[1])
    )
    samples_path = os.path.join(
        analysis_path, "config", args.era, "samples.yaml"
    )
    samples = utilities.get_config(samples_path)
    dataset_cfg = samples.get(args.dataset_name, {})
    is_data = dataset_cfg.get("is_data", False) or "data" in args.dataset_name.lower()

    print(f"[DISCOVERY] Scanning ROOT input: {args.root_input}", flush=True)
    root_files = discover_root_files(args.root_input)
    print(f"[DISCOVERY] Scanning JSON input: {args.json_input}", flush=True)
    json_files = discover_json_files(args.json_input)
    print(
        f"[VALIDATION] {args.dataset_name}: discovered "
        f"{len(root_files)} ROOT and {len(json_files)} JSON files",
        flush=True,
    )
    context = multiprocessing.get_context("spawn")
    root_results = []
    root_valid_count = 0
    root_invalid_count = 0
    root_empty_count = 0
    with context.Pool(args.workers) as pool:
        results = pool.imap_unordered(
            validate_file,
            (
                (path, args.tree, args.retries, args.retry_delay)
                for path in root_files
            ),
            chunksize=1,
        )
        for done, result in enumerate(results, start=1):
            root_results.append(result)
            if result[1]:
                root_valid_count += 1
            elif is_empty_root_result(result):
                root_empty_count += 1
                disposition = (
                    "ignored for data"
                    if is_data
                    else "events skipped; JSON retained for normalization"
                )
                print(f"[ROOT EMPTY] {result[0]}: {disposition}", flush=True)
            else:
                root_invalid_count += 1
                print(f"[ROOT INVALID] {result[0]}: {result[2]}", flush=True)
            if done % args.progress_every == 0 or done == len(root_files):
                print(
                    f"[ROOT PROGRESS] {done}/{len(root_files)} checked; "
                    f"valid={root_valid_count}, empty={root_empty_count}, "
                    f"invalid={root_invalid_count}",
                    flush=True,
                )
    json_results = []
    json_valid_count = 0
    json_invalid_count = 0
    for done, path in enumerate(json_files, start=1):
        result = validate_json(path, args.retries, args.retry_delay)
        json_results.append(result)
        if result[1]:
            json_valid_count += 1
        else:
            json_invalid_count += 1
            print(f"[JSON INVALID] {result[0]}: {result[2]}", flush=True)
        if done % args.progress_every == 0 or done == len(json_files):
            print(
                f"[JSON PROGRESS] {done}/{len(json_files)} checked; "
                f"valid={json_valid_count}, invalid={json_invalid_count}",
                flush=True,
            )
    print("[PAIRING] Matching validated ROOT and JSON files", flush=True)
    (
        valid_roots,
        valid_jsons,
        invalid_roots,
        invalid_jsons,
        ignored_empty_roots,
    ) = pair_results(root_results, json_results)
    print(
        f"[PAIRING] accepted={len(valid_roots)} pairs; "
        f"invalid ROOT={len(invalid_roots)}; invalid JSON={len(invalid_jsons)}",
        flush=True,
    )
    failures = []
    if not valid_roots:
        print(
            "[WARNING] No processable ROOT files remain; downstream histogram "
            "production will write empty histograms.",
            flush=True,
        )
    if not root_results and not json_results:
        print(
            "[INFO] Dataset input is empty; writing a passed manifest with "
            "no valid files.",
            flush=True,
        )
    if is_data and (invalid_roots or invalid_jsons):
        failures.append(
            f"data validation requires zero invalid or unmatched files: "
            f"{len(invalid_roots)} ROOT, "
            f"{len(invalid_jsons)} JSON"
        )
        failures.extend(format_invalid_files(invalid_roots, invalid_jsons))
    elif invalid_roots or invalid_jsons:
        print(
            f"[WARNING] MC validation is skipping {len(invalid_roots)} ROOT "
            f"and {len(invalid_jsons)} JSON invalid or unmatched files",
            flush=True,
        )
    validation_passed = not failures
    write_manifest(
        args.output_manifest,
        "validation",
        era=args.era,
        dataset=args.dataset_name,
        is_data=is_data,
        status="passed" if validation_passed else "failed",
        failures=failures,
        root_input=os.path.abspath(args.root_input),
        json_input=os.path.abspath(args.json_input),
        tree=args.tree,
        valid_root_files=valid_roots,
        valid_json_files=valid_jsons,
        invalid_root_files=invalid_roots,
        ignored_empty_root_files=ignored_empty_roots,
        invalid_json_files=invalid_jsons,
        summary={
            "root_valid": len(valid_roots),
            "root_invalid": len(invalid_roots),
            "root_empty_ignored": len(ignored_empty_roots),
            "json_valid": len(valid_jsons),
            "json_invalid": len(invalid_jsons),
            "json_required": True,
        },
    )
    print(f"[VALIDATION] wrote {args.output_manifest}", flush=True)
    if not validation_passed:
        raise RuntimeError(
            f"Validation failed for {args.dataset_name}: " + "; ".join(failures)
        )


if __name__ == "__main__":
    main()
