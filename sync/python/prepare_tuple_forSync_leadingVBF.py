#!/usr/bin/env python3
"""Sync tuple maker using the pT-leading VBF pair passing mjj/deta cuts."""
import sys

import prepare_tuple_forSync as sync


if __name__ == "__main__":
    if "--vbf-pair-strategy" not in sys.argv:
        sys.argv.extend(["--vbf-pair-strategy", "leading"])
    try:
        sync.main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
