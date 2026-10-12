#!/usr/bin/env python3
"""Compare run:lumi:event lists without relying on ROOT entry order."""
import argparse
from collections import Counter
import json
from pathlib import Path


def compare(left, right):
    a = Counter(line.strip() for line in left.read_text().splitlines() if line.strip())
    b = Counter(line.strip() for line in right.read_text().splitlines() if line.strip())
    return {'identical': a == b, 'left_events': sum(a.values()), 'right_events': sum(b.values()),
            'left_duplicates': sum(a.values()) - len(a), 'right_duplicates': sum(b.values()) - len(b),
            'only_left': dict(a - b), 'only_right': dict(b - a)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('left', type=Path)
    parser.add_argument('right', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = compare(args.left, args.right)
    text = json.dumps(result, indent=2) + '\n'
    print(text, end='')
    if args.output:
        args.output.write_text(text)
    raise SystemExit(0 if result['identical'] else 1)


if __name__ == '__main__':
    main()
