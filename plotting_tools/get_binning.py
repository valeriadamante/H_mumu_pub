if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', required=True, type=float)
    parser.add_argument('--end', required=True, type=float)
    parser.add_argument('--step', required=True, type=float)

    args = parser.parse_args()

    if args.step == 0:
        raise ValueError("Step can't be 0.")

    list = []
    value = args.start
    if args.step > 0:
        while value <= args.end:
            list.append(value)
            value += args.step
    else:
        while value >= args.end:
            list.append(value)
            value += args.step
    rounded_list = [round(num, 6) for num in list]
    print(rounded_list)
