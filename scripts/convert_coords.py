#!/usr/bin/env python3

from argparse import ArgumentParser
from functools import partial


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "ref", type=int, help="reference coordinate of the coordinate system to convert into")
    parser.add_argument("match", type=int,
                        help="matching coordinate from the coordinate system to convert from")
    parser.add_argument("coordinates", type=str, help="coordinate string to convert")
    args = parser.parse_args()

    get_coord = partial(_get_coord, args.ref, args.match)

    converted = []
    parts = args.coordinates.split(",")
    for part in parts:
        start, end = part.split("..")
        converted.append((get_coord(start)-1, get_coord(end)))

    print('"exons": [')
    for start, end in converted:
        print('{"end": ', end, ', "start":', start, '},')
    print("],")


def _get_coord(ref, match, loc):
    diff = match - int(loc)
    return ref - diff


if __name__ == "__main__":
    main()
