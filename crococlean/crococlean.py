"""Command-line interface for CroCoClean."""

import sys
import argparse
import logging
import multiprocessing
from pathlib import Path
import os
import math
from importlib.metadata import version
from crococlean import ab_table_utils
from crococlean.conta_event import ContaminationEventIO
from crococlean.decontaminate import run_decontamination
from crococlean.self_test import self_test


def set_logging() -> None:
    """Configure logging for the command-line application."""
    logging.basicConfig(
        format="%(asctime)s :: %(levelname)s :: %(message)s", level=logging.INFO
    )


def readable_file(fp_str: str) -> Path:
    """Validate that a path points to a readable regular file."""
    fp = Path(fp_str).resolve()

    if not fp.exists():
        raise argparse.ArgumentTypeError(f"{fp} does not exist")
    if not fp.is_file():
        raise argparse.ArgumentTypeError(f"{fp} is not a regular file.")
    if not os.access(fp, os.R_OK):
        raise argparse.ArgumentTypeError(f"{fp} is not readable.")

    return fp


def writable_file(fp_str: str) -> Path:
    """Validate that a path can be used as an output file."""
    fp = Path(fp_str).resolve()

    if fp.exists():
        if fp.is_dir():
            raise argparse.ArgumentTypeError(f"{fp} is a directory, not a file.")
        if not os.access(fp, os.W_OK):
            raise argparse.ArgumentTypeError(f"{fp} is not writable.")
        return fp

    parent_dir = fp.parent or Path(".")
    if not parent_dir.exists():
        raise argparse.ArgumentTypeError(f"directory {parent_dir} does not exist.")
    if not parent_dir.is_dir():
        raise argparse.ArgumentTypeError(f"{parent_dir} is not a directory.")
    if not os.access(parent_dir, os.W_OK):
        raise argparse.ArgumentTypeError(f"directory {parent_dir} is not writable.")

    return fp


def positive_float(value: str) -> float:
    """Validate a strictly positive floating-point value."""
    try:
        fvalue = float(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(f"{value} is not a number") from value_err

    if not math.isfinite(fvalue) or fvalue <= 0:
        raise argparse.ArgumentTypeError("value must be a finite number greater than 0")

    return fvalue


def nproc(value: str) -> int:
    """Validate the number of requested parallel processes."""
    max_nproc = multiprocessing.cpu_count()

    try:
        ivalue = int(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(f"{value} is not an integer") from value_err

    if ivalue <= 0:
        raise argparse.ArgumentTypeError("minimum value is 1")
    if ivalue > max_nproc:
        raise argparse.ArgumentTypeError(f"maximum value is {max_nproc}")

    return ivalue


def get_arguments() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    prog_name = "crococlean"
    prog_version = version(prog_name.lower())
    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=(
            "Decontaminate species abundance profiles affected "
            "by cross-sample contamination."
        ),
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"{prog_name} version {prog_version}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    # Run subcommand
    run_parser = subparsers.add_parser(
        "run",
        help="Decontaminate species abundance profiles.",
    )

    run_parser.add_argument(
        "-s",
        dest="input_table_fp",
        type=readable_file,
        required=True,
        metavar="SPECIES_ABUNDANCE_TABLE",
        help="Input TSV file corresponding to the species abundance table",
    )

    run_parser.add_argument(
        "-c",
        dest="conta_events_fp",
        type=readable_file,
        required=True,
        metavar="CONTAMINATION_EVENTS_FILE",
        help="Input TSV file created by CroCoDeEL listing contamination events.",
    )

    run_parser.add_argument(
        "-o",
        dest="output_table_fp",
        type=writable_file,
        required=True,
        metavar="OUTPUT_TABLE",
        help="Output TSV file containing the original and decontaminated profiles.",
    )

    run_parser.add_argument(
        "--filter-low-ab",
        dest="filtering_ab_thr_factor",
        type=positive_float,
        required=False,
        default=None,
        metavar="AB_THRESHOLD_FACTOR",
        help=(
            "Filter out low-abundance species that may be inaccurately quantified. "
            "In each sample, set the abundance of species to zero if they are up to "
            "%(metavar)s times more abundant than the least abundant species. "
            "Recommended value for MetaPhlAn4: 20 (default: None)"
        ),
    )

    run_parser.add_argument(
        "--nproc",
        dest="nproc",
        type=nproc,
        default=1,
        help="Number of parallel processes performing decontamination (default: %(default)d)",
    )

    # Test subcommand
    subparsers.add_parser(
        "test",
        help="Test the CroCoClean installation.",
    )

    return parser.parse_args(args=sys.argv[1:] or ["--help"])


def run_crococlean(
    input_table_fp: Path,
    conta_events_fp: Path,
    output_table_fp: Path,
    filtering_ab_thr_factor: float | None,
    nproc: int,
) -> None:
    """Run CroCoClean on the specified input files."""
    with open(input_table_fp, "r", encoding="utf8") as input_table_fh:
        input_table = ab_table_utils.read_filter_normalize(
            input_table_fh,
            filtering_ab_thr_factor,
        )

    with open(conta_events_fp, "r", encoding="utf8") as conta_events_fh:
        conta_events = ContaminationEventIO.read_tsv(conta_events_fh)

    corrected_table = run_decontamination(
        input_table,
        conta_events,
        nproc,
    )

    corrected_table.to_csv(
        output_table_fp,
        sep="\t",
        index=True,
    )

    logging.info(
        "Corrected species abundance table saved in %s",
        output_table_fp,
    )


# pylint: disable=missing-function-docstring
def main() -> None:
    set_logging()
    args = get_arguments()

    if args.command == "test":
        self_test(run_crococlean)
        sys.exit(0)

    run_crococlean(
        args.input_table_fp,
        args.conta_events_fp,
        args.output_table_fp,
        args.filtering_ab_thr_factor,
        args.nproc,
    )


if __name__ == "__main__":
    main()
