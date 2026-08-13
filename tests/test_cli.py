"""Unit tests for the CroCoClean command-line interface."""

import argparse
from pathlib import Path

import pytest

from crococlean.crococlean import (
    get_arguments,
    positive_float,
    nproc,
    readable_file,
    writable_file,
)


def test_readable_file(tmp_path):
    """Test that a readable file is accepted."""
    filepath = tmp_path / "input.tsv"
    filepath.write_text("test\n")

    result = readable_file(str(filepath))

    assert result == filepath.resolve()


def test_readable_file_does_not_exist(tmp_path):
    """Test that a missing file is rejected."""
    filepath = tmp_path / "missing.tsv"

    with pytest.raises(argparse.ArgumentTypeError, match="does not exist"):
        readable_file(str(filepath))


def test_readable_file_directory(tmp_path):
    """Test that a directory is rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="not a regular file",
    ):
        readable_file(str(tmp_path))


def test_writable_existing_file(tmp_path):
    """Test that an existing writable file is accepted."""
    filepath = tmp_path / "output.tsv"
    filepath.write_text("test\n")

    result = writable_file(str(filepath))

    assert result == filepath.resolve()


def test_writable_file_directory(tmp_path):
    """Test that a directory is rejected as an output file."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="is a directory",
    ):
        writable_file(str(tmp_path))


def test_writable_file_missing_parent(tmp_path):
    """Test that a missing parent directory is rejected."""
    filepath = tmp_path / "missing" / "output.tsv"

    with pytest.raises(
        argparse.ArgumentTypeError,
        match="does not exist",
    ):
        writable_file(str(filepath))


def test_writable_file_new_file(tmp_path):
    """Test that a new file in a writable directory is accepted."""
    filepath = tmp_path / "output.tsv"

    result = writable_file(str(filepath))

    assert result == filepath.resolve()


@pytest.mark.parametrize(
    "value, expected",
    [
        ("1", 1.0),
        ("0.1", 0.1),
        ("1.5", 1.5),
        ("1e-3", 0.001),
    ],
)
def test_positive_float_valid(value, expected):
    """Test that positive floating-point values are accepted."""
    assert positive_float(value) == pytest.approx(expected)


@pytest.mark.parametrize(
    "value",
    ["0", "-1", "-0.1"],
)
def test_positive_float_non_positive(value):
    """Test that non-positive values are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="finite number greater than 0",
    ):
        positive_float(value)


@pytest.mark.parametrize(
    "value",
    [
        "abc",
        "1.2.3",
        "",
    ],
)
def test_positive_float_invalid(value):
    """Test that non-numeric values are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match=f"{value} is not a number",
    ):
        positive_float(value)


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
def test_positive_float_special_values(value):
    """Test that non-finite floating-point values are rejected."""
    with pytest.raises(argparse.ArgumentTypeError):
        positive_float(value)


@pytest.mark.parametrize(
    "value",
    ["0", "-1", "abc"],
)
def test_nproc_invalid(value):
    """Test that invalid process counts are rejected."""
    with pytest.raises(argparse.ArgumentTypeError):
        nproc(value)


def test_nproc_too_many():
    """Test that more processes than available CPUs are rejected."""
    with pytest.raises(argparse.ArgumentTypeError):
        nproc("999999")


def test_nproc_valid():
    """Test that a valid process count is accepted."""
    assert nproc("1") == 1


def test_get_arguments(monkeypatch, tmp_path):
    """Test parsing the command-line arguments."""
    input_file = tmp_path / "input.tsv"
    conta_file = tmp_path / "contamination.tsv"
    output_file = tmp_path / "output.tsv"

    input_file.write_text("species_name\tsample1\nspecies1\t1\n")
    conta_file.write_text("")

    monkeypatch.setattr(
        "sys.argv",
        [
            "crococlean",
            "-s",
            str(input_file),
            "-c",
            str(conta_file),
            "-o",
            str(output_file),
            "--filter-low-ab",
            "20",
            "--nproc",
            "2",
        ],
    )

    args = get_arguments()

    assert args.input_table_fp == input_file.resolve()
    assert args.conta_events_fp == conta_file.resolve()
    assert args.output_table_fp == output_file.resolve()
    assert args.filtering_ab_thr_factor == 20
    assert args.nproc == 2
