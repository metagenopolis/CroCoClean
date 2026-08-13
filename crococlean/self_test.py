"""End-to-end self-test for CroCoClean."""

import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable

import pandas as pd


def self_test(
    run_crococlean: Callable[
        [Path, Path, Path, float | None, int],
        None,
    ],
) -> None:
    """Run an end-to-end test of the CroCoClean installation."""
    logging.info("Running CroCoClean installation test...")

    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        input_table_fp = tmp_path / "species_abundance.tsv"
        conta_events_fp = tmp_path / "contamination_events.tsv"
        output_table_fp = tmp_path / "output.tsv"

        input_table_fp.write_text(
            "species_name\tsource_1\tsource_2\ttarget\n"
            "species_1\t0.40\t0.05\t0.10\n"
            "species_2\t0.25\t0.30\t0.20\n"
            "species_3\t0.15\t0.25\t0.30\n"
            "species_4\t0.10\t0.20\t0.25\n"
            "species_5\t0.06\t0.12\t0.10\n"
            "species_6\t0.04\t0.08\t0.05\n",
            encoding="utf8",
        )

        conta_events_fp.write_text(
            "source\ttarget\trate\tprobability\t"
            "contamination_specific_species\n"
            "source_1\ttarget\t0.05\t1.0\tspecies_1\n"
            "source_2\ttarget\t0.03\t1.0\tspecies_2\n",
            encoding="utf8",
        )

        # Disable low-abundance filtering. The test should exercise
        # the decontamination algorithm itself.
        run_crococlean(
            input_table_fp,
            conta_events_fp,
            output_table_fp,
            filtering_ab_thr_factor=None,
            nproc=1,
        )

        result = pd.read_csv(
            output_table_fp,
            sep="\t",
            index_col=0,
        )

    # Original target minimum non-zero abundance = 0.05.
    #
    # For source_1:
    #   species_1 = 0.10 - 0.05 * 0.40 = 0.080
    #   species_2 = 0.20 - 0.05 * 0.25 = 0.1875
    #   species_3 = 0.30 - 0.05 * 0.15 = 0.2925
    #   species_4 = 0.25 - 0.05 * 0.10 = 0.245
    #   species_5 = 0.10 - 0.05 * 0.06 = 0.097
    #   species_6 = 0.05 - 0.05 * 0.04 = 0.048
    #
    # species_1 is contamination-specific and species_6 is below
    # the minimum non-zero abundance, so both are set to zero.
    expected_source_1 = pd.Series(
        [
            0.0,
            0.1875,
            0.2925,
            0.245,
            0.097,
            0.0,
        ],
        index=result.index,
        name="target_deconta_source_1",
    )
    expected_source_1 /= expected_source_1.sum()

    # For source_2:
    #   species_1 = 0.10 - 0.03 * 0.05 = 0.0985
    #   species_2 = 0.20 - 0.03 * 0.30 = 0.191
    #   species_3 = 0.30 - 0.03 * 0.25 = 0.2925
    #   species_4 = 0.25 - 0.03 * 0.20 = 0.244
    #   species_5 = 0.10 - 0.03 * 0.12 = 0.0964
    #   species_6 = 0.05 - 0.03 * 0.08 = 0.0476
    #
    # species_2 is contamination-specific and species_6 is below
    # the minimum non-zero abundance, so both are set to zero.
    expected_source_2 = pd.Series(
        [
            0.0985,
            0.0,
            0.2925,
            0.244,
            0.0964,
            0.0,
        ],
        index=result.index,
        name="target_deconta_source_2",
    )
    expected_source_2 /= expected_source_2.sum()

    expected = pd.DataFrame(
        {
            "source_1": [0.40, 0.25, 0.15, 0.10, 0.06, 0.04],
            "source_2": [0.05, 0.30, 0.25, 0.20, 0.12, 0.08],
            "target": [0.10, 0.20, 0.30, 0.25, 0.10, 0.05],
            "target_deconta_source_1": expected_source_1,
            "target_deconta_source_2": expected_source_2,
        },
        index=pd.Index(
            [
                "species_1",
                "species_2",
                "species_3",
                "species_4",
                "species_5",
                "species_6",
            ],
            name="species_name",
        ),
    )

    pd.testing.assert_frame_equal(
        result,
        expected,
        check_exact=False,
        rtol=1e-12,
        atol=1e-12,
    )

    logging.info("CroCoClean installation test passed.")
