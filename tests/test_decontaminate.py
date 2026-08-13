import logging

import pandas as pd
import pytest

from crococlean.conta_event import ContaminationEvent
from crococlean.decontaminate import (
    DecontaminationWorker,
    _warn_high_contamination_rates,
    _warn_multiple_contamination_sources,
    run_decontamination,
)


@pytest.fixture
def species_ab_table():
    """Return a simple species abundance table for testing."""
    return pd.DataFrame(
        {
            "source": [0.4, 0.3, 0.3],
            "target": [0.2, 0.3, 0.5],
        },
        index=["species_1", "species_2", "species_3"],
    )


@pytest.fixture
def worker(species_ab_table):
    """Return a decontamination worker."""
    return DecontaminationWorker(species_ab_table)


def test_decontaminate_subtracts_contamination(worker):
    """Test subtraction of the estimated contamination."""
    event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.1,
        conta_line_species=[],
    )

    corrected = worker.decontaminate(event)

    # Contamination subtraction:
    # species_1: 0.2 - 0.1 * 0.4 = 0.16
    # species_2: 0.3 - 0.1 * 0.3 = 0.27
    # species_3: 0.5 - 0.1 * 0.3 = 0.47
    #
    # Minimum non-zero target abundance = 0.2.
    # Therefore species_1 (0.16) is removed.
    expected = pd.Series(
        [0.0, 0.27, 0.47],
        index=["species_1", "species_2", "species_3"],
        name="target_deconta_source",
    )
    expected /= expected.sum()

    pd.testing.assert_series_equal(corrected, expected)


def test_decontaminate_normalizes_profile(worker):
    """Test that the corrected profile is normalized."""
    event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.01,
        conta_line_species=[],
    )

    corrected = worker.decontaminate(event)

    assert corrected.sum() == pytest.approx(1.0)


def test_decontaminate_removes_low_abundance_species():
    """Test removal of corrected abundances below the detection threshold."""
    table = pd.DataFrame(
        {
            "source": [0.0, 0.0, 1.0],
            "target": [0.01, 0.49, 0.50],
        },
        index=["species_1", "species_2", "species_3"],
    )

    worker = DecontaminationWorker(table)

    event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.495,
        conta_line_species=[],
    )

    corrected = worker.decontaminate(event)

    # Minimum non-zero abundance in the original target = 0.01.
    #
    # species_3:
    # 0.50 - 0.495 * 1.0 = 0.005
    #
    # 0.005 < 0.01, so species_3 is removed.
    assert corrected["species_3"] == 0.0


def test_decontaminate_removes_contamination_specific_species(worker):
    """Test removal of contamination-specific species."""
    event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.01,
        conta_line_species=["species_1"],
    )

    corrected = worker.decontaminate(event)

    assert corrected["species_1"] == 0.0
    assert corrected.sum() == pytest.approx(1.0)


def test_decontaminate_returns_zero_profile_if_everything_removed():
    """Test that an entirely removed profile remains zero."""
    table = pd.DataFrame(
        {
            "source": [0.5, 0.5],
            "target": [0.01, 0.01],
        },
        index=["species_1", "species_2"],
    )

    worker = DecontaminationWorker(table)

    event = ContaminationEvent(
        source="source",
        target="target",
        rate=1.0,
        conta_line_species=["species_1", "species_2"],
    )

    corrected = worker.decontaminate(event)

    assert (corrected == 0).all()
    assert corrected.sum() == 0


def test_decontaminate_sets_profile_name(worker):
    """Test the name assigned to the corrected profile."""
    event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.01,
        conta_line_species=[],
    )

    corrected = worker.decontaminate(event)

    assert corrected.name == "target_deconta_source"


def test_warn_multiple_contamination_sources(caplog):
    """Test warning when a target has multiple contamination sources."""
    events = [
        ContaminationEvent(
            source="source_1",
            target="target",
            rate=0.05,
            conta_line_species=[],
        ),
        ContaminationEvent(
            source="source_2",
            target="target",
            rate=0.03,
            conta_line_species=[],
        ),
    ]

    with caplog.at_level(logging.WARNING):
        _warn_multiple_contamination_sources(events)

    assert (
        "Multiple contamination sources were detected for 1 target sample."
        in caplog.text
    )
    assert (
        "Each decontamination will be performed independently."
        in caplog.text
    )

def test_warn_multiple_contamination_sources_multiple_targets(caplog):
    """Test warning when multiple targets have multiple sources."""
    events = [
        ContaminationEvent(
            source="source_1",
            target="target_1",
            rate=0.05,
            conta_line_species=[],
        ),
        ContaminationEvent(
            source="source_2",
            target="target_1",
            rate=0.03,
            conta_line_species=[],
        ),
        ContaminationEvent(
            source="source_3",
            target="target_2",
            rate=0.04,
            conta_line_species=[],
        ),
        ContaminationEvent(
            source="source_4",
            target="target_2",
            rate=0.02,
            conta_line_species=[],
        ),
    ]

    with caplog.at_level(logging.WARNING):
        _warn_multiple_contamination_sources(events)

    assert (
        "Multiple contamination sources were detected for 2 target samples."
        in caplog.text
    )
    assert (
        "Each decontamination will be performed independently."
        in caplog.text
    )


def test_no_warning_for_single_contamination_source(caplog):
    """Test that no warning is emitted for a single source."""
    events = [
        ContaminationEvent(
            source="source",
            target="target",
            rate=0.05,
            conta_line_species=[],
        )
    ]

    with caplog.at_level(logging.WARNING):
        _warn_multiple_contamination_sources(events)

    assert not caplog.records


def test_warn_high_contamination_rates(caplog):
    """Test warning for contamination rates above the cutoff."""
    events = [
        ContaminationEvent(
            source="source",
            target="target",
            rate=0.10,
            conta_line_species=[],
        )
    ]

    with caplog.at_level(logging.WARNING):
        _warn_high_contamination_rates(events)

    assert "at least 10%" in caplog.text
    assert "may be less accurate" in caplog.text


def test_no_warning_below_high_contamination_cutoff(caplog):
    """Test that no warning is emitted below the cutoff."""
    events = [
        ContaminationEvent(
            source="source",
            target="target",
            rate=0.09,
            conta_line_species=[],
        )
    ]

    with caplog.at_level(logging.WARNING):
        _warn_high_contamination_rates(events)

    assert not caplog.records


def test_high_contamination_rate_custom_cutoff(caplog):
    """Test the configurable high contamination cutoff."""
    events = [
        ContaminationEvent(
            source="source",
            target="target",
            rate=0.15,
            conta_line_species=[],
        )
    ]

    with caplog.at_level(logging.WARNING):
        _warn_high_contamination_rates(
            events,
            high_rate_cutoff=0.20,
        )

    assert not caplog.records


def test_run_decontamination_keeps_original_profiles():
    """Test that original and corrected profiles are returned."""
    table = pd.DataFrame(
        {
            "source": [0.4, 0.3, 0.3],
            "target": [0.2, 0.3, 0.5],
        },
        index=["species_1", "species_2", "species_3"],
    )

    events = [
        ContaminationEvent(
            source="source",
            target="target",
            rate=0.01,
            conta_line_species=[],
        )
    ]

    result = run_decontamination(
        table,
        events,
        nproc=1,
    )

    assert "source" in result.columns
    assert "target" in result.columns
    assert "target_deconta_source" in result.columns

    pd.testing.assert_series_equal(
        result["source"],
        table["source"],
    )

    pd.testing.assert_series_equal(
        result["target"],
        table["target"],
    )


def test_run_decontamination_with_no_events():
    """Test that the original table is returned when there are no events."""
    table = pd.DataFrame(
        {
            "sample_1": [0.2, 0.8],
            "sample_2": [0.5, 0.5],
        },
        index=["species_1", "species_2"],
    )

    result = run_decontamination(
        table,
        [],
        nproc=1,
    )

    pd.testing.assert_frame_equal(result, table)
