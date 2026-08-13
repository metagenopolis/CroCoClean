import logging
from collections import Counter
from multiprocessing import Pool
from time import perf_counter
from tqdm import tqdm
import pandas as pd
from crococlean.conta_event import ContaminationEvent


def run_decontamination(
    species_ab_table: pd.DataFrame, conta_events: list[ContaminationEvent], nproc: int
) -> pd.DataFrame:
    if not conta_events:
        return species_ab_table.copy()

    _warn_high_contamination_rates(conta_events)
    _warn_multiple_contamination_sources(conta_events)

    start = perf_counter()
    logging.info(
        "Performing decontamination using %d process%s...",
        nproc,
        "" if nproc == 1 else "es",
    )

    worker = DecontaminationWorker(species_ab_table)
    corrected_profiles = []

    with Pool(processes=nproc) as pool:
        all_tasks = pool.imap_unordered(
            worker.decontaminate,
            conta_events,
            chunksize=100,
        )
        pbar = tqdm(
            all_tasks,
            total=len(conta_events),
            leave=False,
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} contamination events processed{postfix}",
        )

        for corrected in pbar:
            corrected_profiles.append(corrected)

    logging.info("Decontamination completed in %.1f seconds", perf_counter() - start)

    return pd.concat(
        [species_ab_table, pd.concat(corrected_profiles, axis=1)],
        axis=1,
    )


def _warn_multiple_contamination_sources(
    conta_events: list[ContaminationEvent],
) -> None:
    target_counts = Counter(event.target for event in conta_events)

    n_targets_multiple_sources = sum(count > 1 for count in target_counts.values())

    if n_targets_multiple_sources:
        logging.warning(
            "Multiple contamination sources were detected for %d target sample%s. "
            "Each decontamination will be performed independently.",
            n_targets_multiple_sources,
            "" if n_targets_multiple_sources == 1 else "s",
        )
def _warn_high_contamination_rates(
    conta_events: list[ContaminationEvent],
    high_rate_cutoff: float = 0.10,
) -> None:
    high_rate_events = [
        event for event in conta_events
        if event.rate >= high_rate_cutoff
    ]

    if high_rate_events:
        logging.warning(
            "%d contamination event%s have a contamination rate of at least %.0f%%.",
            len(high_rate_events),
            "" if len(high_rate_events) == 1 else "s",
            high_rate_cutoff * 100,
        )
        logging.warning(
            "Decontamination may be less accurate at high "
            "contamination levels, and some species introduced by "
            "contamination may be missed and therefore not set to zero.",
        )


# pylint: disable=too-few-public-methods
class DecontaminationWorker:
    def __init__(
        self,
        species_ab_table: pd.DataFrame,
    ) -> None:
        self.species_ab_table = species_ab_table

    def decontaminate(self, conta_event: ContaminationEvent) -> pd.Series:
        source = self.species_ab_table[conta_event.source]
        target = self.species_ab_table[conta_event.target]

        # Minimum non-zero abundance in the original target profile
        min_non_zero = target[target > 0].min()

        # Remove the estimated contamination
        corrected = target - conta_event.rate * source

        # Remove abundances below the detection threshold
        corrected[corrected < min_non_zero] = 0.0

        # Species known to be exclusively introduced by contamination
        corrected.loc[conta_event.conta_line_species] = 0.0

        # Renormalize to sum to one
        total = corrected.sum()
        if total > 0:
            corrected /= total

        corrected.name = f"{conta_event.target}_deconta_{conta_event.source}"

        return corrected
