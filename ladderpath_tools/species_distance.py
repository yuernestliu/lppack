"""
Version 3.0.0
Author: codex
Revised for grouped-species distance calculation
Date: 2026.04.09

Multi-sequence species analysis based on Ladderpath:
1. get_ladderpath_species merges targets into species-level targets.
2. get_distances calculates distances from a species-level lpjson.
3. get_distance_table provides the complete operation in one call.

The current implementation validates species_info against target multiplicity.
A target ID may appear across species or repeatedly within one species, but its
total number of occurrences must match the multiplicity implied by lpjson.

Repeated source targets use dictionary input, duplications_info follows the
`{target_id: count}` convention, and species_info uses real target IDs.
"""

from __future__ import annotations

import copy
from collections import Counter

from . import lambda_from_laddergraph as lfl

__all__ = [
    "get_ladderpath_species",
    "get_distances",
    "get_distance_table",
]


def _normalize_calc_json(lpjson_species: dict) -> dict:
    """
    Prepare an object for `lambda_from_laddergraph.get()`.

    Target and ladderon keys are converted to integers, short records are
    padded, and species-level duplications_info is cleared.
    """
    calc_json = copy.deepcopy(lpjson_species)

    if "targets" in calc_json:
        calc_json["targets"] = {int(k): v for k, v in calc_json["targets"].items()}
        for target_info in calc_json["targets"].values():
            if len(target_info) < 3:
                target_info.extend([""] * (3 - len(target_info)))

    if "ladderons" in calc_json:
        calc_json["ladderons"] = {int(k): v for k, v in calc_json["ladderons"].items()}
        for ladderon_info in calc_json["ladderons"].values():
            if len(ladderon_info) < 3:
                ladderon_info.extend([""] * (3 - len(ladderon_info)))

    calc_json["duplications_info"] = {}
    return calc_json


def _resolve_input_ids(input_ids: list[int], targets: dict) -> list[int]:
    """
    Resolve supplied identifiers as real target IDs.

    Original sequence-position identifiers are not supported.
    """
    target_key_set = {int(k) for k in targets.keys()}
    resolved: list[int] = []

    for raw_tid in input_ids:
        tid = int(raw_tid)
        if tid not in target_key_set:
            raise ValueError(
                f"Input ID {tid} is not present in lpjson['targets']. "
                "species_info must use real target IDs; original sequence-position "
                "identifiers are not supported."
            )
        resolved.append(tid)

    return resolved


def _expected_target_counts(target_ids: list[int], duplications_info: dict) -> dict[int, int]:
    """
    Recover the expected source occurrence count for every target.

    A target omitted from duplications_info occurs once; a stored value is the
    target's total multiplicity.
    """
    dup_counts = {int(k): int(v) for k, v in (duplications_info or {}).items()}

    expected_counts: dict[int, int] = {}
    for tid in target_ids:
        count = dup_counts.get(int(tid), 1)
        if count < 1:
            raise ValueError(
                f"Invalid multiplicity for target ID {tid} in duplications_info: {count}"
            )
        expected_counts[int(tid)] = count

    return expected_counts


def _validate_species_partition(
    species_entries: list[dict],
    target_ids: list[int],
    duplications_info: dict,
) -> tuple[dict[int, int], dict[int, int]]:
    """
    Validate species_info occurrence counts against lpjson multiplicities.

    A target ID may occur in multiple species or repeat within one species.
    """
    expected_targets = set(target_ids)
    actual_counts: Counter[int] = Counter()

    for entry in species_entries:
        species_key = entry["key"]
        for tid in entry["resolved_ids"]:
            if tid not in expected_targets:
                raise ValueError(
                    f"{species_key} contains target ID {tid}, which is not present in lpjson['targets']."
                )
            actual_counts[tid] += 1

    expected_counts = _expected_target_counts(target_ids, duplications_info)
    mismatches = []
    for tid in sorted(expected_counts.keys(), key=abs):
        expected = expected_counts[tid]
        actual = actual_counts.get(tid, 0)
        if actual != expected:
            mismatches.append(f"{tid}: expected={expected}, actual={actual}")

    if mismatches:
        raise ValueError(
            "Target occurrence counts in species_info do not match lpjson multiplicities. "
            + "; ".join(mismatches)
        )

    return expected_counts, dict(actual_counts)


def get_ladderpath_species(lpjson: dict, species_info: dict, verbose: bool = False) -> dict:
    """
    Merge targets into species-level targets and update ladderon coordinates.

    Parameters:
        lpjson: Source Ladderpath JSON.
        species_info:
            {
                "species1": ("Morado", [-1, -2]),
                "species2": ("Musa_ornata", [-3, -4])
            }
        verbose: Print intermediate mappings when True.

    Returns:
        A species-level lpjson.

    `species_mapping` stores new species target IDs and names. The result clears
    duplications_info because distances use the merged targets and ladderons.
    species_info must use real target IDs, and its occurrence counts must match
    the source multiplicities.
    """
    processed_json = copy.deepcopy(lpjson)
    targets = processed_json.get("targets", {})
    if not targets:
        raise ValueError("The input lpjson has no targets.")
    if not species_info:
        raise ValueError("species_info must not be empty.")

    target_ids = sorted((int(k) for k in targets.keys()), key=abs)
    duplications_info = processed_json.get("duplications_info", {})

    if verbose:
        expected_counts = _expected_target_counts(target_ids, duplications_info)
        print(f"Unique target count: {len(targets)}")
        print("Partition validation mode: target_multiplicity_check")
        print("Available target IDs:")
        print(f"  {target_ids}")
        print("Expected occurrence count for each target:")
        print(f"  {expected_counts}")

    species_entries = []
    for species_key, payload in species_info.items():
        if not isinstance(payload, (tuple, list)) or len(payload) != 2:
            raise ValueError(
                f"species_info[{species_key!r}] must have the form (name, [ids])."
            )

        species_name, raw_ids = payload
        if not isinstance(raw_ids, (list, tuple)) or len(raw_ids) == 0:
            raise ValueError(
                f"The ID list in species_info[{species_key!r}] must not be empty."
            )

        normalized_raw_ids = [int(x) for x in raw_ids]
        resolved_ids = _resolve_input_ids(normalized_raw_ids, targets=targets)

        species_entries.append(
            {
                "key": species_key,
                "name": str(species_name),
                "raw_ids": normalized_raw_ids,
                "resolved_ids": resolved_ids,
            }
        )

    expected_counts, actual_counts = _validate_species_partition(
        species_entries,
        target_ids=target_ids,
        duplications_info=duplications_info,
    )

    if verbose:
        print("\nResolved species members:")
        for entry in species_entries:
            print(
                f"  {entry['key']} ({entry['name']}): "
                f"input={entry['raw_ids']} -> resolved_target_ids={entry['resolved_ids']}"
            )
        print("\nMembership count validation:")
        print(f"  expected_counts={expected_counts}")
        print(f"  actual_counts={actual_counts}")

    target_abs_max = max(abs(int(k)) for k in targets.keys())
    input_abs_max = max(abs(tid) for entry in species_entries for tid in entry["raw_ids"])
    base_abs_max = max(target_abs_max, input_abs_max)

    new_groups: dict[str, list[int]] = {}
    species_mapping: dict[str, str] = {}

    for idx, entry in enumerate(species_entries, start=1):
        new_id = -(base_abs_max + idx)
        new_id_str = str(new_id)
        new_groups[new_id_str] = entry["resolved_ids"]
        species_mapping[new_id_str] = entry["name"]

    targets_str_map = {str(k): v for k, v in targets.items()}
    new_targets: dict[str, list] = {}
    new_targets_offsets: dict[int, dict[str, list[int]]] = {}

    for new_target_id, old_ids in new_groups.items():
        merged_comp = []
        total_length = 0

        for oid in old_ids:
            key = str(oid)
            if key not in targets_str_map:
                raise ValueError(f"Target ID {oid} is not present in the source targets.")

            new_targets_offsets.setdefault(oid, {}).setdefault(new_target_id, []).append(
                total_length
            )

            old_data = targets_str_map[key]
            merged_comp.extend(old_data[0])
            total_length += old_data[1]

        new_targets[new_target_id] = [merged_comp, total_length, "", 1]

    if "ladderons" in processed_json:
        for _, ladderon_data in processed_json["ladderons"].items():
            if len(ladderon_data) < 4 or not isinstance(ladderon_data[3], dict):
                continue

            original_pos = ladderon_data[3]
            new_posi: dict[str, list[int]] = {}

            for old_target, positions in original_pos.items():
                ot = int(old_target)
                if ot in new_targets_offsets:
                    for new_target, base_offsets in new_targets_offsets[ot].items():
                        for base_offset in base_offsets:
                            shifted = [base_offset + int(pos) for pos in positions]
                            new_posi.setdefault(str(new_target), []).extend(shifted)

            ladderon_data[3] = new_posi

    processed_json["targets"] = new_targets
    processed_json["duplications_info"] = {}
    processed_json["species_mapping"] = species_mapping

    if verbose:
        print("\nNew species mapping:")
        for new_id, name in species_mapping.items():
            print(f"  {new_id} -> {name}")

    return processed_json


def get_distances(
    lpjson_species: dict,
    method: str = "jaccard",
    verbose: bool = False,
) -> tuple[dict[str, float], dict[int, str]]:
    """
    Calculate pairwise distances from a species-level lpjson.

    Parameters:
        lpjson_species: Output from `get_ladderpath_species()`.
        method: Either "jaccard" or "dice".
        verbose: Print intermediate results when True.

    Returns:
        distances_table:
            For example {"1,2": 0.62, "1,3": 0.71, ...}
        name_map:
            For example {1: "Morado", 2: "Musa_ornata"}
    """
    method = method.lower()
    if method not in {"jaccard", "dice"}:
        raise ValueError("method must be either 'jaccard' or 'dice'.")

    species_mapping = lpjson_species.get("species_mapping") or {}
    if not species_mapping:
        raise ValueError("The input has no species_mapping. Run get_ladderpath_species() first.")

    ordered_target_ids = list(species_mapping.keys())
    if len(ordered_target_ids) < 2:
        raise ValueError("At least two species are required to calculate distances.")

    calc_json = _normalize_calc_json(lpjson_species)

    lambda_single: dict[str, float] = {}
    if verbose:
        print(f"Calculating single-species lambda for {len(ordered_target_ids)} species.")

    for tid in ordered_target_ids:
        tid_int = int(tid)
        lam = lfl.get([tid_int], calc_json)
        if lam is None:
            raise ValueError(f"Could not calculate single-species lambda: target_id={tid}")
        lambda_single[tid] = lam

    display_label_map = {
        tid: idx for idx, tid in enumerate(ordered_target_ids, start=1)
    }
    sorted_targets = [(tid, display_label_map[tid]) for tid in ordered_target_ids]

    distances_table: dict[str, float] = {}

    for i in range(len(sorted_targets)):
        for j in range(i + 1, len(sorted_targets)):
            tid1, label1 = sorted_targets[i]
            tid2, label2 = sorted_targets[j]

            L1 = lambda_single[tid1]
            L2 = lambda_single[tid2]
            L12 = lfl.get([int(tid1), int(tid2)], calc_json)

            if L12 is None:
                raise ValueError(
                    f"Could not calculate joint lambda: target_id_pair=({tid1}, {tid2})"
                )

            numerator = L12 - (L1 + L2) / 2
            if method == "dice":
                denominator = (L1 + L2) / 2
            else:
                denominator = L12 / 2

            if denominator == 0:
                raise ValueError(
                    f"Distance denominator is zero: pair=({tid1}, {tid2}), method={method}"
                )

            distances_table[f"{label1},{label2}"] = float(numerator / denominator)

    name_map = {
        display_label_map[tid]: str(species_mapping[tid]) for tid in ordered_target_ids
    }

    if verbose:
        print("Single-species lambda:")
        for tid, label in sorted_targets:
            print(f"  {label}: {name_map[label]} ({tid}) -> {lambda_single[tid]}")
        print("Distance table:")
        for key, value in distances_table.items():
            print(f"  {key}: {value}")

    return distances_table, name_map


def get_distance_table(
    lpjson: dict,
    species_info: dict,
    method: str = "jaccard",
    verbose: bool = False,
) -> tuple[dict[str, float], dict[int, str], dict]:
    """
    Merge species groups and calculate their distance table in one call.

    Parameters:
        lpjson: Source Ladderpath JSON.
        species_info: Species grouping with real target IDs.
        method: "jaccard" or "dice"; defaults to "jaccard".
        verbose: Print intermediate results when True.

    Returns:
        distance_table, name_map, lpjson_merged
    """
    lpjson_merged = get_ladderpath_species(lpjson, species_info, verbose=verbose)
    distance_table, name_map = get_distances(
        lpjson_merged, method=method, verbose=verbose
    )
    return distance_table, name_map, lpjson_merged
