#!/usr/bin/env python3

"""
Version 6.6
Authors: ecsLab (Jing Wang, Xinting Yuan, Yu Liu, Codex, et al.)
Date: 2026.06.12

Public interface:
Use the single-level sort_method argument only:
contribution_folded, contribution_frozen, or multiplicity.
Multiplicity vocabularies always contain exactly vocab_size tokens and reserve
space for every basic building block. Only ladderons may be removed to satisfy
the final vocab_size limit.

Compiled-backend note:
If further acceleration is needed, keep this file as the Python reference
implementation and only offload the dense fastpath hot loops. Prefer trying a
Numba kernel first: flatten children/parents into integer arrays, keep lpjson
loading / readonly checks / fallback / CSV output in Python, and require exact
selection-order equivalence with this file before enabling the backend.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import os
import sys
import types
import warnings
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import Any

from tqdm import tqdm


SORT_METHOD_CONTRIBUTION_FOLDED = "contribution_folded"
SORT_METHOD_CONTRIBUTION_FROZEN = "contribution_frozen"
SORT_METHOD_MULTIPLICITY = "multiplicity"
VALID_SORT_METHODS = {
    SORT_METHOD_CONTRIBUTION_FOLDED,
    SORT_METHOD_CONTRIBUTION_FROZEN,
    SORT_METHOD_MULTIPLICITY,
}
_LAST_VALID_DENSE_FASTPATH_CACHE: tuple[
    object,
    object,
    tuple[list[int], list[list[tuple[int, int]]], list[list[tuple[int, int]]], list[int], list[int], list[int]],
] | None = None
_LAST_READONLY_LPJSON_CACHE: tuple[object, object] | None = None
_LAST_POM_MULTIPLICITY_CACHE: tuple[object, dict[str, int], dict[int, int]] | None = None


class _NullProgress:
    def __enter__(self) -> "_NullProgress":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False

    def update(self, n: int = 1) -> None:
        return None


def make_progress(total: int, desc: str, unit: str, iter_logs: bool):
    if iter_logs:
        return tqdm(total=total, desc=desc, unit=unit, disable=False)
    return _NullProgress()


def warn_iter_logs_may_be_slow() -> None:
    warnings.warn(
        "iter_logs=True enables full per-iteration diagnostics and can be substantially slower. "
        "Use it only for debugging or algorithm inspection.",
        RuntimeWarning,
        stacklevel=3,
    )


def import_ladderpath():
    try:
        from .. import ladderpath as lp  # type: ignore

        return lp
    except Exception:
        sys.modules.setdefault("graphviz", types.SimpleNamespace())
        repo_root = os.environ.get("LP_DEVELOPERS_REPO_ROOT")
        if repo_root and repo_root not in sys.path:
            sys.path.append(repo_root)
        try:
            import lp_developers.ladderpath as lp
        except Exception as exc:
            raise ImportError(
                "Could not import lp_developers.ladderpath. Run this module inside the "
                "lp_developers package, install lp_developers, or set LP_DEVELOPERS_REPO_ROOT "
                "to the repository root before calling sortout."
            ) from exc

        return lp


def _to_jsonable(data: Any) -> Any:
    if isinstance(data, Mapping):
        return {key: _to_jsonable(value) for key, value in data.items()}
    if isinstance(data, tuple):
        return [_to_jsonable(value) for value in data]
    if isinstance(data, list):
        return [_to_jsonable(value) for value in data]
    return data


def save_json(data: Any, file_path: str | Path) -> None:
    with open(file_path, "w", encoding="utf-8") as file:
        json.dump(_to_jsonable(data), file, ensure_ascii=False, indent=4)


def save_csv_rows(rows: list[dict[str, Any]], file_path: str | Path) -> None:
    fieldnames = [
        "token",
        "score",
    ]
    with open(file_path, "w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def freeze_readonly(data: Any) -> Any:
    if isinstance(data, MappingProxyType):
        return data
    if isinstance(data, Mapping):
        return MappingProxyType({key: freeze_readonly(value) for key, value in data.items()})
    if isinstance(data, list):
        return tuple(freeze_readonly(value) for value in data)
    if isinstance(data, tuple):
        return tuple(freeze_readonly(value) for value in data)
    return data


def load_lpjson(
    lpjson_or_path: dict[str, Any] | str | Path,
    readonly: bool = False,
) -> dict[str, Any]:
    lp = import_ladderpath()
    global _LAST_READONLY_LPJSON_CACHE

    if isinstance(lpjson_or_path, Mapping):
        lpjson = lpjson_or_path
        if not isinstance(lpjson_or_path, MappingProxyType):
            lp.fill_lpjson_STR(lpjson)
    else:
        lpjson = lp.load_ladderpath_json(str(lpjson_or_path))
        lp.fill_lpjson_STR(lpjson)
    if readonly:
        if isinstance(lpjson, MappingProxyType):
            return lpjson
        if _LAST_READONLY_LPJSON_CACHE is not None:
            cached_source, cached_readonly = _LAST_READONLY_LPJSON_CACHE
            if cached_source is lpjson:
                return cached_readonly
        readonly_lpjson = freeze_readonly(lpjson)
        _LAST_READONLY_LPJSON_CACHE = (lpjson, readonly_lpjson)
        return readonly_lpjson
    return lpjson


def validate_sort_method(sort_method: str) -> str:
    if sort_method not in VALID_SORT_METHODS:
        raise ValueError(
            f"Unsupported sort_method={sort_method!r}. "
            f"Expected one of {sorted(VALID_SORT_METHODS)!r}."
        )
    return sort_method


def build_static_relations(
    ladderons: dict[int, list[Any]],
) -> tuple[dict[int, set[int]], dict[int, set[int]]]:
    subcomponents: dict[int, set[int]] = {
        ladderon_id: {item for item in info[0] if isinstance(item, int) and item in ladderons}
        for ladderon_id, info in ladderons.items()
    }
    superstructures: dict[int, set[int]] = {
        ladderon_id: {
            int(parent_id)
            for parent_id in info[3]
            if int(parent_id) >= 0 and int(parent_id) in ladderons
        }
        for ladderon_id, info in ladderons.items()
    }

    def closure(
        node: int,
        graph: dict[int, set[int]],
        cache: dict[int, set[int]],
    ) -> set[int]:
        if node in cache:
            return cache[node]

        affected: set[int] = set()
        for nxt in graph.get(node, ()):
            affected.add(nxt)
            affected.update(closure(nxt, graph, cache))
        cache[node] = affected
        return affected

    sub_cache: dict[int, set[int]] = {}
    super_cache: dict[int, set[int]] = {}
    subcomponent_closure = {
        ladderon_id: closure(ladderon_id, subcomponents, sub_cache)
        for ladderon_id in ladderons
    }
    superstructure_closure = {
        ladderon_id: closure(ladderon_id, superstructures, super_cache)
        for ladderon_id in ladderons
    }
    return subcomponent_closure, superstructure_closure


def calculate_importance(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    ladderon_id: int,
    ladderon_importance: dict[int, int],
    visited: set[int],
) -> int:
    if ladderon_id in ladderon_importance:
        return ladderon_importance[ladderon_id]
    if ladderon_id in visited:
        return 0

    visited.add(ladderon_id)
    importance = 0
    for target_id, occurrences in ladderons[ladderon_id][3].items():
        occurrence_count = len(occurrences)
        if int(target_id) < 0:
            importance += occurrence_count * targets[target_id][3]
        else:
            importance += occurrence_count * calculate_importance(
                ladderons,
                targets,
                int(target_id),
                ladderon_importance,
                visited,
            )

    ladderon_importance[ladderon_id] = importance
    return importance


def process_element(
    ladderon_id: int,
    ladderons: dict[int, list[Any]],
    selected_keys: set[int],
    memo: dict[int, int],
) -> int:
    if ladderon_id not in ladderons:
        return 0
    if ladderon_id in memo:
        return memo[ladderon_id]

    total = 0
    for item in ladderons[ladderon_id][0]:
        if isinstance(item, str):
            total += len(item)
        elif item in selected_keys:
            total += 1
        elif item in ladderons:
            total += process_element(item, ladderons, selected_keys, memo)

    memo[ladderon_id] = total
    return total


def build_raw_length_cache(ladderons: dict[int, list[Any]]) -> dict[int, int]:
    raw_lengths: dict[int, int] = {}
    for ladderon_id in ladderons:
        raw_lengths[ladderon_id] = process_element(ladderon_id, ladderons, set(), raw_lengths)
    return raw_lengths


def build_temp_ladderons_light(ladderons: dict[int, list[Any]]) -> dict[int, list[Any]]:
    temp: dict[int, list[Any]] = {}
    for ladderon_id, info in ladderons.items():
        temp[ladderon_id] = [info[0], None, None, dict(info[3])]
    return temp


def try_build_dense_fastpath_state(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    mode_name: str,
) -> tuple[
    str | None,
    tuple[list[int], list[list[tuple[int, int]]], list[list[tuple[int, int]]], list[int], list[int], list[int]] | None,
]:
    n_ladderons = len(ladderons)
    seen = bytearray(n_ladderons)
    candidate_order = [0] * n_ladderons
    literal_costs = [0] * n_ladderons
    target_base_importance = [0] * n_ladderons
    children: list[list[tuple[int, int]]] = [[] for _ in range(n_ladderons)]
    parents: list[list[tuple[int, int]]] = [[] for _ in range(n_ladderons)]

    for order, (ladderon_id, info) in enumerate(ladderons.items()):
        if not isinstance(ladderon_id, int) or ladderon_id < 0 or ladderon_id >= n_ladderons:
            return (
                f"{mode_name} fast path requires ladderon IDs to be dense integers in [0, n-1]; "
                f"got invalid ladderon_id={ladderon_id!r} for n={n_ladderons}. "
                "Falling back to the generic lazy implementation.",
                None,
            )
        if seen[ladderon_id]:
            return (
                f"{mode_name} fast path requires unique dense ladderon IDs, but duplicate IDs were found. "
                "Falling back to the generic lazy implementation.",
                None,
            )
        seen[ladderon_id] = 1
        candidate_order[ladderon_id] = order

        child_counts: dict[int, int] = {}
        literal_cost = 0
        for item in info[0]:
            if isinstance(item, str):
                literal_cost += len(item)
            else:
                if item < 0 or item >= n_ladderons or item <= ladderon_id:
                    return (
                        f"{mode_name} fast path requires the standard ladderpath DAG ordering "
                        "(child ladderon IDs must be greater than parent IDs). "
                        f"Found parent={ladderon_id}, child={item}. "
                        "Falling back to the generic lazy implementation.",
                        None,
                    )
                child_counts[item] = child_counts.get(item, 0) + 1

        literal_costs[ladderon_id] = literal_cost
        child_items = list(child_counts.items())
        children[ladderon_id] = child_items
        for child_id, count in child_items:
            parents[child_id].append((ladderon_id, count))

        target_contribution = 0
        for parent_id, positions in info[3].items():
            if isinstance(parent_id, int) and parent_id >= 0:
                if parent_id >= n_ladderons or parent_id >= ladderon_id:
                    return (
                        f"{mode_name} fast path requires the standard ladderpath DAG ordering "
                        "(stored parent ladderon IDs must be smaller than the current ladderon ID). "
                        f"Found ladderon_id={ladderon_id}, parent_id={parent_id}. "
                        "Falling back to the generic lazy implementation.",
                        None,
                    )
            if int(parent_id) < 0:
                target_contribution += len(positions) * targets[int(parent_id)][3]
        target_base_importance[ladderon_id] = target_contribution

    if not all(seen):
        return (
            f"{mode_name} fast path requires ladderon IDs to be dense integers in [0, n-1], "
            "but at least one ID is missing. Falling back to the generic lazy implementation.",
            None,
        )

    raw_lengths = [0] * n_ladderons
    for ladderon_id in range(n_ladderons - 1, -1, -1):
        total = literal_costs[ladderon_id]
        for child_id, count in children[ladderon_id]:
            total += count * raw_lengths[child_id]
        raw_lengths[ladderon_id] = total

    importances = [0] * n_ladderons
    for ladderon_id in range(n_ladderons):
        total = target_base_importance[ladderon_id]
        for parent_id, count in parents[ladderon_id]:
            total += count * importances[parent_id]
        importances[ladderon_id] = total

    return None, (candidate_order, children, parents, literal_costs, raw_lengths, importances)


def get_cached_dense_fastpath_state(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    mode_name: str,
) -> tuple[
    str | None,
    tuple[list[int], list[list[tuple[int, int]]], list[list[tuple[int, int]]], list[int], list[int], list[int]] | None,
]:
    global _LAST_VALID_DENSE_FASTPATH_CACHE

    if _LAST_VALID_DENSE_FASTPATH_CACHE is not None:
        cached_ladderons, cached_targets, cached_state = _LAST_VALID_DENSE_FASTPATH_CACHE
        if cached_ladderons is ladderons and cached_targets is targets:
            return None, cached_state

    issue, state = try_build_dense_fastpath_state(ladderons, targets, mode_name)
    if issue is None and state is not None:
        _LAST_VALID_DENSE_FASTPATH_CACHE = (ladderons, targets, state)
    return issue, state


def score_candidate_frozen_length_only(expand_cost: int, importance_value: int) -> int:
    return expand_cost * (importance_value - 1) - importance_value


def score_candidate_folded_from_costs(
    expand_cost: int,
    importance_value: int,
    dict_cost: int,
) -> int:
    return (expand_cost - 1) * importance_value - dict_cost


def candidate_dictionary_delta_folded(
    ladderon_id: int,
    selected_keys: set[int],
    d_cache: dict[int, int],
    stable_ladderons: dict[int, list[Any]],
    superstructure_closure: dict[int, set[int]],
) -> int:
    selected_with_candidate = set(selected_keys)
    selected_with_candidate.add(ladderon_id)

    memo: dict[int, int] = {}
    total_d = process_element(ladderon_id, stable_ladderons, selected_with_candidate, memo)
    affected_selected = superstructure_closure.get(ladderon_id, set()) & selected_keys
    for selected_id in affected_selected:
        updated = process_element(selected_id, stable_ladderons, selected_with_candidate, memo)
        total_d += updated - d_cache[selected_id]

    return total_d


def score_candidate_folded(
    ladderon_id: int,
    stable_ladderons: dict[int, list[Any]],
    selected_keys: set[int],
    importance_value: int,
    selected_dictionary_cost_cache: dict[int, int],
    superstructure_closure: dict[int, set[int]],
) -> tuple[int, int, int]:
    expand_cost = process_element(ladderon_id, stable_ladderons, selected_keys, {})
    dict_cost = candidate_dictionary_delta_folded(
        ladderon_id=ladderon_id,
        selected_keys=selected_keys,
        d_cache=selected_dictionary_cost_cache,
        stable_ladderons=stable_ladderons,
        superstructure_closure=superstructure_closure,
    )
    score = (expand_cost - 1) * importance_value - dict_cost
    return expand_cost, dict_cost, score


def build_vocab_rows(
    lpjson: Mapping[str, Any],
    selected_ladderons: dict[int, list[Any]],
    selection_trace: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    trace_by_ladderon_id = {
        int(entry["ladderon_id"]): entry for entry in (selection_trace or [])
    }

    rows: list[dict[str, Any]] = []
    idx = 0
    for ladderon_id, val in selected_ladderons.items():
        trace = trace_by_ladderon_id.get(ladderon_id, {})
        rows.append(
            {
                "token": val[2],
                "score": trace.get("contribution", ""),
                "token_id": idx,
                "ladderon_id": ladderon_id,
                "source": "ladderon",
                "selection_rank": trace.get("selection_rank", idx + 1),
                "contribution": trace.get("contribution", ""),
                "L_prime": trace.get("L_prime", ""),
                "F": trace.get("F", ""),
                "D": trace.get("D", ""),
                "sort_method": trace.get("sort_method", ""),
            }
        )
        idx += 1

    for bbb in lpjson["basic_building_blocks"]:
        rows.append(
            {
                "token": bbb,
                "score": 1,
                "token_id": idx,
                "ladderon_id": -idx,
                "source": "basic_building_block",
                "selection_rank": "",
                "contribution": 1,
                "L_prime": "",
                "F": 1,
                "D": "",
                "sort_method": "",
            }
        )
        idx += 1

    return rows


def build_multiplicity_vocab_rows(
    lpjson: Mapping[str, Any],
    selected_ladderons: dict[int, list[Any]],
    vocab_size: int,
) -> list[dict[str, Any]]:
    basic_building_blocks = list(lpjson["basic_building_blocks"])
    if vocab_size < len(basic_building_blocks):
        raise ValueError(
            f"vocab_size={vocab_size} is smaller than the number of basic building blocks "
            f"({len(basic_building_blocks)}), so a complete fallback vocabulary cannot fit."
        )

    token_multiplicity, ladderon_multiplicity = build_multiplicity_from_pom(lpjson)
    ranked_entries: list[tuple[str, int, str, int | None]] = [
        (info[2], ladderon_multiplicity[int(ladderon_id)], "ladderon", int(ladderon_id))
        for ladderon_id, info in lpjson["ladderons"].items()
    ]
    ranked_entries.extend(
        (token, token_multiplicity[token], "basic_building_block", None)
        for token in basic_building_blocks
    )
    ranked_entries.sort(key=lambda entry: -entry[1])

    selected_ids = {int(ladderon_id) for ladderon_id in selected_ladderons}
    combined_top = ranked_entries[:vocab_size]
    rows: list[dict[str, Any]] = []
    included_basic: set[str] = set()

    for token, score, source, ladderon_id in combined_top:
        # The size limit may remove ladderons only; every basic building block survives.
        if source == "ladderon" and ladderon_id not in selected_ids:
            continue
        if source == "basic_building_block":
            included_basic.add(token)
        rows.append(
            {
                "token": token,
                "score": score,
                "token_id": len(rows),
                "ladderon_id": ladderon_id if ladderon_id is not None else -(len(rows) + 1),
                "source": source,
                "selection_rank": len(rows) + 1,
                "contribution": score,
                "L_prime": "",
                "F": score,
                "D": "",
                "sort_method": SORT_METHOD_MULTIPLICITY,
            }
        )

    for token in basic_building_blocks:
        if token in included_basic:
            continue
        rows.append(
            {
                "token": token,
                "score": 0,
                "token_id": len(rows),
                "ladderon_id": -(len(rows) + 1),
                "source": "basic_building_block",
                "selection_rank": "",
                "contribution": 0,
                "L_prime": "",
                "F": 0,
                "D": "",
                "sort_method": SORT_METHOD_MULTIPLICITY,
            }
        )

    final_basic = {
        row["token"]
        for row in rows
        if row["source"] == "basic_building_block"
    }
    expected_ladderon_count = vocab_size - len(basic_building_blocks)
    actual_ladderon_count = sum(row["source"] == "ladderon" for row in rows)
    if final_basic != set(basic_building_blocks):
        raise RuntimeError("Internal multiplicity vocabulary error: not all basic building blocks were retained.")
    if actual_ladderon_count != expected_ladderon_count:
        raise RuntimeError(
            "Internal multiplicity vocabulary error: expected "
            f"{expected_ladderon_count} ladderons, got {actual_ladderon_count}."
        )
    if len(rows) != vocab_size:
        raise RuntimeError(
            f"Internal multiplicity vocabulary size mismatch: expected {vocab_size}, got {len(rows)}."
        )
    return rows


def update_ladderons(temp_ladderons: dict[int, list[Any]], selected_id: int) -> None:
    if selected_id not in temp_ladderons:
        return

    for item in temp_ladderons[selected_id][0]:
        if isinstance(item, int) and item in temp_ladderons:
            temp_ladderons[item][3].pop(selected_id, None)


def update_branch_importance(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    affected_subcomponents: set[int],
    ladderon_importance: dict[int, int],
) -> None:
    for ladderon_id in affected_subcomponents:
        ladderon_importance.pop(ladderon_id, None)

    for ladderon_id in affected_subcomponents:
        if ladderon_id in ladderons:
            calculate_importance(ladderons, targets, ladderon_id, ladderon_importance, set())


def _normalize_vocab_sizes(vocab_sizes: list[int]) -> list[int]:
    normalized = sorted({int(vocab_size) for vocab_size in vocab_sizes if int(vocab_size) > 0})
    if not normalized:
        raise ValueError("vocab_sizes must contain at least one positive integer.")
    return normalized


def build_multiplicity_from_pom(
    lpjson: Mapping[str, Any],
) -> tuple[dict[str, int], dict[int, int]]:
    global _LAST_POM_MULTIPLICITY_CACHE
    if _LAST_POM_MULTIPLICITY_CACHE is not None:
        cached_source, cached_token_multiplicity, cached_ladderon_multiplicity = _LAST_POM_MULTIPLICITY_CACHE
        if cached_source is lpjson:
            return cached_token_multiplicity, cached_ladderon_multiplicity

    lp = import_ladderpath()
    mutable_lpjson = _to_jsonable(lpjson) if isinstance(lpjson, MappingProxyType) else lpjson
    pom, _pom_str = lp.POM_from_JSON(mutable_lpjson, display_str=False)

    token_multiplicity = {}
    for layer in pom.values():
        for token, count in layer.items():
            token_multiplicity[token] = int(count)

    multiplicity_by_ladderon_id = {
        int(ladderon_id): token_multiplicity[ladderon_info[2]]
        for ladderon_id, ladderon_info in lpjson["ladderons"].items()
    }

    _LAST_POM_MULTIPLICITY_CACHE = (lpjson, token_multiplicity, multiplicity_by_ladderon_id)

    return token_multiplicity, multiplicity_by_ladderon_id


def build_ladderon_multiplicity_from_pom(lpjson: Mapping[str, Any]) -> dict[int, int]:
    return build_multiplicity_from_pom(lpjson)[1]


def select_ladderon_checkpoints_by_multiplicity(
    ladderons: dict[int, list[Any]],
    vocab_sizes: list[int],
    multiplicity_by_ladderon_id: dict[int, int],
    selection_trace: list[dict[str, Any]] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    checkpoints = _normalize_vocab_sizes(vocab_sizes)
    max_vocab_size = checkpoints[-1]
    checkpoint_set = set(checkpoints)

    ranked_ladderons = sorted(
        ladderons.items(),
        key=lambda item: (-multiplicity_by_ladderon_id[int(item[0])], int(item[0])),
    )

    selected: dict[int, list[Any]] = {}
    snapshots: dict[int, dict[int, list[Any]]] = {}
    limit = min(max_vocab_size, len(ranked_ladderons))
    for rank, (ladderon_id, ladderon_info) in enumerate(ranked_ladderons[:limit], start=1):
        score = multiplicity_by_ladderon_id[int(ladderon_id)]
        if selection_trace is not None:
            selection_trace.append(
                {
                    "selection_rank": rank,
                    "ladderon_id": ladderon_id,
                    "token": ladderon_info[2],
                    "contribution": score,
                    "L_prime": "",
                    "F": score,
                    "D": "",
                    "sort_method": SORT_METHOD_MULTIPLICITY,
                }
            )
        selected[ladderon_id] = ladderon_info
        if rank in checkpoint_set:
            snapshots[rank] = dict(selected)

    final_snapshot = dict(selected)
    for checkpoint in checkpoints:
        snapshots.setdefault(checkpoint, final_snapshot)

    return snapshots


def _select_ladderon_checkpoints_by_importance_folded(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_sizes: list[int],
    iter_logs: bool,
    selection_trace: list[dict[str, Any]] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    checkpoints = _normalize_vocab_sizes(vocab_sizes)
    max_vocab_size = checkpoints[-1]
    checkpoint_set = set(checkpoints)

    selected: dict[int, list[Any]] = {}
    snapshots: dict[int, dict[int, list[Any]]] = {}
    temp_ladderons = build_temp_ladderons_light(ladderons)

    subcomponent_closure, superstructure_closure = build_static_relations(ladderons)
    raw_length_cache = build_raw_length_cache(ladderons)
    ladderon_importance: dict[int, int] = {}
    all_iteration_logs: list[dict[str, Any]] = [] if iter_logs else []
    candidate_order = {ladderon_id: index for index, ladderon_id in enumerate(ladderons)}

    for ladderon_id in temp_ladderons:
        calculate_importance(temp_ladderons, targets, ladderon_id, ladderon_importance, set())

    selected_dictionary_cost_cache: dict[int, int] = {}
    expand_cost_cache: dict[int, int] = dict(raw_length_cache)
    dict_cost_cache: dict[int, int] = dict(raw_length_cache)
    score_cache: dict[int, int] = {
        ladderon_id: (raw_length_cache[ladderon_id] - 1) * ladderon_importance[ladderon_id] - raw_length_cache[ladderon_id]
        for ladderon_id in temp_ladderons
    }
    heap_versions: dict[int, int] = {ladderon_id: 0 for ladderon_id in temp_ladderons}
    score_heap: list[tuple[int, int, int, int]] = [
        (-score_cache[ladderon_id], candidate_order[ladderon_id], heap_versions[ladderon_id], ladderon_id)
        for ladderon_id in temp_ladderons
    ]
    heapq.heapify(score_heap)

    def pop_best_candidate() -> int | None:
        while score_heap:
            neg_score, _order, version, ladderon_id = heapq.heappop(score_heap)
            if ladderon_id not in temp_ladderons:
                continue
            if version != heap_versions.get(ladderon_id):
                continue
            return ladderon_id
        return None

    with make_progress(
        total=max_vocab_size,
        desc="Iterate_token_selection",
        unit=" iteration",
        iter_logs=iter_logs,
    ) as pbar:
        while temp_ladderons and len(selected) < max_vocab_size:
            if iter_logs:
                iteration_logs: dict[int, Any] = {}
                for ladderon_id in temp_ladderons:
                    iteration_logs[ladderon_id] = {
                        "L_prime": expand_cost_cache[ladderon_id],
                        "ref_cost": 1,
                        "F": ladderon_importance[ladderon_id],
                        "D": dict_cost_cache[ladderon_id],
                        "C": score_cache[ladderon_id],
                    }
                all_iteration_logs.append({"iteration": len(selected) + 1, "ladderons": iteration_logs})

            selected_id = pop_best_candidate()
            if selected_id is None:
                break

            selected_ancestors_before = superstructure_closure.get(selected_id, set()) & set(selected.keys())
            if selection_trace is not None:
                selection_trace.append(
                    {
                        "selection_rank": len(selected) + 1,
                        "ladderon_id": selected_id,
                        "token": ladderons[selected_id][2],
                        "contribution": score_cache[selected_id],
                        "L_prime": expand_cost_cache[selected_id],
                        "F": ladderon_importance[selected_id],
                        "D": dict_cost_cache[selected_id],
                        "sort_method": SORT_METHOD_CONTRIBUTION_FOLDED,
                    }
                )
            selected[selected_id] = ladderons[selected_id]
            if len(selected) in checkpoint_set:
                snapshots[len(selected)] = dict(selected)

            selected_keys_after = set(selected.keys())

            update_ladderons(temp_ladderons, selected_id)
            affected_subcomponents = subcomponent_closure.get(selected_id, set()) & set(temp_ladderons.keys())
            update_branch_importance(temp_ladderons, targets, affected_subcomponents, ladderon_importance)

            affected_selected_tokens = set(selected_ancestors_before)
            affected_selected_tokens.add(selected_id)
            selected_cost_memo: dict[int, int] = {}
            for ladderon_id in affected_selected_tokens:
                selected_dictionary_cost_cache[ladderon_id] = process_element(
                    ladderon_id,
                    ladderons,
                    selected_keys_after,
                    selected_cost_memo,
                )

            temp_ladderons.pop(selected_id, None)
            expand_cost_cache.pop(selected_id, None)
            dict_cost_cache.pop(selected_id, None)
            score_cache.pop(selected_id, None)
            heap_versions.pop(selected_id, None)

            dirty_candidates = (
                affected_subcomponents | superstructure_closure.get(selected_id, set())
            ) & set(temp_ladderons.keys())
            for selected_ancestor in selected_ancestors_before:
                dirty_candidates |= subcomponent_closure.get(selected_ancestor, set()) & set(temp_ladderons.keys())

            for ladderon_id in dirty_candidates:
                expand_cost, dict_cost, score = score_candidate_folded(
                    ladderon_id=ladderon_id,
                    stable_ladderons=ladderons,
                    selected_keys=selected_keys_after,
                    importance_value=ladderon_importance[ladderon_id],
                    selected_dictionary_cost_cache=selected_dictionary_cost_cache,
                    superstructure_closure=superstructure_closure,
                )
                expand_cost_cache[ladderon_id] = expand_cost
                dict_cost_cache[ladderon_id] = dict_cost
                score_cache[ladderon_id] = score
                heap_versions[ladderon_id] = heap_versions.get(ladderon_id, 0) + 1
                heapq.heappush(
                    score_heap,
                    (-score, candidate_order[ladderon_id], heap_versions[ladderon_id], ladderon_id),
                )

            pbar.update(1)

    if iter_logs:
        save_json(all_iteration_logs, "importance_debug_log.json")

    final_snapshot = dict(selected)
    for checkpoint in checkpoints:
        snapshots.setdefault(checkpoint, final_snapshot)

    return snapshots


def _select_ladderon_checkpoints_by_importance_folded_dense_fastpath(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_sizes: list[int],
    iter_logs: bool,
    selection_trace: list[dict[str, Any]] | None = None,
    dense_state: tuple[list[int], list[list[tuple[int, int]]], list[list[tuple[int, int]]], list[int], list[int], list[int]] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    checkpoints = _normalize_vocab_sizes(vocab_sizes)
    max_vocab_size = checkpoints[-1]
    checkpoint_set = set(checkpoints)
    n_ladderons = len(ladderons)

    selected: dict[int, list[Any]] = {}
    snapshots: dict[int, dict[int, list[Any]]] = {}
    all_iteration_logs: list[dict[str, Any]] = [] if iter_logs else []

    if dense_state is None:
        _issue, dense_state = get_cached_dense_fastpath_state(ladderons, targets, "folded")
        if dense_state is None:
            raise ValueError("folded dense fast path state could not be built.")
    candidate_order, children, parents, _literal_costs, base_expand_cost_cache, base_ladderon_importance = dense_state
    expand_cost_cache = list(base_expand_cost_cache)
    ladderon_importance = list(base_ladderon_importance)
    dictionary_importance = [0] * n_ladderons
    dict_cost_cache = list(expand_cost_cache)
    score_cache = [
        score_candidate_folded_from_costs(
            expand_cost_cache[ladderon_id],
            ladderon_importance[ladderon_id],
            dict_cost_cache[ladderon_id],
        )
        for ladderon_id in range(n_ladderons)
    ]
    heap_versions = [0] * n_ladderons
    score_heap: list[tuple[int, int, int, int]] = [
        (-score_cache[ladderon_id], candidate_order[ladderon_id], 0, ladderon_id)
        for ladderon_id in range(n_ladderons)
    ]
    heapq.heapify(score_heap)

    alive = bytearray(b"\x01") * n_ladderons
    changed_nodes: list[int] = []
    changed_marks = [0] * n_ladderons
    change_epoch = 0

    def pop_best_candidate() -> int | None:
        while score_heap:
            neg_score, _order, version, ladderon_id = heapq.heappop(score_heap)
            if not alive[ladderon_id]:
                continue
            if version != heap_versions[ladderon_id]:
                continue
            return ladderon_id
        return None

    with make_progress(
        total=max_vocab_size,
        desc="Iterate_token_selection",
        unit=" iteration",
        iter_logs=iter_logs,
    ) as pbar:
        while len(selected) < max_vocab_size and len(selected) < n_ladderons:
            if iter_logs:
                iteration_logs: dict[int, Any] = {}
                for ladderon_id in range(n_ladderons):
                    if not alive[ladderon_id]:
                        continue
                    iteration_logs[ladderon_id] = {
                        "L_prime": expand_cost_cache[ladderon_id],
                        "ref_cost": 1,
                        "F": ladderon_importance[ladderon_id],
                        "D": dict_cost_cache[ladderon_id],
                        "C": score_cache[ladderon_id],
                    }
                all_iteration_logs.append({"iteration": len(selected) + 1, "ladderons": iteration_logs})

            selected_id = pop_best_candidate()
            if selected_id is None:
                break

            selected_cost = expand_cost_cache[selected_id]
            selected_importance = ladderon_importance[selected_id]
            selected_dictionary_gain = dictionary_importance[selected_id]
            selected_dict_cost = selected_cost - (selected_cost - 1) * selected_dictionary_gain
            selected_score = score_cache[selected_id]

            if selection_trace is not None:
                selection_trace.append(
                    {
                        "selection_rank": len(selected) + 1,
                        "ladderon_id": selected_id,
                        "token": ladderons[selected_id][2],
                        "contribution": selected_score,
                        "L_prime": selected_cost,
                        "F": selected_importance,
                        "D": selected_dict_cost,
                        "sort_method": SORT_METHOD_CONTRIBUTION_FOLDED,
                    }
                )

            selected[selected_id] = ladderons[selected_id]
            alive[selected_id] = 0
            if len(selected) in checkpoint_set:
                snapshots[len(selected)] = dict(selected)

            length_delta = 1 - selected_cost
            change_epoch += 1
            changed_nodes.clear()

            if length_delta:
                stack: list[tuple[int, int]] = []
                for parent_id, count in parents[selected_id]:
                    if alive[parent_id]:
                        stack.append((parent_id, count * length_delta))
                while stack:
                    ladderon_id, delta = stack.pop()
                    if delta == 0 or not alive[ladderon_id]:
                        continue
                    expand_cost_cache[ladderon_id] += delta
                    if changed_marks[ladderon_id] != change_epoch:
                        changed_marks[ladderon_id] = change_epoch
                        changed_nodes.append(ladderon_id)
                    for parent_id, count in parents[ladderon_id]:
                        if alive[parent_id]:
                            stack.append((parent_id, count * delta))

            target_importance_delta = -selected_importance
            dictionary_importance_delta = 1 - selected_dictionary_gain
            if target_importance_delta or dictionary_importance_delta:
                stack: list[tuple[int, int, int]] = []
                for child_id, count in children[selected_id]:
                    if alive[child_id]:
                        stack.append(
                            (
                                child_id,
                                count * target_importance_delta,
                                count * dictionary_importance_delta,
                            )
                        )
                while stack:
                    ladderon_id, target_delta, dictionary_delta = stack.pop()
                    if (target_delta == 0 and dictionary_delta == 0) or not alive[ladderon_id]:
                        continue
                    if target_delta:
                        ladderon_importance[ladderon_id] += target_delta
                    if dictionary_delta:
                        dictionary_importance[ladderon_id] += dictionary_delta
                    if changed_marks[ladderon_id] != change_epoch:
                        changed_marks[ladderon_id] = change_epoch
                        changed_nodes.append(ladderon_id)
                    for child_id, count in children[ladderon_id]:
                        if alive[child_id]:
                            stack.append(
                                (
                                    child_id,
                                    count * target_delta,
                                    count * dictionary_delta,
                                )
                            )

            dict_cost_cache[selected_id] = selected_dict_cost

            for ladderon_id in changed_nodes:
                dict_cost = (
                    expand_cost_cache[ladderon_id]
                    - (expand_cost_cache[ladderon_id] - 1) * dictionary_importance[ladderon_id]
                )
                dict_cost_cache[ladderon_id] = dict_cost
                new_score = score_candidate_folded_from_costs(
                    expand_cost_cache[ladderon_id],
                    ladderon_importance[ladderon_id],
                    dict_cost,
                )
                if new_score == score_cache[ladderon_id]:
                    continue
                score_cache[ladderon_id] = new_score
                heap_versions[ladderon_id] += 1
                heapq.heappush(
                    score_heap,
                    (-new_score, candidate_order[ladderon_id], heap_versions[ladderon_id], ladderon_id),
                )

            pbar.update(1)

    if iter_logs:
        save_json(all_iteration_logs, "importance_debug_log.json")

    final_snapshot = dict(selected)
    for checkpoint in checkpoints:
        snapshots.setdefault(checkpoint, final_snapshot)

    return snapshots


def _select_ladderon_checkpoints_by_importance_folded_incremental(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_sizes: list[int],
    iter_logs: bool,
    selection_trace: list[dict[str, Any]] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    fastpath_issue, dense_state = get_cached_dense_fastpath_state(ladderons, targets, "folded")
    if fastpath_issue is None and dense_state is not None:
        return _select_ladderon_checkpoints_by_importance_folded_dense_fastpath(
            ladderons=ladderons,
            targets=targets,
            vocab_sizes=vocab_sizes,
            iter_logs=iter_logs,
            selection_trace=selection_trace,
            dense_state=dense_state,
        )
    warnings.warn(fastpath_issue, RuntimeWarning, stacklevel=2)
    return _select_ladderon_checkpoints_by_importance_folded(
        ladderons=ladderons,
        targets=targets,
        vocab_sizes=vocab_sizes,
        iter_logs=iter_logs,
        selection_trace=selection_trace,
    )


def _select_ladderon_checkpoints_by_importance_frozen_lazy(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_sizes: list[int],
    iter_logs: bool,
    selection_trace: list[dict[str, Any]] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    checkpoints = _normalize_vocab_sizes(vocab_sizes)
    max_vocab_size = checkpoints[-1]
    checkpoint_set = set(checkpoints)

    selected: dict[int, list[Any]] = {}
    snapshots: dict[int, dict[int, list[Any]]] = {}
    temp_ladderons = build_temp_ladderons_light(ladderons)
    temp_keys = set(temp_ladderons.keys())
    subcomponent_closure, superstructure_closure = build_static_relations(ladderons)

    raw_length_cache = build_raw_length_cache(ladderons)
    ladderon_importance: dict[int, int] = {}
    all_iteration_logs: list[dict[str, Any]] = [] if iter_logs else []
    candidate_order = {ladderon_id: index for index, ladderon_id in enumerate(ladderons)}
    for ladderon_id in temp_keys:
        calculate_importance(temp_ladderons, targets, ladderon_id, ladderon_importance, set())

    selected_keys: set[int] = set()
    expand_cost_cache: dict[int, int] = dict(raw_length_cache)
    dict_cost_cache: dict[int, int] = dict(raw_length_cache)
    score_cache: dict[int, int] = {
        ladderon_id: score_candidate_frozen_length_only(raw_length_cache[ladderon_id], ladderon_importance[ladderon_id])
        for ladderon_id in temp_keys
    }
    heap_versions: dict[int, int] = {ladderon_id: 0 for ladderon_id in temp_keys}
    dirty_candidates: set[int] = set()
    score_heap: list[tuple[int, int, int, int]] = [
        (-score_cache[ladderon_id], candidate_order[ladderon_id], heap_versions[ladderon_id], ladderon_id)
        for ladderon_id in temp_keys
    ]
    heapq.heapify(score_heap)

    def mark_dirty(candidates: set[int]) -> None:
        for candidate_id in candidates:
            if candidate_id in temp_keys:
                dirty_candidates.add(candidate_id)

    def recompute_candidate(ladderon_id: int) -> None:
        expand_cost = process_element(ladderon_id, ladderons, selected_keys, {})
        expand_cost_cache[ladderon_id] = expand_cost
        dict_cost_cache[ladderon_id] = expand_cost
        score_cache[ladderon_id] = score_candidate_frozen_length_only(expand_cost, ladderon_importance[ladderon_id])
        heap_versions[ladderon_id] = heap_versions.get(ladderon_id, 0) + 1
        heapq.heappush(
            score_heap,
            (-score_cache[ladderon_id], candidate_order[ladderon_id], heap_versions[ladderon_id], ladderon_id),
        )
        dirty_candidates.discard(ladderon_id)

    def pop_best_candidate() -> int | None:
        while score_heap:
            neg_score, _order, version, ladderon_id = heapq.heappop(score_heap)
            if ladderon_id not in temp_keys:
                continue
            if version != heap_versions.get(ladderon_id):
                continue
            if ladderon_id in dirty_candidates:
                recompute_candidate(ladderon_id)
                continue
            return ladderon_id
        return None

    with make_progress(
        total=max_vocab_size,
        desc="Iterate_token_selection",
        unit=" iteration",
        iter_logs=iter_logs,
    ) as pbar:
        while temp_keys and len(selected) < max_vocab_size:
            if iter_logs:
                iteration_logs: dict[int, Any] = {}
                for ladderon_id in temp_keys:
                    iteration_logs[ladderon_id] = {
                        "L_prime": expand_cost_cache[ladderon_id],
                        "ref_cost": 1,
                        "F": ladderon_importance[ladderon_id],
                        "D": dict_cost_cache[ladderon_id],
                        "C": score_cache[ladderon_id],
                    }
                all_iteration_logs.append({"iteration": len(selected) + 1, "ladderons": iteration_logs})

            selected_id = pop_best_candidate()
            if selected_id is None:
                break

            if selection_trace is not None:
                selection_trace.append(
                    {
                        "selection_rank": len(selected) + 1,
                        "ladderon_id": selected_id,
                        "token": ladderons[selected_id][2],
                        "contribution": score_cache[selected_id],
                        "L_prime": expand_cost_cache[selected_id],
                        "F": ladderon_importance[selected_id],
                        "D": expand_cost_cache[selected_id],
                        "sort_method": SORT_METHOD_CONTRIBUTION_FROZEN,
                    }
                )
            selected[selected_id] = ladderons[selected_id]
            selected_keys.add(selected_id)
            if len(selected) in checkpoint_set:
                snapshots[len(selected)] = dict(selected)

            update_ladderons(temp_ladderons, selected_id)
            affected_subcomponents = subcomponent_closure.get(selected_id, set()) & temp_keys
            update_branch_importance(temp_ladderons, targets, affected_subcomponents, ladderon_importance)

            temp_ladderons.pop(selected_id, None)
            temp_keys.discard(selected_id)
            expand_cost_cache.pop(selected_id, None)
            dict_cost_cache.pop(selected_id, None)
            score_cache.pop(selected_id, None)
            heap_versions.pop(selected_id, None)
            dirty_candidates.discard(selected_id)

            mark_dirty((affected_subcomponents | superstructure_closure.get(selected_id, set())) & temp_keys)

            pbar.update(1)

    if iter_logs:
        save_json(all_iteration_logs, "importance_debug_log.json")

    final_snapshot = dict(selected)
    for checkpoint in checkpoints:
        snapshots.setdefault(checkpoint, final_snapshot)

    return snapshots


def _select_ladderon_checkpoints_by_importance_frozen_dense_fastpath(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_sizes: list[int],
    iter_logs: bool,
    selection_trace: list[dict[str, Any]] | None = None,
    dense_state: tuple[list[int], list[list[tuple[int, int]]], list[list[tuple[int, int]]], list[int], list[int], list[int]] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    checkpoints = _normalize_vocab_sizes(vocab_sizes)
    max_vocab_size = checkpoints[-1]
    checkpoint_set = set(checkpoints)
    n_ladderons = len(ladderons)

    selected: dict[int, list[Any]] = {}
    snapshots: dict[int, dict[int, list[Any]]] = {}
    all_iteration_logs: list[dict[str, Any]] = [] if iter_logs else []

    if dense_state is None:
        _issue, dense_state = get_cached_dense_fastpath_state(ladderons, targets, "frozen")
        if dense_state is None:
            raise ValueError("frozen dense fast path state could not be built.")
    candidate_order, children, parents, _literal_costs, base_expand_cost_cache, base_ladderon_importance = dense_state
    expand_cost_cache = list(base_expand_cost_cache)
    ladderon_importance = list(base_ladderon_importance)

    score_cache = [
        score_candidate_frozen_length_only(expand_cost_cache[ladderon_id], ladderon_importance[ladderon_id])
        for ladderon_id in range(n_ladderons)
    ]
    heap_versions = [0] * n_ladderons
    score_heap: list[tuple[int, int, int, int]] = [
        (-score_cache[ladderon_id], candidate_order[ladderon_id], 0, ladderon_id)
        for ladderon_id in range(n_ladderons)
    ]
    heapq.heapify(score_heap)

    alive = bytearray(b"\x01") * n_ladderons
    changed_nodes: list[int] = []
    changed_marks = [0] * n_ladderons
    change_epoch = 0

    def pop_best_candidate() -> int | None:
        while score_heap:
            neg_score, _order, version, ladderon_id = heapq.heappop(score_heap)
            if not alive[ladderon_id]:
                continue
            if version != heap_versions[ladderon_id]:
                continue
            return ladderon_id
        return None

    with make_progress(
        total=max_vocab_size,
        desc="Iterate_token_selection",
        unit=" iteration",
        iter_logs=iter_logs,
    ) as pbar:
        while len(selected) < max_vocab_size and len(selected) < n_ladderons:
            if iter_logs:
                iteration_logs: dict[int, Any] = {}
                for ladderon_id in range(n_ladderons):
                    if not alive[ladderon_id]:
                        continue
                    iteration_logs[ladderon_id] = {
                        "L_prime": expand_cost_cache[ladderon_id],
                        "ref_cost": 1,
                        "F": ladderon_importance[ladderon_id],
                        "D": expand_cost_cache[ladderon_id],
                        "C": score_cache[ladderon_id],
                    }
                all_iteration_logs.append({"iteration": len(selected) + 1, "ladderons": iteration_logs})

            selected_id = pop_best_candidate()
            if selected_id is None:
                break

            selected_cost = expand_cost_cache[selected_id]
            selected_importance = ladderon_importance[selected_id]
            selected_score = score_cache[selected_id]

            if selection_trace is not None:
                selection_trace.append(
                    {
                        "selection_rank": len(selected) + 1,
                        "ladderon_id": selected_id,
                        "token": ladderons[selected_id][2],
                        "contribution": selected_score,
                        "L_prime": selected_cost,
                        "F": selected_importance,
                        "D": selected_cost,
                        "sort_method": SORT_METHOD_CONTRIBUTION_FROZEN,
                    }
                )

            selected[selected_id] = ladderons[selected_id]
            alive[selected_id] = 0
            if len(selected) in checkpoint_set:
                snapshots[len(selected)] = dict(selected)

            change_epoch += 1
            changed_nodes.clear()
            length_delta = 1 - selected_cost
            if length_delta:
                stack: list[tuple[int, int]] = []
                for parent_id, count in parents[selected_id]:
                    if alive[parent_id]:
                        stack.append((parent_id, count * length_delta))
                while stack:
                    ladderon_id, delta = stack.pop()
                    if delta == 0 or not alive[ladderon_id]:
                        continue
                    expand_cost_cache[ladderon_id] += delta
                    if changed_marks[ladderon_id] != change_epoch:
                        changed_marks[ladderon_id] = change_epoch
                        changed_nodes.append(ladderon_id)
                    for parent_id, count in parents[ladderon_id]:
                        if alive[parent_id]:
                            stack.append((parent_id, count * delta))

            importance_delta = -selected_importance
            if importance_delta:
                stack = []
                for child_id, count in children[selected_id]:
                    if alive[child_id]:
                        stack.append((child_id, count * importance_delta))
                while stack:
                    ladderon_id, delta = stack.pop()
                    if delta == 0 or not alive[ladderon_id]:
                        continue
                    ladderon_importance[ladderon_id] += delta
                    if changed_marks[ladderon_id] != change_epoch:
                        changed_marks[ladderon_id] = change_epoch
                        changed_nodes.append(ladderon_id)
                    for child_id, count in children[ladderon_id]:
                        if alive[child_id]:
                            stack.append((child_id, count * delta))

            for ladderon_id in changed_nodes:
                new_score = score_candidate_frozen_length_only(
                    expand_cost_cache[ladderon_id],
                    ladderon_importance[ladderon_id],
                )
                if new_score == score_cache[ladderon_id]:
                    continue
                score_cache[ladderon_id] = new_score
                heap_versions[ladderon_id] += 1
                heapq.heappush(
                    score_heap,
                    (-new_score, candidate_order[ladderon_id], heap_versions[ladderon_id], ladderon_id),
                )

            pbar.update(1)

    if iter_logs:
        save_json(all_iteration_logs, "importance_debug_log.json")

    final_snapshot = dict(selected)
    for checkpoint in checkpoints:
        snapshots.setdefault(checkpoint, final_snapshot)

    return snapshots


def _select_ladderon_checkpoints_by_importance_frozen_incremental(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_sizes: list[int],
    iter_logs: bool,
    selection_trace: list[dict[str, Any]] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    fastpath_issue, dense_state = get_cached_dense_fastpath_state(ladderons, targets, "frozen")
    if fastpath_issue is None and dense_state is not None:
        return _select_ladderon_checkpoints_by_importance_frozen_dense_fastpath(
            ladderons=ladderons,
            targets=targets,
            vocab_sizes=vocab_sizes,
            iter_logs=iter_logs,
            selection_trace=selection_trace,
            dense_state=dense_state,
        )
    warnings.warn(fastpath_issue, RuntimeWarning, stacklevel=2)
    return _select_ladderon_checkpoints_by_importance_frozen_lazy(
        ladderons=ladderons,
        targets=targets,
        vocab_sizes=vocab_sizes,
        iter_logs=iter_logs,
        selection_trace=selection_trace,
    )


def select_ladderon_checkpoints_by_sort_method(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_sizes: list[int],
    iter_logs: bool,
    sort_method: str = SORT_METHOD_CONTRIBUTION_FOLDED,
    selection_trace: list[dict[str, Any]] | None = None,
    lpjson: Mapping[str, Any] | None = None,
) -> dict[int, dict[int, list[Any]]]:
    if iter_logs:
        warn_iter_logs_may_be_slow()
    sort_method = validate_sort_method(sort_method)
    if sort_method == SORT_METHOD_MULTIPLICITY:
        if lpjson is None:
            raise ValueError("lpjson is required when sort_method='multiplicity'.")
        checkpoints = _normalize_vocab_sizes(vocab_sizes)
        basic_count = len(lpjson["basic_building_blocks"])
        if checkpoints[0] < basic_count:
            raise ValueError(
                f"vocab_size={checkpoints[0]} is smaller than the number of basic building blocks "
                f"({basic_count}), so a complete fallback vocabulary cannot fit."
            )
        ladderon_sizes = [checkpoint - basic_count for checkpoint in checkpoints if checkpoint > basic_count]
        if not ladderon_sizes:
            return {checkpoint: {} for checkpoint in checkpoints}
        ladderon_snapshots = select_ladderon_checkpoints_by_multiplicity(
            ladderons=ladderons,
            vocab_sizes=ladderon_sizes,
            multiplicity_by_ladderon_id=build_ladderon_multiplicity_from_pom(lpjson),
            selection_trace=selection_trace,
        )
        return {
            checkpoint: ladderon_snapshots.get(checkpoint - basic_count, {})
            for checkpoint in checkpoints
        }
    if sort_method == SORT_METHOD_CONTRIBUTION_FROZEN:
        return _select_ladderon_checkpoints_by_importance_frozen_incremental(
            ladderons=ladderons,
            targets=targets,
            vocab_sizes=vocab_sizes,
            iter_logs=iter_logs,
            selection_trace=selection_trace,
        )
    if sort_method != SORT_METHOD_CONTRIBUTION_FOLDED:
        raise ValueError(f"Unhandled sort_method={sort_method!r}.")
    return _select_ladderon_checkpoints_by_importance_folded_incremental(
        ladderons=ladderons,
        targets=targets,
        vocab_sizes=vocab_sizes,
        iter_logs=iter_logs,
        selection_trace=selection_trace,
    )


def select_ladderons_by_sort_method(
    ladderons: dict[int, list[Any]],
    targets: dict[int, list[Any]],
    vocab_size: int,
    iter_logs: bool,
    sort_method: str = SORT_METHOD_CONTRIBUTION_FOLDED,
    selection_trace: list[dict[str, Any]] | None = None,
    lpjson: Mapping[str, Any] | None = None,
) -> dict[int, list[Any]]:
    snapshots = select_ladderon_checkpoints_by_sort_method(
        ladderons=ladderons,
        targets=targets,
        vocab_sizes=[vocab_size],
        iter_logs=iter_logs,
        sort_method=sort_method,
        selection_trace=selection_trace,
        lpjson=lpjson,
    )
    return snapshots[vocab_size]


def _select_token_from_loaded_lpjson(
    lpjson: Mapping[str, Any],
    vocab_size: int,
    save_sort_token_file: str | None = None,
    show_selectedladderons: bool = False,
    iter_logs: bool = False,
    sort_method: str = SORT_METHOD_CONTRIBUTION_FOLDED,
    selection_trace: list[dict[str, Any]] | None = None,
) -> dict[int, list[Any]]:
    ladderons, targets = lpjson["ladderons"], lpjson["targets"]
    selected_ladderons = select_ladderons_by_sort_method(
        ladderons=ladderons,
        targets=targets,
        vocab_size=vocab_size,
        iter_logs=iter_logs,
        sort_method=sort_method,
        selection_trace=selection_trace,
        lpjson=lpjson,
    )

    if save_sort_token_file:
        if not save_sort_token_file.endswith(".json"):
            save_sort_token_file = f"{save_sort_token_file}.json"
        save_json(selected_ladderons, save_sort_token_file)

    if show_selectedladderons:
        print("Selected Ladderons:", selected_ladderons.keys())

    return selected_ladderons


def select_token(
    input_lpjson: dict[str, Any] | str | Path,
    vocab_size: int,
    save_sort_token_file: str | None = None,
    show_selectedladderons: bool = False,
    iter_logs: bool = False,
    sort_method: str = SORT_METHOD_CONTRIBUTION_FOLDED,
    selection_trace: list[dict[str, Any]] | None = None,
) -> dict[int, list[Any]]:
    sort_method = validate_sort_method(sort_method)
    lpjson = load_lpjson(input_lpjson, readonly=sort_method != SORT_METHOD_MULTIPLICITY)
    return _select_token_from_loaded_lpjson(
        lpjson=lpjson,
        vocab_size=vocab_size,
        save_sort_token_file=save_sort_token_file,
        show_selectedladderons=show_selectedladderons,
        iter_logs=iter_logs,
        sort_method=sort_method,
        selection_trace=selection_trace,
    )


def select_token_checkpoints(
    input_lpjson: dict[str, Any] | str | Path,
    vocab_sizes: list[int],
    iter_logs: bool = False,
    sort_method: str = SORT_METHOD_CONTRIBUTION_FOLDED,
) -> dict[int, dict[int, list[Any]]]:
    sort_method = validate_sort_method(sort_method)
    lpjson = load_lpjson(input_lpjson, readonly=sort_method != SORT_METHOD_MULTIPLICITY)
    ladderons, targets = lpjson["ladderons"], lpjson["targets"]
    return select_ladderon_checkpoints_by_sort_method(
        ladderons=ladderons,
        targets=targets,
        vocab_sizes=vocab_sizes,
        iter_logs=iter_logs,
        sort_method=sort_method,
        lpjson=lpjson,
    )


def get_vocab_rows(
    lpjson: dict[str, Any] | str | Path,
    vocab_size: int,
    show_selectedladderons: bool = False,
    save_sort_token_lpjson: str | None = None,
    iter_logs: bool = False,
    sort_method: str = SORT_METHOD_CONTRIBUTION_FOLDED,
) -> list[dict[str, Any]]:
    sort_method = validate_sort_method(sort_method)
    normalized_lpjson = load_lpjson(lpjson, readonly=sort_method != SORT_METHOD_MULTIPLICITY)
    selection_trace: list[dict[str, Any]] = []
    selected_ladderons = _select_token_from_loaded_lpjson(
        lpjson=normalized_lpjson,
        vocab_size=vocab_size,
        save_sort_token_file=save_sort_token_lpjson,
        show_selectedladderons=show_selectedladderons,
        iter_logs=iter_logs,
        sort_method=sort_method,
        selection_trace=selection_trace,
    )
    if sort_method == SORT_METHOD_MULTIPLICITY:
        return build_multiplicity_vocab_rows(normalized_lpjson, selected_ladderons, vocab_size)
    return build_vocab_rows(normalized_lpjson, selected_ladderons, selection_trace)


def get(
    lpjson: dict[str, Any] | str | Path,
    vocab_size: int,
    save_vacab_file: str | None = None,
    show_selectedladderons: bool = False,
    save_sort_token_lpjson: str | None = None,
    iter_logs: bool = False,
    sort_method: str = SORT_METHOD_CONTRIBUTION_FOLDED,
) -> dict[str, list[int]]:
    sort_method = validate_sort_method(sort_method)
    normalized_lpjson = load_lpjson(lpjson, readonly=sort_method != SORT_METHOD_MULTIPLICITY)
    selection_trace: list[dict[str, Any]] = []
    sort_lpjson = _select_token_from_loaded_lpjson(
        lpjson=normalized_lpjson,
        vocab_size=vocab_size,
        save_sort_token_file=save_sort_token_lpjson,
        show_selectedladderons=show_selectedladderons,
        iter_logs=iter_logs,
        sort_method=sort_method,
        selection_trace=selection_trace,
    )

    if sort_method == SORT_METHOD_MULTIPLICITY:
        vocab_rows = build_multiplicity_vocab_rows(normalized_lpjson, sort_lpjson, vocab_size)
        output_data = {
            row["token"]: [idx, int(row["ladderon_id"])]
            for idx, row in enumerate(vocab_rows)
        }
    else:
        vocab_rows = build_vocab_rows(normalized_lpjson, sort_lpjson, selection_trace)
        output_data: dict[str, list[int]] = {}
        idx = 0
        for ladderon_id, val in sort_lpjson.items():
            output_data[val[2]] = [idx, ladderon_id]
            idx += 1

        for bbb in normalized_lpjson["basic_building_blocks"]:
            output_data[bbb] = [idx, -idx]
            idx += 1

    if save_vacab_file:
        save_path = Path(save_vacab_file)
        if save_path.suffix.lower() != ".csv":
            save_path = save_path.with_suffix(".csv")
        save_csv_rows(vocab_rows, save_path)

    return output_data


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Ladderpath ladderon selection.")
    parser.add_argument("--input", required=True, help="Path to lpjson JSON file.")
    parser.add_argument("--vocab-size", type=int, required=True, help="Number of ladderons to select.")
    parser.add_argument("--save-vocab-file", help="Path to save final vocab CSV.")
    parser.add_argument("--save-sort-token-lpjson", help="Path to save selected ladderons JSON.")
    parser.add_argument("--show-selectedladderons", action="store_true")
    parser.add_argument(
        "--iter-logs",
        action="store_true",
        help="Enable full per-iteration diagnostics. This is mainly for debugging and can be much slower.",
    )
    parser.add_argument(
        "--sort-method",
        default=SORT_METHOD_CONTRIBUTION_FOLDED,
        choices=sorted(VALID_SORT_METHODS),
        help=(
            "'contribution_folded' uses the folded contribution score; "
            "'contribution_frozen' uses the frozen contribution score; "
            "'multiplicity' ranks ladderons directly by ladderpath multiplicity."
        ),
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    get(
        lpjson=Path(args.input),
        vocab_size=args.vocab_size,
        save_vacab_file=args.save_vocab_file,
        show_selectedladderons=args.show_selectedladderons,
        save_sort_token_lpjson=args.save_sort_token_lpjson,
        iter_logs=args.iter_logs,
        sort_method=args.sort_method,
    )


if __name__ == "__main__":
    main()
