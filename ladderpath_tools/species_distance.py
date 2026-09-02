"""
Version 6.1
Author: codex
Revised for grouped-species distance calculation
Date: 2026.04.09

基于 Ladderpath 理论提供多序列物种分析工具，主要包含：
1. get_ladderpath_species: 按 species_info 将多个 targets 合并为 species-level target。
2. get_distances: 基于 species-level lpjson 计算物种间距离。
3. get_distance_table: 一步式包装接口，同时返回距离表、名称映射和 merged lpjson。

与 v6.0 相比，这一版按当前 duplication 语义修正了 species_info 校验：
- 允许同一个 target ID 出现在不同 species 中；
- 允许同一个 target ID 在同一个 species 中重复出现；
- 校验规则改为：每个 target ID 在 species_info 中的出现次数，
  必须与当前 lpjson 所隐含的重复次数逐个一致。

这与当前 ladderpath.py 的约定保持一致：
- 若输入存在重复，strs 必须是 dict；
- duplications_info 使用 `{target_id: count}` 格式；
- species_info 使用真实 target ID，而不是原始位置编号。
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
    为 lambda 计算准备一个可直接送入 `lambda_from_laddergraph.get()` 的 JSON。

    关键点：
    - targets / ladderons 的 key 转成 int；
    - 补足 target/ladderon 记录长度；
    - species-level duplications_info 直接置空。
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
    将用户传入的编号解析为真实 target ID。

    当前版本只支持真实 target ID，不再支持“原始序列位置编号”。
    """
    target_key_set = {int(k) for k in targets.keys()}
    resolved: list[int] = []

    for raw_tid in input_ids:
        tid = int(raw_tid)
        if tid not in target_key_set:
            raise ValueError(
                f"输入编号 {tid} 不在 lpjson['targets'] 中。"
                "当前版本的 species_info 必须使用真实 target ID，"
                "不再支持原始序列位置编号。"
            )
        resolved.append(tid)

    return resolved


def _expected_target_counts(target_ids: list[int], duplications_info: dict) -> dict[int, int]:
    """
    根据当前 duplications_info 语义，恢复每个 target 在原始输入中的应有出现次数。

    当前 ladderpath 约定：
    - 若某个 target 没出现在 duplications_info 中，则该 target 只出现 1 次；
    - 若某个 target 出现在 duplications_info 中，则其值就是该 target 的总重复次数。
    """
    dup_counts = {int(k): int(v) for k, v in (duplications_info or {}).items()}

    expected_counts: dict[int, int] = {}
    for tid in target_ids:
        count = dup_counts.get(int(tid), 1)
        if count < 1:
            raise ValueError(
                f"duplications_info 中 target ID {tid} 的重复次数非法: {count}"
            )
        expected_counts[int(tid)] = count

    return expected_counts


def _validate_species_partition(
    species_entries: list[dict],
    target_ids: list[int],
    duplications_info: dict,
) -> tuple[dict[int, int], dict[int, int]]:
    """
    检查 species_info 中每个 target ID 的出现次数，是否与 lpjson 隐含的重复次数一致。

    注意：
    - 当前版本允许同一个 target ID 出现在多个 species 中；
    - 也允许同一个 target ID 在同一个 species 中重复出现；
    - 这里按每个 target 的 multiplicity 做逐个校验。
    """
    expected_targets = set(target_ids)
    actual_counts: Counter[int] = Counter()

    for entry in species_entries:
        species_key = entry["key"]
        for tid in entry["resolved_ids"]:
            if tid not in expected_targets:
                raise ValueError(
                    f"{species_key} 中的 target ID {tid} 不在当前 lpjson['targets'] 中。"
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
            "species_info 中各 target ID 的出现次数与 lpjson 隐含的重复次数不一致。 "
            + "; ".join(mismatches)
        )

    return expected_counts, dict(actual_counts)


def get_ladderpath_species(lpjson: dict, species_info: dict, verbose: bool = False) -> dict:
    """
    根据 species_info 把多个 target 合并为 species-level target，并同步更新 ladderons 坐标。

    参数:
        lpjson: 原始 Ladderpath JSON。
        species_info:
            {
                "species1": ("Morado", [-1, -2]),
                "species2": ("Musa_ornata", [-3, -4])
            }
        verbose: 是否打印中间映射信息。

    返回:
        species-level lpjson。

    注意:
    - `species_mapping` 只保存 `新物种 target ID -> 物种名`
    - `duplications_info` 会被显式置空，因为 species-level distance 只基于合并后的 targets/ladderons
    - 当前版本要求 `species_info` 使用真实 target ID
    - 当前版本允许重复 target ID，但要求每个 target ID 的出现次数都与原始输入一致
    """
    processed_json = copy.deepcopy(lpjson)
    targets = processed_json.get("targets", {})
    if not targets:
        raise ValueError("输入 lpjson 缺少 targets。")
    if not species_info:
        raise ValueError("species_info 不能为空。")

    target_ids = sorted((int(k) for k in targets.keys()), key=abs)
    duplications_info = processed_json.get("duplications_info", {})

    if verbose:
        expected_counts = _expected_target_counts(target_ids, duplications_info)
        print(f"唯一 targets 数: {len(targets)}")
        print("分组校验模式: target_multiplicity_check")
        print("可用 target IDs:")
        print(f"  {target_ids}")
        print("各 target 应有出现次数:")
        print(f"  {expected_counts}")

    species_entries = []
    for species_key, payload in species_info.items():
        if not isinstance(payload, (tuple, list)) or len(payload) != 2:
            raise ValueError(
                f"species_info[{species_key!r}] 格式错误，应为 (name, [ids])。"
            )

        species_name, raw_ids = payload
        if not isinstance(raw_ids, (list, tuple)) or len(raw_ids) == 0:
            raise ValueError(
                f"species_info[{species_key!r}] 的序列列表不能为空。"
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
        print("\n物种成员解析:")
        for entry in species_entries:
            print(
                f"  {entry['key']} ({entry['name']}): "
                f"input={entry['raw_ids']} -> resolved_target_ids={entry['resolved_ids']}"
            )
        print("\n成员数校验:")
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
                raise ValueError(f"Target ID {oid} 在原始 targets 中不存在。")

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
        print("\n新物种映射:")
        for new_id, name in species_mapping.items():
            print(f"  {new_id} -> {name}")

    return processed_json


def get_distances(
    lpjson_species: dict,
    method: str = "jaccard",
    verbose: bool = False,
) -> tuple[dict[str, float], dict[int, str]]:
    """
    对 species-level lpjson 计算两两物种距离。

    参数:
        lpjson_species: `get_ladderpath_species()` 的输出。
        method: 只允许 "jaccard" 或 "dice"。
        verbose: 是否打印中间结果。

    返回:
        distances_table:
            形如 {"1,2": 0.62, "1,3": 0.71, ...}
        name_map:
            形如 {1: "Morado", 2: "Musa_ornata"}
    """
    method = method.lower()
    if method not in {"jaccard", "dice"}:
        raise ValueError("method 只允许是 'jaccard' 或 'dice'。")

    species_mapping = lpjson_species.get("species_mapping") or {}
    if not species_mapping:
        raise ValueError("输入 JSON 缺少 species_mapping。请先运行 get_ladderpath_species。")

    ordered_target_ids = list(species_mapping.keys())
    if len(ordered_target_ids) < 2:
        raise ValueError("至少需要两个物种才能计算距离。")

    calc_json = _normalize_calc_json(lpjson_species)

    lambda_single: dict[str, float] = {}
    if verbose:
        print(f"开始计算单物种 Lambda，共 {len(ordered_target_ids)} 个物种。")

    for tid in ordered_target_ids:
        tid_int = int(tid)
        lam = lfl.get([tid_int], calc_json)
        if lam is None:
            raise ValueError(f"无法计算单物种 Lambda: target_id={tid}")
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
                    f"无法计算联合 Lambda: target_id_pair=({tid1}, {tid2})"
                )

            numerator = L12 - (L1 + L2) / 2
            if method == "dice":
                denominator = (L1 + L2) / 2
            else:
                denominator = L12 / 2

            if denominator == 0:
                raise ValueError(
                    f"距离计算分母为 0: pair=({tid1}, {tid2}), method={method}"
                )

            distances_table[f"{label1},{label2}"] = float(numerator / denominator)

    name_map = {
        display_label_map[tid]: str(species_mapping[tid]) for tid in ordered_target_ids
    }

    if verbose:
        print("单物种 Lambda:")
        for tid, label in sorted_targets:
            print(f"  {label}: {name_map[label]} ({tid}) -> {lambda_single[tid]}")
        print("距离表:")
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
    一步式包装接口：
    先按物种分组合并 targets，再输出距离表，并返回 merged lpjson。

    参数:
        lpjson: 原始 Ladderpath JSON。
        species_info: 物种分组信息。当前必须使用真实 target IDs。
        method: 只允许 "jaccard" 或 "dice"，默认是 "jaccard"。
        verbose: 是否打印中间结果。

    返回:
        distance_table, name_map, lpjson_merged
    """
    lpjson_merged = get_ladderpath_species(lpjson, species_info, verbose=verbose)
    distance_table, name_map = get_distances(
        lpjson_merged, method=method, verbose=verbose
    )
    return distance_table, name_map, lpjson_merged
