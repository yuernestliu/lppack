"""
Version 3.0.0
Authors: ecsLab (Jingwen Zhang et al.)
    revised by Xiaojun Hu
    updated by codex
Date: 2026.03.31

Calculate lambda for selected targets (leaves) in an existing laddergraph.

Public constraints:
1. `lpjson` must come from `lp.get_ladderpath()` or
   `lp.load_ladderpath_json()`, not directly from `json.load()`.
2. The current implementation does not support an `lpjson` with duplicates;
   `duplications_info` must be empty.
3. `leaves` must be a non-empty list containing either all integers or all strings.
"""

from .. import ladderpath as lp


def get_target_ID(two_leaf_strs0, lpjson):

    strs_subset = two_leaf_strs0.copy()
    ids = []
    for ID, info in lpjson["targets"].items():
        if info[2] in strs_subset:
            ids.append(ID)
            strs_subset.remove(info[2])
    if len(ids) != len(two_leaf_strs0):
        raise ValueError(
            "Not all target strings were found in lpjson['targets']."
        )
    return ids

def separate_basic_blocks(COMP):

    new_COMP=[]
    block_count=0
    for i in COMP:
        if type(i)==str:
            block_count+=len(i)
        else:
            new_COMP.append(i)
    return new_COMP,block_count

def if_id_in_COMP(COMP):

    for i in COMP:
        if type(i) == int:
            return True
    return False


def _validate_lpjson_ready(lpjson):
    if not isinstance(lpjson, dict):
        raise ValueError("lpjson must be a dictionary.")

    if "targets" not in lpjson or "ladderons" not in lpjson:
        raise ValueError("lpjson is missing the 'targets' or 'ladderons' field.")

    if not lpjson["targets"]:
        raise ValueError("lpjson['targets'] must not be empty.")

    target_key = next(iter(lpjson["targets"].keys()))
    ladderon_key = next(iter(lpjson["ladderons"].keys())) if lpjson["ladderons"] else None

    if not isinstance(target_key, int) or (
        ladderon_key is not None and not isinstance(ladderon_key, int)
    ):
        raise ValueError(
            "lpjson contains keys of the wrong type. "
            "Use lp.get_ladderpath(...) or lp.load_ladderpath_json(...) to obtain lpjson; "
            "do not use json.load() directly."
        )

    if lpjson.get("duplications_info"):
        raise ValueError(
            "lambda_from_laddergraph.get() does not currently support an lpjson "
            "with duplicates. Deduplicate the input or use an object with empty "
            "duplications_info."
        )


def _validate_leaves(leaves):
    if not isinstance(leaves, list) or len(leaves) == 0:
        raise ValueError("leaves must be a non-empty list.")

    first = leaves[0]
    if isinstance(first, int):
        if not all(isinstance(x, int) for x in leaves):
            raise ValueError("All leaves must have the same type; do not mix integers and strings.")
    elif isinstance(first, str):
        if not all(isinstance(x, str) for x in leaves):
            raise ValueError("All leaves must have the same type; do not mix integers and strings.")
    else:
        raise ValueError("leaves must contain only integers or only strings.")

    if len(set(leaves)) != len(leaves):
        raise ValueError(
            "lambda_from_laddergraph.get() does not support repeated leaves. "
            "Deduplicate leaves before calculating lambda."
        )

def get(leaves, lpjson):
    _validate_leaves(leaves)
    _validate_lpjson_ready(lpjson)

    lp.fill_lpjson_STR(lpjson)


    if isinstance(leaves[0], int):
        for leaf in leaves:
            if leaf not in lpjson['targets']:
                raise ValueError(f"Target ID {leaf} is not present in lpjson['targets'].")
        ids = leaves
    else:
        ids = get_target_ID(leaves.copy(), lpjson)

    basic_block_num = 0
    comp = []
    ladderon_dic = {}
    for i in ids:
        comp += lpjson["targets"][i][0]
    pool = set()

    while if_id_in_COMP(comp):
        ids_temp, basic_blocks = separate_basic_blocks(comp)
        comp = []
        basic_block_num += basic_blocks
        id_set = set()
        for i in ids_temp:
            if i in ladderon_dic:
                ladderon_dic[i] += 1
            else:
                ladderon_dic[i] = 1
            if i not in id_set and i not in pool:
                comp += lpjson["ladderons"][i][0]
                pool.add(i)
                id_set.add(i)

    for i in comp:
        basic_block_num += len(i)
    count = basic_block_num

    for ladderon, times in ladderon_dic.items():
        count += (times - 1)
    return count
