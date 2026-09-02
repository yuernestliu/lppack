"""
Version 1.4
Author: codex
Date: 2026.04.09

Standalone phylogeny builder.

This module only does four things:
- accept a precomputed `distance_table`
- accept a precomputed `name_map`
- build a phylogenetic tree with UPGMA or NJ
- optionally reroot the tree with an outgroup
- optionally save a `.newick` / `.nwk` file
- optionally save a tree figure file

Example:

    import ladderpath_tools.phylo as lp_phylo

    lp_phylo.phylo(
        distance_table,
        name_map,
        method="nj",
        newick_file_name="my_tree.nwk",
        fig_file_name="my_tree.png",
    )
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from io import StringIO
from math import isfinite
import os

__all__ = ["phylo"]


def _normalize_name_map(name_map: Mapping) -> dict[int, str]:
    if not isinstance(name_map, Mapping) or not name_map:
        raise ValueError("name_map must be a non-empty mapping.")

    normalized: dict[int, str] = {}
    seen_names: set[str] = set()

    for raw_label, raw_name in name_map.items():
        try:
            label = int(raw_label)
        except Exception as exc:
            raise ValueError(f"name_map key {raw_label!r} cannot be converted to an integer label.") from exc

        name = str(raw_name).strip()
        if not name:
            raise ValueError(f"name_map[{raw_label!r}] is empty.")
        if label in normalized:
            raise ValueError(f"Duplicate label found in name_map: {label}")
        if name in seen_names:
            raise ValueError(f"Duplicate species name found in name_map: {name!r}")

        normalized[label] = name
        seen_names.add(name)

    if len(normalized) < 2:
        raise ValueError("At least two species are required to build a phylogenetic tree.")

    return normalized


def _parse_pair_key(key) -> tuple[int, int]:
    if isinstance(key, str):
        parts = [part.strip() for part in key.split(",")]
        if len(parts) != 2:
            raise ValueError(f"Invalid distance_table key {key!r}; expected '1,2'.")
        left, right = parts
    elif isinstance(key, Sequence) and not isinstance(key, (str, bytes)) and len(key) == 2:
        left, right = key
    else:
        raise ValueError(
            f"Invalid distance_table key {key!r}; expected '1,2' or a tuple/list of length 2."
        )

    try:
        a = int(left)
        b = int(right)
    except Exception as exc:
        raise ValueError(f"Labels in distance_table key {key!r} cannot be converted to integers.") from exc

    if a == b:
        raise ValueError(f"Invalid distance_table key {key!r}; self-distances are not allowed.")

    return (a, b) if a < b else (b, a)


def _normalize_distance_table(
    distance_table: Mapping,
    normalized_name_map: dict[int, str],
) -> tuple[list[int], dict[tuple[int, int], float]]:
    if not isinstance(distance_table, Mapping) or not distance_table:
        raise ValueError("distance_table must be a non-empty mapping.")

    sorted_labels = sorted(normalized_name_map.keys())
    label_set = set(sorted_labels)
    distance_lookup: dict[tuple[int, int], float] = {}
    labels_in_table: set[int] = set()

    for raw_key, raw_value in distance_table.items():
        pair = _parse_pair_key(raw_key)
        if pair in distance_lookup:
            raise ValueError(f"Duplicate species pair found in distance_table: {pair}")

        try:
            distance = float(raw_value)
        except Exception as exc:
            raise ValueError(
                f"distance_table[{raw_key!r}]={raw_value!r} cannot be converted to float."
            ) from exc

        if not isfinite(distance):
            raise ValueError(f"distance_table[{raw_key!r}] is not a finite number.")
        if distance < 0:
            raise ValueError(f"distance_table[{raw_key!r}] is negative; distances must be >= 0.")

        if pair[0] not in label_set or pair[1] not in label_set:
            raise ValueError(
                f"Species pair {pair} in distance_table is not fully defined in name_map."
            )

        distance_lookup[pair] = distance
        labels_in_table.update(pair)

    missing_pairs = []
    for i in range(len(sorted_labels)):
        for j in range(i + 1, len(sorted_labels)):
            pair = (sorted_labels[i], sorted_labels[j])
            if pair not in distance_lookup:
                missing_pairs.append(f"{pair[0]},{pair[1]}")

    if missing_pairs:
        raise ValueError(
            "distance_table does not cover all pairwise species distances. Missing: "
            + ", ".join(missing_pairs)
        )

    missing_labels = sorted(label_set - labels_in_table)
    if missing_labels:
        raise ValueError(
            "Some species never appear in distance_table, so the tree cannot be built. Missing labels: "
            + ", ".join(str(label) for label in missing_labels)
        )

    return sorted_labels, distance_lookup


def _build_condensed_distance_matrix(
    sorted_labels: list[int],
    distance_lookup: dict[tuple[int, int], float],
) -> list[float]:
    condensed = []
    for i in range(len(sorted_labels)):
        for j in range(i + 1, len(sorted_labels)):
            condensed.append(distance_lookup[(sorted_labels[i], sorted_labels[j])])
    return condensed


def _build_lower_triangular_distance_matrix(
    sorted_labels: list[int],
    distance_lookup: dict[tuple[int, int], float],
) -> list[list[float]]:
    lower_triangular = []
    for i in range(len(sorted_labels)):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                row.append(distance_lookup[(sorted_labels[j], sorted_labels[i])])
        lower_triangular.append(row)
    return lower_triangular


def _scipy_tree_to_newick(tree_node, leaf_names: list[str], parent_dist: float) -> str:
    branch_length = max(0.0, parent_dist - tree_node.dist)
    if tree_node.is_leaf():
        return f"{leaf_names[tree_node.id]}:{branch_length}"

    left = _scipy_tree_to_newick(tree_node.get_left(), leaf_names, tree_node.dist)
    right = _scipy_tree_to_newick(tree_node.get_right(), leaf_names, tree_node.dist)
    return f"({left},{right}):{branch_length}"


def _normalize_newick_file_name(newick_file_name: str) -> str:
    file_name = str(newick_file_name).strip()
    if not file_name:
        raise ValueError("newick_file_name must not be an empty string.")
    if file_name.endswith((".newick", ".nwk")):
        return file_name
    return f"{file_name}.nwk"


def _normalize_figure_file_name(fig_file_name: str) -> str:
    file_name = str(fig_file_name).strip()
    if not file_name:
        raise ValueError("fig_file_name must not be an empty string.")

    root, ext = os.path.splitext(file_name)
    if ext.lower() in {".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tif", ".tiff"}:
        return file_name
    return f"{file_name}.png"


def _clean_tiny_branch_lengths(tree, atol: float = 1e-12) -> None:
    # 清理浮点误差导致的 -0.0 / 极小负数分支长度
    for clade in tree.find_clades():
        branch_length = clade.branch_length
        if branch_length is not None and abs(branch_length) < atol:
            clade.branch_length = 0.0


def _resolve_outgroup_terminal_name(
    tree,
    outgroup_name,
    normalized_name_map: dict[int, str],
) -> str:
    candidate = str(outgroup_name).strip()
    if not candidate:
        raise ValueError("outgroup_name must not be an empty string.")

    terminal_names = [
        clade.name for clade in tree.find_clades() if clade.is_terminal() and clade.name
    ]
    if candidate in terminal_names:
        return candidate

    lower_candidate = candidate.lower()
    if lower_candidate.startswith("species"):
        suffix = lower_candidate[len("species") :]
        if suffix.isdigit():
            label = int(suffix)
            if label in normalized_name_map:
                return normalized_name_map[label]

    if candidate.isdigit():
        label = int(candidate)
        if label in normalized_name_map:
            return normalized_name_map[label]

    raise ValueError(
        f"Outgroup {candidate!r} was not found in the tree. "
        f"You can use a terminal name like {terminal_names} "
        f"or a species alias like {[f'species{label}' for label in sorted(normalized_name_map)]}."
    )


def _reroot_tree_with_outgroup(
    tree,
    outgroup_name,
    normalized_name_map: dict[int, str],
) -> str:
    resolved_name = _resolve_outgroup_terminal_name(
        tree, outgroup_name, normalized_name_map
    )
    outgroup = tree.find_any(name=resolved_name)
    if outgroup is None:
        raise ValueError(f"Resolved outgroup {resolved_name!r} was not found in the tree.")

    tree.root_with_outgroup(outgroup)
    tree.rooted = True
    return resolved_name


def phylo(
    distance_table: Mapping,
    name_map: Mapping,
    method: str = "nj",
    outgroup_name: str | None = None,
    newick_file_name: str | None = None,
    fig_file_name: str | None = None,
    show_plot: bool = True,
    verbose: bool = True,
):
    """
    Build and draw a phylogenetic tree from a precomputed distance table.

    Parameters:
        distance_table:
            For example {"1,2": 0.62, "1,3": 0.71, ...}
            Tuple/list keys are also supported: {(1, 2): 0.62, (1, 3): 0.71, ...}
        name_map:
            For example {1: "Morado", 2: "Musa_ornata"}
        method:
            Must be either "upgma" or "nj"
        outgroup_name:
            If provided, reroots the tree by the given outgroup name.
            This is especially useful for rooting an NJ tree, which is unrooted by default.
        newick_file_name:
            If provided, saves a Newick file; supports .newick / .nwk
        fig_file_name:
            If provided, saves a figure file; supports .png / .pdf / .svg / .jpg / .jpeg / .tif / .tiff
            If no suffix is given, ".png" is appended automatically.
        show_plot:
            Whether to display the tree immediately. In Jupyter Notebook, keep this as True
            if you want to see the tree inline right after running the cell.
        verbose:
            Whether to print progress messages

    Returns:
        A Biopython tree object
    """
    method = str(method).lower().strip()
    if method not in {"upgma", "nj"}:
        raise ValueError("method must be either 'upgma' or 'nj'.")

    normalized_name_map = _normalize_name_map(name_map)
    sorted_labels, distance_lookup = _normalize_distance_table(
        distance_table, normalized_name_map
    )
    species_names = [normalized_name_map[label] for label in sorted_labels]

    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import scipy.cluster.hierarchy as sch
        from Bio import Phylo
        from Bio.Phylo import TreeConstruction
    except ImportError as exc:
        raise ImportError(
            "phylo() requires matplotlib, scipy, and Biopython (`Bio`). "
            "Install it with: pip install biopython"
        ) from exc

    if verbose:
        print(f"Building phylogenetic tree with method={method.upper()} for {len(species_names)} species.")

    if method == "upgma":
        condensed_mat = _build_condensed_distance_matrix(sorted_labels, distance_lookup)
        linkage_matrix = sch.linkage(condensed_mat, method="average", optimal_ordering=True)
        scipy_tree, _ = sch.to_tree(linkage_matrix, rd=True)
        newick_str = _scipy_tree_to_newick(scipy_tree, species_names, scipy_tree.dist) + ";"
        tree = Phylo.read(StringIO(newick_str), "newick")
        tree.rooted = True
    else:
        lower_tri_mat = _build_lower_triangular_distance_matrix(
            sorted_labels, distance_lookup
        )
        distance_matrix = TreeConstruction.DistanceMatrix(species_names, lower_tri_mat)
        tree = TreeConstruction.DistanceTreeConstructor().nj(distance_matrix)
        for clade in tree.find_clades():
            if not clade.is_terminal():
                clade.name = None
        tree.rooted = False

    resolved_outgroup_name = None
    if outgroup_name is not None:
        resolved_outgroup_name = _reroot_tree_with_outgroup(
            tree, outgroup_name, normalized_name_map
        )

    if verbose:
        if outgroup_name is not None:
            print(
                "Root status: rerooted with "
                f"outgroup={outgroup_name!r} (resolved to {resolved_outgroup_name!r})."
            )
        elif method == "upgma":
            print("Root status: rooted by default (UPGMA).")
        else:
            print("Root status: unrooted by default (NJ).")

    _clean_tiny_branch_lengths(tree)
    newick_buffer = StringIO()
    Phylo.write(tree, newick_buffer, "newick", format_branch_length="%s")
    newick_str = newick_buffer.getvalue()

    if newick_file_name is not None:
        normalized_file_name = _normalize_newick_file_name(newick_file_name)
        with open(normalized_file_name, "w", encoding="utf-8") as file:
            file.write(newick_str)
        if verbose:
            print(f"Newick file saved to: {normalized_file_name}")

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    Phylo.draw(tree, do_show=False, axes=ax)
    tree_state = "Rooted" if tree.rooted else "Unrooted"
    ax.set_title(f"Phylogeny Tree ({method.upper()}, {tree_state})")
    plt.tight_layout()

    if fig_file_name is not None:
        normalized_fig_file_name = _normalize_figure_file_name(fig_file_name)
        fig.savefig(normalized_fig_file_name, bbox_inches="tight")
        if verbose:
            print(f"Figure file saved to: {normalized_fig_file_name}")

    if show_plot:
        # 在纯 Agg 后端下 plt.show() 会产生无意义警告；Notebook inline/backend 正常显示
        backend = matplotlib.get_backend().lower()
        if backend != "agg":
            plt.show()
    plt.close(fig)

    return tree
