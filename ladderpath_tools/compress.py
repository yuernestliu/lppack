#!/usr/bin/env python3

"""
Version 3.0.0
Authors: ecsLab & Codex
Date: 26y0508

Compression entrypoint for `lppack.ladderpath_tools.compress`.

Main policy:
    - Keep only the latest mainline compression path here.
    - Keep the public API short and readable.
    - Keep experimental model families outside the public entry point.

Public API:
    - encoding2streams(lpjson, ...) -> encoded_streams
    - entropy_coding(encoded_streams, ...) -> bits(bytes)
    - compress2bit(lpjson, ...) -> bits(bytes)
"""

from __future__ import annotations

from collections import Counter
import heapq
import os
import sys
import types
from pathlib import Path
from typing import Any


DEFAULT_SEP = "@"
DEFAULT_REINDEX = True

SORT_METHOD_CONTRIBUTION_FOLDED = "contribution_folded"
SORT_METHOD_CONTRIBUTION_FROZEN = "contribution_frozen"
SORT_METHOD_MULTIPLICITY = "multiplicity"
VALID_SORT_METHODS = {
    SORT_METHOD_CONTRIBUTION_FOLDED,
    SORT_METHOD_CONTRIBUTION_FROZEN,
    SORT_METHOD_MULTIPLICITY,
}
DEFAULT_SORT_METHOD = SORT_METHOD_CONTRIBUTION_FOLDED

TOKEN_TYPE_MODE_BITPACK = "bitpack"
TOKEN_TYPE_MODE_PREV1_CONTEXT = "prev1_context"
TOKEN_TYPE_MODE_AUTO_BEST = "auto_best"
VALID_TOKEN_TYPE_MODES = {
    TOKEN_TYPE_MODE_BITPACK,
    TOKEN_TYPE_MODE_PREV1_CONTEXT,
    TOKEN_TYPE_MODE_AUTO_BEST,
}
DEFAULT_TOKEN_TYPE_MODE = TOKEN_TYPE_MODE_PREV1_CONTEXT

MAGIC = b"LPE1"
LPE_FORMAT_VERSION = 1
STREAM_BITPACK = 0
STREAM_HUFFMAN = 1
STREAM_DENSE_HUFFMAN = 2
STREAM_REF_BUCKET_POWER2_DENSE = 5
STREAM_ARITHMETIC = 6
STREAM_TOKEN_PREV1_CONTEXT = 13

ARITH_MAX = (1 << 32) - 1
ARITH_HALF = 1 << 31
ARITH_QUARTER = 1 << 30
ARITH_THREE_QUARTER = 3 << 30


def _append_optional_repo_root() -> None:
    repo_root = os.environ.get("LPPACK_REPO_ROOT")
    if repo_root and repo_root not in sys.path:
        sys.path.append(repo_root)


def import_ladderpath():
    try:
        from .. import ladderpath as lp  # type: ignore

        return lp
    except Exception:
        _append_optional_repo_root()
        sys.modules.setdefault("graphviz", types.SimpleNamespace())
        try:
            import lppack.ladderpath as lp
        except Exception as exc:
            raise ImportError(
                "Could not import lppack.ladderpath. Install lppack, "
                "run this module inside the lppack package, or set "
                "LPPACK_REPO_ROOT to the repository parent directory."
            ) from exc

        return lp


def import_sortout():
    try:
        from . import sortout as lp_sortout  # type: ignore

        return lp_sortout
    except Exception:
        _append_optional_repo_root()
        sys.modules.setdefault("graphviz", types.SimpleNamespace())
        try:
            import lppack.ladderpath_tools.sortout as lp_sortout
        except Exception as exc:
            raise ImportError(
                "Could not import lppack.ladderpath_tools.sortout. Install "
                "lppack, run this module inside the lppack package, "
                "or set LPPACK_REPO_ROOT to the repository parent directory."
            ) from exc

        return lp_sortout


def load_lpjson(lpjson_or_path: dict[str, Any] | str | Path) -> dict[str, Any]:
    lp = import_ladderpath()

    if isinstance(lpjson_or_path, dict):
        lpjson = lpjson_or_path
        lp.fill_lpjson_STR(lpjson)
        return lpjson

    lpjson = lp.load_ladderpath_json(str(lpjson_or_path))
    lp.fill_lpjson_STR(lpjson)
    return lpjson


def _write_payload_if_requested(payload: bytes, output_path: str | Path | None) -> None:
    if output_path is None:
        return
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _validate_token_type_mode(token_type_mode: str) -> str:
    if token_type_mode not in VALID_TOKEN_TYPE_MODES:
        raise ValueError(
            f"Unsupported token_type_mode={token_type_mode!r}. "
            f"Expected one of {sorted(VALID_TOKEN_TYPE_MODES)!r}."
        )
    return token_type_mode


def _validate_sort_method(sort_method: str) -> str:
    if sort_method not in VALID_SORT_METHODS:
        raise ValueError(
            f"Unsupported sort_method={sort_method!r}. "
            f"Expected one of {sorted(VALID_SORT_METHODS)!r}."
        )
    return sort_method


def _resolve_selected_ladderon_ids(
    lpjson: dict[str, Any],
    *,
    top_n: int | None,
    sort_method: str,
    selected_ladderon_ids: list[int] | None,
) -> list[int] | None:
    if top_n is not None and selected_ladderon_ids is not None:
        raise ValueError("Pass either top_n or selected_ladderon_ids, not both.")

    if selected_ladderon_ids is not None:
        return list(selected_ladderon_ids)

    if top_n is None:
        return None

    lp_sortout = import_sortout()
    selected = lp_sortout.select_token(
        input_lpjson=lpjson,
        vocab_size=int(top_n),
        show_selectedladderons=False,
        iter_logs=False,
        sort_method=_validate_sort_method(sort_method),
    )
    return list(selected.keys())


def _append_token(out: list[Any], token: Any, sep: str) -> None:
    if isinstance(token, str) and token == sep:
        raise ValueError(
            f"Literal token conflicts with separator {sep!r}. "
            "Choose a different separator."
        )
    if isinstance(token, str) and out and isinstance(out[-1], str) and out[-1] != sep:
        out[-1] += token
        return
    out.append(token)


def _expand_item_to_tokens(
    item: Any,
    ladderons: dict[int, list[Any]],
    selected_set: set[int],
    local_index: dict[int, int],
    sep: str,
    memo: dict[int, tuple[Any, ...]],
) -> list[Any]:
    if isinstance(item, str):
        return [item]

    if item in selected_set:
        return [local_index[item]]

    if item in memo:
        return list(memo[item])

    if item not in ladderons:
        raise KeyError(f"Missing ladderon_id={item} while expanding tokens.")

    tokens: list[Any] = []
    for child in ladderons[item][0]:
        for token in _expand_item_to_tokens(child, ladderons, selected_set, local_index, sep, memo):
            _append_token(tokens, token, sep)

    memo[item] = tuple(tokens)
    return list(tokens)


def _encode_component_sequence(
    components: list[Any],
    ladderons: dict[int, list[Any]],
    selected_set: set[int],
    local_index: dict[int, int],
    sep: str,
    memo: dict[int, tuple[Any, ...]],
) -> list[Any]:
    tokens: list[Any] = []
    for item in components:
        for token in _expand_item_to_tokens(item, ladderons, selected_set, local_index, sep, memo):
            _append_token(tokens, token, sep)
    return tokens


def _compress_lpjson(
    lpjson: dict[str, Any],
    *,
    sep: str,
    selected_ladderon_ids: list[int] | None,
) -> list[Any]:
    input_type = lpjson.get("input_type")
    if input_type is not None and input_type != "list":
        raise ValueError("The input_type for targets must be 'list' to use structural compression.")

    lp = import_ladderpath()
    targets = lpjson["targets"]
    ladderons = lpjson["ladderons"]

    if selected_ladderon_ids is None:
        included_ladderons: list[int] = []
        for info_list in targets.values():
            for item in info_list[0]:
                if isinstance(item, int) and item not in included_ladderons:
                    included_ladderons.append(item)
        i = 0
        while i < len(included_ladderons):
            for item in ladderons[included_ladderons[i]][0]:
                if isinstance(item, int) and item not in included_ladderons:
                    included_ladderons.append(item)
            i += 1
        included_ladderons.sort(reverse=True)
    else:
        included_ladderons = [ladderon_id for ladderon_id in selected_ladderon_ids if ladderon_id in ladderons]

    local_index = {ladderon_id: i for i, ladderon_id in enumerate(included_ladderons)}
    selected_set = set(included_ladderons)
    memo: dict[int, tuple[Any, ...]] = {}

    target_payload: list[Any] = []
    targets_list = lp.reconstruct_targets_list(lpjson)
    for target_id in targets_list:
        target_payload.append(sep)
        encoded = _encode_component_sequence(
            targets[target_id][0],
            ladderons,
            selected_set,
            local_index,
            sep,
            memo,
        )
        for token in encoded:
            _append_token(target_payload, token, sep)

    ladderon_payload: list[Any] = []
    for ladderon_id in included_ladderons:
        ladderon_payload.append(sep)
        encoded = _encode_component_sequence(
            ladderons[ladderon_id][0],
            ladderons,
            selected_set,
            local_index,
            sep,
            memo,
        )
        for token in encoded:
            _append_token(ladderon_payload, token, sep)

    return [len(targets_list), *target_payload, *ladderon_payload]


def _split_segments(compressed_list: list[Any]) -> tuple[int, str, list[list[Any]], list[list[Any]]]:
    if len(compressed_list) < 2:
        raise ValueError("Compressed list is too short.")

    n_targets = compressed_list[0]
    sep = compressed_list[1]
    if not isinstance(n_targets, int):
        raise TypeError("The first element of compressed_list must be an int target count.")
    if not isinstance(sep, str):
        raise TypeError("The second element of compressed_list must be a string separator.")

    segments: list[list[Any]] = []
    current: list[Any] = []
    for item in compressed_list[1:]:
        if item == sep:
            segments.append(current)
            current = []
        else:
            current.append(item)
    segments.append(current)
    segments = segments[1:]

    if len(segments) < n_targets:
        raise ValueError("Compressed list does not contain enough target segments.")

    return n_targets, sep, segments[:n_targets], segments[n_targets:]


def _build_compressed_list(
    n_targets: int,
    sep: str,
    target_segments: list[list[Any]],
    ladderon_segments: list[list[Any]],
) -> list[Any]:
    compressed = [n_targets]
    for segment in [*target_segments, *ladderon_segments]:
        compressed.append(sep)
        compressed.extend(segment)
    return compressed


def _flatten_segment_family_for_entropy(segments: list[list[Any]]) -> dict[str, Any]:
    token_counts: list[int] = []
    ref_counts: list[int] = []
    literal_counts: list[int] = []
    literal_byte_counts: list[int] = []
    token_types: list[int] = []
    ref_ids: list[int] = []
    literal_lengths: list[int] = []
    literal_bytes: list[int] = []

    for segment in segments:
        seg_ref_count = 0
        seg_literal_count = 0
        seg_literal_bytes = 0
        token_counts.append(len(segment))

        for token in segment:
            if isinstance(token, int):
                token_types.append(1)
                ref_ids.append(token)
                seg_ref_count += 1
            elif isinstance(token, str):
                token_types.append(0)
                token_bytes = token.encode("utf-8")
                literal_lengths.append(len(token_bytes))
                literal_bytes.extend(token_bytes)
                seg_literal_count += 1
                seg_literal_bytes += len(token_bytes)
            else:
                raise TypeError(f"Unsupported token type in segment: {type(token)!r}")

        ref_counts.append(seg_ref_count)
        literal_counts.append(seg_literal_count)
        literal_byte_counts.append(seg_literal_bytes)

    return {
        "segment_count": len(segments),
        "segment_token_counts": token_counts,
        "segment_ref_counts": ref_counts,
        "segment_literal_counts": literal_counts,
        "segment_literal_byte_counts": literal_byte_counts,
        "token_types": token_types,
        "ref_ids": ref_ids,
        "literal_lengths": literal_lengths,
        "literal_bytes": literal_bytes,
    }


def _build_entropy_ready_json(compressed_list: list[Any]) -> dict[str, Any]:
    n_targets, sep, target_segments, ladderon_segments = _split_segments(compressed_list)
    target_streams = _flatten_segment_family_for_entropy(target_segments)
    dictionary_streams = _flatten_segment_family_for_entropy(ladderon_segments)

    return {
        "format": "lp-entropy-input-v1",
        "description": "Lightweight entropy-ready export for the mainline compress.py.",
        "n_targets": n_targets,
        "n_dictionary_rules": len(ladderon_segments),
        "sep": sep,
        "target": target_streams,
        "dictionary": {
            "rule_lengths": dictionary_streams["segment_token_counts"],
            "rule_ref_counts": dictionary_streams["segment_ref_counts"],
            "rule_literal_counts": dictionary_streams["segment_literal_counts"],
            "rule_literal_byte_counts": dictionary_streams["segment_literal_byte_counts"],
            "token_types": dictionary_streams["token_types"],
            "ref_ids": dictionary_streams["ref_ids"],
            "literal_lengths": dictionary_streams["literal_lengths"],
            "literal_bytes": dictionary_streams["literal_bytes"],
        },
    }


def _collect_direct_ref_counts(
    target_segments: list[list[Any]],
    ladderon_segments: list[list[Any]],
) -> tuple[list[int], list[int]]:
    n_rules = len(ladderon_segments)
    target_direct_counts = [0] * n_rules
    dict_direct_counts = [0] * n_rules

    for segment in target_segments:
        for token in segment:
            if isinstance(token, int):
                target_direct_counts[token] += 1

    for segment in ladderon_segments:
        for token in segment:
            if isinstance(token, int):
                dict_direct_counts[token] += 1

    return target_direct_counts, dict_direct_counts


def _reindex_compressed_list_latest(compressed_list: list[Any]) -> tuple[list[Any], dict[str, Any]]:
    n_targets, sep, target_segments, ladderon_segments = _split_segments(compressed_list)
    n_rules = len(ladderon_segments)
    target_direct_counts, dict_direct_counts = _collect_direct_ref_counts(target_segments, ladderon_segments)

    rule_order = sorted(
        range(n_rules),
        key=lambda old_id: (
            -(target_direct_counts[old_id] + 2.0 * dict_direct_counts[old_id]),
            -target_direct_counts[old_id],
            -dict_direct_counts[old_id],
            old_id,
        ),
    )
    remap = {old_id: new_id for new_id, old_id in enumerate(rule_order)}

    def remap_segment(segment: list[Any]) -> list[Any]:
        out: list[Any] = []
        for token in segment:
            out.append(remap[token] if isinstance(token, int) else token)
        return out

    reordered_ladderon_segments = [remap_segment(ladderon_segments[old_id]) for old_id in rule_order]
    remapped_target_segments = [remap_segment(segment) for segment in target_segments]
    reindexed = _build_compressed_list(n_targets, sep, remapped_target_segments, reordered_ladderon_segments)

    metadata = {
        "name": "ref_usage_bucket_power2_v2",
        "description": "Latest mainline ref-id reindex for bucket_power2_dense coding.",
        "n_rules": n_rules,
        "rule_order_old_to_new": rule_order,
        "target_direct_counts_old_order": target_direct_counts,
        "dict_direct_counts_old_order": dict_direct_counts,
    }
    return reindexed, metadata


def encode_varint(value: int) -> bytes:
    if value < 0:
        raise ValueError("Varint only supports non-negative integers.")

    out = bytearray()
    while True:
        to_write = value & 0x7F
        value >>= 7
        if value:
            out.append(to_write | 0x80)
        else:
            out.append(to_write)
            return bytes(out)


def build_huffman_code_lengths(freqs: dict[int, int]) -> dict[int, int]:
    if not freqs:
        return {}
    if len(freqs) == 1:
        symbol = next(iter(freqs))
        return {symbol: 1}

    heap: list[tuple[int, int, int, Any]] = []
    next_id = 0
    for symbol, freq in sorted(freqs.items()):
        heapq.heappush(heap, (freq, symbol, next_id, ("leaf", symbol)))
        next_id += 1

    while len(heap) > 1:
        freq_a, min_a, _id_a, node_a = heapq.heappop(heap)
        freq_b, min_b, _id_b, node_b = heapq.heappop(heap)
        merged = ("node", node_a, node_b)
        heapq.heappush(heap, (freq_a + freq_b, min(min_a, min_b), next_id, merged))
        next_id += 1

    root = heap[0][3]
    lengths: dict[int, int] = {}

    def walk(node: Any, depth: int) -> None:
        if node[0] == "leaf":
            lengths[node[1]] = max(depth, 1)
            return
        walk(node[1], depth + 1)
        walk(node[2], depth + 1)

    walk(root, 0)
    return lengths


def build_canonical_huffman_codes(lengths: dict[int, int]) -> dict[int, tuple[int, int]]:
    codes: dict[int, tuple[int, int]] = {}
    code = 0
    prev_len = 0
    for bit_length, symbol in sorted((length, symbol) for symbol, length in lengths.items()):
        code <<= bit_length - prev_len
        codes[symbol] = (code, bit_length)
        code += 1
        prev_len = bit_length
    return codes


codec = types.SimpleNamespace(
    Counter=Counter,
    encode_varint=encode_varint,
    build_huffman_code_lengths=build_huffman_code_lengths,
    build_canonical_huffman_codes=build_canonical_huffman_codes,
)


def pack_bits(values: list[int]) -> bytes:
    out = bytearray((len(values) + 7) // 8)
    for idx, value in enumerate(values):
        if value not in {0, 1}:
            raise ValueError(f"Bit-packed stream only supports 0/1 values, got {value!r}")
        if value:
            out[idx // 8] |= 1 << (7 - (idx % 8))
    return bytes(out)


def pack_nibbles(values: list[int]) -> bytes:
    out = bytearray((len(values) + 1) // 2)
    for idx, value in enumerate(values):
        if not 0 <= value <= 0x0F:
            raise ValueError(f"Nibble-packed stream only supports 0..15 values, got {value!r}")
        byte_idx = idx // 2
        if idx % 2 == 0:
            out[byte_idx] = value << 4
        else:
            out[byte_idx] |= value
    return bytes(out)


def pack_variable_width_bits(chunks: list[tuple[int, int]]) -> bytes:
    bit_buffer = 0
    bit_count = 0
    out = bytearray()
    for value, width in chunks:
        if width < 0:
            raise ValueError("Bit width must be non-negative.")
        if width == 0:
            continue
        if value < 0 or value >= (1 << width):
            raise ValueError(f"Value {value} does not fit into {width} bits.")
        bit_buffer = (bit_buffer << width) | value
        bit_count += width
        while bit_count >= 8:
            shift = bit_count - 8
            out.append((bit_buffer >> shift) & 0xFF)
            bit_buffer &= (1 << shift) - 1 if shift else 0
            bit_count = shift
    if bit_count:
        out.append((bit_buffer << (8 - bit_count)) & 0xFF)
    return bytes(out)


class BitWriter:
    def __init__(self) -> None:
        self._buffer = 0
        self._count = 0
        self._out = bytearray()

    def write(self, bit: int) -> None:
        self._buffer = (self._buffer << 1) | (bit & 1)
        self._count += 1
        if self._count == 8:
            self._out.append(self._buffer & 0xFF)
            self._buffer = 0
            self._count = 0

    def to_bytes(self) -> bytes:
        if self._count:
            self._out.append((self._buffer << (8 - self._count)) & 0xFF)
            self._buffer = 0
            self._count = 0
        return bytes(self._out)


def build_arithmetic_model(freqs: dict[int, int]) -> tuple[list[int], list[int], int]:
    if any(freq <= 0 for freq in freqs.values()):
        raise ValueError("Arithmetic model frequencies must be positive.")
    symbols = sorted(freqs)
    cumulative = [0]
    total = 0
    for symbol in symbols:
        total += freqs[symbol]
        cumulative.append(total)
    return symbols, cumulative, total


def encode_arithmetic_values(values: list[int], freqs: dict[int, int]) -> bytes:
    if not values:
        return b""

    symbols, cumulative, total = build_arithmetic_model(freqs)
    bounds = {symbol: (cumulative[idx], cumulative[idx + 1]) for idx, symbol in enumerate(symbols)}

    low = 0
    high = ARITH_MAX
    pending_bits = 0
    writer = BitWriter()

    def emit_with_follow(bit: int) -> None:
        nonlocal pending_bits
        writer.write(bit)
        follow = 1 - bit
        for _ in range(pending_bits):
            writer.write(follow)
        pending_bits = 0

    for value in values:
        symbol_low, symbol_high = bounds[value]
        span = high - low + 1
        high = low + (span * symbol_high // total) - 1
        low = low + (span * symbol_low // total)

        while True:
            if high < ARITH_HALF:
                emit_with_follow(0)
            elif low >= ARITH_HALF:
                emit_with_follow(1)
                low -= ARITH_HALF
                high -= ARITH_HALF
            elif low >= ARITH_QUARTER and high < ARITH_THREE_QUARTER:
                pending_bits += 1
                low -= ARITH_QUARTER
                high -= ARITH_QUARTER
            else:
                break
            low = (low << 1) & ARITH_MAX
            high = ((high << 1) & ARITH_MAX) | 1

    pending_bits += 1
    if low < ARITH_QUARTER:
        emit_with_follow(0)
    else:
        emit_with_follow(1)

    return writer.to_bytes()


def encode_arithmetic_stream(values: list[int]) -> bytes:
    out = bytearray()
    out.append(STREAM_ARITHMETIC)
    out.extend(codec.encode_varint(len(values)))

    if not values:
        out.extend(codec.encode_varint(0))
        out.extend(codec.encode_varint(0))
        return bytes(out)

    freqs = dict(codec.Counter(values))
    out.extend(codec.encode_varint(len(freqs)))
    for symbol, freq in sorted(freqs.items()):
        out.extend(codec.encode_varint(symbol))
        out.extend(codec.encode_varint(freq))

    bitstream = encode_arithmetic_values(values, freqs)
    out.extend(codec.encode_varint(len(bitstream)))
    out.extend(bitstream)
    return bytes(out)


def build_power2_buckets(universe_size: int) -> list[tuple[int, int, int]]:
    buckets: list[tuple[int, int, int]] = []
    start = 0
    width = 1
    while start < universe_size:
        end = min(universe_size - 1, start + width - 1)
        bucket_width = end - start + 1
        extra_bits = (bucket_width - 1).bit_length()
        buckets.append((start, end, extra_bits))
        start = end + 1
        width *= 2
    return buckets


def encode_symbol_huffman_stream(values: list[int]) -> bytes:
    out = bytearray()
    out.append(STREAM_HUFFMAN)
    out.extend(codec.encode_varint(len(values)))

    if not values:
        out.extend(codec.encode_varint(0))
        out.extend(codec.encode_varint(0))
        return bytes(out)

    freqs = dict(codec.Counter(values))
    lengths = codec.build_huffman_code_lengths(freqs)
    codes = codec.build_canonical_huffman_codes(lengths)

    out.extend(codec.encode_varint(len(freqs)))
    for symbol, freq in sorted(freqs.items()):
        out.extend(codec.encode_varint(symbol))
        out.extend(codec.encode_varint(freq))

    bit_buffer = 0
    bit_count = 0
    bitstream = bytearray()
    for value in values:
        code_value, code_len = codes[value]
        bit_buffer = (bit_buffer << code_len) | code_value
        bit_count += code_len
        while bit_count >= 8:
            shift = bit_count - 8
            bitstream.append((bit_buffer >> shift) & 0xFF)
            bit_buffer &= (1 << shift) - 1 if shift else 0
            bit_count = shift

    if bit_count:
        bitstream.append((bit_buffer << (8 - bit_count)) & 0xFF)

    out.extend(codec.encode_varint(len(bitstream)))
    out.extend(bitstream)
    return bytes(out)


def encode_dense_huffman_stream(values: list[int], universe_size: int) -> bytes:
    out = bytearray()
    out.append(STREAM_DENSE_HUFFMAN)
    out.extend(codec.encode_varint(len(values)))
    out.extend(codec.encode_varint(universe_size))

    if universe_size < 0:
        raise ValueError("universe_size must be non-negative.")
    if any(value < 0 or value >= universe_size for value in values):
        raise ValueError("Value exceeds dense-huffman universe size.")

    freqs = dict(codec.Counter(values))
    lengths_map = codec.build_huffman_code_lengths(freqs) if freqs else {}
    dense_lengths = [int(lengths_map.get(symbol, 0)) for symbol in range(universe_size)]
    max_len = max(dense_lengths, default=0)
    bits_per_length = 4 if max_len <= 0x0F else 8
    out.append(bits_per_length)

    if bits_per_length == 4:
        packed_lengths = pack_nibbles(dense_lengths)
    else:
        if max_len > 0xFF:
            raise ValueError(f"Huffman code length {max_len} exceeds 8-bit dense table support.")
        packed_lengths = bytes(dense_lengths)

    out.extend(codec.encode_varint(len(packed_lengths)))
    out.extend(packed_lengths)

    codes = codec.build_canonical_huffman_codes(lengths_map) if lengths_map else {}
    bit_buffer = 0
    bit_count = 0
    bitstream = bytearray()
    for value in values:
        code_value, code_len = codes[value]
        bit_buffer = (bit_buffer << code_len) | code_value
        bit_count += code_len
        while bit_count >= 8:
            shift = bit_count - 8
            bitstream.append((bit_buffer >> shift) & 0xFF)
            bit_buffer &= (1 << shift) - 1 if shift else 0
            bit_count = shift

    if bit_count:
        bitstream.append((bit_buffer << (8 - bit_count)) & 0xFF)

    out.extend(codec.encode_varint(len(bitstream)))
    out.extend(bitstream)
    return bytes(out)


def encode_bitpacked_stream(values: list[int]) -> bytes:
    packed = pack_bits(values)
    out = bytearray()
    out.append(STREAM_BITPACK)
    out.extend(codec.encode_varint(len(values)))
    out.extend(codec.encode_varint(len(packed)))
    out.extend(packed)
    return bytes(out)


def _encode_token_prev1_context_stream(
    values: list[int],
    segment_counts: list[int],
) -> bytes:
    out = bytearray()
    out.append(STREAM_TOKEN_PREV1_CONTEXT)
    out.extend(codec.encode_varint(len(values)))

    if sum(segment_counts) != len(values):
        raise ValueError("segment_counts do not sum to the token_types length.")

    prefix_bits: list[int] = []
    ctx0_values: list[int] = []
    ctx1_values: list[int] = []
    idx = 0

    for count in segment_counts:
        if count == 0:
            continue
        local_head = min(count, 2)
        prefix = values[idx : idx + local_head]
        if any(bit not in {0, 1} for bit in prefix):
            raise ValueError("Token context-range stream only supports 0/1 values.")
        prefix_bits.extend(prefix)
        idx += local_head

        if count <= 2:
            continue

        prev = values[idx - 1]
        for _ in range(count - 2):
            value = values[idx]
            if value not in {0, 1}:
                raise ValueError("Token context-range stream only supports 0/1 values.")
            if prev == 0:
                ctx0_values.append(value)
            else:
                ctx1_values.append(value)
            prev = value
            idx += 1

    out.extend(encode_bitpacked_stream(prefix_bits))
    out.extend(encode_arithmetic_stream(ctx0_values))
    out.extend(encode_arithmetic_stream(ctx1_values))
    return bytes(out)


def _encode_token_type_stream(values: list[int], token_type_mode: str, segment_counts: list[int]) -> bytes:
    mode = _validate_token_type_mode(token_type_mode)
    if mode == TOKEN_TYPE_MODE_BITPACK:
        return encode_bitpacked_stream(values)
    if mode == TOKEN_TYPE_MODE_PREV1_CONTEXT:
        return _encode_token_prev1_context_stream(values, segment_counts)

    candidates = [
        encode_bitpacked_stream(values),
        _encode_token_prev1_context_stream(values, segment_counts),
    ]
    return min(candidates, key=len)


def _encode_ref_bucket_power2_dense_stream(values: list[int], universe_size: int) -> bytes:
    if universe_size < 0:
        raise ValueError("universe_size must be non-negative.")
    if any(value < 0 or value >= universe_size for value in values):
        raise ValueError("ref id exceeds bucket-power2 universe size.")

    buckets = build_power2_buckets(universe_size)
    bucket_meta_by_value: list[tuple[int, int, int]] = [(0, 0, 0)] * universe_size
    for bucket_id, (start, end, extra_bits) in enumerate(buckets):
        for value in range(start, end + 1):
            bucket_meta_by_value[value] = (bucket_id, start, extra_bits)

    bucket_ids: list[int] = []
    extra_chunks: list[tuple[int, int]] = []
    for value in values:
        bucket_id, bucket_start, extra_bits = bucket_meta_by_value[value]
        bucket_ids.append(bucket_id)
        if extra_bits:
            extra_chunks.append((value - bucket_start, extra_bits))

    bucket_stream = encode_dense_huffman_stream(bucket_ids, len(buckets))
    extra_bytes = pack_variable_width_bits(extra_chunks)

    out = bytearray()
    out.append(STREAM_REF_BUCKET_POWER2_DENSE)
    out.extend(codec.encode_varint(universe_size))
    out.extend(bucket_stream)
    out.extend(codec.encode_varint(len(extra_bytes)))
    out.extend(extra_bytes)
    return bytes(out)


def _encode_entropy_payload(
    entropy_json: dict[str, Any],
    *,
    token_type_mode: str,
) -> bytes:
    out = bytearray()
    out.extend(MAGIC)
    out.append(LPE_FORMAT_VERSION)

    sep_bytes = entropy_json["sep"].encode("utf-8")
    out.extend(codec.encode_varint(len(sep_bytes)))
    out.extend(sep_bytes)

    n_targets = int(entropy_json["n_targets"])
    n_rules = int(entropy_json["n_dictionary_rules"])
    out.extend(codec.encode_varint(n_targets))
    out.extend(codec.encode_varint(n_rules))

    target = entropy_json["target"]
    dictionary = entropy_json["dictionary"]

    target_segment_counts = [int(x) for x in target["segment_token_counts"]]
    rule_lengths = [int(x) for x in dictionary["rule_lengths"]]
    target_token_types = [int(x) for x in target["token_types"]]
    dict_token_types = [int(x) for x in dictionary["token_types"]]

    out.extend(encode_symbol_huffman_stream(target_segment_counts))
    out.extend(encode_symbol_huffman_stream(rule_lengths))
    out.extend(_encode_token_type_stream(target_token_types, token_type_mode, target_segment_counts))
    out.extend(_encode_ref_bucket_power2_dense_stream([int(x) for x in target["ref_ids"]], n_rules))
    out.extend(encode_dense_huffman_stream([int(x) for x in target["literal_lengths"]], (max([int(x) for x in target["literal_lengths"]], default=-1) + 1)))
    out.extend(encode_dense_huffman_stream([int(x) for x in target["literal_bytes"]], 256))
    out.extend(_encode_token_type_stream(dict_token_types, token_type_mode, rule_lengths))
    out.extend(_encode_ref_bucket_power2_dense_stream([int(x) for x in dictionary["ref_ids"]], n_rules))
    out.extend(encode_dense_huffman_stream([int(x) for x in dictionary["literal_lengths"]], (max([int(x) for x in dictionary["literal_lengths"]], default=-1) + 1)))
    out.extend(encode_dense_huffman_stream([int(x) for x in dictionary["literal_bytes"]], 256))
    return bytes(out)


def encoding2streams(
    lpjson: dict[str, Any] | str | Path,
    *,
    top_n: int | None = None,
    sort_method: str = DEFAULT_SORT_METHOD,
    selected_ladderon_ids: list[int] | None = None,
    sep: str = DEFAULT_SEP,
    reindex: bool = DEFAULT_REINDEX,
) -> dict[str, Any]:
    normalized_lpjson = load_lpjson(lpjson)
    resolved_sort_method = _validate_sort_method(sort_method)
    resolved_ids = _resolve_selected_ladderon_ids(
        normalized_lpjson,
        top_n=top_n,
        sort_method=resolved_sort_method,
        selected_ladderon_ids=selected_ladderon_ids,
    )

    compressed_list = _compress_lpjson(
        normalized_lpjson,
        sep=sep,
        selected_ladderon_ids=resolved_ids,
    )

    reindex_metadata = None
    if reindex:
        compressed_list, reindex_metadata = _reindex_compressed_list_latest(compressed_list)

    encoded_streams = _build_entropy_ready_json(compressed_list)
    encoded_streams["compress_py_profile"] = "lppack_3.0.0"
    encoded_streams["sort_method"] = resolved_sort_method
    encoded_streams["reindex"] = reindex_metadata if reindex_metadata is not None else {"name": "none"}
    if top_n is not None:
        encoded_streams["top_n"] = int(top_n)
    if resolved_ids is not None:
        encoded_streams["selected_ladderon_ids"] = list(resolved_ids)
    return encoded_streams


def entropy_coding(
    encoded_streams: dict[str, Any],
    *,
    token_type_mode: str = DEFAULT_TOKEN_TYPE_MODE,
    output_path: str | Path | None = None,
) -> bytes:
    payload = _encode_entropy_payload(
        encoded_streams,
        token_type_mode=token_type_mode,
    )
    _write_payload_if_requested(payload, output_path)
    return payload


def compress2bit(
    lpjson: dict[str, Any] | str | Path,
    *,
    top_n: int | None = None,
    sort_method: str = DEFAULT_SORT_METHOD,
    selected_ladderon_ids: list[int] | None = None,
    sep: str = DEFAULT_SEP,
    reindex: bool = DEFAULT_REINDEX,
    token_type_mode: str = DEFAULT_TOKEN_TYPE_MODE,
    output_path: str | Path | None = None,
) -> bytes:
    encoded_streams = encoding2streams(
        lpjson,
        top_n=top_n,
        sort_method=sort_method,
        selected_ladderon_ids=selected_ladderon_ids,
        sep=sep,
        reindex=reindex,
    )
    return entropy_coding(
        encoded_streams,
        token_type_mode=token_type_mode,
        output_path=output_path,
    )


compress = encoding2streams
encode = entropy_coding
compress_to_bytes = compress2bit


__all__ = [
    "encoding2streams",
    "entropy_coding",
    "compress2bit",
    "compress",
    "encode",
    "compress_to_bytes",
    "SORT_METHOD_CONTRIBUTION_FOLDED",
    "SORT_METHOD_CONTRIBUTION_FROZEN",
    "SORT_METHOD_MULTIPLICITY",
    "TOKEN_TYPE_MODE_BITPACK",
    "TOKEN_TYPE_MODE_PREV1_CONTEXT",
    "TOKEN_TYPE_MODE_AUTO_BEST",
]
