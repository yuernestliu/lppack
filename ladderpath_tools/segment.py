"""Unified segmentation utilities for vocab-based tokenization.

Version 1.2.0
Date: 2026.04.01

This module provides a small public API:

    segment(vocab, text) -> list[str] or list[int]

`vocab` can come from:
1. A CSV / TSV / JSON file path.
2. An in-memory list of tokens.
3. An in-memory list of tuples such as:
   - [("token", importance), ...]
   - [("token", token_id, importance), ...]
4. An in-memory dict such as:
   - {"token": importance}
   - {"token": [token_id, importance]}
   - {"token": (token_id, importance)}

Supported decoding methods:
1. `max_matching`: greedy longest-match.
2. `pathpiece`: global minimum-token decoding with longest tie-break.
3. `dag_dp_Importance`: DAG + dynamic programming scored by importance.
4. `dag_dp_logImportance`: DAG + dynamic programming scored by log-importance.
5. `jieba_dp`: jieba precise mode with HMM disabled.
6. `jieba_dp_hmm`: jieba precise mode with HMM enabled.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union


KNOWN_TOKEN_COLS = {"token", "piece", "word", "term", "vocab", "text"}
KNOWN_ID_COLS = {"id", "token_id", "index", "idx"}
KNOWN_IMPORTANCE_COLS = {
    "importance",
    "weight",
    "score",
    "priority",
    "freq",
    "frequency",
    "rank_score",
}
TRIE_END = None
VOCAB_CACHE: Dict[Tuple[Any, ...], "SegmentVocab"] = {}


@dataclass(frozen=True)
class VocabEntry:
    token: str
    token_id: int
    importance: Optional[float] = None


@dataclass(frozen=True)
class SegmentResult:
    tokens: List[str]
    token_ids: List[int]
    score: Tuple[float, ...]
    method: str


class SegmentVocab:
    """Vocabulary container with trie index and flexible loaders."""

    def __init__(self, entries: Sequence[VocabEntry]):
        if not entries:
            raise ValueError("vocab is empty")

        self.entries: List[VocabEntry] = []
        self.token_to_id: Dict[str, int] = {}
        self.id_to_token: Dict[int, str] = {}
        self.token_to_importance: Dict[str, float] = {}
        self.has_explicit_importance = False
        self.max_token_len = 1
        self.trie: Dict[Any, Dict[Any, Any]] = {}

        for entry in entries:
            token = str(entry.token)
            if not token:
                continue
            if token in self.token_to_id:
                raise ValueError(f"duplicate token in vocab: {token!r}")
            if entry.token_id in self.id_to_token:
                raise ValueError(f"duplicate token id in vocab: {entry.token_id}")

            self.entries.append(entry)
            self.token_to_id[token] = entry.token_id
            self.id_to_token[entry.token_id] = token
            if entry.importance is not None:
                self.token_to_importance[token] = float(entry.importance)
                self.has_explicit_importance = True
            self.max_token_len = max(self.max_token_len, len(token))
            self._insert_token(token)

        if not self.entries:
            raise ValueError("vocab does not contain any valid token")

    @classmethod
    def from_source(
        cls,
        vocab: Union["SegmentVocab", str, Path, Mapping[str, Any], Sequence[Any]],
        *,
        token_col: Optional[Union[str, int]] = None,
        id_col: Optional[Union[str, int]] = None,
        importance_col: Optional[Union[str, int]] = None,
        has_header: Optional[bool] = None,
        delimiter: Optional[str] = None,
        encoding: str = "utf-8-sig",
    ) -> "SegmentVocab":
        if isinstance(vocab, SegmentVocab):
            return vocab
        if isinstance(vocab, (str, Path)):
            cache_key = _make_vocab_cache_key(
                Path(vocab),
                token_col=token_col,
                id_col=id_col,
                importance_col=importance_col,
                has_header=has_header,
                delimiter=delimiter,
                encoding=encoding,
            )
            cached = VOCAB_CACHE.get(cache_key)
            if cached is not None:
                return cached
        entries = _load_vocab_entries(
            vocab,
            token_col=token_col,
            id_col=id_col,
            importance_col=importance_col,
            has_header=has_header,
            delimiter=delimiter,
            encoding=encoding,
        )
        loaded = cls(entries)
        if isinstance(vocab, (str, Path)):
            VOCAB_CACHE[cache_key] = loaded
        return loaded

    def match_spans(self, text: str, start: int) -> List[Tuple[int, str]]:
        node = self.trie
        spans: List[Tuple[int, str]] = []
        for idx in range(start, min(len(text), start + self.max_token_len)):
            node = node.get(text[idx])
            if node is None:
                break
            token = node.get(TRIE_END)
            if token is not None:
                spans.append((idx + 1, token))
        return spans

    def resolve_id(self, token: str, unknown_id: int = -1) -> int:
        if token in self.token_to_id:
            return self.token_to_id[token]
        if "[UNK]" in self.token_to_id:
            return self.token_to_id["[UNK]"]
        return unknown_id

    def importance_of(self, token: str, default: float = 1.0) -> float:
        return self.token_to_importance.get(token, default)

    def _insert_token(self, token: str) -> None:
        node = self.trie
        for ch in token:
            node = node.setdefault(ch, {})
        node[TRIE_END] = token


class Segmenter:
    """Segmenter with multiple vocab-driven decoding methods."""

    def __init__(self, vocab: Union[SegmentVocab, str, Path, Mapping[str, Any], Sequence[Any]], **kwargs: Any):
        self.vocab = SegmentVocab.from_source(vocab, **kwargs)

    def segment(
        self,
        text: Union[str, Sequence[str]],
        *,
        method: str = "dag_dp_logImportance",
        output: str = "tokens",
        unknown_id: int = -1,
    ) -> Union[List[str], List[int], List[List[str]], List[List[int]]]:
        if isinstance(text, str):
            result = self.segment_one(text, method=method, output=output, unknown_id=unknown_id)
            return result

        outputs = []
        for one_text in text:
            outputs.append(self.segment_one(one_text, method=method, output=output, unknown_id=unknown_id))
        return outputs

    def segment_one(
        self,
        text: str,
        *,
        method: str = "dag_dp_logImportance",
        output: str = "tokens",
        unknown_id: int = -1,
    ) -> Union[List[str], List[int]]:
        result = self.segment_with_details(text, method=method, unknown_id=unknown_id)
        if output == "tokens":
            return result.tokens
        if output == "ids":
            return result.token_ids
        raise ValueError("output must be 'tokens' or 'ids'")

    def segment_with_details(
        self,
        text: str,
        *,
        method: str = "dag_dp_logImportance",
        unknown_id: int = -1,
    ) -> SegmentResult:
        method = method.lower()
        if method == "dag_dp_importance":
            tokens, score = self._segment_dag_dp_importance(text)
            token_ids = [self.vocab.resolve_id(token, unknown_id=unknown_id) for token in tokens]
            return SegmentResult(tokens=tokens, token_ids=token_ids, score=score, method="dag_dp_Importance")
        if method == "dag_dp_logimportance":
            tokens, score = self._segment_dag_dp_log_importance(text)
            token_ids = [self.vocab.resolve_id(token, unknown_id=unknown_id) for token in tokens]
            return SegmentResult(tokens=tokens, token_ids=token_ids, score=score, method="dag_dp_logImportance")
        if method == "pathpiece":
            tokens = self._segment_pathpiece(text)
            token_ids = [self.vocab.resolve_id(token, unknown_id=unknown_id) for token in tokens]
            return SegmentResult(tokens=tokens, token_ids=token_ids, score=tuple(), method="pathpiece")
        if method == "max_matching":
            tokens = self._segment_max_matching(text)
            token_ids = [self.vocab.resolve_id(token, unknown_id=unknown_id) for token in tokens]
            return SegmentResult(tokens=tokens, token_ids=token_ids, score=tuple(), method="max_matching")
        if method == "jieba_dp":
            tokens = self._segment_jieba(text, hmm=False)
            token_ids = [self.vocab.resolve_id(token, unknown_id=unknown_id) for token in tokens]
            return SegmentResult(tokens=tokens, token_ids=token_ids, score=tuple(), method="jieba_dp")
        if method == "jieba_dp_hmm":
            tokens = self._segment_jieba(text, hmm=True)
            token_ids = [self.vocab.resolve_id(token, unknown_id=unknown_id) for token in tokens]
            return SegmentResult(tokens=tokens, token_ids=token_ids, score=tuple(), method="jieba_dp_hmm")
        raise ValueError(
            "method must be 'max_matching', 'pathpiece', 'dag_dp_Importance', "
            "'dag_dp_logImportance', 'jieba_dp', or 'jieba_dp_hmm'"
        )

    def _segment_max_matching(self, text: str) -> List[str]:
        if not text:
            return []

        tokens: List[str] = []
        start = 0
        while start < len(text):
            matches = self.vocab.match_spans(text, start)
            if matches:
                end, token = max(matches, key=lambda item: item[0] - start)
                tokens.append(token)
                start = end
                continue
            tokens.append(text[start : start + 1])
            start += 1
        return tokens

    def _build_dag(self, text: str) -> Dict[int, List[Tuple[int, str]]]:
        dag: Dict[int, List[Tuple[int, str]]] = {}
        for idx in range(len(text)):
            matches = self.vocab.match_spans(text, idx)
            dag[idx] = matches if matches else [(idx + 1, text[idx : idx + 1])]
        return dag

    def _segment_dag_dp_importance(self, text: str) -> Tuple[List[str], Tuple[float, int]]:
        if not text:
            return [], tuple()

        n = len(text)
        dag = self._build_dag(text)
        route: Dict[int, Tuple[Tuple[float, int], int, str]] = {n: ((0.0, n), n, "")}

        for idx in range(n - 1, -1, -1):
            best_choice: Optional[Tuple[Tuple[float, int], int, str]] = None
            for end, token in dag[idx]:
                local_score = self._importance_score(token)
                next_score = route[end][0]
                score = (local_score + next_score[0], end)
                choice = (score, end, token)
                if best_choice is None or choice[0] > best_choice[0]:
                    best_choice = choice
            if best_choice is None:
                raise RuntimeError("failed to decode DAG path")
            route[idx] = best_choice

        tokens: List[str] = []
        pos = 0
        while pos < n:
            _, next_pos, token = route[pos]
            tokens.append(token)
            pos = next_pos
        return tokens, route[0][0]

    def _segment_dag_dp_log_importance(self, text: str) -> Tuple[List[str], Tuple[float, int]]:
        if not text:
            return [], tuple()

        n = len(text)
        dag = self._build_dag(text)
        route: Dict[int, Tuple[Tuple[float, int], int, str]] = {n: ((0.0, n), n, "")}

        for idx in range(n - 1, -1, -1):
            best_choice: Optional[Tuple[Tuple[float, int], int, str]] = None
            for end, token in dag[idx]:
                next_score = route[end][0]
                score = (self._log_importance_score(token) + next_score[0], end)
                choice = (score, end, token)
                if best_choice is None or choice[0] > best_choice[0]:
                    best_choice = choice
            if best_choice is None:
                raise RuntimeError("failed to decode segmentation path")
            route[idx] = best_choice

        tokens: List[str] = []
        pos = 0
        while pos < n:
            _, next_pos, token = route[pos]
            tokens.append(token)
            pos = next_pos
        return tokens, route[0][0]

    def _segment_pathpiece(self, text: str) -> List[str]:
        if not text:
            return []

        n = len(text)
        widths = [0] * n
        path_length = [float("inf")] * n

        for end in range(n):
            max_width = min(self.vocab.max_token_len, end + 1)
            for width in range(1, max_width + 1):
                start = end - width + 1
                token = text[start : end + 1]
                if token not in self.vocab.token_to_id:
                    continue

                candidate_length = 1 if start == 0 else path_length[start - 1] + 1
                if candidate_length < path_length[end]:
                    widths[end] = width
                    path_length[end] = candidate_length
                elif candidate_length == path_length[end]:
                    # Match the provided pathpiece implementation: when the
                    # optimal path length is tied, keep the longer current token.
                    widths[end] = width

        tokens: List[str] = []
        end = n - 1
        while end >= 0:
            width = widths[end]
            if width == 0:
                tokens.append(text[end : end + 1])
                end -= 1
                continue
            start = end - width + 1
            tokens.append(text[start : end + 1])
            end -= width

        tokens.reverse()
        return tokens

    def _segment_jieba(self, text: str, *, hmm: bool) -> List[str]:
        if not text:
            return []

        try:
            from jieba import Tokenizer as JiebaTokenizer
        except ImportError as exc:
            raise ImportError("jieba methods require the jieba package") from exc

        tokenizer = JiebaTokenizer(dictionary=None)
        tokenizer.initialized = True
        for token in self.vocab.token_to_id:
            freq = max(1, int(round(self._raw_weight_of(token))))
            tokenizer.add_word(token, freq=freq)
        return list(tokenizer.cut(text, cut_all=False, HMM=hmm))

    def _importance_score(self, token: str) -> float:
        # If importance is globally missing, every token is treated as equally
        # important with no score preference between paths.
        if not self.vocab.has_explicit_importance:
            return 0.0
        if token in self.vocab.token_to_id:
            return self.vocab.importance_of(token, default=1.0)
        return 1.0

    def _log_importance_score(self, token: str) -> float:
        # Missing importance is treated as equal weight 1, so log-score is 0.
        return math.log(self._raw_weight_of(token))

    def _raw_weight_of(self, token: str) -> float:
        if token not in self.vocab.token_to_id:
            return 1e-8

        if self.vocab.has_explicit_importance:
            value = self.vocab.importance_of(token, default=1.0)
            return max(float(value), 1e-8)

        return 1.0


def load_vocab(
    vocab: Union[SegmentVocab, str, Path, Mapping[str, Any], Sequence[Any]],
    **kwargs: Any,
) -> SegmentVocab:
    return SegmentVocab.from_source(vocab, **kwargs)


def segment(
    vocab: Union[SegmentVocab, str, Path, Mapping[str, Any], Sequence[Any]],
    text: Union[str, Sequence[str]],
    *,
    method: str = "dag_dp_logImportance",
    output: str = "tokens",
    unknown_id: int = -1,
    **kwargs: Any,
) -> Union[List[str], List[int], List[List[str]], List[List[int]]]:
    """Public convenience wrapper.

    Examples:
        segment(vocab, "南京市长江大桥")
        segment(vocab, "南京市长江大桥", output="ids")
        segment(vocab, ["文本1", "文本2"], method="max_matching")
    """

    segmenter = Segmenter(vocab, **kwargs)
    return segmenter.segment(text, method=method, output=output, unknown_id=unknown_id)


def _load_vocab_entries(
    vocab: Union[str, Path, Mapping[str, Any], Sequence[Any]],
    *,
    token_col: Optional[Union[str, int]],
    id_col: Optional[Union[str, int]],
    importance_col: Optional[Union[str, int]],
    has_header: Optional[bool],
    delimiter: Optional[str],
    encoding: str,
) -> List[VocabEntry]:
    if isinstance(vocab, (str, Path)):
        return _load_vocab_entries_from_path(
            Path(vocab),
            token_col=token_col,
            id_col=id_col,
            importance_col=importance_col,
            has_header=has_header,
            delimiter=delimiter,
            encoding=encoding,
        )
    if isinstance(vocab, Mapping):
        return _load_vocab_entries_from_mapping(vocab)
    if isinstance(vocab, Sequence):
        return _load_vocab_entries_from_sequence(vocab)
    raise TypeError("unsupported vocab type")


def _load_vocab_entries_from_path(
    path: Path,
    *,
    token_col: Optional[Union[str, int]],
    id_col: Optional[Union[str, int]],
    importance_col: Optional[Union[str, int]],
    has_header: Optional[bool],
    delimiter: Optional[str],
    encoding: str,
) -> List[VocabEntry]:
    suffix = path.suffix.lower()
    if suffix in {".csv", ".tsv"}:
        return _load_vocab_entries_from_table(
            path,
            token_col=token_col,
            id_col=id_col,
            importance_col=importance_col,
            has_header=has_header,
            delimiter=delimiter,
            encoding=encoding,
        )
    if suffix == ".json":
        with path.open("r", encoding=encoding) as f:
            data = json.load(f)
        if isinstance(data, Mapping):
            return _load_vocab_entries_from_mapping(data)
        if isinstance(data, list):
            return _load_vocab_entries_from_sequence(data)
        raise ValueError("unsupported json vocab format")
    raise ValueError("vocab file must be .csv, .tsv, or .json")


def _load_vocab_entries_from_table(
    path: Path,
    *,
    token_col: Optional[Union[str, int]],
    id_col: Optional[Union[str, int]],
    importance_col: Optional[Union[str, int]],
    has_header: Optional[bool],
    delimiter: Optional[str],
    encoding: str,
) -> List[VocabEntry]:
    if delimiter is None:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","

    with path.open("r", encoding=encoding, newline="") as f:
        rows = list(csv.reader(f, delimiter=delimiter))

    rows = [row for row in rows if _row_has_content(row)]
    if not rows:
        raise ValueError("vocab file is empty")

    if has_header is None:
        first_row_labels = {cell.strip().lower() for cell in rows[0]}
        has_header = bool(
            first_row_labels
            & (KNOWN_TOKEN_COLS | KNOWN_ID_COLS | KNOWN_IMPORTANCE_COLS)
        )

    entries: List[VocabEntry] = []
    next_id = 0

    if has_header:
        header = [cell.strip() for cell in rows[0]]
        body = rows[1:]
        token_idx = _resolve_named_or_index(token_col, header, default_by_name=KNOWN_TOKEN_COLS, default_index=0)
        id_idx = _resolve_named_or_index(id_col, header, default_by_name=KNOWN_ID_COLS, default_index=None)
        if importance_col is None:
            importance_idx = token_idx + 1 if token_idx + 1 < len(header) else None
        else:
            importance_idx = _resolve_named_or_index(
                importance_col,
                header,
                default_by_name=KNOWN_IMPORTANCE_COLS,
                default_index=None,
            )
        for row in body:
            token = _get_cell_preserve_token(row, token_idx)
            if not token:
                continue
            row_id = _parse_optional_int(_get_cell(row, id_idx))
            row_importance = _parse_optional_float(_get_cell(row, importance_idx))
            if row_id is None:
                row_id = next_id
            entries.append(VocabEntry(token=token, token_id=row_id, importance=row_importance))
            next_id = max(next_id, row_id + 1)
        return entries

    col_count = max(len(row) for row in rows)
    token_idx = int(token_col) if isinstance(token_col, int) else 0
    if isinstance(id_col, int):
        id_idx = id_col
    else:
        id_idx = None

    if isinstance(importance_col, int):
        importance_idx = importance_col
    elif importance_col is None:
        importance_idx = token_idx + 1 if token_idx + 1 < col_count else None
    else:
        importance_idx = None

    for row in rows:
        token = _get_cell_preserve_token(row, token_idx)
        if not token:
            continue
        row_id = _parse_optional_int(_get_cell(row, id_idx))
        row_importance = _parse_optional_float(_get_cell(row, importance_idx))
        if row_id is None:
            row_id = next_id
        entries.append(VocabEntry(token=token, token_id=row_id, importance=row_importance))
        next_id = max(next_id, row_id + 1)
    return entries


def _load_vocab_entries_from_mapping(data: Mapping[str, Any]) -> List[VocabEntry]:
    entries: List[VocabEntry] = []
    next_id = 0
    for token, value in data.items():
        if isinstance(value, Mapping):
            token_id = value.get("id")
            importance = value.get("importance")
        elif isinstance(value, (list, tuple)):
            if len(value) == 0:
                token_id = None
                importance = None
            elif len(value) == 1:
                token_id = None
                importance = value[0]
            else:
                token_id = value[0] if isinstance(value[0], int) else None
                importance = value[1]
        elif isinstance(value, (int, float)):
            token_id = None
            importance = value
        else:
            token_id = None
            importance = None

        if token_id is None:
            token_id = next_id
        entries.append(
            VocabEntry(
                token=str(token),
                token_id=int(token_id),
                importance=_parse_optional_float(importance),
            )
        )
        next_id = max(next_id, int(token_id) + 1)
    return entries


def _load_vocab_entries_from_sequence(data: Sequence[Any]) -> List[VocabEntry]:
    entries: List[VocabEntry] = []
    next_id = 0
    for item in data:
        if isinstance(item, str):
            token = item
            token_id = next_id
            importance = None
        elif isinstance(item, Mapping):
            token = str(item.get("token", item.get("piece", item.get("word", ""))))
            if not token:
                continue
            token_id = item.get("id", next_id)
            importance = item.get("importance", item.get("weight"))
        elif isinstance(item, (list, tuple)):
            if len(item) == 0:
                continue
            if len(item) == 1:
                token = str(item[0])
                token_id = next_id
                importance = None
            elif len(item) == 2:
                token = str(item[0])
                token_id = next_id
                importance = item[1]
            else:
                token = str(item[0])
                token_id = item[1]
                importance = item[2]
        else:
            raise TypeError(f"unsupported vocab row: {item!r}")

        entries.append(
            VocabEntry(
                token=token,
                token_id=int(token_id),
                importance=_parse_optional_float(importance),
            )
        )
        next_id = max(next_id, int(token_id) + 1)
    return entries


def _resolve_named_or_index(
    selector: Optional[Union[str, int]],
    header: Sequence[str],
    *,
    default_by_name: Iterable[str],
    default_index: Optional[int],
) -> Optional[int]:
    if isinstance(selector, int):
        return selector
    if isinstance(selector, str):
        lowered = [name.strip().lower() for name in header]
        try:
            return lowered.index(selector.strip().lower())
        except ValueError as exc:
            raise ValueError(f"column not found in header: {selector!r}") from exc
    lowered = [name.strip().lower() for name in header]
    for candidate in default_by_name:
        if candidate in lowered:
            return lowered.index(candidate)
    return default_index


def _get_cell(row: Sequence[Any], idx: Optional[int]) -> str:
    if idx is None or idx >= len(row):
        return ""
    return str(row[idx]).strip()


def _get_cell_preserve_token(row: Sequence[Any], idx: Optional[int]) -> str:
    if idx is None or idx >= len(row):
        return ""
    return str(row[idx])


def _parse_optional_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return int(text)


def _parse_optional_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return float(text)


def _looks_like_id_column(values: Sequence[Any]) -> bool:
    ints: List[int] = []
    for value in values:
        text = str(value).strip()
        if not text:
            return False
        try:
            ints.append(int(text))
        except ValueError:
            return False
    unique = sorted(set(ints))
    if len(unique) != len(ints):
        return False
    return unique == list(range(unique[0], unique[0] + len(unique)))


def _row_has_content(row: Sequence[Any]) -> bool:
    if not row:
        return False
    return any(str(cell) != "" for cell in row)


def _make_vocab_cache_key(
    path: Path,
    *,
    token_col: Optional[Union[str, int]],
    id_col: Optional[Union[str, int]],
    importance_col: Optional[Union[str, int]],
    has_header: Optional[bool],
    delimiter: Optional[str],
    encoding: str,
) -> Tuple[Any, ...]:
    resolved = path.expanduser().resolve()
    stat = resolved.stat()
    return (
        str(resolved),
        stat.st_mtime_ns,
        stat.st_size,
        token_col,
        id_col,
        importance_col,
        has_header,
        delimiter,
        encoding,
    )
