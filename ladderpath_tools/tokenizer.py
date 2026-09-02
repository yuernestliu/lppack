"""Tokenizer wrapper built on top of segment.py.

Version 1.4.0
Date: 2026.04.01
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

try:
    from segment import SegmentResult, SegmentVocab, Segmenter, load_vocab
except ImportError:  # pragma: no cover
    from .segment import SegmentResult, SegmentVocab, Segmenter, load_vocab


VocabSource = Union[SegmentVocab, str, Path, Mapping[str, Any], Sequence[Any]]


@dataclass
class TokenizerConfig:
    method: str = "dag_dp_logImportance"
    unknown_id: int = -1
    token_col: Optional[Union[str, int]] = None
    id_col: Optional[Union[str, int]] = None
    importance_col: Optional[Union[str, int]] = None
    has_header: Optional[bool] = None
    delimiter: Optional[str] = None
    encoding: str = "utf-8-sig"


class LadderpathTokenizer:
    """High-level tokenizer built on the segmentation methods.

    Public API preference:
    - `encode(...)` is the primary entry point
    - `decode(...)` restores text from token ids
    """

    def __init__(
        self,
        vocab: Optional[VocabSource] = None,
        *,
        method: str = "dag_dp_logImportance",
        unknown_id: int = -1,
        token_col: Optional[Union[str, int]] = None,
        id_col: Optional[Union[str, int]] = None,
        importance_col: Optional[Union[str, int]] = None,
        has_header: Optional[bool] = None,
        delimiter: Optional[str] = None,
        encoding: str = "utf-8-sig",
    ) -> None:
        self.config = TokenizerConfig(
            method=method,
            unknown_id=unknown_id,
            token_col=token_col,
            id_col=id_col,
            importance_col=importance_col,
            has_header=has_header,
            delimiter=delimiter,
            encoding=encoding,
        )
        self.vocab: Optional[SegmentVocab] = None
        self.segmenter: Optional[Segmenter] = None
        if vocab is not None:
            self.load_vocab(vocab)

    @property
    def token2id(self) -> Dict[str, int]:
        self._require_vocab()
        return self.vocab.token_to_id

    @property
    def id2token(self) -> Dict[int, str]:
        self._require_vocab()
        return self.vocab.id_to_token

    @property
    def vocab_size(self) -> int:
        self._require_vocab()
        return len(self.vocab.token_to_id)

    @property
    def max_token_len(self) -> int:
        self._require_vocab()
        return self.vocab.max_token_len

    def load_vocab(self, vocab: VocabSource, **overrides: Any) -> "LadderpathTokenizer":
        loader_kwargs = self._loader_kwargs(**overrides)
        self.vocab = load_vocab(vocab, **loader_kwargs)
        self.segmenter = Segmenter(self.vocab)
        return self

    def set_method(self, method: str) -> "LadderpathTokenizer":
        self.config.method = method
        return self

    def encode(
        self,
        texts: Union[str, Sequence[str]],
        *,
        method: Optional[str] = None,
        output_mode: str = "ids",
        str_sep: str = "|",
        unknown_id: Optional[int] = None,
    ) -> Union[
        List[int],
        List[str],
        str,
        SegmentResult,
        List[List[int]],
        List[List[str]],
        List[str],
        List[SegmentResult],
    ]:
        self._require_segmenter()
        normalized = output_mode.lower()
        resolved_method = method or self.config.method
        resolved_unknown_id = self.config.unknown_id if unknown_id is None else unknown_id

        if normalized in {"id", "ids"}:
            return self.segmenter.segment(
                texts,
                method=resolved_method,
                output="ids",
                unknown_id=resolved_unknown_id,
            )

        if normalized in {"token", "tokens"}:
            return self.segmenter.segment(
                texts,
                method=resolved_method,
                output="tokens",
                unknown_id=resolved_unknown_id,
            )

        if normalized in {"detail", "details"}:
            if isinstance(texts, str):
                return self.segmenter.segment_with_details(
                    texts,
                    method=resolved_method,
                    unknown_id=resolved_unknown_id,
                )
            return [
                self.segmenter.segment_with_details(
                    one_text,
                    method=resolved_method,
                    unknown_id=resolved_unknown_id,
                )
                for one_text in texts
            ]

        if normalized != "str":
            raise ValueError("output_mode must be 'ids', 'tokens', 'str', or 'details'")

        tokenized = self.segmenter.segment(
            texts,
            method=resolved_method,
            output="tokens",
            unknown_id=resolved_unknown_id,
        )
        if isinstance(texts, str):
            return str_sep.join(tokenized)
        return [str_sep.join(tokens) for tokens in tokenized]

    def decode(
        self,
        token_ids: Union[Sequence[int], Sequence[Sequence[int]]],
        *,
        unknown_token: str = "[UNK]",
    ) -> Union[str, List[str]]:
        self._require_vocab()
        is_batch = bool(token_ids) and isinstance(token_ids[0], (list, tuple))

        sequences: Sequence[Sequence[int]]
        if is_batch:
            sequences = token_ids  # type: ignore[assignment]
        else:
            sequences = [token_ids]  # type: ignore[list-item]

        decoded: List[str] = []
        for seq in sequences:
            decoded.append("".join(self.id2token.get(token_id, unknown_token) for token_id in seq))

        if is_batch:
            return decoded
        return decoded[0]

    def vocab_summary(self) -> Dict[str, Any]:
        self._require_vocab()
        return {
            "vocab_size": self.vocab_size,
            "max_token_len": self.max_token_len,
            "has_explicit_importance": self.vocab.has_explicit_importance,
            "default_method": self.config.method,
            "unknown_id": self.config.unknown_id,
        }

    def _loader_kwargs(self, **overrides: Any) -> Dict[str, Any]:
        return {
            "token_col": overrides.get("token_col", self.config.token_col),
            "id_col": overrides.get("id_col", self.config.id_col),
            "importance_col": overrides.get("importance_col", self.config.importance_col),
            "has_header": overrides.get("has_header", self.config.has_header),
            "delimiter": overrides.get("delimiter", self.config.delimiter),
            "encoding": overrides.get("encoding", self.config.encoding),
        }

    def _require_vocab(self) -> None:
        if self.vocab is None:
            raise RuntimeError("vocab is not loaded; call load_vocab(...) first")

    def _require_segmenter(self) -> None:
        if self.segmenter is None:
            raise RuntimeError("tokenizer is not ready; call load_vocab(...) first")
