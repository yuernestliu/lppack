# Changelog

All notable public changes to `lppack` are documented here.

## 3.0.0 — 2026-09-02

### Added

- Standard install support through `pyproject.toml`.
- Public top-level API through `import lppack as lp`.
- NVIDIA CUDA backend selectable with `mode="gpu-cuda"`.
- Apple Silicon / MLX Metal backend selectable with `mode="macm"`.
- Unicode character-coordinate handling in accelerated backends.
- Entropy-coded Ladderpath compression output.
- Ladderon ranking and vocabulary construction.
- Vocabulary-based segmentation and `LadderpathTokenizer`.
- Species-distance and phylogenetic-tree utilities.
- English and Chinese public documentation.
- Machine-readable software citation metadata.

### Changed

- Unified the package, all Python modules, all computation backends, and the `lpjson` format under version `3.0.0`.
- Standardized the public package name as `lppack`.
- Made one normal installation include all Python dependencies used by the package.
- Replaced manual `sys.path` setup with normal `pip` installation.
- Standardized backend selection as `cpu`, `gpu-cuda`, or `macm`.
- Standardized repeated input targets as insertion-ordered `dict[str, int]` values.
- Updated `reconstruct_targets_list()` for the current multiplicity representation.

### Removed

- Legacy development-only import names and manual path setup.

### Compatibility notes

- Python 3.10 or newer is required.
- `lpjson` 3.0.0 should be treated as a new public format identifier; older saved objects are not guaranteed to be byte-compatible.
- The `.lpe1` binary compression format has an independent internal format identifier.
