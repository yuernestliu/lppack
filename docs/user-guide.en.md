# lppack 3.0.1 User Guide

[中文操作手册](user-guide.zh-CN.md) · [Project README](../README.md)

This guide documents the public interfaces in `lppack` 3.0.1 and the `lpjson` 3.0.0 format.

## Contents

- [1. Installation](#1-installation)
- [2. Core concepts](#2-core-concepts)
- [3. Computing a Ladderpath](#3-computing-a-ladderpath)
- [4. The lpjson 3.0.0 format](#4-the-lpjson-300-format)
- [5. Inspecting, saving, and loading results](#5-inspecting-saving-and-loading-results)
- [6. Drawing laddergraphs](#6-drawing-laddergraphs)
- [7. CUDA and Apple Silicon backends](#7-cuda-and-apple-silicon-backends)
- [8. Compression](#8-compression)
- [9. Lambda from an existing laddergraph](#9-lambda-from-an-existing-laddergraph)
- [10. Ladderon ranking and vocabulary construction](#10-ladderon-ranking-and-vocabulary-construction)
- [11. Segmentation and tokenization](#11-segmentation-and-tokenization)
- [12. Species distances and phylogeny](#12-species-distances-and-phylogeny)
- [13. Reproducibility and backend comparisons](#13-reproducibility-and-backend-comparisons)
- [14. Migrating from 2.x](#14-migrating-from-2x)
- [15. Troubleshooting](#15-troubleshooting)

## 1. Installation

### 1.1 Requirements

- Python 3.10 or newer.
- Installation includes the Python dependencies used by every feature.
- Accelerated backends still require compatible hardware and system software.

### 1.2 Install from GitHub

```bash
python -m pip install git+https://github.com/yuernestliu/lppack.git
```

### 1.3 Clone and install in editable mode

```bash
git clone https://github.com/yuernestliu/lppack.git
cd lppack
python -m pip install -e .
```

Run Python from any directory after installation; do not add a personal checkout path to `sys.path`.

### 1.4 Installed dependencies

A normal installation includes NumPy, Numba, pydivsufsort, Graphviz's Python binding, tqdm, Jieba, SciPy, Matplotlib, and Biopython. MLX is installed automatically on compatible macOS arm64 environments.

### 1.5 Graphviz system program

The `graphviz` package installed by `pip` is the Python interface used by the drawing function. It does not install the Graphviz system application or its `dot` executable.

The core Ladderpath computations do not use Graphviz. If the system application is absent or its installation is unsuccessful, computing, saving, loading, and inspecting Ladderpaths still work normally. Install the system application only when you need to render a laddergraph as PDF, PNG, or another image format:

```bash
# macOS with Homebrew
brew install graphviz

# Ubuntu / Debian
sudo apt-get install graphviz

# Conda
conda install -c conda-forge graphviz python-graphviz
```

On Windows, the simplest installation method is Windows Package Manager:

```powershell
winget install graphviz
```

Alternatively, download a 64-bit EXE installer from the [official Graphviz download page](https://graphviz.org/download/). After installation, open a new Command Prompt or PowerShell window and run:

```powershell
dot -V
```

If Windows reports that `dot` is not recognized, add Graphviz to `PATH` manually:

1. Find Graphviz's `bin` directory. The usual location is `C:\Program Files\Graphviz\bin`; use the actual installation directory if yours is different.
2. Open the Start menu and search for **Edit the system environment variables**.
3. Select **Environment Variables**. Under **User variables**, select **Path**, choose **Edit**, and then choose **New**.
4. Add the full path to Graphviz's `bin` directory and confirm all dialogs.
5. Close and reopen the terminal, IDE, or notebook application so that it receives the updated `PATH`.
6. Run `dot -V` again to verify the configuration.

### 1.6 Verify the installation

```bash
python -c "import lppack; print(lppack.__version__, lppack.LPJSON_VERSION)"
```

Expected output:

```text
3.0.1 3.0.0
```

## 2. Core concepts

- **Target**: an input sequence. A computation can include one or many targets.
- **Basic building block**: one character from the input alphabet.
- **Ladderon**: a reusable sequence component represented by literals and/or lower-level ladderons.
- **Ladderpath index**: the size of the resulting Ladderpath representation.
- **Order index**: `size-index - ladderpath-index`.
- **Size index**: total target length, including target multiplicities.
- **POM**: a partially ordered multiset view of the hierarchy.
- **Laddergraph**: the graph representation of targets, ladderons, and their reuse relations.
- **lpjson**: the standard dictionary representation returned by the package.

## 3. Computing a Ladderpath

The main entry point is `lppack.get_ladderpath()`:

```python
import lppack as lp

targets = [
    "ABCDBCDBCDCDEFEF",
    "AA",
    "DBCDBCA",
    "BCDCDAEF",
]

lpjson = lp.get_ladderpath(
    targets,
    fill_ladderons_STR=False,
    save_file_name=None,
    show_version=True,
    mode="cpu",
)
```

### 3.1 Parameters

- `targets0`: a `list[str]` of unique targets or an insertion-ordered `dict[str, int]` mapping each target to its multiplicity.
- `fill_ladderons_STR`: populate the cached explicit sequence in every ladderon and target when `True`. The default is `False` to keep the result smaller.
- `save_file_name`: JSON output path. No file is written when it is `None`.
- `show_version`: print package and lpjson format versions.
- `mode`: exactly one of `"cpu"`, `"gpu-cuda"`, or `"macm"`.

### 3.2 Input rules

- Every target must be a string with at least two characters.
- A `list` must not contain duplicate strings.
- Use `dict[str, int]` for repeated targets.
- Dictionary insertion order defines target ID order.

Repeated-target example:

```python
targets = {
    "ABCDBCDBCDCDEFEF": 1,
    "AA": 3,
    "DBCDBCA": 2,
    "BCDCDAEF": 1,
}

lpjson = lp.get_ladderpath(targets, show_version=False, mode="cpu")
```

### 3.3 Return value

`get_ladderpath()` returns an in-memory Python dictionary. It writes the same object as JSON when `save_file_name` is provided.

```python
print(lpjson["ladderpath-index"])
print(lpjson["order-index"])
print(lpjson["size-index"])
```

## 4. The lpjson 3.0.0 format

A compact illustrative object looks like this:

```jsonc
{
  "info": "3.0.0",
  "ladderpath-index": 20,
  "order-index": 13,
  "size-index": 33,
  "basic_building_blocks": ["A", "D", "C", "E", "F", "B"],
  "ladderons": {
    "0": [[1, 1], 6, "", {"-1": [3], "-3": [0]}],
    "1": [["D", 4], 3, "", {"0": [0, 3]}]
  },
  "targets": {
    "-1": [["A", 4, 0, 2, 3, 3], 16, "", 1],
    "-2": [["AA"], 2, "", 1]
  },
  "duplications_info": {},
  "input_type": "list"
}
```

JSON object keys are strings on disk. `load_ladderpath_json()` restores numeric IDs to Python integers.

### 4.1 Top-level fields

- `info`: lpjson format version, currently `"3.0.0"`.
- `ladderpath-index`: Ladderpath index.
- `order-index`: size index minus Ladderpath index.
- `size-index`: total target length including multiplicity.
- `basic_building_blocks`: basic characters in the representation.
- `ladderons`: mapping from ladderon ID `0, 1, 2, ...` to a ladderon record.
- `targets`: mapping from target ID `-1, -2, -3, ...` to a target record.
- `duplications_info`: multiplicity mapping for repeated dictionary inputs; otherwise empty.
- `input_type`: `"list"` or `"dict"`.

### 4.2 Ladderon records

Each ladderon is `[COMP, LEN, STR, POS]`:

- `COMP`: ordered literals and lower-level ladderon IDs.
- `LEN`: full character length.
- `STR`: cached full sequence, or an empty string when the cache has not been filled.
- `POS`: parent ID to occurrence-position list.

### 4.3 Target records

Each target is `[COMP, LEN, STR, REP]`:

- `COMP`: ordered literals and ladderon IDs.
- `LEN`: target length.
- `STR`: cached full sequence, or an empty string.
- `REP`: target multiplicity.

## 5. Inspecting, saving, and loading results

### 5.1 Display the three indices

```python
indices = lp.disp3index(lpjson)
# (ladderpath-index, order-index, size-index)
```

### 5.2 Obtain the POM

```python
pom, pom_text = lp.POM_from_JSON(lpjson, display_str=True)
```

`pom` is a level-indexed dictionary. `pom_text` is a compact string representation when `display_str=True`.

### 5.3 Fill or clear string caches

```python
lp.fill_lpjson_STR(lpjson)
lp.clear_lpjson_STR(lpjson)
```

Both functions modify the object in place.

### 5.4 Save and load

```python
lp.save_ladderpath_json(lpjson, "result.json")
loaded = lp.load_ladderpath_json("result.json")
```

Use the package loader before calling tools that require integer target and ladderon IDs.

### 5.5 Other helpers

```python
target_id_list = lp.reconstruct_targets_list(lpjson)
basic_counts = lp.bbb_multiplicity(lpjson)
```

For dictionary input, reconstructed repeated IDs are grouped in target insertion order.

## 6. Drawing laddergraphs

The core computations in the preceding sections do not require the Graphviz system application. Install it only to render a laddergraph; the Python binding is already included in the normal installation.

```python
graph = lp.draw_laddergraph(
    lpjson,
    show_longer_than=0,
    style="ellipse",
    TargetNames_UserDefined={-1: "sample-1", -2: "sample-2"},
    warning_n_ladderons_to_show=500,
    rankdir="BT",
    color="grey",
    figsize="8,8",
    save_fig_name="figures/laddergraph",
    figformat="pdf",
    cleanGVfile=True,
    render_fig=True,
    show_node_size=True,
)
```

Important parameters:

- `show_longer_than`: include ladderons longer than this threshold.
- `style`: `ellipse`, `ellipse-OnlyShowTargetID`, `box`, or `box-OnlyShowTargetID`.
- `TargetNames_UserDefined`: optional target-ID-to-label mapping.
- `rankdir`: `BT`, `TB`, `LR`, or `RL`.
- `render_fig=False`: save the Graphviz source without rendering an image.
- `show_node_size`: append sequence lengths to node labels.

The function fills the `STR` cache in the supplied object.

## 7. CUDA and Apple Silicon backends

### 7.1 NVIDIA CUDA

The environment must provide a working NVIDIA driver and CUDA runtime compatible with Numba.

```python
lpjson = lp.get_ladderpath(targets, mode="gpu-cuda")
```

### 7.2 Apple Silicon / MLX

Use native arm64 Python on an Apple Silicon Mac:

```python
lpjson = lp.get_ladderpath(targets, mode="macm")
```

The MLX backend requires macOS, native arm64 Python, and an available Metal backend. It does not silently fall back to CPU.

### 7.3 Non-ASCII input

The accelerated backends automatically switch to Unicode character coordinates when any input contains non-ASCII characters. This avoids mixing character positions with UTF-8 byte offsets.

## 8. Compression

The compression module converts `lpjson` to structured streams and then to entropy-coded bytes.

```python
from lppack.ladderpath_tools import compress as lp_compress

encoded_streams = lp_compress.encoding2streams(
    lpjson,
    top_n=3000,
    sort_method="contribution_folded",
    selected_ladderon_ids=None,
    sep="@",
    reindex=True,
)

payload = lp_compress.entropy_coding(
    encoded_streams,
    token_type_mode="prev1_context",
    output_path="result.lpe1",
)
```

The combined entry point is:

```python
payload = lp_compress.compress2bit(
    lpjson,
    top_n=3000,
    sort_method="contribution_folded",
    token_type_mode="prev1_context",
    output_path="result.lpe1",
)
```

Valid sorting methods are `contribution_folded`, `contribution_frozen`, and `multiplicity`. Valid token-type modes are `prev1_context`, `bitpack`, and `auto_best`.

The `.lpe1` serialization format has its own format identifier; it is independent of the `lppack` and `lpjson` version numbers.

## 9. Lambda from an existing laddergraph

```python
from lppack.ladderpath_tools import lambda_from_laddergraph as lp_lambda

lambda_one = lp_lambda.get([-1], lpjson)
lambda_group = lp_lambda.get([-1, -3], lpjson)
```

Target strings may be used instead of IDs:

```python
lambda_group = lp_lambda.get(
    ["ABCDBCDBCDCDEFEF", "DBCDBCA"],
    lpjson,
)
```

Constraints:

- `leaves` must be a non-empty list containing all integers or all strings.
- The selected targets must exist in the supplied object.
- The current implementation does not accept an `lpjson` with duplicates.
- Load saved data with `load_ladderpath_json()`, not plain `json.load()`.

## 10. Ladderon ranking and vocabulary construction

```python
from lppack.ladderpath_tools import sortout

vocab = sortout.get(
    lpjson,
    vocab_size=100,
    save_vacab_file="vocab.csv",
    show_selectedladderons=False,
    save_sort_token_lpjson=None,
    iter_logs=False,
    sort_method="contribution_folded",
)
```

The returned mapping has the form:

```python
{
    "the": [0, 4049],
    " and ": [1, 2663],
    "in": [2, 4682],
}
```

`iter_logs=True` enables detailed progress diagnostics and can be substantially slower. The public sorting methods are `contribution_folded`, `contribution_frozen`, and `multiplicity`.

## 11. Segmentation and tokenization

### 11.1 Segment text directly

```python
from lppack.ladderpath_tools.segment import segment

vocab = {
    "Chapter": 100.0,
    " ": 50.0,
    "1": 80.0,
    "Loomings": 120.0,
}

tokens = segment(vocab, "Chapter 1 Loomings", output="tokens")
ids = segment(vocab, "Chapter 1 Loomings", output="ids", unknown_id=-1)
```

`vocab` may be a CSV, TSV, or JSON path, or an in-memory list/dictionary. `text` may be one string or a sequence of strings.

Methods:

- `max_matching`
- `pathpiece`
- `dag_dp_Importance`
- `dag_dp_logImportance`
- `jieba_dp`
- `jieba_dp_hmm`

Jieba is installed as part of the normal `lppack` installation.

Table vocabularies can specify `token_col`, `id_col`, `importance_col`, `has_header`, `delimiter`, and `encoding`.

### 11.2 Use `LadderpathTokenizer`

```python
from lppack.ladderpath_tools.tokenizer import LadderpathTokenizer

tokenizer = LadderpathTokenizer(
    "vocab.csv",
    method="dag_dp_logImportance",
    token_col="token",
    has_header=True,
)

ids = tokenizer.encode("Chapter 1 Loomings")
tokens = tokenizer.encode("Chapter 1 Loomings", output_mode="tokens")
details = tokenizer.encode("Chapter 1 Loomings", output_mode="details")
text = tokenizer.decode(ids)
```

`encode()` supports `ids`, `tokens`, `str`, and `details` output modes, plus batch input. `decode()` accepts one ID sequence or a batch of ID sequences.

## 12. Species distances and phylogeny

### 12.1 Build a distance table

```python
import lppack as lp
from lppack.ladderpath_tools import species_distance

targets = {
    "ATGGAGAAATCTCAA": 1,
    "TTTCGAATTCT": 1,
    "AGGATATTTATAAAT": 1,
    "ATCGATTTTGTT": 2,
    "CGGATGACTTCCT": 1,
}

lpjson = lp.get_ladderpath(targets, show_version=False)

species_info = {
    "species1": ("Morado", [-1, -2]),
    "species2": ("Musa_ornata", [-3, -4]),
    "species3": ("Kluai_Khai", [-5, -4]),
}

distance_table, name_map, merged = species_distance.get_distance_table(
    lpjson,
    species_info,
    method="jaccard",
    verbose=True,
)
```

`species_info` must use real target IDs. Across the full mapping, occurrences of each target ID must match the multiplicity represented by the source object. The supported distance methods are `jaccard` and `dice`.

### 12.2 Build a phylogenetic tree

The required SciPy, Matplotlib, and Biopython packages are installed with `lppack`:

```python
from lppack.ladderpath_tools import phylo

tree = phylo.phylo(
    distance_table,
    name_map,
    method="nj",
    outgroup_name="species1",
    newick_file_name="my_tree.nwk",
    fig_file_name="my_tree.png",
    show_plot=True,
)
```

- `method` is `nj` or `upgma`.
- NJ is unrooted by default; provide `outgroup_name` to reroot it.
- UPGMA is rooted by construction.
- Newick output supports `.nwk` and `.newick`.
- Figure output supports PNG, PDF, SVG, JPEG, and TIFF.

## 13. Reproducibility and backend comparisons

- The CPU backend may choose different equivalent candidates when internal ties occur. Byte-identical output is not guaranteed across repeated runs.
- The CUDA and Apple Silicon backends are designed to be deterministic for the same input and backend.
- Exact equality between CPU and accelerated output is not guaranteed for every tie case.
- If exact identity matters, compare the complete `lpjson` objects on representative inputs in the intended environment.
- CUDA simulation can check code paths but does not represent real GPU performance.

## 14. Migrating from 2.x

Use the installed public package:

```python
import lppack as lp
```

Do not add a private checkout path to `sys.path`. Select the backend through the `mode` parameter. Duplicate list entries are no longer accepted; use a multiplicity dictionary. The software version is `3.0.1`; the `lpjson` format identifier remains `3.0.0`.

Older saved objects may use a different schema identifier and should not be assumed byte-compatible with 3.0.0. Validate required fields before using them in a 3.0.0 workflow.

## 15. Troubleshooting

### `ModuleNotFoundError: No module named 'lppack'`

Install the repository with `python -m pip install -e .` from its root, or install it directly from GitHub. Do not run by adding the repository itself to `sys.path`.

### Graphviz Python package is installed but rendering fails

Install the Graphviz system executable and make sure `dot` is on `PATH`:

```bash
dot -V
```

### CUDA mode fails at import or device initialization

Confirm that Numba imports successfully and can access a compatible NVIDIA CUDA runtime.

### `macm` reports that MLX or Metal is unavailable

Use native arm64 Python on Apple Silicon macOS. Confirm that `platform.machine()` returns `arm64` and that MLX can access Metal.

### Duplicate list input raises an error

Convert the input to `dict[str, int]`, where each value is the total multiplicity of that target.

### A saved JSON has string IDs

Load it with `lppack.load_ladderpath_json()` rather than `json.load()`.
