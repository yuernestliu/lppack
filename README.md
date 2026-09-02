# lppack

**Ladderpath Python Toolkit — version 3.0.0**

[中文 README](README.zh-CN.md) · [English user guide](docs/user-guide.en.md) · [中文操作手册](docs/user-guide.zh-CN.md) · [Changelog](docs/CHANGELOG.md) · [Citation](CITATION.cff)

`lppack` discovers reusable, nested structure in one or more sequences and represents the result as a standard `lpjson` object. The package includes the Ladderpath algorithm, laddergraph visualization, partially ordered multiset (POM) analysis, compression, ladderon ranking, vocabulary-based segmentation and tokenization, inter-species distance calculation, and phylogenetic tree construction.

## Requirements

- Python 3.10 or newer.
- Installing `lppack` installs all Python dependencies for every included feature.
- CUDA and Apple Silicon acceleration still require compatible hardware and system software.

## Installation

Install the current release directly from GitHub:

```bash
python -m pip install git+https://github.com/yuernestliu/lppack.git
```

For development or editable use:

```bash
git clone https://github.com/yuernestliu/lppack.git
cd lppack
python -m pip install -e .
```

Both installation methods install NumPy, Numba, pydivsufsort, MLX on compatible Apple Silicon systems, Graphviz's Python binding, tqdm, Jieba, SciPy, Matplotlib, and Biopython automatically.

The `graphviz` package installed by `pip` is a Python interface. It does not install the Graphviz system application or its `dot` executable. The core Ladderpath computations do not use `dot`: if the Graphviz system application is not installed, or its installation is unsuccessful, you can still compute and inspect Ladderpaths normally. The system application is needed only when rendering laddergraphs as PDF, PNG, or other image files.

Install the Graphviz system application only if you need image rendering:

```bash
# macOS with Homebrew
brew install graphviz

# Ubuntu / Debian
sudo apt-get install graphviz

# Conda
conda install -c conda-forge graphviz python-graphviz
```

On Windows, install Graphviz with Windows Package Manager:

```powershell
winget install graphviz
```

Alternatively, download a 64-bit EXE installer from the [official Graphviz download page](https://graphviz.org/download/). Open a new Command Prompt or PowerShell window after installation and run `dot -V`. If Windows reports that `dot` is not recognized, add Graphviz to `PATH` manually:

1. Find Graphviz's `bin` directory, usually `C:\Program Files\Graphviz\bin`.
2. Open the Start menu, search for **Edit the system environment variables**, and select **Environment Variables**.
3. Under **User variables**, select **Path**, choose **Edit** and then **New**, and add the `bin` directory.
4. Confirm all dialogs, close and reopen your terminal or IDE, and run `dot -V` again.

## Quick start

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

ladderpath_index, order_index, size_index = lp.disp3index(lpjson)
```

The three returned measures are:

- **Ladderpath index**: the size of the Ladderpath representation.
- **Order index**: `size-index - ladderpath-index`.
- **Size index**: the total length of all targets, including multiplicity.

### Repeated targets

A `list[str]` must not contain duplicate strings. Represent repeated targets with an insertion-ordered `dict[str, int]`:

```python
targets = {
    "ABCDBCDBCDCDEFEF": 1,
    "AA": 3,
    "DBCDBCA": 2,
    "BCDCDAEF": 1,
}

lpjson = lp.get_ladderpath(targets, show_version=False, mode="cpu")
```

Every target must contain at least two characters.

## Computation backends

| Mode | Platform | Notes |
| --- | --- | --- |
| `cpu` | Any supported Python platform | Default; standard-library implementation |
| `gpu-cuda` | NVIDIA GPU with a working CUDA environment | Uses Numba, NumPy, and pydivsufsort |
| `macm` | Native arm64 Python on Apple Silicon macOS | Uses MLX Metal, NumPy, and pydivsufsort |

The public call is the same for every backend:

```python
lpjson_cuda = lp.get_ladderpath(targets, mode="gpu-cuda")
lpjson_mac = lp.get_ladderpath(targets, mode="macm")
```

Hardware acceleration is used only when requested. Selecting an unavailable backend raises an error instead of silently falling back to CPU.

## Common operations

### Save and load `lpjson`

```python
lp.save_ladderpath_json(lpjson, "result.json")
lpjson = lp.load_ladderpath_json("result.json")
```

Use `load_ladderpath_json()` rather than plain `json.load()` when the integer IDs in `targets`, `ladderons`, and position maps must be restored.

### Obtain a POM representation

```python
pom, pom_text = lp.POM_from_JSON(lpjson, display_str=True)
```

### Draw a laddergraph

```python
graph = lp.draw_laddergraph(
    lpjson,
    save_fig_name="laddergraph",
    figformat="pdf",
    render_fig=True,
    show_node_size=True,
)
```

Only image rendering requires the Graphviz system executable; all core Ladderpath computations work without it. Its Python binding is installed with `lppack`.

## Included tools

| Module | Purpose |
| --- | --- |
| `ladderpath_tools.compress` | Encode `lpjson` as streams and entropy-coded `.lpe1` bytes |
| `ladderpath_tools.lambda_from_laddergraph` | Recalculate lambda for selected targets in an existing laddergraph |
| `ladderpath_tools.sortout` | Rank ladderons and build vocabularies |
| `ladderpath_tools.segment` | Segment text with a vocabulary |
| `ladderpath_tools.tokenizer` | High-level encode/decode wrapper |
| `ladderpath_tools.species_distance` | Build Jaccard- or Dice-style inter-species distance tables |
| `ladderpath_tools.phylo` | Build NJ or UPGMA phylogenetic trees |

## Documentation

- [English user guide](docs/user-guide.en.md)
- [中文操作手册](docs/user-guide.zh-CN.md)
- [Changelog](docs/CHANGELOG.md)
- [Machine-readable citation metadata](CITATION.cff)

## `lpjson` and reproducibility

The software release and `lpjson` format are both version `3.0.0`. Saved objects record the format in `lpjson["info"]`.

The CPU backend may select different equivalent candidates when internal ties occur, so repeated runs are not guaranteed to produce byte-identical `lpjson` objects. The CUDA and Apple Silicon backends are designed to be deterministic for the same input and backend. Validate exact cross-backend equality for any workflow that depends on byte-for-byte identity.

## Migration from 2.x

- Use `import lppack as lp` instead of importing `ladderpath` or a private development package.
- Install the repository with `pip` instead of modifying `sys.path`.
- Select acceleration through `mode="cpu"`, `mode="gpu-cuda"`, or `mode="macm"`.
- Duplicate strings are represented with `dict[str, int]`, not repeated entries in a list.
- The `lpjson` format identifier is now `3.0.0`.

## Related research

Ladderpath is associated with the following publications:

1. **[Method][Linguistics]** [*Measuring language complexity from hierarchical reuse of recurring patterns*](https://arxiv.org/abs/2606.11531) (2026).
2. **[LLM][Biology]** [*Compression-Based Tokenization Improves Language Modeling of Hierarchical Genomic Structure*](https://langtaosha.org.cn/lts/zh_Hans/preprint/view/51) (2025).
3. **[Biology][Phylogenetics]** [*Structural Compression for Phylogenetic Inference under Alignment Instability and Indel-Rich Evolution*](https://arxiv.org/abs/2607.15308) (2026).
4. **[Method][Compression]** [*Text Distance from Nested and Hierarchical Repetitions: A Compression-Based Perspective*](https://arxiv.org/abs/2607.05416) (2026).
5. **[Biology][Development]** [*Diversity-repeatability trade-off governs hierarchical cell fate coding in multicellular development*](https://www.chinaxiv.org/abs/202512.00064) (2025).
6. **[Method]** [*Ladderpath: An Efficient Algorithm for Revealing Nested Hierarchy in Sequences*](https://linkinghub.elsevier.com/retrieve/pii/S2590005625002905) (2025).
7. **[Network][AI]** [*Correlating Measures of Hierarchical Structures in Artificial Neural Networks with Their Performance*](https://www.nature.com/articles/s44260-024-00015-x) (2024).
8. **[Biology][Drug]** [*Discovery of Highly Bioactive Peptides through Hierarchical Structural Information and Molecular Dynamics Simulations*](https://pubs.acs.org/jcisd8/article-abstract/64/21/8164/852816/Discovery-of-Highly-Bioactive-Peptides-through) (2024).
9. **[Biology]** [*Evolutionary Tinkering Enriches the Hierarchical and Nested Structures in Amino Acid Sequences*](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.023215) (2024).
10. **[Method]** [*Ladderpath Approach: How Tinkering and Reuse Increase Complexity and Information*](https://www.mdpi.com/1099-4300/24/8/1082) (2022).

## Citation

If you use this software, cite the relevant Ladderpath paper and the software release. Machine-readable software citation metadata is provided in [CITATION.cff](CITATION.cff). A commonly used BibTeX citation for this release is:

```bibtex
@software{lppack_3_0_0,
  author  = {{ecsLab}},
  title   = {lppack: Ladderpath Python Toolkit},
  year    = {2026},
  version = {3.0.0},
  url     = {https://github.com/yuernestliu/lppack}
}
```

## License

Licensed under the [Apache License 2.0](LICENSE).

## Highlights in 3.0.0

- One public package and import path: `import lppack as lp`.
- CPU, NVIDIA CUDA, and Apple Silicon/MLX computation backends.
- Unicode-aware accelerated processing for non-ASCII sequences.
- Standard `lpjson` format version `3.0.0`.
- Tools for compression, ladderon ranking, segmentation, tokenization, species distances, and phylogeny.
- Installable directly from GitHub through standard Python packaging.
- English and Chinese user guides.
