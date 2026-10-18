# lppack

**梯径 Ladderpath Python 工具包 — 版本 3.0.0**

[English README](README.md) · [中文操作手册](docs/user-guide.zh-CN.md) · [English user guide](docs/user-guide.en.md) · [版本记录](docs/CHANGELOG.md) · [引用信息](CITATION.cff)

`lppack` 用于发现一个或多个序列中可复用的嵌套结构，并将结果表示为标准 `lpjson`。本项目包括梯径主算法、梯图可视化、偏序多重集 POM、压缩、梯元排序、基于词表的分词与 tokenizer、物种距离和系统发生树等工具。

## 环境要求

- Python 3.10 或更高版本。
- 安装 `lppack` 时会一次安装全部功能所需的 Python 依赖。
- CUDA 和 Apple Silicon 加速仍要求相应的硬件及系统环境。

## 安装

直接从 GitHub 安装当前版本：

```bash
python -m pip install git+https://github.com/yuernestliu/lppack.git
```

如需修改源码或使用 editable 安装：

```bash
git clone https://github.com/yuernestliu/lppack.git
cd lppack
python -m pip install -e .
```

上述两种安装方式都会自动安装 NumPy、Numba、pydivsufsort、兼容 Apple Silicon 环境中的 MLX、Graphviz Python binding、tqdm、Jieba、SciPy、Matplotlib 和 Biopython。

`pip` 安装的 `graphviz` 是 Python 接口，它不会安装 Graphviz 系统程序，也不包含其中的 `dot` 可执行文件。梯径的主体计算不使用 `dot`：即使没有安装 Graphviz 系统程序，或者系统程序没有安装成功，仍然可以正常计算和查看梯径结果。只有把梯图渲染为 PDF、PNG 等图片时，才需要 Graphviz 系统程序。

如果需要渲染图片，请另外安装 Graphviz 系统程序：

```bash
# macOS + Homebrew
brew install graphviz

# Ubuntu / Debian
sudo apt-get install graphviz

# Conda
conda install -c conda-forge graphviz python-graphviz
```

Windows 用户可以先通过 Windows Package Manager 安装：

```powershell
winget install graphviz
```

也可以从 [Graphviz 官方下载页面](https://graphviz.org/download/) 下载 64 位 EXE 安装程序。安装完成后，打开一个新的命令提示符或 PowerShell 窗口，运行 `dot -V`。如果 Windows 提示无法识别 `dot` 命令，请按下面的步骤手动把 Graphviz 加入 `PATH`：

1. 找到 Graphviz 的 `bin` 文件夹，通常是 `C:\Program Files\Graphviz\bin`。
2. 打开开始菜单，搜索“编辑系统环境变量”（英文系统中为 **Edit the system environment variables**），然后点击“环境变量”。
3. 在“用户变量”中选择 `Path`，依次点击“编辑”和“新建”，填入上述 `bin` 文件夹路径。
4. 确认并关闭所有窗口，重新打开终端或 IDE，然后再次运行 `dot -V`。

## 快速开始

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

三个指标分别为：

- **梯径度（Ladderpath index）**：梯径表示的规模。
- **有序度（Order index）**：`size-index - ladderpath-index`。
- **规模度（Size index）**：所有 targets 的总长度，包含重复次数。

### 重复 targets

`list[str]` 中不能包含重复字符串。重复 targets 要用保持插入顺序的 `dict[str, int]` 表示：

```python
targets = {
    "ABCDBCDBCDCDEFEF": 1,
    "AA": 3,
    "DBCDBCA": 2,
    "BCDCDAEF": 1,
}

lpjson = lp.get_ladderpath(targets, show_version=False, mode="cpu")
```

每个 target 至少包含两个字符。

## 计算后端

| Mode | 平台 | 说明 |
| --- | --- | --- |
| `cpu` | 支持 Python 的通用平台 | 默认后端，只使用标准库 |
| `gpu-cuda` | 已正确配置 CUDA 的 NVIDIA GPU | 使用 Numba、NumPy 和 pydivsufsort |
| `macm` | Apple Silicon macOS 原生 arm64 Python | 使用 MLX Metal、NumPy 和 pydivsufsort |

三种后端使用同一个公开入口：

```python
lpjson_cuda = lp.get_ladderpath(targets, mode="gpu-cuda")
lpjson_mac = lp.get_ladderpath(targets, mode="macm")
```

只有明确指定时才会使用硬件加速。如果指定的后端不可用，程序会报错，不会静默退回 CPU。

## 常用操作

### 保存和读取 `lpjson`

```python
lp.save_ladderpath_json(lpjson, "result.json")
lpjson = lp.load_ladderpath_json("result.json")
```

需要恢复 `targets`、`ladderons` 和位置映射中的整数 ID 时，应使用 `load_ladderpath_json()`，不要直接使用 `json.load()`。

### 获取 POM 表示

```python
pom, pom_text = lp.POM_from_JSON(lpjson, display_str=True)
```

### 绘制梯图

```python
graph = lp.draw_laddergraph(
    lpjson,
    save_fig_name="laddergraph",
    figformat="pdf",
    render_fig=True,
    show_node_size=True,
)
```

只有渲染图片时才需要 Graphviz 系统程序；所有梯径主体计算在没有该系统程序时也能正常运行。它的 Python binding 已随 `lppack` 安装。

## 工具模块

| 模块 | 功能 |
| --- | --- |
| `ladderpath_tools.compress` | 将 `lpjson` 编码为 streams 和熵编码 `.lpe1` 字节 |
| `ladderpath_tools.lambda_from_laddergraph` | 从已有梯图重新计算指定 targets 的 lambda |
| `ladderpath_tools.sortout` | 梯元排序并构建词表 |
| `ladderpath_tools.segment` | 基于词表分词 |
| `ladderpath_tools.tokenizer` | 高层 encode/decode 封装 |
| `ladderpath_tools.species_distance` | 计算 Jaccard 或 Dice 形式的物种距离 |
| `ladderpath_tools.phylo` | 构建 NJ 或 UPGMA 系统发生树 |

## 文档

- [中文操作手册](docs/user-guide.zh-CN.md)
- [English user guide](docs/user-guide.en.md)
- [版本记录](docs/CHANGELOG.md)
- [机器可读引用信息](CITATION.cff)

## `lpjson` 与可复现性

软件发行版本和 `lpjson` 格式版本均为 `3.0.0`。保存的对象会在 `lpjson["info"]` 中记录格式版本。

CPU 后端在内部出现等价候选时可能选择不同结果，因此重复运行不保证得到逐字节完全一致的 `lpjson`。CUDA 和 Apple Silicon 后端的设计目标是同一输入、同一后端产生确定性结果。如果工作流要求跨后端逐字节一致，应针对实际数据进行验证。

## 从 2.x 迁移

- 使用 `import lppack as lp`，不再直接导入 `ladderpath` 或私有开发包。
- 通过 `pip` 安装仓库，不再手动修改 `sys.path`。
- 通过 `mode="cpu"`、`mode="gpu-cuda"` 或 `mode="macm"` 选择后端。
- 重复字符串使用 `dict[str, int]` 表示，不能在 list 中重复出现。
- `lpjson` 格式标识更新为 `3.0.0`。

## 相关论文

Ladderpath 与以下论文相关：

1. **[Method][Linguistics]** [*Measuring language complexity from hierarchical reuse of recurring patterns*](https://arxiv.org/abs/2606.11531)（2026）。
2. **[LLM][Biology]** [*Compression-Based Tokenization Improves Language Modeling of Hierarchical Genomic Structure*](https://langtaosha.org.cn/lts/zh_Hans/preprint/view/51)（2025）。
3. **[Biology][Phylogenetics]** [*Structural Compression for Phylogenetic Inference under Alignment Instability and Indel-Rich Evolution*](https://arxiv.org/abs/2607.15308)（2026）。
4. **[Method][Compression]** [*Text Distance from Nested and Hierarchical Repetitions: A Compression-Based Perspective*](https://arxiv.org/abs/2607.05416)（2026）。
5. **[Biology][Development]** [*Diversity-repeatability trade-off governs hierarchical cell fate coding in multicellular development*](https://www.chinaxiv.org/abs/202512.00064)（2025）。
6. **[Method]** [*Ladderpath: An Efficient Algorithm for Revealing Nested Hierarchy in Sequences*](https://linkinghub.elsevier.com/retrieve/pii/S2590005625002905)（2025）。
7. **[Network][AI]** [*Correlating Measures of Hierarchical Structures in Artificial Neural Networks with Their Performance*](https://www.nature.com/articles/s44260-024-00015-x)（2024）。
8. **[Biology][Drug]** [*Discovery of Highly Bioactive Peptides through Hierarchical Structural Information and Molecular Dynamics Simulations*](https://pubs.acs.org/jcisd8/article-abstract/64/21/8164/852816/Discovery-of-Highly-Bioactive-Peptides-through)（2024）。
9. **[Biology]** [*Evolutionary Tinkering Enriches the Hierarchical and Nested Structures in Amino Acid Sequences*](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.023215)（2024）。
10. **[Method]** [*Ladderpath Approach: How Tinkering and Reuse Increase Complexity and Information*](https://www.mdpi.com/1099-4300/24/8/1082)（2022）。

## 引用

如果在研究中使用本软件，请引用软件版本，并根据实际使用的方法或应用，从上方[相关论文](#相关论文)中选择相应论文引用。机器可读的软件引用信息见 [CITATION.cff](CITATION.cff)。本版本常用的 BibTeX 引用格式如下：

```bibtex
@software{lppack_3_0_0,
  author  = {Liu, Yu and Zhang, Zecheng and Zhu, Hu and Deng, Xiaodong and Zhang, Jingwen and Que, Linjie and {ecsLab}},
  title   = {lppack: Ladderpath Python Toolkit},
  year    = {2026},
  version = {3.0.0},
  url     = {https://github.com/yuernestliu/lppack}
}
```

## 许可证

本项目使用 [Apache License 2.0](LICENSE)。

## 3.0.0 主要内容

- 统一的公开包名与导入方式：`import lppack as lp`。
- 支持 CPU、NVIDIA CUDA 和 Apple Silicon/MLX 三种计算后端。
- 加速后端支持中文等非 ASCII 序列的 Unicode 字符坐标。
- 标准 `lpjson` 格式版本为 `3.0.0`。
- 包含压缩、梯元排序、分词、tokenizer、物种距离和系统发生树工具。
- 支持通过标准 Python 打包方式直接从 GitHub 安装。
- 提供中文和英文操作手册。
