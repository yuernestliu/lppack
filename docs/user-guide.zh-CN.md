# lppack 3.0.1 中文操作手册

[English user guide](user-guide.en.md) · [项目 README](../README.zh-CN.md)

本文档说明 `lppack` 3.0.1 的公开接口和 `lpjson` 3.0.0 标准格式。

## 目录

- [1. 安装](#1-安装)
- [2. 核心概念](#2-核心概念)
- [3. 计算梯径](#3-计算梯径)
- [4. lpjson 3.0.0 格式](#4-lpjson-300-格式)
- [5. 查看、保存和读取结果](#5-查看保存和读取结果)
- [6. 绘制梯图](#6-绘制梯图)
- [7. CUDA 和 Apple Silicon 后端](#7-cuda-和-apple-silicon-后端)
- [8. 压缩](#8-压缩)
- [9. 从已有梯图计算 lambda](#9-从已有梯图计算-lambda)
- [10. 梯元排序与词表构建](#10-梯元排序与词表构建)
- [11. 分词与 tokenizer](#11-分词与-tokenizer)
- [12. 物种距离与系统发生树](#12-物种距离与系统发生树)
- [13. 可复现性与后端比较](#13-可复现性与后端比较)
- [14. 从 2.x 迁移](#14-从-2x-迁移)
- [15. 常见问题](#15-常见问题)

## 1. 安装

### 1.1 环境要求

- Python 3.10 或更高版本。
- 安装时会包含全部功能使用的 Python 依赖。
- 加速后端仍要求兼容的硬件和系统环境。

### 1.2 从 GitHub 安装

```bash
python -m pip install git+https://github.com/yuernestliu/lppack.git
```

### 1.3 克隆并以 editable 模式安装

```bash
git clone https://github.com/yuernestliu/lppack.git
cd lppack
python -m pip install -e .
```

安装后可以从任意目录启动 Python，不需要再把个人电脑上的仓库路径加入 `sys.path`。

### 1.4 自动安装的依赖

普通安装会一次安装 NumPy、Numba、pydivsufsort、Graphviz Python binding、tqdm、Jieba、SciPy、Matplotlib 和 Biopython。MLX 会在兼容的 macOS arm64 环境中自动安装。

### 1.5 Graphviz 系统程序

`pip` 安装的 `graphviz` 是绘图功能使用的 Python 接口，它不会安装 Graphviz 系统程序，也不包含其中的 `dot` 可执行文件。

梯径的主体计算不使用 Graphviz。即使没有安装 Graphviz 系统程序，或者该系统程序没有安装成功，梯径的计算、保存、读取和结果查看仍然可以正常运行。只有需要把梯图渲染为 PDF、PNG 等图片时，才需要安装系统程序：

```bash
# macOS + Homebrew
brew install graphviz

# Ubuntu / Debian
sudo apt-get install graphviz

# Conda
conda install -c conda-forge graphviz python-graphviz
```

Windows 用户最简单的安装方式是使用 Windows Package Manager：

```powershell
winget install graphviz
```

也可以从 [Graphviz 官方下载页面](https://graphviz.org/download/) 下载 64 位 EXE 安装程序。安装完成后，打开一个新的命令提示符或 PowerShell 窗口并运行：

```powershell
dot -V
```

如果 Windows 提示无法识别 `dot` 命令，请按下面的步骤手动把 Graphviz 加入 `PATH`：

1. 找到 Graphviz 的 `bin` 文件夹。通常是 `C:\Program Files\Graphviz\bin`；如果安装在其他位置，请使用实际路径。
2. 打开开始菜单，搜索“编辑系统环境变量”（英文系统中为 **Edit the system environment variables**）。
3. 点击“环境变量”，在“用户变量”中选择 `Path`，依次点击“编辑”和“新建”。
4. 填入 Graphviz `bin` 文件夹的完整路径，然后确认并关闭所有窗口。
5. 关闭并重新打开终端、IDE 或 notebook 程序，使其读取更新后的 `PATH`。
6. 再次运行 `dot -V`，确认配置成功。

### 1.6 验证安装

```bash
python -c "import lppack; print(lppack.__version__, lppack.LPJSON_VERSION)"
```

预期输出：

```text
3.0.1 3.0.0
```

## 2. 核心概念

- **Target**：输入序列；一次计算可以包含一个或多个 targets。
- **Basic building block**：输入字母表中的单个字符，是最基本单元。
- **梯元（Ladderon）**：由字符串片段和/或低层梯元构成的可复用序列结构。
- **梯径度（Ladderpath index）**：梯径表示的规模。
- **有序度（Order index）**：`size-index - ladderpath-index`。
- **规模度（Size index）**：所有 targets 的总长度，包含重复次数。
- **POM**：梯径层级结构的偏序多重集表示。
- **梯图（Laddergraph）**：targets、梯元及复用关系的图表示。
- **lpjson**：程序返回的标准 Python 字典表示。

## 3. 计算梯径

主要入口是 `lppack.get_ladderpath()`：

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

### 3.1 参数

- `targets0`：不含重复项的 `list[str]`，或者将 target 映射到重复次数的保持插入顺序的 `dict[str, int]`。
- `fill_ladderons_STR`：为 `True` 时填充每个梯元和 target 的完整字符串缓存；默认为 `False`，以减小结果体积。
- `save_file_name`：JSON 输出路径；为 `None` 时不写文件。
- `show_version`：是否显示软件版本和 lpjson 格式版本。
- `mode`：只能是 `"cpu"`、`"gpu-cuda"` 或 `"macm"`。

### 3.2 输入规则

- 每个 target 必须是长度至少为 2 的字符串。
- `list` 中不能包含重复字符串。
- 有重复 targets 时使用 `dict[str, int]`。
- 字典的插入顺序决定 target ID 顺序。

重复 targets 示例：

```python
targets = {
    "ABCDBCDBCDCDEFEF": 1,
    "AA": 3,
    "DBCDBCA": 2,
    "BCDCDAEF": 1,
}

lpjson = lp.get_ladderpath(targets, show_version=False, mode="cpu")
```

### 3.3 返回值

`get_ladderpath()` 返回内存中的 Python 字典。如果设置 `save_file_name`，同时将相同对象写为 JSON。

```python
print(lpjson["ladderpath-index"])
print(lpjson["order-index"])
print(lpjson["size-index"])
```

## 4. lpjson 3.0.0 格式

下面是一个精简的示意对象：

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

JSON 文件中的对象 key 是字符串。`load_ladderpath_json()` 会把数字 ID 恢复成 Python 整数。

### 4.1 顶层字段

- `info`：lpjson 格式版本，当前为 `"3.0.0"`。
- `ladderpath-index`：梯径度。
- `order-index`：规模度减去梯径度。
- `size-index`：包含重复次数的 targets 总长度。
- `basic_building_blocks`：表示中出现的最基本字符。
- `ladderons`：梯元 ID `0, 1, 2, ...` 到梯元记录的映射。
- `targets`：target ID `-1, -2, -3, ...` 到 target 记录的映射。
- `duplications_info`：字典输入中重复 target 的重数映射；没有重复时为空。
- `input_type`：`"list"` 或 `"dict"`。

### 4.2 梯元记录

每个梯元采用 `[COMP, LEN, STR, POS]`：

- `COMP`：按顺序排列的字符串片段和低层梯元 ID。
- `LEN`：完整字符长度。
- `STR`：完整字符串缓存；未填充时为空字符串。
- `POS`：上家 ID 到出现位置列表的映射。

### 4.3 Target 记录

每个 target 采用 `[COMP, LEN, STR, REP]`：

- `COMP`：按顺序排列的字符串片段和梯元 ID。
- `LEN`：target 长度。
- `STR`：完整字符串缓存；未填充时为空字符串。
- `REP`：该 target 的重复次数。

## 5. 查看、保存和读取结果

### 5.1 显示三个指标

```python
indices = lp.disp3index(lpjson)
# (ladderpath-index, order-index, size-index)
```

### 5.2 获取 POM

```python
pom, pom_text = lp.POM_from_JSON(lpjson, display_str=True)
```

`pom` 是按层级编号的字典。`display_str=True` 时，`pom_text` 是紧凑的字符串表示。

### 5.3 填充或清空字符串缓存

```python
lp.fill_lpjson_STR(lpjson)
lp.clear_lpjson_STR(lpjson)
```

两个函数都会直接修改传入对象。

### 5.4 保存与读取

```python
lp.save_ladderpath_json(lpjson, "result.json")
loaded = lp.load_ladderpath_json("result.json")
```

调用需要整数 target/梯元 ID 的工具前，应使用本项目的读取函数。

### 5.5 其他辅助函数

```python
target_id_list = lp.reconstruct_targets_list(lpjson)
basic_counts = lp.bbb_multiplicity(lpjson)
```

对于字典输入，重构的重复 target ID 会按照 target 插入顺序分组。

## 6. 绘制梯图

前面章节中的梯径主体计算不需要 Graphviz 系统程序。只有渲染梯图时才需要安装它；Python binding 已包含在普通安装中。

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

主要参数：

- `show_longer_than`：只显示长度大于此阈值的梯元。
- `style`：`ellipse`、`ellipse-OnlyShowTargetID`、`box` 或 `box-OnlyShowTargetID`。
- `TargetNames_UserDefined`：可选的 target ID 到显示名称映射。
- `rankdir`：`BT`、`TB`、`LR` 或 `RL`。
- `render_fig=False`：只保存 Graphviz 源文件，不渲染图片。
- `show_node_size`：在节点标签中附加序列长度。

该函数会填充传入对象中的 `STR` 缓存。

## 7. CUDA 和 Apple Silicon 后端

### 7.1 NVIDIA CUDA

系统必须提供与 Numba 兼容、可正常使用的 NVIDIA 驱动和 CUDA runtime。

```python
lpjson = lp.get_ladderpath(targets, mode="gpu-cuda")
```

### 7.2 Apple Silicon / MLX

在 Apple Silicon Mac 上使用原生 arm64 Python：

```python
lpjson = lp.get_ladderpath(targets, mode="macm")
```

MLX 后端要求 macOS、原生 arm64 Python 和可用的 Metal 后端；它不会静默退回 CPU。

### 7.3 非 ASCII 输入

只要任一输入包含非 ASCII 字符，加速后端就会自动使用 Unicode 字符坐标，避免把字符位置与 UTF-8 字节偏移混用。

## 8. 压缩

压缩模块先把 `lpjson` 转换为结构化 streams，再生成熵编码字节。

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

也可以一步完成：

```python
payload = lp_compress.compress2bit(
    lpjson,
    top_n=3000,
    sort_method="contribution_folded",
    token_type_mode="prev1_context",
    output_path="result.lpe1",
)
```

可用的排序方法是 `contribution_folded`、`contribution_frozen` 和 `multiplicity`。可用的 token type mode 是 `prev1_context`、`bitpack` 和 `auto_best`。

`.lpe1` 序列化格式有独立的格式标识，不等同于 `lppack` 或 `lpjson` 的版本号。

## 9. 从已有梯图计算 lambda

```python
from lppack.ladderpath_tools import lambda_from_laddergraph as lp_lambda

lambda_one = lp_lambda.get([-1], lpjson)
lambda_group = lp_lambda.get([-1, -3], lpjson)
```

也可以使用 target 字符串：

```python
lambda_group = lp_lambda.get(
    ["ABCDBCDBCDCDEFEF", "DBCDBCA"],
    lpjson,
)
```

限制：

- `leaves` 必须是非空 list，并且元素要么全是整数，要么全是字符串。
- 指定的 targets 必须存在于传入对象中。
- 当前实现不接受带 duplicates 的 `lpjson`。
- 从文件读取时必须使用 `load_ladderpath_json()`，不能直接使用 `json.load()`。

## 10. 梯元排序与词表构建

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

返回结果形式如下：

```python
{
    "the": [0, 4049],
    " and ": [1, 2663],
    "in": [2, 4682],
}
```

`iter_logs=True` 会开启详细进度诊断，速度可能明显变慢。公开排序方法是 `contribution_folded`、`contribution_frozen` 和 `multiplicity`。

## 11. 分词与 tokenizer

### 11.1 直接分词

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

`vocab` 可以是 CSV、TSV、JSON 文件路径，也可以是内存中的 list 或 dict。`text` 可以是单个字符串或字符串序列。

支持的方法：

- `max_matching`
- `pathpiece`
- `dag_dp_Importance`
- `dag_dp_logImportance`
- `jieba_dp`
- `jieba_dp_hmm`

Jieba 已包含在 `lppack` 的普通安装依赖中。

使用表格词表时，可以设置 `token_col`、`id_col`、`importance_col`、`has_header`、`delimiter` 和 `encoding`。

### 11.2 使用 `LadderpathTokenizer`

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

`encode()` 支持 `ids`、`tokens`、`str` 和 `details` 输出模式，也支持批量输入。`decode()` 可以接收一组 ID 或多组 ID。

## 12. 物种距离与系统发生树

### 12.1 构建距离表

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

`species_info` 必须使用真实 target ID。在整个 mapping 中，每个 target ID 的出现次数必须与源对象记录的重数一致。支持的距离方法是 `jaccard` 和 `dice`。

### 12.2 构建系统发生树

所需的 SciPy、Matplotlib 和 Biopython 会随 `lppack` 安装：

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

- `method` 为 `nj` 或 `upgma`。
- NJ 默认是无根树；可以通过 `outgroup_name` 重新定根。
- UPGMA 构建时即为有根树。
- Newick 输出支持 `.nwk` 和 `.newick`。
- 图片输出支持 PNG、PDF、SVG、JPEG 和 TIFF。

## 13. 可复现性与后端比较

- CPU 后端在内部出现等价候选时可能选择不同结果，重复运行不保证逐字节一致。
- CUDA 和 Apple Silicon 后端的设计目标是同一输入、同一后端产生确定性结果。
- 在某些 tie 情况下，CPU 与加速后端不保证完全相同。
- 如果工作流要求严格一致，应在目标环境中使用代表性输入比较完整 `lpjson`。
- CUDA simulator 可用于检查代码路径，但不能代表真实 GPU 性能。

## 14. 从 2.x 迁移

统一使用安装后的公开包：

```python
import lppack as lp
```

不再把私有仓库路径加入 `sys.path`。通过 `mode` 参数选择后端。list 中不再允许重复 target，应改用重数字典。软件版本为 `3.0.1`，`lpjson` 格式标识仍为 `3.0.0`。

旧版本保存的对象可能采用不同格式标识，不能默认与 3.0.0 逐字节兼容。在 3.0.0 工作流中使用前，应检查所需字段。

## 15. 常见问题

### `ModuleNotFoundError: No module named 'lppack'`

在仓库根目录执行 `python -m pip install -e .`，或者直接从 GitHub 安装。不要只把仓库自身加入 `sys.path`。

### 已安装 Python graphviz，但无法渲染

安装 Graphviz 系统程序，并确认 `dot` 在 `PATH` 中：

```bash
dot -V
```

### CUDA 模式在导入或初始化设备时失败

确认 Numba 可以正常导入，并且能够访问兼容的 NVIDIA CUDA runtime。

### `macm` 提示 MLX 或 Metal 不可用

请在 Apple Silicon macOS 上使用原生 arm64 Python。确认 `platform.machine()` 返回 `arm64`，并且 MLX 可以访问 Metal。

### list 中的重复 target 引发错误

把输入转换成 `dict[str, int]`，其中 value 是该 target 的总重复次数。

### 读取的 JSON 中 ID 是字符串

使用 `lppack.load_ladderpath_json()`，不要直接使用 `json.load()`。
