# lp_developers

`lp_developers` 是梯径 Ladderpath 的 Python 开发仓库，用于把一个或多个字符串序列计算为标准 `lpjson`，并围绕 `lpjson` 提供梯图可视化、偏序多重集 POM、压缩、梯元排序、分词、距离矩阵和系统发生树等工具。

当前主要入口是：

```python
import lp_developers.ladderpath as lp

lpjson = lp.get_ladderpath(
    strs,
    fill_ladderons_STR=False,
    save_file_name=None,
    show_version=True,
    mode='cpu',
)
```

`mode` 只接受 `'cpu'`、`'gpu-cuda'`、`'macm'`。

## 文档

- 详细操作说明见 [操作手册.md](操作手册.md)。
- 完整版本记录见 [_backup_readme_files/CHANGELOG.md](_backup_readme_files/CHANGELOG.md)。
- 当前 `lpjson` 标准格式版本为 `v1.1`。

## 安装

目前推荐直接从 GitHub 克隆或下载源码后，将仓库所在目录加入 Python path。

```bash
git clone https://github.com/yuernestliu/lp_developers.git
```

如果 `lp_developers` 放在 `/Users/ecslab/codes/` 下，可这样导入：

```python
import sys

sys.path.append("/Users/ecslab/codes")

import lp_developers.ladderpath as lp
```

基础 CPU 计算主要使用 Python 标准库。按功能可选择安装以下依赖：

```bash
pip install graphviz tqdm
```

如果需要 CUDA GPU 后端：

```bash
pip install numpy numba pydivsufsort
```

如果需要 Mac M 系列 Apple GPU 后端，需要使用原生 arm64 Python 环境：

```bash
pip install numpy pydivsufsort mlx
```

如果需要系统发生树绘制：

```bash
pip install scipy matplotlib biopython
```

如果需要渲染梯图图片，除 Python 包 `graphviz` 外，还需要安装系统 Graphviz 程序。`brew install graphviz` 只是 macOS + Homebrew 用户的安装示例，不是必须使用 `brew`。

```bash
# macOS, Homebrew
brew install graphviz

# Ubuntu / Debian
sudo apt-get install graphviz

# Conda
conda install -c conda-forge graphviz python-graphviz
```

## 快速开始

### CPU 计算

```python
import lp_developers.ladderpath as lp

strs = ['ABCDBCDBCDCDEFEF', 'AA', 'DBCDBCA', 'BCDCDAEF']

lpjson = lp.get_ladderpath(
    strs,
    fill_ladderons_STR=False,
    save_file_name=None,
    show_version=True,
    mode='cpu',
)
```

如果输入中有重复 target，请使用 `dict[str, int]` 表示重复次数：

```python
strs = {
    'ABCDBCDBCDCDEFEF': 1,
    'AA': 3,
    'DBCDBCA': 2,
    'BCDCDAEF': 1,
}

lpjson = lp.get_ladderpath(strs, show_version=False, mode='cpu')
```

输入规则：

- `list[str]` 中不能有重复字符串。
- 如果有重复字符串，必须用 `dict[str, int]`。
- target 不能是单字符字符串。

### CUDA GPU 计算

CUDA 后端文件是 `ladderpath_cuda.py`，公开入口仍然通过 `ladderpath.py` 调用：

```python
lpjson = lp.get_ladderpath(strs, show_version=True, mode='gpu-cuda')
```

该模式需要可用的 NVIDIA CUDA 环境，以及 `numba`、`numpy`、`pydivsufsort`。
对于中文等非 ASCII 输入，程序会自动使用 Unicode 字符坐标，避免 UTF-8 字节坐标与字符坐标混用。

### Mac M 系列 GPU 计算

Mac M 系列后端文件是 `ladderpath_macm.py`，使用 Apple Silicon / MLX Metal：

```python
lpjson = lp.get_ladderpath(strs, show_version=True, mode='macm')
```

该模式需要 macOS、Apple Silicon 原生 arm64 Python 环境，以及 `mlx`、`numpy`、`pydivsufsort`。
对于中文等非 ASCII 输入，程序会自动使用 Unicode 字符坐标。

## 常用功能

### 查看三个指标

```python
index3 = lp.disp3index(lpjson)
```

返回或显示梯径度、有序度和规模度。

### POM 表示

```python
pom, pom_str = lp.POM_from_JSON(lpjson, display_str=False)
```

`pom` 是偏序多重集的字典表示，`pom_str` 是可读字符串表示。

### 画梯图

```python
g = lp.draw_laddergraph(
    lpjson,
    show_longer_than=0,
    style="ellipse",
    save_fig_name="laddergraph",
    figformat="pdf",
    render_fig=True,
    show_node_size=True,
)
```

如果只想输出 `.gv` 文件而不渲染图片：

```python
g = lp.draw_laddergraph(
    lpjson,
    save_fig_name="laddergraph",
    render_fig=False,
)
```

### 保存和读取 lpjson

```python
lp.save_ladderpath_json(lpjson, "result.json")
lpjson = lp.load_ladderpath_json("result.json")
```

`load_ladderpath_json()` 会自动把 JSON 中的字典 key 转回整数。

## `ladderpath_tools`

本仓库还包含若干围绕 `lpjson` 的工具模块：

- `compress.py`: 将 `lpjson` 编码为 streams，并进一步做熵编码压缩。
- `lambda_from_laddergraph.py`: 从已有 laddergraph / lpjson 中计算指定 leaves 的 lambda。
- `sortout.py`: 对梯元进行排序，支持 `contribution_folded`、`contribution_frozen`、`multiplicity`。
- `segment.py`: 基于词表进行分词切割。
- `tokenizer.py`: `LadderpathTokenizer`，封装 encode / decode。
- `species_distance.py`: 基于 `lpjson` 计算 species distance table。
- `phylo.py`: 基于距离矩阵生成系统发生树。
- `orderrate.py`: 计算 eta 的开发工具。

详细示例见 [操作手册.md](操作手册.md)。

## 版本、兼容性与确定性

当前主要版本：

- `ladderpath.py`: `v2.0.21`
- `helpers_laddergraph.py`: `v2.0.21`
- `ladderpath_cuda.py`: `v1.4`
- `ladderpath_macm.py`: `v1.1`
- `lpjson`: `v1.1`

确定性说明：

- CPU 后端目前不承诺完全确定性；在存在等价候选或内部 tie-breaking 的情况下，运行结果允许出现小差异，通常表现为 `lpjson` 结构或梯径度的轻微差别。
- CUDA 后端和 macm 后端当前目标是确定性的：同一输入在同一后端重复运行，应输出一致的完整 `lpjson`。
- CUDA 与 macm 后端目标是保持完整 `lpjson` 一致。
- CPU 后端和 GPU 后端在部分样本上可能存在很小差异；通常梯径度数值接近。

需要严格对比后端输出时，可使用 `benchmarks/cuda_macm_parity.py`。

## Benchmark 与测试

CUDA simulator 与 macm 的一致性检查示例：

```bash
NUMBA_ENABLE_CUDASIM=1 python benchmarks/cuda_macm_parity.py emit --backend cuda_sim
python benchmarks/cuda_macm_parity.py emit --backend macm
python benchmarks/cuda_macm_parity.py emit --backend cpu
python benchmarks/cuda_macm_parity.py compare
```

注意：CUDA simulator 只适合验证结果，不代表真实 CUDA GPU 速度。

## License

见 [LICENSE](LICENSE)。
