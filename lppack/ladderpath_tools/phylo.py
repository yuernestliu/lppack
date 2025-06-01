"""
Version 2.0.4
Authors: ecsLab (Yu Liu et al.)
Data: 2024.12.15

利用梯径计算距离矩阵，然后画 系统发生树 和 最小生成树
"""

import ladderpath as lp
import ladderpath_tools.lp_dis_laddergraph2points as lp_dis_laddergraph2points

import matplotlib.pyplot as plt
from typing import List
import networkx as nx

from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch

from Bio import Phylo
from Bio.Phylo import TreeConstruction
from io import StringIO
import os



class LADDERPATH_DISTANCE_MATRIX: #必须先排除重复序列、重复命名
    def __init__(self, strs:List[str], strsNames:List[str]=None,
        method='lp_dis_laddergraph2points',
        lpjson=None):

        if len(strs) != len(set(strs)):
            print('Wrong & Abort: There are duplications in strs!!!')
            return
        self.strs = strs

        if strsNames:
            if len(strsNames) != len(set(strsNames)):
                print('Wrong & Abort: There are duplications in strsNames!!!')
                return
            self.strsNames = strsNames
        else:
            self.strsNames = list(range(len(strs)))
        self.mat = self._cal_mat_lp22(method, lpjson)  # 在初始化时计算距离矩阵 List[float]


    def _cal_mat_lp22(self, method, lpjson): # 计算两两序列之间的梯径距离
        # method 只能是 lp_dis_laddergraph2points-absolute/relative, independent2points-absolute/relative
        dist_mat_lp22 = [] # 距离矩阵的1维压缩形式，即对称距离矩阵的上三角
        n = len(self.strs)
        if method.endswith('-absolute'): # 相对距离还是绝对距离
            absolute = True
        elif method.endswith('-relative'):
            absolute = False
        else:
            return

        if method.startswith('lp_dis_laddergraph2points-'):
            assert lpjson is not None
            L_cache = [lp_dis_laddergraph2points.lp_dis_laddergraph2points([s], lpjson)-1 for s in self.strs] # L=lambda-1 对单个target来说
            for i in range(n):
                L1 = L_cache[i]
                for j in range(i+1, n):
                    L2 = L_cache[j]
                    two_leaf = [self.strs[i], self.strs[j]]
                    L12 = lp_dis_laddergraph2points.lp_dis_laddergraph2points(two_leaf, lpjson)-2 # L=lambda-2 对2个target来说
                    dis = 2*L12 - (L1+L2) # 绝对距离
                    if absolute:
                        dist_mat_lp22.append(dis)
                    else:
                        dist_mat_lp22.append( dis/(L1+L2) )

        elif method.startswith('independent2points-'):
            # 缓存每个单个元素的 ladderpath 计算结果，避免重复调用
            L_cache = [lp.find_ladderpath([s])[0]-1 for s in self.strs] # L=lambda-1 对单个target来说
            for i in range(n):
                L1 = L_cache[i]
                for j in range(i+1, n):
                    L2 = L_cache[j]
                    L12 = lp.find_ladderpath( [self.strs[i], self.strs[j]] )[0] -2 # L=lambda-2 对2个target来说
                    dis = 2*L12 - (L1+L2)
                    if absolute:
                        dist_mat_lp22.append(dis)
                    else:
                        dist_mat_lp22.append( dis/(L1+L2) )
                        # dis = 1 - 2*(lambda1+lambda2-lambda12)/(lambda1+lambda2-2)

        else:
            print('!!!Wrong: method can only be: lp_dis_laddergraph2points-absolute (or -relative), independent2points-.')

        return dist_mat_lp22



    # 根据距离矩阵画系统发生树：Unweighted Pair-Group Method with Arithmetic Mean（UPGMA）
    def phylotree_upgma_newick(self, outgroupName=None, save_newick_file=None, save_fig_pdf_name=None):
        if save_newick_file is None:
            save_newick_file = '___temp_upgma.newick'
        if os.path.exists(save_newick_file):
            print('!!!Wrong: save_newick_file already exists.')
            return
        linkage_matrix = sch.linkage(self.mat, method='average', optimal_ordering=True)
        tree, _ = sch.to_tree(linkage_matrix, rd=True)  # 生成树对象
        tree.rooted = True  # 确保树是有根的
        newick_texts = _to_newick(tree, self.strsNames, tree.dist) + ";"
        with open(save_newick_file, 'w') as file:
            file.write(newick_texts)

        if outgroupName is not None:
            tree = Phylo.read(save_newick_file, "newick")
            outgroup = tree.find_any(name=outgroupName)  # 替换为你外群节点的名称
            tree.root_with_outgroup(outgroup)
            tree.rooted = True  # 确保树是有根的
            newick_io = StringIO() # 使用 StringIO 将树输出为 Newick 格式
            Phylo.write(tree, newick_io, "newick")
            newick_str = newick_io.getvalue() # 获取 Newick 字符串
            with open(save_newick_file, 'w') as file:
                file.write(newick_str)

        ax = draw_tree_from_newick(save_newick_file)
        plt.tight_layout()
        if save_fig_pdf_name:
            plt.savefig(save_fig_pdf_name+'.pdf', format='pdf', bbox_inches='tight')
        plt.show()

        if save_newick_file == '___temp_upgma.newick':
            os.remove(save_newick_file)


    # 根据距离矩阵画系统发生树：Neighbor-Joining（nj）
    def phylotree_nj_newick(self, outgroupName=None, save_newick_file=None, save_fig_pdf_name=None):
        if save_newick_file is None:
            save_newick_file = '___temp_nj.newick'
        if os.path.exists(save_newick_file):
            print('!!!Wrong: save_newick_file already exists.')
            return

        full_matrix = squareform(self.mat) # 转换为完整的矩阵
        lower_triangular_matrix = [list(full_matrix[i, :i+1]) for i in range(full_matrix.shape[0])] # 提取下三角部分
        dm = TreeConstruction.DistanceMatrix(self.strsNames, lower_triangular_matrix) # 创建DistanceMatrix 对象
        constructor = TreeConstruction.DistanceTreeConstructor()
        nj_tree = constructor.nj(dm)
        delete_inner_names(nj_tree) # 把内部节点的名称都删除

        if outgroupName is not None:
            outgroup = nj_tree.find_any(name=outgroupName)  # 替换为你外群节点的名称
            nj_tree.root_with_outgroup(outgroup)
        nj_tree.rooted = True  # 确保树是有根的

        newick_io = StringIO() # 使用 StringIO 将树输出为 Newick 格式
        Phylo.write(nj_tree, newick_io, "newick")
        newick_str = newick_io.getvalue() # 获取 Newick 字符串
        with open(save_newick_file, 'w') as file:
            file.write(newick_str)

        ax = draw_tree_from_newick(save_newick_file)
        plt.tight_layout()
        if save_fig_pdf_name:
            plt.savefig(save_fig_pdf_name+'.pdf', format='pdf', bbox_inches='tight')
        plt.show()

        if save_newick_file == '___temp_nj.newick':
            os.remove(save_newick_file)


    # 画最小生成树
    def minimum_spanning_tree(self, save_fig_pdf_name=None, show_weights=False):
        """
        计算并绘制最小生成树。距离矩阵是n*n的上三角矩阵, self.mat[i,j]表示分子i和j的距离。
        - show_weights: 是否显示边权重
        - save_fig_name: 是否保存图像, 若为None则不保存
        """
        # Kruskal算法
        def kruskal(nodes, edges):
            parent = {node: node for node in nodes}
            def find(node):
                if parent[node] != node:
                    parent[node] = find(parent[node])
                return parent[node]
            def union(node1, node2):
                parent[find(node2)] = find(node1)
            
            mst = []
            for node1, node2, weight in sorted(edges, key=lambda x: x[2]):
                if find(node1) != find(node2):
                    mst.append((node1, node2, weight))
                    union(node1, node2)
                    if len(mst) == len(nodes) - 1:
                        break
            return mst

        # 转换距离矩阵
        dis_mat = squareform(self.mat)
        if dis_mat.shape[0] != dis_mat.shape[1]:
            print('The input matrix is not a square matrix!')
            return

        nodes = self.strsNames  # 使用 strsNames 作为节点
        edges = [(nodes[i], nodes[j], dis_mat[i, j]) for i in range(len(nodes)) for j in range(i+1, len(nodes))]

        # 计算最小生成树
        mst_edges = kruskal(nodes, edges)
        mst_edges_with_weights = [(n1, n2, {'weight': round(w, 4)}) for n1, n2, w in mst_edges]

        # 创建无向图并绘制
        G = nx.Graph()
        G.add_edges_from(mst_edges_with_weights)

        pos = nx.kamada_kawai_layout(G)  # 使用Kamada-Kawai布局
        nx.draw(G, pos, node_color='green', node_size=300, font_size=13, font_color='black', 
                edge_color='grey', width=2, alpha=1, with_labels=True)

        if show_weights:
            nx.draw_networkx_edge_labels(G, pos, edge_labels=nx.get_edge_attributes(G, 'weight'))  # 显示权重

        if save_fig_pdf_name:
            plt.savefig(f'{save_fig_pdf_name}.pdf', format='pdf')
        plt.show()



# ============= associated functions =============
def _to_newick(tree_node, leaf_names, parent_dist, newick=''):
    if tree_node.is_leaf():
        return f"{leaf_names[tree_node.id]}:{parent_dist - tree_node.dist}{newick}"
    else:
        left = _to_newick(tree_node.get_left(), leaf_names, tree_node.dist)
        right = _to_newick(tree_node.get_right(), leaf_names, tree_node.dist)
        return f"({left},{right}):{parent_dist - tree_node.dist}{newick}"


def delete_inner_names(tree): # 把nj树的内部节点的名称都删除
    for clade in tree.find_clades():
        if not clade.is_terminal(): 
            clade.name = None  # 将内部节点的名称设为 None



# ============= 其他独立函数 =============
def draw_tree_from_newick(newickFile, outgroupName=None):  # newickFile必须带有.newick
    tree = Phylo.read(newickFile, "newick")
    if outgroupName is not None:
        outgroup = tree.find_any(name=outgroupName)  # 替换为你外群节点的名称
        tree.root_with_outgroup(outgroup)
    tree.rooted = True  # 确保树是有根的
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300) # figsize=(8, 6)
    Phylo.draw(tree, do_show=False, axes=ax)
    return ax
