"""
Version 2.0
Authors: ecsLab (Yu Liu et al.)
Data: 2024.09.13

利用梯径计算距离矩阵，然后画 系统发生树 和 最小生成树
"""

import ladderpath as lp
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from typing import List
import networkx as nx
from scipy.spatial.distance import squareform



class LADDERPATH_DISTANCE_MATRIX: #必须先排除重复序列、重复命名
    def __init__(self, strs:List[str], strsNames:List[str]=None): 
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
        self.mat = self._cal_mat_lp22()  # 在初始化时计算距离矩阵 List[float]


    def _cal_mat_lp22(self): # 计算两两序列之间的梯径距离
        dist_mat_lp22 = [] # 距离矩阵的1维压缩形式，即对称距离矩阵的上三角
        n = len(self.strs)
        # 缓存每个单个元素的 ladderpath 计算结果，避免重复调用
        lambda_cache = [lp.find_ladderpath([s])[0] for s in self.strs]
        for i in range(n):
            lambda1 = lambda_cache[i]
            for j in range(i+1, n):
                lambda2 = lambda_cache[j]
                lambda12 = lp.find_ladderpath( [self.strs[i], self.strs[j]] )[0]
                dis = 1 - 2*(lambda1+lambda2-lambda12)/(lambda1+lambda2-2)
                dist_mat_lp22.append(dis)
        return dist_mat_lp22



    # 画系统发生树
    def phylotree(self, save_fig_name=None, figformat='pdf'):
        sch.dendrogram(sch.linkage(self.mat, method='average', optimal_ordering=True),
                       leaf_rotation=0, 
                       orientation='left', 
                       labels=self.strsNames
                       )
        plt.tight_layout()
        if save_fig_name:
            plt.savefig(f'{save_fig_name}.{figformat}', format=figformat)
        plt.show()



    # 画最小生成树
    def minimum_spanning_tree(self, save_fig_name=None, figformat='pdf', show_weights=False):
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

        if save_fig_name:
            plt.savefig(f'{save_fig_name}.{figformat}', format=figformat)
        plt.show()