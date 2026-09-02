"""
Version 2.0.21
Authors: ecsLab
Data: 2026.06.12
"""

from graphviz import Digraph
from .ladderpath import POM_from_JSON


VALID_STYLES = ("ellipse", "ellipse-OnlyShowTargetID", "box", "box-OnlyShowTargetID")


def ellipse_len(seq):  # 画梯图的associated函数
    # length = np.sqrt(len(seq))/2
    length = (len(seq) ** (1 / 3)) / 2
    return length


def _warn_if_too_many_ladderons(ladderons, show_longer_than, warning_n_ladderons_to_show):
    n_ladderons_to_show = sum(
        1 for infoList in ladderons.values()
        if len(infoList[2]) > show_longer_than
    )
    if n_ladderons_to_show > warning_n_ladderons_to_show:
        print(f'Warning: too many ladderons (>{n_ladderons_to_show}) to draw!!!')


def _normalize_style(style):
    if style not in VALID_STYLES:
        print('! Wrong: the parameter \'style\' can only be either \'ellipse\', \'ellipse-OnlyShowTargetID\', \'box\', \'box-OnlyShowTargetID\'')
        return None, False

    if style == "box-OnlyShowTargetID":
        return "box", True
    if style == "ellipse-OnlyShowTargetID":
        return "ellipse", True
    return style, False


def _build_lookup_tables(lpjson):
    targets = lpjson["targets"]
    ladderons = lpjson["ladderons"]
    targets_2ID = {val[2]: ID for ID, val in targets.items()}  # targets_strs按输入来，不一定是-1，-2，-3...
    ladderons_2ID = {val[2]: ID for ID, val in ladderons.items()}
    bbbs_2ID = {val: f"b{i}" for i, val in enumerate(lpjson["basic_building_blocks"])}  # 按bbb List中的顺序来，赋予其ID b1, b2,...

    multi_from_pom, _ = POM_from_JSON(lpjson)
    multi_bbb = multi_from_pom[0]
    multi_ladderons = {}
    for level, sub_dict in multi_from_pom.items():
        if level == 0:
            continue
        multi_ladderons.update(sub_dict)
    return targets, ladderons, targets_2ID, ladderons_2ID, bbbs_2ID, multi_bbb, multi_ladderons


def _target_names(TargetNames_UserDefined, targets_2ID):
    if TargetNames_UserDefined is not None:
        return TargetNames_UserDefined

    names = {}
    for target, ID in targets_2ID.items():
        names[ID] = str(ID)
    return names


def _new_graph(rankdir, figsize):
    g = Digraph()
    if figsize is None:
        g.attr(rankdir=rankdir)
    else:
        g.attr(rankdir=rankdir, size=figsize)
    return g


def _ellipse_attrs(seq):
    length = ellipse_len(seq)
    return {"width": str(length), "height": str(length / 2)}


def _node_size_suffix(seq, show_node_size):
    return f"[{len(seq)}]" if show_node_size else ""


def _multiplicity_suffix(multiplicity):
    return f" ({multiplicity})" if multiplicity > 1 else ""


def _target_label(target, ID, ladderons_2ID, TargetNames_UserDefined, onlyshowtargetid, show_node_size):
    if target in ladderons_2ID:  # 说明这个target也是梯元
        return f"& {ladderons_2ID[target]}"

    target_name = TargetNames_UserDefined[ID] if onlyshowtargetid else target
    return f"{target_name}{_node_size_suffix(target, show_node_size)}"


def _basic_block_label(bbb, multi_bbb):
    return f"{bbb}{_multiplicity_suffix(multi_bbb[bbb])}"


def _ladderon_label(ldID, ladderon_self, multiplicity, style, show_node_size):
    base_label = ladderon_self if style == "box" else str(ldID)
    return f"{base_label}{_node_size_suffix(ladderon_self, show_node_size)}{_multiplicity_suffix(multiplicity)}"


def _draw_box_targets(g, targets_2ID, ladderons_2ID, TargetNames_UserDefined, onlyshowtargetid, show_node_size):
    g.attr("node", shape="box")
    for target, ID in targets_2ID.items():  # draw all target sequences
        label = _target_label(target, ID, ladderons_2ID, TargetNames_UserDefined, onlyshowtargetid, show_node_size)
        g.node(str(ID), label=label)


def _draw_ellipse_targets(g, targets_2ID, ladderons_2ID, TargetNames_UserDefined, onlyshowtargetid, show_node_size):
    g.attr("node", shape="ellipse")
    for target, ID in targets_2ID.items():  # draw all target sequences
        label = _target_label(target, ID, ladderons_2ID, TargetNames_UserDefined, onlyshowtargetid, show_node_size)
        if onlyshowtargetid:
            g.node(str(ID), label=label)
        else:
            g.node(str(ID), label=label, **_ellipse_attrs(target))


def _draw_basic_block_nodes(g, style, bbbs_2ID, multi_bbb, color):
    for bbb, IDstr in bbbs_2ID.items():
        g.node(IDstr, label=_basic_block_label(bbb, multi_bbb), shape="hexagon", color=color)


def _draw_targets(g, style, targets_2ID, ladderons_2ID, bbbs_2ID, multi_bbb, 
    TargetNames_UserDefined, onlyshowtargetid, show_longer_than, color, show_node_size):
    if style == "box":  # a detailed version, showing ladderons length
        _draw_box_targets(g, targets_2ID, ladderons_2ID, TargetNames_UserDefined, onlyshowtargetid, show_node_size)
    else:  # style == "ellipse"
        _draw_ellipse_targets(g, targets_2ID, ladderons_2ID, TargetNames_UserDefined, onlyshowtargetid, show_node_size)

    if show_longer_than == 0:  # display basic building blocks or not
        _draw_basic_block_nodes(g, style, bbbs_2ID, multi_bbb, color)


def _draw_ladderons(g, ladderons, multi_ladderons, style, show_longer_than, color, show_node_size):
    for ldID, infoList in ladderons.items():
        if len(infoList[2]) <= show_longer_than:
            continue

        edge_counts = {}
        for linkedTo, positions in infoList[3].items():  # 遍历POS
            edge_counts[str(ldID), str(linkedTo)] = len(positions)

        if not edge_counts:
            continue

        ladderonSelf = infoList[2]
        multi = multi_ladderons[ladderonSelf]
        if style == "box":
            g.node(str(ldID), label=_ladderon_label(ldID, ladderonSelf, multi, style, show_node_size), 
                style="filled", color=color)
        else:  # style == "ellipse"
            lenk = ellipse_len(ladderonSelf)
            g.node(
                str(ldID),
                label=_ladderon_label(ldID, ladderonSelf, multi, style, show_node_size),
                width=str(lenk),
                height=str(lenk / 2),
                style="filled",
                color=color,
            )

        for id_linkedto, times in edge_counts.items():
            for _ in range(times):
                g.edge(id_linkedto[0], id_linkedto[1], color=color)


def _draw_duplication_edges(g, lpjson, color):
    for ID, pos in lpjson["duplications_info"].items():  # 处理targets中有重复的情况
        if isinstance(pos, int):
            thismulti = pos - 1  # 字典input形式。duplications_info记录的就是重复次数
        else:
            thismulti = len(pos) - 1  # list input形式。duplications_info记录的是重复出现的位置
        for _ in range(thismulti):
            g.edge(str(lpjson["targets"][ID][0][0]), str(ID), color=color)


def _draw_basic_block_edges(g, targets, ladderons, bbbs_2ID, color):
    for ldID, infoList in ladderons.items():
        for comp0 in infoList[0]:  # COMP
            if isinstance(comp0, str):
                for bbb in comp0:
                    g.edge(bbbs_2ID[bbb], str(ldID), color=color)

    for ldID, infoList in targets.items():
        if infoList[2] not in ladderons:  # 如果既是target又是ladderon，这里就不画了，因为上面画过了
            for comp0 in infoList[0]:  # COMP
                if isinstance(comp0, str):
                    for bbb in comp0:
                        g.edge(bbbs_2ID[bbb], str(ldID), color=color)


def _render_if_requested(g, save_fig_name, figformat, cleanGVfile, render_fig):
    if save_fig_name:
        if not save_fig_name.endswith(".gv"):
            save_fig_name += ".gv"
        if render_fig:
            g.render(filename=save_fig_name, format=figformat, cleanup=cleanGVfile)
        else:
            g.save(filename=save_fig_name)


# 画梯图
def draw_laddergraph(lpjson, show_longer_than=0, 
    style="ellipse",
    TargetNames_UserDefined=None,
    warning_n_ladderons_to_show=500,
    rankdir="BT", 
    color="grey", 
    figsize=None,
    save_fig_name=None, 
    figformat="pdf", 
    cleanGVfile=True,
    render_fig=True,
    show_node_size=True):
    # Draw the laddergraph.
    # "show_longer_than": When the length of the ladderon > show_longer_than, this ladderon will be displayed.
    #     Note that "show_longer_than" should always be >= 1, and the basic building blocks are also omitted.
    # "style" dictates how the laddergraph is displayed. It can either be
    #     "ellipse" (the sequence won't be displayed, but the size of the ellipse is positively related to the length of the sequence),
    #     "ellipse-OnlyShowTargetID" (avoid the ellipse of the target being too large),
    #     or "box" (the sequence will be displayed),
    #     or "box-OnlyShowTargetID" (the sequence will be displayed, but only show the target's ID).
    # "TargetNames_UserDefined"：字典，给出每个target代号的对应名字
    # "warning_n_ladderons_to_show": The largest number of ladderons allowed to draw.
    # "rankdir" is the order of the nodes, should be BT (from bottom to top), TB, LR, RL.
    # "color" can be "grey", "red", "#808080", etc.
    # "figsize" is size of graph, format: "6,6", default could be "8,8"
    # "save_fig_name" is the file name of the figure. Defalut with '.gv'
    # "figformat" is the format of the exported figure, which could be "pdf", "png", etc.
    # "cleanGVfile" 是否删除render图片的.gv格式文件
    # "render_fig" 是否调用Graphviz渲染图片；若为False，只保存.gv源文件
    # "show_node_size" 是否在target和ladderon的label中显示节点长度

    ladderons = lpjson["ladderons"]
    _warn_if_too_many_ladderons(ladderons, show_longer_than, warning_n_ladderons_to_show)

    style, onlyshowtargetid = _normalize_style(style)
    if style is None:
        return

    targets, ladderons, targets_2ID, ladderons_2ID, bbbs_2ID, multi_bbb, multi_ladderons = _build_lookup_tables(lpjson)
    TargetNames_UserDefined = _target_names(TargetNames_UserDefined, targets_2ID)

    g = _new_graph(rankdir, figsize)
    _draw_targets(
        g,
        style,
        targets_2ID,
        ladderons_2ID,
        bbbs_2ID,
        multi_bbb,
        TargetNames_UserDefined,
        onlyshowtargetid,
        show_longer_than,
        color,
        show_node_size,
    )
    _draw_ladderons(g, ladderons, multi_ladderons, style, show_longer_than, color, show_node_size)
    _draw_duplication_edges(g, lpjson, color)

    if show_longer_than == 0:  # display the links from basic building blocks
        _draw_basic_block_edges(g, targets, ladderons, bbbs_2ID, color)

    _render_if_requested(g, save_fig_name, figformat, cleanGVfile, render_fig)
    if render_fig:
        return g
    else:
        return None
