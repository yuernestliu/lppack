# Version 1.0
# 25y0424

# ladderpaht-index可能还有问题？------------------


import json
import copy

def ladderpath_merge_group(lpjson, group_info, save_file_name=None):
    if isinstance(lpjson, str):
        with open(lpjson, 'r', encoding='utf-8') as f:
            ladderpath_data = json.load(f)
    else:
        ladderpath_data = copy.deepcopy(lpjson)

    ## 先处理new_groups
    renew_dic={}
    for id, positions in ladderpath_data['duplications_info'].items():
        for p in positions:
            renew_dic[-p-1]=id
    def replace_elements(l,d):
        l_new=[]
        for i in l:
            if i in d:
                l_new.append(d[i])
            else:
                l_new.append(i)
        return l_new

    for new_id,old_ids in group_info.items():
        updated_old_ids=replace_elements(old_ids,renew_dic)
        group_info[new_id]=updated_old_ids
    targets = ladderpath_data["targets"]

    # 创建新的 targets 数据结构
    new_targets = {}

    # 遍历每个新分组规则，合并旧的 targets
    for new_target, old_targets in group_info.items():
        merged_comp = []  # 用来存合并后的 COMP
        total_length = 0  # 记录新 target 的总长度

        # 合并每个旧 target
        for old_target in old_targets:
            if old_target not in targets:
                print(f"Key {old_target} not found in targets!")  # 如果旧 target 找不到，就跳过
                continue

            old_target_data = targets[old_target]
            merged_comp.extend(old_target_data[0])  # 把旧的 COMP 加进来
            total_length += old_target_data[1]  # 累加长度

        # 把新的 target 信息写进 new_targets
        new_targets[new_target] = [
            merged_comp,  # 合并后的 COMP
            total_length,  # 新的总长度
            "",  # 还是空字符串占位符
            1  # 重复次数保持不变
        ]

    # 替换掉原来的 targets，更新成新的 targets
    ladderpath_data["targets"] = new_targets

    # 计算新的 targets 的偏移量，方便更新 ladderons
    ## 这里的偏移量相当于一种字符计数，用来记录各个旧target在新target里面的起始位置。
    new_targets_offsets = {}
    for new_target, old_targets in group_info.items():
        offset = 0  # 累计偏移量
        for old_target in old_targets:
            if old_target not in new_targets_offsets:
                new_targets_offsets[old_target] = {}
            if old_target in targets:  # 仅处理存在的旧 target
                new_targets_offsets[old_target][new_target] = offset
                offset += targets[old_target][1]  # 累加长度
                # 形成一个字典：{旧target的ID:{新target的ID：起始位置}}

    # 定义一个函数来更新 ladderons 的 POS 数据
    def update_pos(ladderon_pos, new_targets_offsets):
        # 参数：单个梯元的POS，新的target偏移量
        new_posi = {}
        for old_target, positions in ladderon_pos.items():
            if old_target in new_targets_offsets:
                for new_target, base_offset in new_targets_offsets[old_target].items():
                    updated_positions = [base_offset + pos for pos in positions]  # 更新位点
                    if new_target not in new_posi:
                        new_posi[new_target]=[]
                    new_posi[new_target].extend(updated_positions)  # 合并更新后的位点
            else:
                # 如果找不到对应的新 target，就保留原来的 POS
                new_posi[old_target] = positions
        return new_posi

    # 遍历所有 ladderons，更新 POS 数据
    for ladderon_id, ladderon_data in ladderpath_data["ladderons"].items():
        original_pos = ladderon_data[3]  # 获取原始 POS 数据
        ladderon_data[3] = update_pos(original_pos, new_targets_offsets)  # 更新 POS
        ladderon_data[2]=''

    ladderpath_data['duplications_info']={}
    try:
        ladderpath_data['file_info']['group_merged'] = 1
    except:
        ladderpath_data['file_info'] = {}
        ladderpath_data['file_info']['group_merged'] = 1

    if save_file_name:
        with open(save_file_name+'.json', "w", encoding="utf-8", newline='') as f:
            json.dump(ladderpath_data, f, ensure_ascii=False, indent=4)

    return ladderpath_data
