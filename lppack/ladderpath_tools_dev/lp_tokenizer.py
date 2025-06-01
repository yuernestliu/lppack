"""
Version 1.0.1
Authors: ecsLab (Jing Wang, Yu Liu, et al.)
Data: 2025.03.24

# Notes (to-do?):
# 优化D值的计算思路，只更新词典中的受影响梯元
# 1000梯元20s算完
# 设定词表大小k，提前结束输出
"""



import sys
import json
from tqdm import tqdm


# ============================================================
# ====================第1部分： 排序的主程序 ====================
# ============================================================

def load_json_ladderons(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        lpjson = json.load(file)
        return lpjson['ladderons']

def save_json(data, file_path):
    with open(file_path, "w", encoding="utf-8") as file:
        json.dump(data, file, ensure_ascii=False, indent=4)

# 计算约化次数K
def calculate_importance(ladderons, ladderon_id, ladderon_importance, visited, selected):
    if ladderon_id in ladderon_importance:
        return ladderon_importance[ladderon_id]
    if ladderon_id in visited:
        return 0
    visited.add(ladderon_id)

    ladderon_info = ladderons[ladderon_id]
    pos = ladderon_info[3]
    importance = 0

    for target_id, occurrences in pos.items():
        occurrence_count = len(occurrences)

        if int(target_id) < 0:
            importance += occurrence_count
        else:
            target_importance = calculate_importance(ladderons, target_id, ladderon_importance, visited, selected)
            importance += occurrence_count * target_importance

    ladderon_importance[ladderon_id] = importance
    return importance

def process_element(selected_id, ladderons, selected_keys):
    total = 0
    if selected_id not in ladderons:
        return total

    components = ladderons[selected_id][0]

    for item in components:
        if isinstance(item, str):
            total += len(item)
        elif isinstance(item, int):
            item_str = str(item)
            if item_str in selected_keys:
                total += 1
            elif item_str in ladderons:
                total += process_element(item_str, ladderons, selected_keys)
    return total

def get_affected_ladderons(ladderons, selected_id):
    def collect_parents(ladderon_id, affected):
        if ladderon_id in affected:
            return
        affected.add(ladderon_id)
        parents = ladderons[ladderon_id][0]
        for parent_id in parents:
            parent_id = str(parent_id)
            if parent_id in ladderons:
                collect_parents(parent_id, affected)

    def collect_children(ladderon_id, affected):
        if ladderon_id in affected:
            return
        affected.add(ladderon_id)
        children = ladderons[ladderon_id][3]
        for child_id in children:
            if int(child_id) >= 0:
                child_id = str(child_id)
                if child_id in ladderons:
                    collect_children(child_id, affected)

    affected_parents = set()
    affected_children = set()

    # 收集子节点和父节点
    collect_children(str(selected_id), affected_children)
    collect_parents(str(selected_id), affected_parents)

    # 从受影响的集合中移除 selected_id
    affected_children.discard(selected_id)
    affected_parents.discard(selected_id)
    return affected_children, affected_parents

def update_ladderons(ladderons, selected_id):
    # 删除所有梯元指向已选中梯元的引用
    for ladderon_id, ladderon_info in ladderons.items():
        pos = ladderon_info[3]
        if selected_id in pos:
            del pos[selected_id]


def update_branch_importance(ladderons, affected_parents, ladderon_importance, selected):

    # 清空受影响梯元的 F 值缓存
    for ladderon_id in affected_parents:
        if ladderon_id in ladderon_importance:
            del ladderon_importance[ladderon_id]

    # 重新计算受影响梯元的 F 值
    for ladderon_id in affected_parents:
        if ladderon_id in ladderons:
            calculate_importance(ladderons, ladderon_id, ladderon_importance, set(), selected)


def select_ladderons_by_importance_optimized_new(ladderons, input_lpjson, vocab_size):
    selected = {}
    ladderon_importance = {}
    l_prime_cache = {}
    stable_ladderons = load_json_ladderons(input_lpjson)
    first_iteration = True  # 标记第一轮迭代

    # 预先计算所有梯元的 F 值并存入缓存
    for ladderon_id in ladderons.keys():
        calculate_importance(ladderons, ladderon_id, ladderon_importance, set(), selected)

    iteration_count = 0
    with tqdm(total=vocab_size, desc="Iterate_token_selection", unit=" iteration") as pbar:
        while ladderons:
            iteration_count += 1
            # print(f"Iteration: {iteration_count}")
            
            # 检查是否达到 vocab_size 的迭代次数
            if iteration_count > vocab_size:
                print(f"达到预设词表大小 ({vocab_size}), 停止排序")
                return selected  # 达到迭代次数，直接返回 selected

            # 计算当前所有梯元的 D 值
            d_prime_cache, fold_dict = calculate_d_value(selected, stable_ladderons)
            c_values = {}

            for ladderon_id, ladderon_info in ladderons.items():
                L_prime = l_prime_cache.get(ladderon_id, process_element(ladderon_id, ladderons, selected.keys()))
                F = ladderon_importance[ladderon_id]

                if first_iteration:
                    D = L_prime  # 第一轮，D = L'
                else:
                    D = cal_cache_d(ladderon_id, ladderons, selected, d_prime_cache, stable_ladderons)  # 第二轮开始用缓存计算
                
                # 计算 C 值
                C = (L_prime - 1) * F - D
                c_values[ladderon_id] = C

            first_iteration = False  # 进入第二轮

            if not c_values:
                break

            # 选择 C 值最大的梯元
            selected_id = max(c_values, key=c_values.get)
            selected[selected_id] = stable_ladderons[selected_id]
            # 找到选中梯元的受影响梯元
            affected_children, affected_parents = get_affected_ladderons(ladderons, selected_id)
            affected_ladderons = affected_children | affected_parents

            # **更新受影响梯元的 L' 和 F 的缓存**
            for affected_ladderon in affected_ladderons:
                if affected_ladderon in ladderons:
                    l_prime_cache[affected_ladderon] = process_element(affected_ladderon, ladderons, selected.keys())

            # **更新梯元的层级结构**
            update_ladderons(ladderons, selected_id)
            update_branch_importance(ladderons, affected_parents, ladderon_importance, selected)

            # **删除已选中的梯元**
            del ladderons[selected_id]
            pbar.update(1)

    return selected




def calculate_d_value(selected_for_fold, ladderons):
    """
    计算 selected_for_fold 中所有梯元的 D' 及其之和。
    仅返回：
    - d_prime_cache: 仅包含 selected_for_fold 中的梯元 D' 值 {梯元id: D' 值}
    - total_d_prime: selected_for_fold 中所有梯元的 D' 之和
    """
    d_prime_cache = {}  # 仅存储 selected_for_fold 内梯元的 D' 值
    total_d_prime = 0  # D' 之和
    selected_keys = set(selected_for_fold.keys())  # 选中的梯元集合
    temp_cache = {}  # 用于存储所有计算过的梯元 D'，但最终不会返回

    def process_element_d_value(item):
        """ 计算单个元素的 D' 值 """
        if isinstance(item, str):
            return len(item)
        elif isinstance(item, int):
            item_str = str(item)
            if item_str in selected_keys:  # 选中的梯元 D' 直接为 1
                return 1
            if item_str in temp_cache:  # 已计算，直接返回
                return temp_cache[item_str]
            if item_str in ladderons:  # 计算 D' 并缓存
                d_value = sum(process_element_d_value(sub_item) for sub_item in ladderons[item_str][0])
                temp_cache[item_str] = d_value  # 存入临时缓存
                return d_value
        return 0

    # 遍历 selected_for_fold 计算 D' 并存入 d_prime_cache
    for key in selected_keys:
        d_value = sum(process_element_d_value(element) for element in ladderons[key][0])
        d_prime_cache[key] = d_value  # 仅缓存 selected_for_fold 中梯元的 D'
        total_d_prime += d_value  # 累加 D'

    return d_prime_cache, total_d_prime  # 仅返回 selected_for_fold 的 D' 值及其总和



# 计算各种指针数
def cal_n_pointers(top_n_token, lpjson_file, lp_token_file):
    sys.setrecursionlimit(5000)

    # 步骤1: 从Moby.json中读取数据
    with open(lpjson_file, 'r', encoding='utf-8') as file:
        data = json.load(file)
    targets = data['targets']

    numLetters = len(data['basic_building_blocks']) #可以直接计算出 numLetters
    # print(numLetters)

    # 步骤2: 从MobyLadderons_sort.json文件中读取字典
    with open(lp_token_file, 'r', encoding='utf-8') as file:
        ladder_data = json.load(file)

    ladder_keys = list(ladder_data.keys())[:top_n_token]
    total_sum = sum(value[1] for key, value in ladder_data.items() if key in ladder_keys)
    # print("词典指针:", total_sum)
    Npointers_dict = total_sum + numLetters # 加上了最基本字符

    ladder_dict = {key: value[0] for key, value in ladder_data.items()}

    total = 0

    def process_element_pointer(item):
        nonlocal total
        if isinstance(item, str):
            total += len(item)
        elif isinstance(item, int):
            item_str = str(item)
            if item_str in ladder_keys:
                total += 1
            elif item_str in ladder_dict:
                for sub_item in ladder_dict[item_str]:
                    process_element_pointer(sub_item)

    for key, value in targets.items():
        first_element = value[0]
        for item in first_element:
            process_element_pointer(item)
    Npointers_target = total
    # print("原文梯径指针:", total)
    # print('总指针:', total_sum + total)

    Npointers_all = Npointers_dict + Npointers_target

    return Npointers_dict, Npointers_target, Npointers_all


# 折叠词表指针=sum(L')
def calculate_top_l_pri_sum(lp_token_file, lpjson_file, top_n_token):
    with open(lp_token_file, 'r', encoding='utf-8') as file:
        selected_ladderons = json.load(file)

    # 获取前 top_n_token 个梯元的键
    ladder_keys = list(selected_ladderons.keys())[:top_n_token]
    
    with open(lpjson_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    numLetters = len(data['basic_building_blocks'])
    
    # 初始化总 L_pri 值
    total_l_pri = 0
    cache = {}  # 创建缓存字典

    def process_element_L_pri(item, ladder_keys, ladder_dict):
        """ 计算单个元素的 L' 值，使用缓存避免重复计算 """
        if isinstance(item, str):
            return len(item)
        elif isinstance(item, int):
            item_str = str(item)
            if item_str in cache:  
                return cache[item_str]  # 直接返回缓存值，避免重复计算
            total = 0
            if item_str in ladder_keys:
                total += 1
            elif item_str in ladder_dict:
                for sub_item in ladder_dict[item_str][0]:  # 获取梯元的第一个元素（列表）
                    total += process_element_L_pri(sub_item, ladder_keys, ladder_dict)  # 递归处理
            cache[item_str] = total  # 存入缓存
            return total

        return 0  

    # 遍历选定的 top_n_token 梯元，计算 L'
    for key in ladder_keys:
        elements = selected_ladderons[key][0]
        total_l_pri += sum(process_element_L_pri(element, ladder_keys, selected_ladderons) for element in elements)
    
    return total_l_pri + numLetters

# 计算加入一个梯元后词典指针变化，即公式中的D值
def cal_cache_d(ladderon,ladderons,selected,d_cache,stable_ladderons):
    affected_children, affected_parents = get_affected_ladderons(stable_ladderons, ladderon)
    sum_d=0
    selected_stable=selected.copy()
    lapped = {k: v for k, v in selected_stable.items() if k in affected_children}
    selected[ladderon]=ladderons[ladderon]  
    ladderon_d=process_element(ladderon,stable_ladderons,selected.keys()) 
    for lap_id,lap_info in lapped.items():
        
        d2=process_element(lap_id,stable_ladderons,selected.keys())
        d1=d_cache[lap_id] 
        d_difference=d2-d1   # 单个梯元的D值变化    对于1号梯元 现为0-11 应为9-11 d2计算有误
        sum_d+=d_difference
    del selected[ladderon]
    return sum_d+ladderon_d



def select_token(input_lpjson,filename,vocab_size,lp_tokenjson=True ,show_selectedladderons=False):
    ladderons = load_json_ladderons(input_lpjson)
    selected_ladderons = select_ladderons_by_importance_optimized_new(ladderons,input_lpjson,vocab_size)
    if lp_tokenjson :
        if show_selectedladderons:
            print(selected_ladderons.keys())
        if filename[-5:] == '.json':
            save_json(selected_ladderons, filename)
        else:
            save_json(selected_ladderons, f'{filename}.json')
    else:
        if show_selectedladderons:
            print(selected_ladderons.keys())
    return selected_ladderons




# ============================================================
# ====================第2部分： Tokenizer类 ===================
# ============================================================
# 这是外部程序调用的入口

class Tokenizer:
    # 初始化分词器
    def __init__(self):
        self.token_dict = {}       # 初始化一个空词典

    
    def train(self, input_lpjson, vocab_size, save_json_file):
        """
        能够用我的排序方法输出筛选后的梯元,输入的lpjson应该是显式缓存了str的
        梯元排序得到一个json文件,是排序后的ladderons的信息
        同时把单字符写入token_dict中
        """
        with open(input_lpjson, "r", encoding="utf-8") as f:
            data = json.load(f)

        # 提取 basic_building_blocks 并写入 self.token_dict（从1开始编号）
        basic_blocks_list = data.get("basic_building_blocks", [])
        self.token_dict = {i + 1: token for i, token in enumerate(basic_blocks_list)}

        basic_blocks = len(basic_blocks_list)
        if vocab_size < basic_blocks:
            raise ValueError("vocab_size must be greater than or equal to the number of basic building blocks")

        # 调用原始排序+保存函数
        select_token(
            input_lpjson=input_lpjson,
            filename=save_json_file,
            vocab_size=vocab_size - basic_blocks,
            lp_tokenjson=True,
            show_selectedladderons=False
        )




    def save_token_dict(self, load_modelname, save_json_file):
        """
        从 load_modelname 中读取 token，在已有的 self.token_dict 基础上继续编号（保持 int -> str），然后保存为 JSON。
        """

        # 获取当前最大编号
        if self.token_dict:
            current_id = max(self.token_dict.keys()) + 1  # int keys
        else:
            self.token_dict = {}
            current_id = 1

        existing_values = set(self.token_dict.values())  # 防止重复 token

        # 读取已有模型中的 token
        with open(load_modelname, "r", encoding="utf-8") as f:
            full_token_data = json.load(f)

        for token_info in full_token_data.values():
            token_str = token_info[2]  # 假设 token_str 是 ladderon 字符串
            if token_str not in existing_values:
                self.token_dict[current_id] = token_str  # ✅ key 保持为 int
                existing_values.add(token_str)
                current_id += 1

        # 保存为 JSON 时，key 必须转为 str
        with open(save_json_file, "w", encoding="utf-8") as f:
            json.dump({str(k): v for k, v in self.token_dict.items()}, f, ensure_ascii=False, indent=4)



    # 读取保存的梯元Token
    def load_from_token_dict_file(self):
        print(self.token_dict)
        


    # 给出文本，返回编码序列
    def encode(self, text):
        """
        对输入文本进行编码，输出 token ID 和单字符列表。
        匹配失败时保留原字符作为 token。
        """
        token_str_to_id = {v: int(k) for k, v in self.token_dict.items()}
        sorted_tokens = sorted(token_str_to_id.keys(), key=len, reverse=True)

        encoded_result = []
        i = 0
        while i < len(text):
            matched = False
            for token in sorted_tokens:
                if text.startswith(token, i):
                    encoded_result.append(token_str_to_id[token])
                    i += len(token)
                    matched = True
                    break
            if not matched:
                # 没有匹配成功，保留当前字符
                encoded_result.append(text[i])
                i += 1

        return encoded_result
