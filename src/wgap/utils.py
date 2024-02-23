import os
import sys
import re

def extract_chromosome_coordinates(input_str):
    # 使用正则表达式匹配字符串
    # 正则表达式解释：
    # (.*?) 匹配任何字符序列，直到遇到下一个模式的开始（非贪婪模式）
    # (\d+)_(\d+)$ 匹配字符串的最后两个由下划线分隔的数字部分，\d+ 匹配一个或多个数字
    match = re.match(r"(.*?)_(\d+)_(\d+)$", input_str)
    if match:
        # 从匹配对象中提取染色体编号（可能包含下划线）和起始结束位置
        prefix, start, end = match.groups()
        # 格式化字符串
        return prefix, int(start), int(end)
    else:
        return None