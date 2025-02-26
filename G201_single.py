import os
import pandas as pd
import numpy as np
import warnings
from Excel_Data_Loader import ExcelDataProcessor

# 禁用所有警告
warnings.filterwarnings('ignore')
target_path = input("请输入excel文件：")
target_path = target_path.strip("'").strip('"')
classify_by = input("请输入分类统计依据（project/tablet），留空表示不分类：").strip().lower()
cut_num=input("请输入过滤低质量数据的STR_AVG阈值")
sutter_cut_num=input("请输入过滤stutter高占比数阈值")
doublet_precent=input("是否计算STR双峰比？ y/n").strip().lower()
if cut_num == '':
    cut_num = 0
if sutter_cut_num=='':
    sutter_cut_num=0

processor = ExcelDataProcessor(target_path)
df, df_geno, df_genoDp, df_noDP = processor.batch_read_excel()

def calculate_statistics(subdf,project_tablet_name):
    """
    计算统计数据，并根据用户输入按 project 或 tablet 分类统计（如果用户输入有效且对应列存在）
    如果分类，则返回一个字典列表，每个字典对应一个组；否则返回一个单一的统计结果字典。
    """
    doublet_avg=0
    if classify_by in ["project","tablet"]:
        subdf =  df[df[classify_by] == project_tablet_name]
    else:
        subdf=df #考虑阈值为空情况

    if subdf.empty or subdf[subdf['STR_AVG'] >= int(cut_num)].empty or subdf[subdf['stutter高占比数'] <= int(sutter_cut_num)].empty:
        return {f'{classify_by}': project_tablet_name}
    subdf = subdf[subdf['STR_AVG'] >= int(cut_num)]
    subdf = subdf[subdf['stutter高占比数'] <= int(sutter_cut_num)]
    doublet_df=subdf[subdf['name']]
    if doublet_precent=="y":
        doublet_index = subdf.index
        real_index=len(doublet_index)
        doublet_avg=0
        for idx in doublet_index:
            avg = processor.process_doublet_precent(idx)
            if isinstance(avg, (float, np.float64)):
                if np.isnan(avg):  # 检查 avg 是否为 NaN
                    real_index -= 1  # 当 avg 为 NaN 时，减少有效的 real_index
                    doublet_df["doublet_avg"].append(avg)
                    continue
                else:
                    doublet_avg += avg  # 累加有效的 avg
                    doublet_df["doublet_avg"].append(round(avg,2))
            else:
                real_index -= 1  # 如果 avg 不是浮动类型，认为其无效，减少 real_index
                doublet_df["doublet_avg"].append(avg)
                continue
        if real_index > 0:
            doublet_avg=round(doublet_avg/real_index,2)
        else:
            doublet_avg = 0
        doublet_avg_str = str(doublet_avg) + '%'
    else:
        doublet_avg_str = 'N/A'
    # A_Typed / auto_loci_typed 平均值
    if 'A_Typed' in subdf.columns:
        normal_avg = round(subdf['A_Typed'].mean(), 2)
    elif 'auto_loci_typed' in subdf.columns:
        normal_avg = round(subdf['auto_loci_typed'].mean(), 2)
    else:
        normal_avg = None

    # X_typed
    if 'X_typed' in subdf.columns:
        x_ac_avg = round(subdf['X_typed'].mean(), 2)
    else:
        x_ac_avg = None

    # Y_Typed / y_loci_typed
    if 'Y_Typed' in subdf.columns:
        y_ac_avg = round(subdf['Y_Typed'].mean(), 2)
    elif 'y_loci_typed' in subdf.columns:
        y_ac_avg = round(subdf['y_loci_typed'].mean(), 2)
    else:
        y_ac_avg = None

    lane_num = len(subdf)
    total_reads_M = round(subdf['总reads'].sum() / 1_000_000)

    if '常核心基因座_Typed' in subdf.columns:
        core_gene_avg = round(subdf['常核心基因座_Typed'].mean(), 2)
    elif 'Auto_AlleleCount' in subdf.columns:
        core_gene_avg = round(subdf['Auto_AlleleCount'].mean(), 2)
    else:
        print("常核心基因座_Typed 列不存在！")
        core_gene_avg = None

    if 'Y核心基因座_Typed' in subdf.columns:
        y_core_avg = round(subdf['Y核心基因座_Typed'].mean(), 2)
    elif 'Y_AlleleCount' in subdf.columns:
        y_core_avg = round(subdf['Y_AlleleCount'].mean(), 2)
    else:
        print("常核心基因座_Typed 列不存在！")
        y_core_avg = None

    sample_reads_avg = round(subdf['总reads'].mean() / 1_000_000, 2)
    effect_reads_avg = round(subdf['有效reads比'].mean(), 2)
    effect_reads_avg_percentage = "{:.2f}%".format(effect_reads_avg * 100)
    # 计算 AVG 和 STD 指标（以指定区间列为例）
    average_value_A = subdf.loc[:, 'Amelogenin':'D21S1270'].mean(axis=1)
    average_value_X = subdf.loc[:, 'DXS10148':'HPRTB'].mean(axis=1)
    average_value_Y = subdf.loc[:, 'Y-indel':'Y-GATA-H4'].mean(axis=1)
    subdf['A.AVG'] = round(average_value_A)
    subdf['X.AVG'] = round(average_value_X)
    subdf['Y.AVG'] = round(average_value_Y)
    subdf['A.STD'] = round(subdf.loc[:, 'Amelogenin':'D21S1270'].std(axis=1) / average_value_A, 2)
    subdf['X.STD'] = round(subdf.loc[:, 'DXS10148':'HPRTB'].std(axis=1) / average_value_X, 2)
    subdf['Y.STD'] = round(subdf.loc[:, 'Y-indel':'Y-GATA-H4'].std(axis=1) / average_value_Y, 2)
    #reads = subdf.loc[:, 'Amelogenin':'Y-GATA-H4'].mean(axis=1)
    reads_df = subdf.loc[:, 'Amelogenin':'Y-GATA-H4']

    A_avg = round(subdf['A.AVG'].mean())
    X_avg = round(subdf['X.AVG'].mean())
    Y_avg = round(subdf['Y.AVG'].mean())
    A_std = round(subdf['A.STD'].mean(), 2)
    X_std = round(subdf['X.STD'].mean(), 2)
    Y_std = round(subdf['Y.STD'].mean(), 2)

    # 新增 STR 指标
    if 'STR_AVG' in subdf.columns:
        str_avg = round(subdf['STR_AVG'].mean())
    elif 'STR均值' in subdf.columns:
        str_avg = round(subdf['STR均值'].mean())
    else:
        str_avg = None

    if 'STR_STD' in subdf.columns:
        str_std = round(subdf['STR_STD'].mean(), 2)
    elif 'STR标准化STD' in subdf.columns:
        str_std = round(subdf['STR标准化STD'].mean(), 2)
    else:
        str_std = None

    # 构造结果字典
    result = {
        '统计个数': lane_num,
        f"{classify_by}": project_tablet_name,
        '总reads（M）': total_reads_M,
        '常位点（平均）': normal_avg,
        'X位点（平均）': x_ac_avg,
        'Y位点（平均）': y_ac_avg,
        '核心常（平均）': core_gene_avg,
        '核心Y（平均）': y_core_avg,
        '单个样本reads(M)': sample_reads_avg,
        'STR reads占比': effect_reads_avg_percentage,
        '常STR双峰比':doublet_avg_str,
        'STR.AVG': str_avg,
        'STR.STD': str_std,
        'A.AVG': A_avg,
        'X.AVG': X_avg,
        'Y.AVG': Y_avg,
        'A.STD': A_std,
        'X.STD': X_std,
        'Y.STD': Y_std
    }
    return result,doublet_df

# 目标路径（使用原始字符串）
# 提示用户输入目标路径
file_dir = os.path.dirname(target_path)
try:
    # 切换到文件所在目录
    os.chdir(file_dir)
    print(f"成功切换到工作目录: {os.getcwd()}")
except FileNotFoundError:
    print(f"指定的目录 {file_dir} 不存在。")
except PermissionError:
    print(f"没有权限访问目录 {file_dir}。")
# 调整后的列名（移除了不需要的列）
columns = [f"{classify_by}", '统计个数', '总reads（M）',
    '常位点（平均）', 'X位点（平均）', 'Y位点（平均）',
    '核心常（平均）', '核心Y（平均）', '单个样本reads(M)',
    'STR reads占比', 'STR.AVG', 'STR.STD',
    'A.AVG', 'X.AVG', 'Y.AVG', 'A.STD', 'X.STD', 'Y.STD','常STR双峰比'
]

columns1=['基因组','基因位点检测率%','测序深度','STR占比']
# 创建结果DataFrame
df_result = pd.DataFrame(columns=columns)
all_doublet_df = pd.DataFrame()
# 批量读取文件
if classify_by == "project":
    # 获取所有唯一的 project 值
    unique_projects = processor.df["project"].unique()

    for proj in unique_projects:
        # 对 excel_dfs[0] 进行切片，筛选出当前 project 的数据
        df_slice = processor.df[processor.df["project"] == proj]
        df_slice = df_slice[pd.to_numeric(df_slice.index, errors="coerce") < 1000]
        # 调用 calculate_statistics 对当前切片进行统计计算
        stats_result,doublet_df_temp = calculate_statistics(df_slice,proj)
        df_result_sub = pd.DataFrame([stats_result])
        if df_result_sub.empty:
            continue
        # 拼接结果
        df_result = pd.concat([df_result, df_result_sub], axis=0, ignore_index=True)
        if doublet_df_temp is not None:
            all_doublet_df = pd.concat([all_doublet_df, doublet_df_temp], axis=0, ignore_index=True)

elif classify_by == "tablet":
    # 获取所有唯一的 project 值
    unique_projects = processor.df["tablet"].unique()

    for proj in unique_projects:
        # 对 excel_dfs[0] 进行切片，筛选出当前 project 的数据
        df_slice = processor.df[processor.df["tablet"] == proj]
        df_slice = df_slice[pd.to_numeric(df_slice.index, errors="coerce") < 1000]
        # 调用 calculate_statistics 对当前切片进行统计计算
        stats_result = calculate_statistics(df_slice,proj)
        df_result_sub = pd.DataFrame([stats_result])
        if df_result_sub.empty:
            continue
        # 拼接结果
        df_result = pd.concat([df_result, df_result_sub], axis=0, ignore_index=True)

else:
        df_slice = processor.df[pd.to_numeric(processor.df.index, errors="coerce") < 1000]
        stats_result =pd.DataFrame([calculate_statistics(df_slice,None)])  # 调用已定义的统计函数
        df_result = pd.concat([df_result, stats_result], axis=0, ignore_index=True)

#shared_columns = df_genoDp.columns[1:].tolist()
excel_name = os.path.splitext(os.path.basename(target_path))[0]
df_results = processor.process_data_by_gene(classify_by,df.loc[:, 'Amelogenin':'Y-GATA-H4'])
output_path = os.path.join(file_dir, f"{excel_name}_统计结果.xlsx")
df_result.to_excel(output_path, index=False)
output_path_split_table = os.path.join(file_dir, f"{excel_name}_统计结果分基因.xlsx")
df_results.to_excel(output_path_split_table, index=False)
print(f"G201统计结果已保存到 {output_path}")
input("按任意键退出...")