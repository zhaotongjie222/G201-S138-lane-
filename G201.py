import os
import pandas as pd
import numpy as np
import warnings

# 禁用所有警告
warnings.filterwarnings('ignore')

def batch_read_excel(target_path):
    """
    批量读取指定路径下所有以 0x.xlsx 结尾的 Excel 文件
    :param target_path: 目标文件夹路径
    :return: 包含所有 DataFrame 的字典 {文件名: DataFrame}
    """
    excel_files = {}
    for root, dirs, files in os.walk(target_path):
        for file in files:
            if file.endswith('0x.xlsx') : #or file.startswith('（G201')
                file_path = os.path.join(root, file)
                try:
                    df = pd.read_excel(file_path,sheet_name="Lane", header=3,index_col="index")  # 修复路径错误
                    excel_files[file] = df
                except Exception as e:
                    print(f"读取文件 {file_path} 时出错: {e}")
    return excel_files


def split_read_excel(target_path):
    """
    批量读取指定路径下所有以 0x.xlsx 结尾的 Excel 文件，获取所有工作表名称，
    header=0，保存为 df 和 df_分表名称，并针对特定工作表进行数据清洗。
    :param target_path: 目标文件夹路径
    :return: 包含所有 DataFrame 的字典 {文件名: {分表名称: DataFrame}}
    """
    excel_files = {}

    for root, dir, files in os.walk(target_path):
        for file in files:
            # 检查文件名格式：以中文括号G201开头，0x.xlsx结尾
            if file.endswith('0x.xlsx'):  #file.startswith('（G201')
                file_path = os.path.join(root, file)
                try:
                    # 一次性读取所有工作表（第一行为列头）
                    sheets_dict = pd.read_excel(
                        file_path,
                        sheet_name=None,
                        header=0,
                        engine='openpyxl'
                    )

                    file_dfs = {}
                    for idx, (sheet_name, sheet_df) in enumerate(sheets_dict.items()):
                        if idx == 0:  # 跳过第一个工作表
                            continue

                        df_key = f'df_{sheet_name}'
                        processed_df = sheet_df.copy()  # 使用已读取的数据副本

                        # 数据清洗逻辑
                        if sheet_name in ['genoDp', 'noDP']:
                            # 删除第二列起 0值占比≥95% 的行
                            cols_to_check = processed_df.iloc[:, 1:]
                            zero_mask = (cols_to_check == 0).mean(axis=1) >= 0.95
                            processed_df = processed_df[~zero_mask]
                        elif sheet_name == 'geno':
                            # 删除第五列起全为空的行
                            cols_to_check = processed_df.iloc[:, 1:]
                            null_mask = cols_to_check.isnull().mean(axis=1)>=0.9
                            processed_df = processed_df[~null_mask]

                        file_dfs[df_key] = processed_df

                    excel_files[file] = file_dfs

                except Exception as e:
                    print(f"处理文件 {file_path} 失败: {str(e)}")

    return excel_files

def calculate_statistics(df):
    df = df[df['STR_AVG'] >= 0]
    global reads, reads_df
    reads_df = []  # 后续会重新赋值

    # 检查是否存在 A_Typed 列，如果存在则使用它，否则使用 auto_loci_typed 列
    if 'A_Typed' in df.columns:
        normal_avg = round(df['A_Typed'].mean(), 2)
    elif 'auto_loci_typed' in df.columns:
        normal_avg = round(df['auto_loci_typed'].mean(), 2)
    else:
        normal_avg = None

    # 检查是否存在 X_typed 列
    if 'X_typed' in df.columns:
        x_ac_avg = round(df['X_typed'].mean(), 2)
    else:
        x_ac_avg = None

    # 检查是否存在 Y_Typed 列
    if 'Y_Typed' in df.columns:
        y_ac_avg = round(df['Y_Typed'].mean(), 2)
    elif 'y_loci_typed' in df.columns:
        y_ac_avg = round(df['y_loci_typed'].mean(), 2)
    else:
        y_ac_avg = None

    lane_num = len(df) - 1
    total_reads_M = round(df['总reads'].sum() / 1_000_000)

    if '常核心基因座_Typed' in df.columns:
        core_gene_avg = round(df['常核心基因座_Typed'].mean(), 2)
    elif 'Auto_AlleleCount' in df.columns:
        core_gene_avg = round(df['Auto_AlleleCount'].mean(), 2)
    else:
        print("常核心基因座_Typed 列不存在！")
        core_gene_avg = None

    if 'Y核心基因座_Typed' in df.columns:
        y_core_avg = round(df['Y核心基因座_Typed'].mean(), 2)
    elif 'Y_AlleleCount' in df.columns:
        y_core_avg = round(df['Y_AlleleCount'].mean(), 2)
    else:
        print("常核心基因座_Typed 列不存在！")
        y_core_avg = None

    sample_reads_avg = round(df['总reads'].mean()/ 1_000_000, 2)
    effect_reads_avg = round(df['有效reads比'].mean(), 2)

    # 将平均值转换为百分比形式（保留两位小数）
    effect_reads_avg_percentage = "{:.2f}%".format(effect_reads_avg * 100)

    # 计算AVG和STD指标
    average_value_A = df.loc[:, 'Amelogenin':'D21S1270'].mean(axis=1)
    average_value_X = df.loc[:, 'DXS10148':'HPRTB'].mean(axis=1)
    average_value_Y = df.loc[:, 'Y-indel':'Y-GATA-H4'].mean(axis=1)
    # 将计算的平均值添加到A.AVG, X.AVG, Y.AVG列（取整）
    df['A.AVG'] = round(average_value_A)
    df['X.AVG'] = round(average_value_X)
    df['Y.AVG'] = round(average_value_Y)
    # STD指标除以对应的平均值，保留两位小数
    df['A.STD'] = round(df.loc[:, 'Amelogenin':'D21S1270'].std(axis=1) / average_value_A, 2)
    df['X.STD'] = round(df.loc[:, 'DXS10148':'HPRTB'].std(axis=1) / average_value_X, 2)
    df['Y.STD'] = round(df.loc[:, 'Y-indel':'Y-GATA-H4'].std(axis=1) / average_value_Y, 2)

    reads = df.loc[:, 'Amelogenin':'Y-GATA-H4'].mean(axis=1)
    reads_df = df.loc[:, 'Amelogenin':'Y-GATA-H4']

    A_avg = round(df['A.AVG'].mean())
    X_avg = round(df['X.AVG'].mean())
    Y_avg = round(df['Y.AVG'].mean())
    A_std = round(df['A.STD'].mean(), 2)
    X_std = round(df['X.STD'].mean(), 2)
    Y_std = round(df['Y.STD'].mean(), 2)

    # 新增STR指标
    if 'STR_AVG' in df.columns:
        str_avg = round(df['STR_AVG'].mean())
    elif 'STR均值' in df.columns:
        str_avg = round(df['STR均值'].mean())
    else:
        str_avg = None

    if 'STR_STD' in df.columns:
        str_std = round(df['STR_STD'].mean(), 2)
    elif 'STR标准化STD' in df.columns:
        str_std = round(df['STR标准化STD'].mean(), 2)
    else:
        str_std = None

    # 返回结果字典
    return {
        '表格名称': '',
        '统计个数': lane_num,
        '总reads（M）': total_reads_M,
        '常位点（平均）': normal_avg,
        'X位点（平均）': x_ac_avg,
        'Y位点（平均）': y_ac_avg,
        '核心常（平均）': core_gene_avg,
        '核心Y（平均）': y_core_avg,
        '单个样本reads(M)': sample_reads_avg,
        'STR reads占比': effect_reads_avg_percentage,
        'STR.AVG': str_avg,
        'STR.STD': str_std,
        'A.AVG': A_avg,
        'X.AVG': X_avg,
        'Y.AVG': Y_avg,
        'A.STD': A_std,
        'X.STD': X_std,
        'Y.STD': Y_std
    }

def calculate_missing_rate(df_geno, start_col):
    target_columns = df_geno.columns[start_col:]

    missing_rates = []
    for col in target_columns:
        valid_count = df_geno[col].count()  # 非缺失值数量
        total_count = len(df_geno)
        missing_rates.append(valid_count / total_count * 100)

    # 返回平均值或按需调整
    return missing_rates
# 目标路径（使用原始字符串）
# 提示用户输入目标路径
target_path = input("请输入目标路径：")
classiy_by=bool(input("按照project还是tablet分类统计？"))
try:
    # 切换工作目录
    os.chdir(target_path)
    print(f"成功切换到工作目录: {os.getcwd()}")
except FileNotFoundError:
    print(f"指定的目录 {target_path} 不存在。")
except PermissionError:
    print(f"没有权限访问目录 {target_path}。")
# 调整后的列名（移除了不需要的列）
columns = [
    '表格名称', '统计个数', '总reads（M）',
    '常位点（平均）', 'X位点（平均）', 'Y位点（平均）',
    '核心常（平均）', '核心Y（平均）', '单个样本reads(M)',
    'STR reads占比', 'STR.AVG', 'STR.STD',
    'A.AVG', 'X.AVG', 'Y.AVG', 'A.STD', 'X.STD', 'Y.STD'
]

columns1=['表格名称','基因组','基因位点检测率','测序深度','STR%']
# 创建结果DataFrame
df_result = pd.DataFrame(columns=columns)

# 批量读取文件
excel_dfs = batch_read_excel(target_path)

# 处理每个文件
for file_name, df in excel_dfs.items():
    stats_result = calculate_statistics(df)  # 使用修正后的函数名
    stats_result['表格名称'] = file_name
    df_result = pd.concat([df_result, pd.DataFrame([stats_result])], ignore_index=True)

df_results = pd.DataFrame(columns=columns1)
excel_dfs = split_read_excel(target_path)

for file_name, df_dict in excel_dfs.items():
    try:
        # 提取必要的数据表
        df_genoDp = df_dict['df_genoDp']
        df_noDP = df_dict['df_noDP']
        df_geno = df_dict['df_geno']

        # 确保列对齐：使用 df_genoDp 的列作为基准
        shared_columns = df_genoDp.columns[1:].tolist()  # 假设第一列为基因名，后续为样本

        # 计算 df_Y_n_sum
        df_Y_n_sum = df_genoDp[shared_columns] + df_noDP[shared_columns]

        # 将 reads_df 中的零值替换为 NaN，避免除零错误
        reads_df_nonzero = reads_df.replace(0, np.nan)

        # 计算 str_df，避免除以零
        str_df = df_Y_n_sum.div(reads_df_nonzero)

        # 计算基因位点检测率
        geno_sample_columns = df_geno.columns[5:5 + len(shared_columns)]  # 动态切片
        missing_rate = (1 - df_geno[geno_sample_columns].isnull().mean(axis=0)).values

        # 构建当前文件的数据字典
        file_data = {
            '表格名称': [file_name] * len(shared_columns),
            '基因组': shared_columns,
            '测序深度': np.round(reads_df.mean(axis=0), 2),
            '基因位点检测率': np.round(missing_rate * 100, 2),  # 转为百分比
            'STR%': np.round(str_df.mean(axis=0), 4)
        }

        # 创建临时 DataFrame 并合并
        df_results1 = pd.DataFrame(file_data)
        df_results = pd.concat([df_results, df_results1], axis=0, ignore_index=True)

    except KeyError as e:
        print(f"文件 {file_name} 缺少必要列: {str(e)}")
    except Exception as e:
        print(f"处理文件 {file_name} 时出错: {str(e)}")

output_path = os.path.join(target_path, "G201统计结果.xlsx")
df_result.to_excel(output_path, index=False)
output_path_split_table = os.path.join(target_path, "G201统计结果分基因.xlsx")
df_results.to_excel(output_path_split_table, index=False)
print(f"G201统计结果已保存到 {output_path}")