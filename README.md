#1. 概述
本程序用于从指定的 Excel 文件中读取多个工作表的数据，依据用户输入的分类依据（"project" 或 "tablet"）和过滤阈值，对数据进行统计计算，并生成两份统计结果文件。

主要功能：
批量读取 Excel 文件中的 "Lane"、"geno"、"genoDp"、"noDP" 工作表数据。
根据用户输入的分类依据对数据进行分组（或整体处理）。
对分组后的数据进行过滤（基于 STR_AVG 阈值）及统计计算（如各类位点的平均值、标准差、总reads、STR指标等）。
计算基因位点检测率和测序深度等指标。
将计算结果分别保存到两个 Excel 文件中：
G201统计结果.xlsx：整体统计结果。
G201统计结果分基因.xlsx：基于基因数据的详细统计结果。
2. 环境依赖
Python 版本：Python 3.x
第三方库：
pandas
numpy
标准库：
os
warnings
安装依赖（如未安装时）：

bash
复制
编辑
pip install pandas numpy
3. 程序结构
程序主要包含以下部分：

3.1 导入模块与配置
导入 os、pandas、numpy、warnings 模块。
禁用所有警告输出，以便控制台输出更简洁。
3.2 数据读取
函数：batch_read_excel(file_path)

功能：
批量读取指定路径下 Excel 文件中四个工作表的数据，分别为 "Lane"、"geno"、"genoDp" 和 "noDP"。

参数：

file_path：Excel 文件的完整路径。
返回值：
返回一个元组 (df, df_geno, df_genoDp, df_noDP)，对应不同工作表的数据，其中：

df：来自 "Lane" 工作表，读取时跳过前 3 行，并以 "index" 列作为索引。
df_geno、df_genoDp、df_noDP：分别对应其他三个工作表，均以第一列作为索引。
异常处理：
如果文件读取过程中出现异常，将捕获并打印错误信息。

3.3 数据统计计算
函数：calculate_statistics(df, project_tablet_name)

功能：
根据用户输入的分类依据，对传入的 DataFrame（通常为 "Lane" 工作表数据的一个切片）进行数据过滤和统计计算。

先根据全局变量 classify_by 判断是否按 "project" 或 "tablet" 进行分组筛选。
对数据进行过滤（保留 STR_AVG 大于或等于用户输入阈值的数据）。
计算各类统计指标，如：
均值计算：常位点（A_Typed 或 auto_loci_typed）、X位点、Y位点。
核心基因：常核心基因和 Y 核心基因的平均值。
reads 统计：总reads（以百万为单位）、单个样本的平均reads以及有效reads比（转换为百分比）。
区间统计：对指定基因区域（例如从 "Amelogenin" 到 "D21S1270"）计算平均值（AVG）和相对标准差（STD）。
同时，还计算 STR 指标（如 STR_AVG 或 STR均值、STR_STD 或 STR标准化STD）。
参数：

df：待处理的 DataFrame（从 "Lane" 工作表切片得到）。
project_tablet_name：当前分组的名称（例如某个 project 的具体名称）。
返回值：
返回一个字典，包含所有统计指标及分组标签，供后续拼接到结果 DataFrame 中。

3.4 缺失率计算
函数：calculate_missing_rate(df_geno, start_col)

功能：
计算传入的 df_geno DataFrame 中，从指定起始列开始的各列缺失值比例，并转换为百分比。

参数：

df_geno：基因检测数据的 DataFrame（来自 "geno" 工作表）。
start_col：计算缺失率的起始列索引。
返回值：
返回一个列表，列表中的每个值对应目标列的缺失率（百分比）。

3.5 主程序逻辑
用户输入：

分类依据（classify_by）：提示用户输入 "project" 或 "tablet"，留空则不进行分组统计。
STR_AVG 过滤阈值（cut_num）：用户输入过滤低质量数据的阈值。
Excel 文件路径（target_path）：提示用户输入 Excel 文件的完整路径。
工作目录设置：
根据文件路径切换到文件所在目录，确保后续文件读写在正确目录下进行。

数据读取与统计：

调用 batch_read_excel 读取 Excel 文件中的数据，并将结果分别赋值给 df0、df_geno、df_genoDp、df_noDP。
根据 classify_by 的取值，对 df0 进行分组（或整体处理），调用 calculate_statistics 对每个分组进行统计计算，并将结果拼接到 df_result DataFrame 中。
另外，通过对 df_genoDp 和 df_noDP 进行数据合并、除零处理（将 0 替换为 NaN）、及计算 STR 指标，进一步生成包含基因组、测序深度、基因位点检测率、STR% 等指标的 DataFrame df_results。
结果输出：
将最终生成的统计结果 DataFrame 分别保存为 Excel 文件：

G201统计结果.xlsx：包含总体统计数据。
G201统计结果分基因.xlsx：包含基于基因的详细统计数据。
4. 使用说明
运行程序：
在命令行中运行该 Python 脚本，例如：

bash
复制
编辑
python your_script.py
输入提示：

分类统计依据：输入 “project” 或 “tablet”（不区分大小写），留空表示不按分组统计。
过滤阈值：输入一个整数，用于过滤 STR_AVG 值低于该阈值的数据。
Excel 文件路径：输入包含数据的 Excel 文件的完整路径，确保文件中包含以下工作表：
"Lane"（含有统计数据，需包含 'STR_AVG'、'总reads' 等列）
"geno"、"genoDp"、"noDP"（用于后续的基因检测率、测序深度等计算）。
输出文件：
程序执行完毕后，将在 Excel 文件所在目录下生成两个文件：

G201统计结果.xlsx
G201统计结果分基因.xlsx
5. 注意事项
数据格式：
请确保 Excel 文件中各工作表的列名和数据格式与程序要求一致，否则可能导致读取或计算错误。

路径问题：
使用程序前请确认输入的 Excel 文件路径正确，并且当前用户具有该目录的读写权限。

警告信息：
程序默认禁用所有警告，若需调试或检查潜在问题，可取消警告禁用。

扩展和维护：
根据实际需求可以进一步扩展统计指标或调整数据处理逻辑，建议在修改前备份原始代码。

6. 版权与作者信息
作者：赵通介
用途：本程序主要用于内部数据统计与分析，未经许可不得用于商业分发。
更新日期：2025/2/20
