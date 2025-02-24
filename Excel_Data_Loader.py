import pandas as pd
import numpy as np


class ExcelDataProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.df, self.df_geno, self.df_genoDp, self.df_noDP = self.batch_read_excel()

    def batch_read_excel(self):
        """
        批量读取指定路径下所有以 0x.xlsx 结尾的 Excel 文件
        :param target_path: 目标文件夹路径
        :return: 包含所有 DataFrame 的元组 (df, df_geno, df_genoDp, df_noDP)
        """
        try:
            df = pd.read_excel(self.file_path, sheet_name="Lane", header=3, index_col="index")
            df_geno = pd.read_excel(self.file_path, sheet_name="geno", header=0, index_col=0)
            df_genoDp = pd.read_excel(self.file_path, sheet_name="genoDp", header=0, index_col=0)
            df_noDP = pd.read_excel(self.file_path, sheet_name="noDP", header=0, index_col=0)

            return df, df_geno, df_genoDp, df_noDP
        except Exception as e:
            print(f"读取文件 {self.file_path} 时出错: {e}")
            return None, None, None, None

    def calculate_missing_rate(self, start_col):
        """
        计算指定起始列之后的列的缺失率
        :param start_col: 起始列索引
        :return: 包含每列缺失率的列表
        """
        if self.df_geno is None:
            print("未成功读取 df_geno 数据。")
            return []

        target_columns = self.df_geno.columns[start_col:]

        missing_rates = []
        for col in target_columns:
            valid_count = self.df_geno[col].count()  # 非缺失值数量
            total_count = len(self.df_geno)
            missing_rates.append(valid_count / total_count * 100)

        # 返回平均值或按需调整
        return missing_rates

    def process_data_by_gene(self, classify_by,reads_df):
        """
        计算测序深度、基因位点检测率和 STR% 并返回 DataFrame
        """

        #shared_columns = self.df_genoDp.columns[0:].tolist()
        shared_columns = [str(x) for x in self.df_genoDp.columns]
        df_results = pd.DataFrame()

        if classify_by in ["project", "tablet"]:
            unique_values = self.df[classify_by].unique()
            for val in unique_values:
                idx = reads_df[self.df[classify_by] == val].index
                print(f'indx是+{idx}')
                print(f'reads_df是{reads_df}')
                df_geno_slice = self.df_geno.loc[idx]
                df_genoDp_slice = self.df_genoDp.loc[idx]
                df_noDP_slice = self.df_noDP.loc[idx]
                reads_df_tmp = reads_df.loc[idx]

                df_Y_n_sum = df_genoDp_slice[shared_columns] + df_noDP_slice[shared_columns]
                #reads_df = reads_df[shared_columns]
                reads_df_nonzero = reads_df_tmp.replace(0, np.nan)
                str_df = df_Y_n_sum.div(reads_df_nonzero)

                geno_sample_columns = df_geno_slice.columns[4:4 + len(shared_columns)]
                missing_rate = (1 - df_geno_slice[geno_sample_columns].isnull().mean()).values
                if missing_rate.all(None) or missing_rate.all(0):
                    continue
                file_data = {
                    '基因组': shared_columns,
                    '测序深度': np.round(reads_df.mean(axis=0)),
                    '基因位点检测率%': np.round(missing_rate * 100),
                    'STR%': np.round(str_df.mean(axis=0) * 100, 2).astype(str) + '%'
                }
                df_results1 = pd.DataFrame(file_data)
                df_results1[classify_by] = val
                df_results = pd.concat([df_results, df_results1], ignore_index=True)
        else:
            df_Y_n_sum = self.df_genoDp[shared_columns] + self.df_noDP[shared_columns]
            #reads_df = self.df_genoDp[shared_columns]
            reads_df_nonzero = reads_df.replace(0, np.nan)
            str_df = df_Y_n_sum.div(reads_df_nonzero)

            geno_sample_columns = self.df_geno.columns[4:4 + len(shared_columns)]
            missing_rate = (1 - self.df_geno[geno_sample_columns].isnull().mean()).values

            file_data = {
                '基因组': shared_columns,
                '测序深度': np.round(reads_df.mean(axis=0)),
                '基因位点检测率%': np.round(missing_rate * 100),
                'STR%': np.round(str_df.mean(axis=0) * 100, 2).astype(str) + '%'

            }
            df_results = pd.DataFrame(file_data)

        return df_results