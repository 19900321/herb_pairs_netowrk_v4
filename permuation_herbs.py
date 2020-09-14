import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.gridspec as gridspec
# some outsides functions need

# the permutation related ones
def permu_test(pooled, sizeZ, sizeY):
    np.random.shuffle(pooled)
    starZ = pooled[:sizeZ]
    starY = pooled[-sizeY:]
    return starZ.mean() - starY.mean()

def permu_p_value( z, y, numSamples):
    pooled = np.hstack([z, y])
    delta = np.mean(z) - np.mean(y)
    estimates = np.array([permu_test(pooled, len(z), len(y)) for i in range(numSamples)])
    diffCount = len(np.where(estimates <= delta)[0])
    hat_asl_perm = 1.0 - (float(diffCount) / float(numSamples))
    return hat_asl_perm

def cal_permu_p_va_repeat(top_col, random_col, numSamples, repeat):

    return [permu_p_value(top_col, random_list, numSamples)
            for i, random_list in enumerate(random_col) if i < repeat]

# the p value related ones
def cal_p_va_repeat(top_col, random_col, repeat):
    return [stats.ttest_ind(top_col, random_list, equal_var=False, nan_policy='omit')[1]
            for i, random_list in enumerate(random_col) if i < repeat]

# the auc related ones
from sklearn.metrics import roc_curve, auc
from sklearn import metrics
from sklearn.metrics import confusion_matrix, accuracy_score, roc_auc_score, f1_score, \
    roc_curve, auc, precision_recall_curve

def roc_cal(top_v, random_v):
    label = [1] * len(top_v) + [2] * len(random_v)  # creat lable by top as 1, random as 2
    values_y = top_v + random_v
    fpr, tpr, _ = metrics.roc_curve(label, values_y, pos_label=2)
    roc_auc = auc(fpr, tpr)
    return {'fpr':fpr, 'tpr':tpr, 'roc':roc_auc}

from collections import defaultdict
def cal_roc_repeat(top_col, random_col, repeat):
    fpr_list = defaultdict(list)
    all_result = [roc_cal(top_col, random_list)
            for i, random_list in enumerate(random_col) if i < repeat]

    return [i['roc'] for i in all_result]


def prc_cal(top_v, random_v):
    label = [1] * len(top_v) + [2] * len(random_v)  # creat lable by top as 1, random as 2
    values_y = top_v + random_v
    precision, recall, _ = precision_recall_curve(label, values_y, pos_label=2)
    prc_auc = auc(recall, precision)
    return {'precision':precision, 'recall':recall, 'roc':prc_auc}

def cal_prc_repeat(top_col, random_col, repeat):
    fpr_list = defaultdict(list)
    all_result = [prc_cal(top_col, random_list)
            for i, random_list in enumerate(random_col) if i < repeat]

    return [i['roc'] for i in all_result]

# this is used for changing data type
def expand_list(df, list_column, new_column):
    lens_of_lists = df[list_column].apply(len)
    origin_rows = range(df.shape[0])
    destination_rows = np.repeat(origin_rows, lens_of_lists)
    non_list_cols = (
      [idx for idx, col in enumerate(df.columns)
       if col != list_column]
    )
    expanded_df = df.iloc[destination_rows, non_list_cols].copy()
    expanded_df[new_column] = (
      [item for items in df[list_column] for item in items]
      )
    expanded_df.reset_index(inplace=True, drop=True)
    return expanded_df


class Result_distance:

    def __init__(self, filename_top, filename_random):
        self.filename_top = filename_top
        self.filename_random = filename_random
        self.data = self.get_combined_data()
        self.pd_melt = self.get_melt()
        self.grouped_data = self.group_list()

    def get_combined_data(self):
        data_top = pd.read_csv(self.filename_top)
        data_random = pd.read_csv(self.filename_random)
        data = pd.concat([data_top,  data_random], keys=['top', 'random'], sort=False).reset_index()
        data = data.rename(index=str, columns={'level_0': 'class'}).drop(['level_1'], axis=1)
        data = data.sort_values(by=['frequency'], ascending=False)
        return data

    def get_melt(self):
        value_vars_use = ['separation', 'closest', 'shortest', 'kernel', 'center']
        id_vars_use = ['herb1', 'herb1_name', 'herb2', 'herb2_name', 'frequency', 'class', 'Ingredient-level distance type']

        pd_melt = pd.melt(self.data, id_vars=id_vars_use,
                          value_vars=value_vars_use,
                          value_name='Distance')
        pd_melt = pd_melt.rename(index=str, columns={'variable': 'Herb-level distance type'})
        pd_melt['herb_ingre_method'] = pd_melt['Herb-level distance type'] + '_' + pd_melt['Ingredient-level distance type']

        return pd_melt


    def get_mean(self):
        pd_mean = pd.pivot_table(self.pd_melt, values='Distance', index=['Herb-level distance type', 'Ingredient-level distance type'], columns = ['class'], aggfunc=np.mean).reset_index()
        pd_mean = pd_mean.rename(columns = {'top':'Distance for top herb pairs', 'random':'Distance for random herb pairs'})
        self.pd_mean = pd_mean

    #  calculated correlation value, frequency_number is we only choose the herb pairs appear larger than frequency_number time
    def get_cor(self):
        pd_cor_pre = self.pd_melt[['herb1', 'herb2', 'frequency', 'herb_ingre_method', 'Distance']]
        pd_cor_pre = pd_cor_pre.pivot_table('Distance',
                                            ['herb1', 'herb2','frequency'],
                                            'herb_ingre_method')
        pd_cor_pre = pd_cor_pre.reset_index('frequency')
        self.correlation = pd_cor_pre.corr(method='spearman')


    # group data by herb method and ingredients method, as well as top or random. Grouped data used for further analysis
    def group_list(self):
        grouped = self.pd_melt.groupby(['Herb-level distance type', 'Ingredient-level distance type', 'class'])['Distance'].apply(list)
        return grouped

    # repeat the top and random comparison, top number is how many top pairs you want to use
    def distance_split(self):
        pf_new = self.grouped_data.unstack()
        top_number = len(list(pf_new['top'])[1])
        random_values = pf_new['random']
        random_list = [[i[x:x + top_number] for x in range(0, len(i), top_number)] for i in random_values]
        pf_new['random'] = random_list
        self.new_grouped_list = pf_new

    # get normal p values, repeat is how many time we random select the random distance to get p value
    def get_p_value(self):
        repeat = len(list(self.new_grouped_list['random'])[1]) - 1
        print('repeat {} time'.format(repeat))
        self.new_grouped_list['p-value_one'] = self.new_grouped_list.apply(
            lambda x: cal_p_va_repeat(x.top, x.random, repeat), axis=1)
        self.new_grouped_list['p-value'] = self.new_grouped_list['p-value_one'].apply(np.mean)
        self.pd_mean['p-value'] = list(self.new_grouped_list['p-value'])


    # get permutation p values, repeat is how many time we random select the random distance,
    # numSamples is the permutation times for each permutation
    def get_permu_p_value(self, numSamples):
        repeat = len(list(self.new_grouped_list['random'])[1]) - 1
        self.new_grouped_list['permu_p-value_one'] = self.new_grouped_list.apply(
            lambda x: cal_permu_p_va_repeat(x.top, x.random, numSamples, repeat), axis=1)
        self.new_grouped_list['permu_p-value'] = self.new_grouped_list['permu_p_values'].apply(np.mean)


    def pre_plot_data_p_value(self):
        pd_new = self.new_grouped_list.reset_index()
        pd_new = pd.melt(pd_new, id_vars=['Herb-level distance type', 'Ingredient-level distance type'],
                          value_vars=['p-value_one', 'permu_p-value_one'],
                          value_name='p-value',
                         var_name='test_method')
        self.plot_data_p_value = expand_list(pd_new, 'p_value', 'p_values')


    def plot_p_value(self, how, file_folder):
        plt.figure(figsize=(27, 15))
        sns.set_context("paper", font_scale=2)
        sns.catplot(x='Herb-level distance type',
                        y='p-value',
                        hue='Ingredient-level distance type',
                        row='test_method',
                        data=self.plot_data_p_value,
                        kind="box",
                        height=4,
                        aspect=4,
                        legend=False)
        plt.legend(loc='upper left', ncol=3)

        if how == 'save':
            plt.savefig('result/{}/p_values.png'.format(file_folder))
        elif how == 'show':
            plt.show()

    def get_auc(self):
        repeat = len(list(self.new_grouped_list['random'])[1]) - 1
        self.new_grouped_list['roc'] = self.new_grouped_list.apply(
            lambda x: cal_roc_repeat(x.top, x.random, repeat), axis=1)
        self.new_grouped_list['AUROC'] = self.new_grouped_list['roc'].apply(np.mean)
        self.pd_mean['AUROC'] = list(self.new_grouped_list['AUROC'])

    def get_prc(self):
        repeat = len(list(self.new_grouped_list['random'])[1]) - 1
        self.new_grouped_list['prc'] = self.new_grouped_list.apply(
            lambda x: cal_prc_repeat(x.top, x.random, repeat), axis=1)
        self.new_grouped_list['AUPRC'] = self.new_grouped_list['prc'].apply(np.mean)
        self.pd_mean['AUPRC'] = list(self.new_grouped_list['AUPRC'])


    # plot the density of random and tops, top_number is how many of top herb pairs to use
    def plot_density(self, how, file_folder):

        import seaborn as sns
        import matplotlib
        import matplotlib.pyplot as plt

        sns.set_context('paper', font_scale=1)
        fig = plt.figure(figsize=(15, 10))
        n = 1
        for herb_m in self.pd_melt['Herb-level distance type'].unique():
            for ingre_m in self.pd_melt['Ingredient-level distance type'].unique():
                ax = fig.add_subplot(5, 5, n)
                n += 1
                #ax.set_ylabel('density', fontsize=25)
                ax.set_alpha(.6)
                ax.set_title('herb:{} ingre:{}'.format(herb_m,ingre_m))
                for source in self.pd_melt['class'].unique():
                    data = self.grouped_data[herb_m][ingre_m][source]
                    sns.distplot(data, hist=False, kde=True, norm_hist=True,
                                 kde_kws={'shade': True, 'linewidth': 8},
                                 label=source)
                if n == 25:
                    ax.legend(loc='upper left', fontsize='large')
        plt.tight_layout()
        if how == 'save':
            fig.savefig('result/{}/density.png'.format(file_folder))
        elif how == 'show':
            plt.show()


    # one way to plot correlation based on the correlation value we calculated
    def plot_correlation_1(self, how, file_folder):
        import seaborn as sns
        # Generate a mask for the upper triangle
        mask = np.triu(np.ones_like(self.correlation, dtype=np.bool))
        sns.set_context("paper", font_scale=1)
        f, ax = plt.subplots(figsize=(11, 9))

        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        sns.heatmap(self.correlation, mask=mask, cmap=cmap, vmax=.3, center=0,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5})

        if how == 'save':
            plt.savefig('result/{}/correlation.png'.format(file_folder))
        elif how == 'show':
            plt.show()


    # another way to plot correlation based on the correlation value we calculated
    def plot_correlation_2(self, how, file_folder):
        import seaborn as sns
        f, ax = plt.subplots(figsize=(11, 9))
        sns.set_context('paper', font_scale=1)
        ax = sns.heatmap(
            self.correlation,
            vmin=-1, vmax=1, center=0,
            cmap=sns.diverging_palette(20, 220, n=200),
            square=True
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            horizontalalignment='right'
        )
        if how == 'save':
            plt.savefig('result/{}/correlation.png'.format(file_folder))
        elif how == 'show':
            plt.show()


    def plot_auc(self, how, file_folder):
        import matplotlib
        # font = {'size': 20}
        # matplotlib.rc('font', **font)
        repeat = len(list(self.new_grouped_list['random'])[1]) - 1
        fig = plt.figure(figsize=(10, 10))
        fig = plt.figure(figsize=(15, 15))
        sns.set_context('paper', font_scale=1)
        n = 1
        for herb_m in self.pd_melt['Herb-level distance type'].unique():
            for ingre_m in self.pd_melt['Ingredient-level distance type'].unique():
                ax = fig.add_subplot(5, 5, n)
                n += 1
                # ax.set_ylabel('auc'); ax.set_alpha(.6)
                ax.set_title('herb:{} ingre:{}'.format(herb_m, ingre_m))
                top = self.new_grouped_list.loc[(herb_m, ingre_m), 'top']
                randoms = self.new_grouped_list.loc[(herb_m, ingre_m), 'random']
                auc_list = []; fpr_list=[];tpr_list=[]
                for i, random_list in enumerate(randoms):
                    if i < repeat:
                        out = roc_cal(top, random_list)
                        fpr, tpr, roc_auc = out['fpr'], out['tpr'], out['roc']
                        fpr_list.append(fpr); tpr_list.append(tpr); auc_list.append(roc_auc)
                        ax.plot(fpr, tpr, lw=2, alpha=0.3)
                ax.plot(np.mean(pd.DataFrame(fpr_list)), np.mean(pd.DataFrame(tpr_list)), '',
                        lw=2, alpha=0.8, color='blue',
                        label='AUC = {:.3f}'.format(np.mean(auc_list)))
                ax.plot([0, 1], [0, 1], '--')
                ax.set_xlim(-0.1, 1.1); ax.set_ylim(-0.1, 1.1)
                ax.set_ylabel('True Positive Rate'); ax.set_xlabel('False Positive Rate')
                ax.legend(loc='lower right', fontsize='large')
            plt.subplots_adjust(hspace=1.2, wspace=0.50); plt.tight_layout()
        if how == 'save':
            plt.savefig('result/{}/roc_auc.png'.format(file_folder))
        elif how == 'show':
            plt.show()

    def plot_prc(self, how, file_folder):
        import matplotlib
        repeat = len(list(self.new_grouped_list['random'])[1]) - 1
        fig = plt.figure(figsize=(15, 15))
        sns.set_context('paper', font_scale=1)
        n = 1
        for herb_m in self.pd_melt['Herb-level distance type'].unique():
            for ingre_m in self.pd_melt['Ingredient-level distance type'].unique():
                ax = fig.add_subplot(5, 5, n)
                n += 1
                ax.set_ylabel('prc'); ax.set_alpha(.6)
                ax.set_title('herb:{} ingre:{}'.format(herb_m, ingre_m))
                top = self.new_grouped_list.loc[(herb_m, ingre_m), 'top']
                randoms = self.new_grouped_list.loc[(herb_m, ingre_m), 'random']
                auc_list = []
                precision_list = []
                recall_list = []
                for i, random_list in enumerate(randoms):
                    if i < repeat:
                        out = prc_cal(top, random_list)
                        precision, recall, prc_auc = out['precision'], out['recall'], out['roc']
                        precision_list.append(precision)
                        recall_list.append(recall)
                        auc_list.append(prc_auc)
                        ax.plot(recall, precision, lw=2, alpha=0.3)

                ax.plot(np.mean(pd.DataFrame(recall_list)), np.mean(pd.DataFrame(precision_list)), '',
                        lw=2, alpha=0.8, color='blue',
                        label='AUC = {:.3f}'.format(np.mean(auc_list)))
                #ax.plot([0, 1], [0, 1], '--')
                ax.set_xlim(-0.1, 1.1)
                ax.set_ylim(-0.1, 1.1)
                ax.set_ylabel('Precision')
                ax.set_xlabel('Recall')
                ax.legend(loc='lower right', fontsize='large')
            plt.subplots_adjust(hspace=1.2, wspace=0.50); plt.tight_layout()
        if how == 'save':
            plt.savefig('result/{}/prc_auc.png'.format(file_folder))
        elif how == 'show':
            plt.show()

    def prepare_raw_result(self):
        new_grouped_list = self.new_grouped_list[['roc', 'prc', 'p-value_one']].reset_index()
        zip_r = list(zip(new_grouped_list['roc'], new_grouped_list['prc'], new_grouped_list['p-value_one']))
        extend_r_list = []
        for i in list(range(len(zip_r))):
            method = list(
                new_grouped_list.iloc[i, :][['Herb-level distance type', 'Ingredient-level distance type']])
            n = len(zip_r[i][0])
            methods_1 = [method[0]] * n
            methods_2 = [method[1]] * n
            extend_r = list(zip(methods_1, methods_2, *zip_r[i]))
            extend_r_list += extend_r

        extend_raw_result = pd.DataFrame(extend_r_list,
                                         columns=['Herb-level distance type',
                                                  'Ingredient-level distance type',
                                                  'AUROC',
                                                  'AUPRC',
                                                  'p-value_one'])
        return extend_raw_result

    def save_analysis(self, file_folder, forbid_include=True):
        self.get_mean()
        import scipy
        #mean_p = scipy.stats.ttest_rel(self.pd_mean['top'], self.pd_mean[
            #'random'])
        #print('the p value about top and random mean value is {}'.format(mean_p))
        if forbid_include == False:
             self.get_cor()
             self.plot_correlation_2('save', file_folder)
        self.group_list()
        self.plot_density('save', file_folder)
        self.distance_split()
        self.get_p_value()
        self.get_auc()
        self.get_prc()
        self.pd_mean.reset_index().to_csv('result/{}/Table 1.csv'.format(file_folder))

        new_grouped_list = self.prepare_raw_result()
        new_grouped_list.to_csv('result/{}/all_records.csv'.format(file_folder))
        #self.get_permu_p_value(10000, repeat_time)  # change smaller one as example
        #self.pre_plot_data_p_value()
        #self.plot_p_value('save', file_folder)
        self.plot_auc('save', file_folder)
        self.plot_prc('save', file_folder)

def main():
    filename_top = 'result/top_200.csv'
    filename_random = 'result/result_random_0_10000.csv'
    result = Result_distance(filename_top, filename_random)
    file_folder = 'top_random_new'
    result.get_mean()
    result.plot_density('save', file_folder)
    result.distance_split()
    result.get_auc()
    result.get_prc()
    result.plot_auc('save', file_folder)
    result.plot_prc('save', file_folder)
    result.save_analysis(file_folder,forbid_include=False )




