import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def change_label(dataset: pd.DataFrame, label_encoding: dict):
    df = dataset.copy()
    df.Label = df.Label.map(label_encoding)
    return df

def get_dict_group_ID(df_intensity, df_id):
    key_groups = list(df_id.ID_group.unique())

    # RT_metabolites : ID Name -> dict_id

    _features = df_id.Features
    _id = df_id.ID

    dict_id = {_features[i] : _id[i] for i in range(len(_features))}

    # ID_group : features -> dict_id_group

    dict_group = {}
    for g in key_groups:
        dict_group[g] = list(df_id[df_id.ID_group == g].Features)
    return dict_group, dict_id, key_groups


def get_mean_dataset(dataset: pd.DataFrame, label_list: list, is_save_file=False):

    '''
    dataset     : pd.DataFrame, column orders = [Name, Label, feature1, feature2, feature3...], features should be int or float type
    label_list  : Group label list eg, Control_group, Treatment_group...
    '''

    results = {}
    for label in label_list:
        df_rev = dataset[dataset.Label == label]
        m = df_rev.mean(numeric_only=True).values
        results[label] = m

    df_mean = pd.DataFrame(results, index=[f for f in dataset.columns if f not in ['Name', 'Label']])

    if is_save_file:
        df_mean.to_csv('mean.csv')
        print('mean.csv files are saved')

    return df_mean


# Figure

def get_group_figure(
        dataset: pd.DataFrame,
        dict_group: dict,
        dict_id: dict,
        group: str,
        order: list,
        kind='barplot', palette='Pastel1', figsize=(20, 18), top_alpha=1.0, is_save_file=False, dpi=600, x_label_rotation='45'
    ):

    '''
    dataset (pd.DataFrame)  : dataset containing metabolite intensity. The dataset should include Label column
    dict_group (dict)       : ID_groups:features # Key metabolties... (interesting features)
    dict_id (dict)          : features:ID # Get metabolties ID
    top_alpha               : hyperparameters of y-axis max value
    '''

    plt.figure(figsize=figsize)

    for i, feat in enumerate(dict_group[group]):

        plt.subplot(5, 6, i + 1)


        if kind == 'barplot':
            sns.barplot(data=dataset, x='Label', y=feat, order=order, palette=palette, errwidth=1.5, errcolor='k', edgecolor='blacK')           

        elif kind == 'stripplot':    
            sns.stripplot(data=dataset, x='Label', y=feat, order=order, jitter=True, palette=palette)
            

            sns.boxplot(data=dataset, x='Label', y=feat, order=order,
                        showmeans=True,
                        showbox=False,
                        medianprops={'visible' : False},
                        whiskerprops={'visible' : False},
                        showcaps=False,
                        meanline=True,
                        showfliers=False,
                        meanprops={'color' : 'k', 'lw': 1})
        
        elif kind == 'boxplot':
            sns.boxplot(data=dataset, x='Label', y=feat, order=order, palette=palette,)

        else:
            raise ValueError('Invalid kinds of plot. [barplot, stripplot, boxplot]')           


        plt.xlabel('')
        plt.ylabel('')
        plt.xticks(fontsize=12)
        plt.title(f'{dict_id[feat]}\n{feat}', fontdict={'fontsize' : 15, 'fontweight' : 'bold'})
        plt.xticks(rotation=x_label_rotation, ha='right')
        bottom, top = plt.ylim()
        plt.ylim((bottom, top*top_alpha))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.7)


    if is_save_file:
        plt.savefig('figure.png', dpi=dpi)
        print('figure.png file is saved')

    plt.show()



def mapping_new_features_from_reference(
        reference_dataset: pd.DataFrame,
        matching_dataset: pd.DataFrame,
        windows_ppm=20, is_save_file=False
    ):

    '''
    reference_dataset should include the Features(name), ID, ID_group, mz columns
    matching_dataset should include the Features(name) column
    
    '''
    windows_1ppm = reference_dataset.mz / 1000000
    reference_dataset['mz-ppm'] = reference_dataset.mz - (windows_1ppm * windows_ppm)
    reference_dataset['mz+ppm'] = reference_dataset.mz + (windows_1ppm * windows_ppm)


    mz_list = list(matching_dataset.mz)
    features_list = list(matching_dataset.Features) # features = names, eg... RT_m/z
    

    results = {}
    for mz, feat in zip(mz_list, features_list):
        dataset = reference_dataset[(reference_dataset['mz-ppm'] <= mz) & (reference_dataset['mz+ppm'] >= mz)]

        if not dataset.empty :
            results[mz] = [feat] + list(dataset.Features) + list(dataset.ID) + list(dataset.ID_group)


    dataset = pd.DataFrame(results).T.rename(columns={0: 'mapped_features',1:'reference_features', 2:'ID', 3:'ID_group'})
    if is_save_file:
        dataset.to_csv('mapping.csv')
        print('mapping.csv files are saved')

    return dataset

