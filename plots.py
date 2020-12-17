
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import copy
from matplotlib import rc


species = [r"$\it{S. thermophilus}$ LMD-9", "$\it{L. acidophilus}$ La-14", "$\it{L. rhamnosus}$ GG", "$\it{L. helveticus}$ CNRZ32"]
codes = ["Sth", "La", "Lr", "Lh"]

def read_results(filepath):
    xls = pd.ExcelWriter(filepath)
    dataframes = []
    try:
        dataframes.append(pd.read_excel(xls, "iCC431"))
    except Exception as e:
        print(e)
    try:
        dataframes.append(pd.read_excel(xls, "iCC470"))
    except Exception as e:
        print(e)
    try:
        dataframes.append(pd.read_excel(xls, "iCC651"))
    except Exception as e:
        print(e)
    try:
        dataframes.append(pd.read_excel(xls, "iCC389"))
    except Exception as e:
        print(e)
    return dataframes


def dataframes_for_robustness(dataframes, compounds):
    i=0
    dataframes_res = []
    for df in dataframes:
        df = df.dropna(how='all')
        carbon_uptake = ['carbon_uptake']
        for column in df.columns:
            temp = column.split('_')[-1]
            if temp.isnumeric():
                carbon_uptake.append(temp)
                carbon_uptake.append(temp)
                carbon_uptake.append(temp)
        df.iloc[0,0] = 'Oxygen'
        df.index = df.iloc[:, 0].to_list()
        df = df.drop("Unnamed: 0", axis=1)
        oxygen = set(df.iloc[0,:].to_list())
        dataframes_dict = {}
        df = df.T
        for cp in compounds:
            if cp in df.columns:
                df[cp].loc[df["status"] == 'infeasible'] = 0
        for unique in oxygen:
            dataframes_dict[unique] = df.loc[df['Oxygen'] == unique]
            dataframes_dict[unique] = dataframes_dict[unique].drop("Oxygen", axis=1)
            dataframes_dict[unique].index = list(range(2, 42, 2))
            for column in dataframes_dict[unique].columns:
                dataframes_dict[unique] = dataframes_dict[unique].rename(columns = {column: column + "__" +unique+"__" + codes[i]})
        df = pd.DataFrame(index =list(range(2, 42, 2)))
        for data in dataframes_dict.keys():
            df = df.merge(dataframes_dict[data], left_index=True, right_index=True)
        df['carbon_uptake'] = list(range(2, 42, 2))
        dataframes_res.append(df)
        i+=1
    dataframes_compounds = {}
    for compound in compounds:
        dataframes_compounds[compound] = pd.DataFrame(index =list(range(2, 42, 2)))
        for data in dataframes_res:
            regex = compound + '.*'
            temp_df = data.filter(regex=regex)
            dataframes_compounds[compound] = dataframes_compounds[compound].merge(temp_df, left_index=True, right_index=True)
        dataframes_compounds[compound]['carbon_uptake'] =  dataframes_compounds[compound].index
    dataframes_compounds = list(dataframes_compounds.values())
    return dataframes_compounds

def dataframes_for_fva(dataframes):
    res_dataframes = []
    for dataframe in dataframes:
        dataframe = dataframe.dropna(how='all')
        columns_index = dataframe.shape[1]
        for i in range(1, columns_index):
            column = dataframe.columns[i]
            if 'Unnamed' not in column:
                old_col = column
                new_column = column
            else:
                new_column = old_col
            s =  dataframe.iloc[0,i]
            if type(s)!=str:
                s=''
            suffix = new_column + "__" + s
            dataframe = dataframe.rename(columns={column: suffix})
        dataframe= dataframe.rename(columns={'Unnamed: 0': "Metabolite"})
        dataframe.set_index('Metabolite',inplace=True)
        res_dataframes.append(dataframe.iloc[1:].T)
    return res_dataframes
def get_groupby(list_of_groups):
    pass

def apply_symetric(dataframes):
    for i in range(len(dataframes)):
        dataframes[i].fillna(0)
        for column in dataframes[i].columns:
            if dataframes[i][column].dtype != str:
                dataframes[i][column]= dataframes[i][column].apply(lambda x: -x if type(x) != str else x)
    return dataframes


def line_plot(directory, filename,dataframes, column_x, column_y, x_label, ylabel,axis = 0, column_y_axis_2=None, ylabel_axis2=None,grid_kws = {"hspace": 0.2, "wspace":0.3},figsize = (12,12), subplots_number=(2,2),
              title=None, titles=species,group_by='-',personalized_labels=None,color_map=None,linestyle_map=None, ylim=None):
    fig, axes = plt.subplots(subplots_number[0],subplots_number[1], figsize = figsize, gridspec_kw = grid_kws )
    fig.suptitle(title, fontsize = 14, y=0.95)
    secondary_axis = []
    k=0
    df = 0
    labels = []
    plots = []
    for i in range(subplots_number[0]):
        for j in range(subplots_number[1]):
            dataframe = dataframes[df]
            if axis==1:
                dataframe.index = dataframe.iloc[:,0].to_list()
                dataframe = dataframe.iloc[:,1:]
                dataframe = dataframe.T
            dataframe = rename_columns(dataframe)
            columny = get_existing_columns(dataframe,column_y )
            columnyaxis2 = get_existing_columns(dataframe,column_y_axis_2)
            dataframe[column_x] = dataframe[column_x].apply(np.abs)
            if group_by!='-': axes[i,j].set_prop_cycle( linestyle=group_by[1]) #color=group_by[0],
            for col in columny:
                if color_map:
                    plots += axes[i,j].plot(dataframe[column_x], dataframe[col], color = color_map[col.split('__')[-1]], linestyle = linestyle_map[col.split('__')[-2].split('_')[-1]])
                else:
                    plots += axes[i, j].plot(dataframe[column_x], dataframe[col])
            if ylim and codes[df] in ylim.keys():
                axes[i, j].set_ylim(ylim[codes[df]][0],ylim[codes[df]][1])
            [labels.append(x) for x in columny if x not in labels]
            axes[i,j].set_title(titles[df], fontstyle='italic')
            if i==subplots_number[0]-1:
                axes[i, j].set_xlabel(x_label, rotation=0, fontsize=10, color="k")
            if j == 0:
                axes[i, j].set_ylabel(ylabel, fontsize=10, color="k")
            if columnyaxis2:
                secondary_axis.append(axes[i,j].twinx())
                plots += secondary_axis[k].plot(dataframe[column_x], dataframe[columnyaxis2], color='r')
                [labels.append(x) for x in columnyaxis2 if x not in labels]
                if j == subplots_number[1]-1:
                    secondary_axis[k].set_ylabel(ylabel_axis2, fontsize=10, color="k")
                k += 1
            df+=1
    if personalized_labels:
        if type(personalized_labels)==dict:
            labels = get_labels_by_order(labels, personalized_labels)
        else:
            labels = personalized_labels
    fig.legend(plots, labels, loc='lower center', ncol=2)
    fig.show()
    fig.savefig(directory + filename + '.png')



def fva_plot(directory, filename,dataframes, column_x = None, column_y= None, x_label= None, y_label= None,axis = 0, column_y_axis_2=None, ylabel_axis2=None,figsize = (12,12), subplots_number=(2,2),title = None,
             titles=species,group_by='-',personalized_labels=None,columns=True, ylim=None):
    grid_kws = {"hspace": 0.25, "wspace": 0.3}  # , gridspec_kw=grid_kws,figsize = (7,18)
    fig, axes = plt.subplots(subplots_number[0], subplots_number[1], figsize=figsize, gridspec_kw=grid_kws)
    fig.suptitle(title, fontsize = 14, y=0.95)
    df=0
    plots = []
    labels = []
    for i in range(subplots_number[0]):
        for j in range(subplots_number[1]):
            dataframe = dataframes[df]
            if column_y[0] in dataframe.columns:
                columnx = column_x[df]
                xlabel = columnx.replace("PKETX","PKT") + x_label
                dataframe[columnx] = dataframe[columnx].apply(np.abs)
                flux, minimum, maximum = dataframe.filter(regex = 'flux', axis=0),dataframe.filter(regex = 'minimum', axis=0),dataframe.filter(regex = 'maximum', axis=0)
                flux[columnx], minimum[columnx], maximum[columnx] = np.arange(0,110,10),np.arange(0,110,10),np.arange(0,110,10)
                if columns:
                    plots += axes[i, j].bar(maximum[columnx], [x[0] for x in maximum[column_y].values.tolist()], width=2.5, label='Maximum')
                    plots += axes[i, j].bar(minimum[columnx], [x[0] for x in minimum[column_y].values.tolist()],
                                            width=2.5, color='orange', label='Minimum')
                if ylim:
                    plots+= axes[i, j].set_ylim(flux[column_y].min().min() - ylim[0], flux[column_y].max().max() + ylim[1])
                axes[i,j].plot(flux[columnx], flux[column_y], label = 'Flux')
                axes[i,j].set_title(titles[df])
                axes[i, j].set_xlabel(xlabel)
                axes[i, j].set_ylabel(y_label)
            df+=1
    if columns:
        h, l = axes[i,j].get_legend_handles_labels()
        fig.legend(h, l , loc='lower center', ncol=3)
    fig.show()
    fig.savefig(directory + filename + '.png')

def phenotypic_phase_plane(directory, filename,dataframes, column_x = None, column_y= None,column_z=None, x_label= None, y_label= None,z_label= None,figsize = (12,12), subplots_number=(2,2),title = None,
             titles=species,group_by='-',personalized_labels=None):
    grid_kws = {"hspace": 0.5, "wspace": 0.5}
    fig, axes= plt.subplots(figsize=figsize, gridspec_kw=grid_kws) #
    fig.suptitle(title, fontsize=14, y=0.95)
    df=0
    plt.axis('off')
    for i in range(subplots_number[0]):
        for j in range(subplots_number[1]):
            dataframe = dataframes[df]
            dataframe.index = dataframe['Unnamed: 0']
            dataframe = dataframe.drop(['Unnamed: 0'], axis=1)
            for column in dataframe.columns:
                new_column = str(-float(column))
                dataframe = dataframe.rename(columns={column:new_column})
            dataframe = dataframe.unstack().reset_index()
            dataframe.columns = [column_x,column_y,column_z]
            ax = fig.add_subplot(2,2,df+1, projection='3d')
            surf  = ax.plot_trisurf(dataframe[column_y], dataframe[column_x],   dataframe[column_z], cmap=plt.cm.coolwarm, linewidth=0.2, antialiased=False)
            ax.set_xlabel(x_label, fontsize = 9)
            ax.set_ylabel(y_label, fontsize = 9)
            ax.set_zlabel(z_label, fontsize = 9)
            ax.set_title(titles[df])
            # fig.colorbar(surf, shrink=0.5, aspect=15)
            df += 1
    fig.show()
    fig.savefig(directory + filename + '.png')



def pie(dataframes, column_x, column_y,directory, filename, title):
    fig, axes = plt.subplots(2,2,figsize=(12,12), subplot_kw=dict(aspect="equal"))
    fig.suptitle(title, fontsize=14, y=0.95)
    axe = axes.ravel()
    i=0
    for dataframe in dataframes:
        wedges, texts = axe[i].pie(dataframe[column_y].head(16), wedgeprops=dict(width=0.3), startangle=-40)
        axe[i].set_xlabel(species[i], rotation=0, fontsize=10, color="k")

        i+=1
    if column_x != 'index':
        fig.legend(dataframes[0][column_x].head(16).to_list(), fontsize=8,loc='lower center' )
    else:
        fig.legend(dataframes[0].head(16).index.to_list(), fontsize=8,loc='lower center' )
    # plt.show()
    fig.savefig(directory + filename + '.png')

def connectivity_pie(directory,directory_to_save, column_x, column_y, filename, title):
    xls = pd.ExcelWriter(directory + 'results_connectivity.xlsx')
    df1 = pd.read_excel(xls, "iCC431")
    df2 = pd.read_excel(xls, "iCC470")
    df3 = pd.read_excel(xls, "iCC651")
    df4 = pd.read_excel(xls, "iCC389")
    pie(dataframes=[df1, df2, df3, df4],
        column_x = column_x,
        column_y = column_y,
        directory = directory_to_save,
        filename = filename,
        title= title)



def get_labels_by_order(old_labels, labels):
    res = []
    for lab in old_labels:
        res.append(labels[lab])
    return res


def get_label(column):
    mapper = {'Sth': 'S.thermophilus', 'La': 'L. acidophilus', 'Lh': 'L. helveticus', 'Lr': 'L. rhamnosus'}
    temp = column.split('__')
    o2=temp[-2].split("_")[-1]
    sp=temp[-1]
    res = 'O2 = ' + o2 +", "+mapper[sp]
    return res



def rename_columns(dataframe):
    columns = ['EX_lcts_e',"EX_glc__aD_e","Lactose","Glucose", 'Alpha D-glucose', 'lcts_e',"glc__aD_e" ]
    for col in columns:
        if col in dataframe.columns:
            dataframe = dataframe.rename(columns={col: "carbon_uptake"})
    return dataframe


def get_existing_columns(dataframe, column_y):
    if type(column_y) != list: return
    column = copy.deepcopy(column_y)
    for element in column_y:
        if element not in dataframe.columns:
            column.remove(element)
    return column

def get_column_name(list_of_columns):
    res = []
    for compound in list_of_columns[0]:
        for oxygen in list_of_columns[1]:
            for species in list_of_columns[2]:
                res.append(compound +  "__EX_o2_e_" + oxygen + '__' + species)
    return res




def robustness_plot(directory):
    xls = pd.ExcelWriter(directory + 'results_robustness_analysis.xlsx')
    df1 = pd.read_excel(xls, "iCC431")
    df2 = pd.read_excel(xls, "iCC464")
    df3 = pd.read_excel(xls, "iCC390")
    df4 = pd.read_excel(xls, "iCC644")
    line_plot(dataframes=[df1,df2,df3,df4],
              column_x=None,
              column_y='flux_maximum',
              ylabel='Growth Rate ($h^{-1}$)',
              directory=directory,
              filename="Robustness")

