from plots import *





if __name__ == "__main__":
    _results_directory = r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Python scripts\LAB__git\lab_models\results/"
    _plots_directory = r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Python scripts\LAB__git\lab_models\plots_plt/"

#### Connectivity
    connectivity_pie(directory=_results_directory,
                     directory_to_save = _plots_directory,
                     column_x='Unnamed: 0',
                     column_y="sum",
                     filename="Connectivity",
                     title='Connectivity Analysis'
                     )

    ##### Robustness
    # dataframes = read_results(_results_directory + "results_production_robustness_analysis" + ".xlsx")
    # compounds =  ['EX_lac__L_e', 'EX_ac_e', 'EX_etoh_e', 'EX_for_e','EX_fum_e', 'EX_succ_e','EX_acald_e','EX_alac__S_e', 'EX_diact_e', 'EX_actn__R_e']
    # dataframes = dataframes_for_robustness(dataframes,compounds)
    # cols = [ (compounds), ('0.0','5.0'), ('Sth', 'La', 'Lr', 'Lh')]
    # cols = get_column_name(cols)
    # line_plot(directory= _plots_directory,
    #           filename='results_production_robustness_analysis',
    #           dataframes= dataframes,
    #           column_x='carbon_uptake',
    #           column_y=cols,
    #           x_label='Carbon Source Uptake ($mmol.gDW^{-1}.h^{-1}$)',
    #           ylabel='Production Rate ($h^{-1}$)',
    #           grid_kws = {"hspace": 0.3, "wspace":0.3},
    #           figsize=(12,16),
    #           subplots_number = (5,2) ,
    #           title='Robustness Analysis',
    #           titles = ['Lactate', 'Acetate', "Ethanol", 'Formate', 'Fumarate', 'Succinate', 'Acetaldehyde','Acetolactate', 'Diacetyl', 'Acetoin'],
    #           group_by = (['blue','green','red','black','blue','green','red','black'],['-','-','-','-','--','--','--','--']),
    #           personalized_labels =  ['O2=0, ' +  r"$\it{S. thermophilus}$ LMD-9", 'O2=0, ' + "$\it{L. acidophilus}$ La-14", 'O2=0, ' + "$\it{L. rhamnosus}$ GG",'O2=0, ' +  "$\it{L. helveticus}$ CNRZ32",
    #                                   'O2=5, ' +  r"$\it{S. thermophilus}$ LMD-9", 'O2=5, ' + "$\it{L. acidophilus}$ La-14", 'O2=5, ' + "$\it{L. rhamnosus}$ GG",'O2=5, ' +  "$\it{L. helveticus}$ CNRZ32"],
    #           color_map={"Sth": 'blue', 'La': 'green', 'Lh': 'red', 'Lr': 'black'},
    #           linestyle_map = {'0.0': '-', '5.0':'--'}
    #           )


    ##### Maximizing for growth
    dataframes = read_results(_results_directory + "results_carbon_source_analysis" + ".xlsx")
    dataframes = [append_growth_rate(dataframes)]
    plot_gr(directory=_plots_directory,
              filename='growth_rate',
              dataframes=dataframes,
              column_x="carbon_uptake",
              column_y=["Sth", "La", "Lr","Lh"],
              x_label='Percentage of Carbon Uptake (%)',
              ylabel='Growth Rate ($h^{-1}$)',
              subplots_number=1,
              title='Maximizing for Growth',
              titles=None,
              figsize=(12,4),
              personalized_labels={'Sth':  r"$\it{S. thermophilus}$ LMD-9", 'La':"$\it{L. acidophilus}$ La-14", 'Lr': "$\it{L. rhamnosus}$ GG",'Lh':  "$\it{L. helveticus}$ CNRZ32" },
              )

    dataframes = read_results(_results_directory + "results_carbon_source_analysis" + ".xlsx")
    dataframes = change_df(dataframes, [1,1,1,1])
    # dataframes=apply_symetric(dataframes)
    line_plot(directory=_plots_directory,
              filename='results_carbon_source_analysis',
              dataframes=dataframes,
              column_x= "carbon_uptake",
              column_y=["ac_e", "etoh_e","for_e"],
              x_label='Carbon Source Uptake ($mmol.gDW^{-1}.h^{-1}$)',
              ylabel='Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
              subplots_number= (2,2),
              axis=0,
              # title= 'Maximizing for Growth',
              title=None,
              figsize=(12, 8),
              column_y_axis_2= ["lac__L_e"],
              ylabel_axis2= "Lactate production ($mmol.gDW^{-1}.h^{-1}$)",
              personalized_labels={'ac_e': 'Acetate', 'etoh_e': 'Ethanol', 'for_e' : "Formate", "lac__L_e":  'Lactate', 'Sth':''},
              ylim = {"La": (0,0.7)}
              )

    combine_images(directory=_plots_directory,
                   filename = "maximizing_for_growth",
                   filenames = ["growth_rate.png", "results_carbon_source_analysis.png"])



    # ##### Minimizing substrate uptake
    # dataframes = read_results(_results_directory + "results_minimum_substrate_analysis" + ".xlsx")
    # # dataframes=apply_symetric(dataframes)
    # line_plot(directory=_plots_directory,
    #           filename='results_minimum_substrate_analysis',
    #           dataframes=dataframes,
    #           column_x="e_Biomass",
    #           column_y=["ac_e", "etoh_e","for_e"],
    #           x_label='Growth Rate ($h^{-1}$)',
    #           ylabel='Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #           axis=1,
    #           title='Minimizing Substrate Uptake',
    #           column_y_axis_2=["lac__L_e"],
    #           ylabel_axis2="Lactate production ($mmol.gDW^{-1}.h^{-1}$)",
    #           personalized_labels={'ac_e': 'Acetate', 'etoh_e': 'Ethanol', 'for_e': "Formate", "lac__L_e": 'Lactate'},
    #           ylim={"La": (0, 1)}
    #           )

    # # # #####
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_no_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # dataframes = apply_symetric(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_no_growth_NOGrowth',
    #          dataframes=dataframes,
    #          column_x= ['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["e_Biomass"],
    #          x_label=' expression',
    #          y_label='Growth Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title = 'Growth Rate',
    #          columns=False,
    #          )
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # dataframes = apply_symetric(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_no_growth_Growth',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["e_Biomass"],
    #          x_label=' expression',
    #          y_label='Growth Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title='Growth Rate',
    #          columns=False,
    #          ylim = [0.1,0.1]
    #          )
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_no_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_no_growth_Lactate',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["lac__L_e"],
    #          x_label=' expression',
    #          y_label='Lactate Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title = 'L-Lactate',
    #          ylim= {"Sth": [2,2], "Lr": [2,2],"La": [2,2], "Lh": [2,2]}
    #          )
    #
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_no_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_no_growth_Acetate',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["ac_e"],
    #          x_label=' expression',
    #          y_label='Acetate Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title = 'Acetate',
    #          ylim = {"Lr": [0.5, 3]}
    #          )
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_no_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_no_growth_Formate',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["for_e"],
    #          x_label=' expression',
    #          y_label='Formate Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title='Formate'
    #          )
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_no_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_no_growth_Ethanol',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["etoh_e"],
    #          x_label=' expression',
    #          y_label='Ethanol Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title='Ethanol'
    #          )
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_growth_Lactate',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["lac__L_e"],
    #          x_label=' expression',
    #          y_label='Lactate Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title='L-Lactate',
    #          ylim= {"Sth": [2,2], "Lr": [2,2],"La": [2,2], "Lh": [2,2]}
    #          )
    #
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_growth_Acetate',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["ac_e"],
    #          x_label=' expression',
    #          y_label='Acetate Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title='Acetate'
    #          )
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_growth_Formate',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["for_e"],
    #          x_label=' expression',
    #          y_label='Formate Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title='Formate'
    #          )
    # dataframes = read_results(_results_directory + "results_enzyme_fva_analysis_growth" + ".xlsx")
    # dataframes = dataframes_for_fva(dataframes)
    # fva_plot(directory=_plots_directory,
    #          filename='results_enzyme_fva_analysis_growth_Ethanol',
    #          dataframes=dataframes,
    #          column_x=['PFL', 'PKETX', 'PFL', 'PKETX'],
    #          column_y=["etoh_e"],
    #          x_label=' expression',
    #          y_label='Ethanol Production Rate ($mmol.gDW^{-1}.h^{-1}$)',
    #          title='Ethanol',
    #          ylim= {"Sth": [0,0.1], "Lr": [0.85,0.4]}
    #          )

    # dataframes = read_results(_results_directory + "results_phenotypic_phase_plane_analysis" + ".xlsx")
    # phenotypic_phase_plane(directory=_plots_directory,
    #          filename='results_phenotypic_phase_plane_analysis',
    #          dataframes=dataframes,
    #          column_x= 'carbon uptake',
    #          column_y='Second axis',
    #          column_z='e_Biomass',
    #          figsize = (12,12),
    #          x_label='Acetate',
    #          y_label='Carbon Source',
    #          z_label='Growth rate ($h^{-1}$)',
    #          title='Phenotypic Phase Plane'
    #          )