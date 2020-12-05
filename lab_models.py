from collections import namedtuple

from cobra.flux_analysis import pfba
from numpy import arange
from pandas import DataFrame

from flux_analysis import ModelAnalysis


def analysis_pipeline(models_analysis, analysis):
    Solution = namedtuple(typename='Solution', field_names=[model_analysis.model.id
                                                            for model_analysis in models_analysis])

    _results = {}

    for model_analysis, _analysis in zip(models_analysis, analysis):
        _results[model_analysis.model.id] = run_analysis(model_analysis=model_analysis,
                                                         analysis=_analysis)

    # noinspection PyArgumentList
    return Solution(**_results)


def run_analysis(model_analysis, analysis):
    if not analysis:
        analysis = {}

    field_names = list(analysis.keys())

    field_names.insert(0, 'id')
    field_names.insert(1, 'model_analysis')
    ModelSolution = namedtuple(typename='ModelSolution', field_names=field_names)

    _results = {'id': model_analysis.model.id, 'model_analysis': model_analysis}

    for _analysis, configuration in analysis.items():
        print(f'Running {_analysis} for {model_analysis.model.id}')
        method = getattr(model_analysis, _analysis, lambda **kwargs: DataFrame())

        _results[_analysis] = method(**configuration)

    # noinspection PyArgumentList
    return ModelSolution(**_results)


def build_analysis(model_analysis,
                   model_id,
                   conditions=None,
                   atp_tuning_conditions=None,
                   carbon_sources_conditions=None,
                   analysis_to_drop=None,
                   atp='atp_c',
                   h2o='h2o_c',
                   adp='adp_c',
                   h='h_c',
                   pi='pi_c',
                   growth_atp_linspace=None,
                   atp_m='ATP_Maintenance',
                   atp_m_linspace=None,
                   atp_hydrolysis='ATP_Maintenance',
                   main_metabolites=None,
                   carbon_sources=None,
                   amino_acids=None,
                   carbon_source='EX_lcts_e',
                   carbon_source_linspace=None,
                   substrates=None,
                   objective=None,
                   objective_linspace=None,
                   enzyme='PFLh',
                   rxns_to_track=None,
                   minimum_growth=None,
                   ):

    if not conditions:
        conditions = 'wild_type_conditions.xlsx'

    if not atp_tuning_conditions:
        atp_tuning_conditions = conditions

    if not carbon_sources_conditions:
        carbon_sources_conditions = conditions

    if not analysis_to_drop:
        analysis_to_drop = []

    if not growth_atp_linspace:
        growth_atp_linspace = arange(0, 52, 2)

    if not atp_m_linspace:
        atp_m_linspace = arange(0.0, 9.5, 0.5)

    if not main_metabolites:
        main_metabolites = ['h',
                            'h2o',
                            'atp',
                            'pi',
                            'adp',
                            'ppi',
                            'nad',
                            'nadh',
                            'glu__L_',
                            'nadp',
                            'nadph',
                            'NH4',
                            'amp',
                            'co2',
                            'pyr',
                            'ACP',
                            'CoA']

    if not carbon_sources:
        carbon_sources = {
            'EX_lcts_e': (-27.624, 999999),
            'EX_sucr_e': (-13.812, 999999),
            'EX_cellb_e': (-13.812, 999999),
            'EX_xyl__D_e': (-27.624, 999999),
            'EX_arab__D_e': (-27.624, 999999),
            'EX_glc__aD_e': (-27.624, 999999),
            'EX_glc_D_B_e': (-27.624, 999999),
            'EX_fru_B_e': (-27.624, 999999),
            'EX_man_e': (-27.624, 999999),
            'EX_gal_e': (-27.624, 999999),
            'EX_rib__D_e': (-27.624, 999999)}

    if not amino_acids:
        amino_acids = ['EX_leu__L_e',
                       'EX_glu__L__e',
                       'EX_phe__L_e',
                       'EX_val__L_e',
                       'EX_asn__L_e',
                       'EX_trp__L_e',
                       'EX_tyr__L_e',
                       'EX_met__L_e',
                       'EX_ile__L_e',
                       'EX_gln__L_e',
                       'EX_gly_e',
                       'EX_lys__L_e',
                       'EX_his__L_e',
                       'EX_asp__L_e',
                       'EX_arg__L_e',
                       'EX_cys__L_e',
                       'EX_pro__L_e',
                       'EX_ser__L_e',
                       'EX_ala__L_e',
                       'EX_thr__L_e', ]

    if not carbon_source_linspace:
        carbon_source_linspace = list(range(2, 42, 2))

    if not substrates:
        substrates = [
            'EX_lcts_e',
            'EX_leu__L_e',
            'EX_glu__L__e',
            'EX_phe__L_e',
            'EX_val__L_e',
            'EX_asn__L_e',
            'EX_trp__L_e',
            'EX_tyr__L_e',
            'EX_met__L_e',
            'EX_ile__L_e',
            'EX_gln__L_e',
            'EX_gly_e',
            'EX_lys__L_e',
            'EX_his__L_e',
            'EX_asp__L_e',
            'EX_arg__L_e',
            'EX_cys__L_e',
            'EX_pro__L_e',
            'EX_ser__L_e',
            'EX_ala__L_e',
            'EX_thr__L_e',
        ]

    if not objective:
        objective = model_analysis.biomass_reaction.id

    if not objective_linspace:
        objective_linspace = arange(0.1, 1.6, 0.1)

    if not rxns_to_track:
        rxns_to_track = ['PDHm']

    if minimum_growth is None:
        minimum_growth = pfba(model_analysis.model)[model_analysis.biomass_reaction.id] * 0.1

    conditions_sheet = model_id
    output_sheet = model_id

    analysis_configuration = {}

    output = 'growth_atp_tuning_results.xlsx'
    analysis_configuration['growth_atp_tuning'] = dict(conditions=conditions,
                                                       sheet=conditions_sheet,
                                                       atp=atp,
                                                       h2o=h2o,
                                                       adp=adp,
                                                       h=h,
                                                       pi=pi,
                                                       growth_atp_linspace=growth_atp_linspace,
                                                       output=output,
                                                       output_sheet=output_sheet)

    output = 'maintenance_atp_tuning_results.xlsx'
    analysis_configuration['maintenance_atp_tuning'] = dict(conditions=conditions,
                                                            sheet=conditions_sheet,
                                                            atp_m=atp_m,
                                                            atp_m_linspace=atp_m_linspace,
                                                            output=output,
                                                            output_sheet=output_sheet)

    output = 'atp_tuning.xlsx'
    analysis_configuration['atp_tuning'] = dict(conditions=atp_tuning_conditions,
                                                sheet=conditions_sheet,
                                                atp_hydrolysis=atp_hydrolysis,
                                                output=output,
                                                output_sheet=output_sheet)

    output = 'wild_type_simulation_results.xlsx'
    analysis_configuration['summary'] = dict(conditions=conditions, sheet=conditions_sheet, output=output)

    output = 'connectivity_results.xlsx'
    analysis_configuration['connectivity'] = dict(main_metabolites=main_metabolites,
                                                  output=output,
                                                  output_sheet=output_sheet)

    output = 'topological_analysis_results.xlsx'
    analysis_configuration['topological_analysis'] = dict(output=output,
                                                          output_sheet=output_sheet)

    output = 'carbon_sources_results.xlsx'
    analysis_configuration['carbon_sources'] = dict(conditions=carbon_sources_conditions,
                                                    sheet=conditions_sheet,
                                                    carbon_sources=carbon_sources,
                                                    output=output,
                                                    output_sheet=output_sheet,
                                                    minimum_growth=minimum_growth)

    output = 'amino_acids_results.xlsx'
    analysis_configuration['amino_acids'] = dict(conditions=conditions,
                                                 sheet=conditions_sheet,
                                                 amino_acids=amino_acids,
                                                 output=output,
                                                 output_sheet=output_sheet,
                                                 minimum_growth=minimum_growth)

    output = 'minimal_requirements_results.xlsx'
    analysis_configuration['minimal_requirements'] = dict(conditions=conditions,
                                                          sheet=conditions_sheet,
                                                          output=output,
                                                          output_sheet=output_sheet,
                                                          minimum_growth=minimum_growth)

    output = 'carbon_source_results.xlsx'
    analysis_configuration['carbon_source_analysis'] = dict(conditions=conditions,
                                                            sheet=conditions_sheet,
                                                            carbon_source=carbon_source,
                                                            carbon_source_linspace=carbon_source_linspace,
                                                            output=output,
                                                            output_sheet=output_sheet,
                                                            minimum_growth=minimum_growth)

    output = 'minimum_substrate_results.xlsx'
    analysis_configuration['minimum_substrate_analysis'] = dict(conditions=conditions,
                                                                sheet=conditions_sheet,
                                                                substrates=substrates,
                                                                objective=objective,
                                                                objective_linspace=objective_linspace,
                                                                output=output,
                                                                output_sheet=output_sheet,
                                                                minimum_growth=minimum_growth)

    output = 'enzyme_fva_growth.xlsx'
    analysis_configuration['enzyme_fva_analysis_growth'] = dict(conditions=conditions,
                                                                sheet=conditions_sheet,
                                                                enzyme=enzyme,
                                                                rxns_to_track=rxns_to_track,
                                                                output=output,
                                                                output_sheet=output_sheet,
                                                                minimum_growth=minimum_growth)

    output = 'enzyme_fva_no_growth.xlsx'
    analysis_configuration['enzyme_fva_analysis_no_growth'] = dict(conditions=conditions,
                                                                   sheet=output_sheet,
                                                                   enzyme=enzyme,
                                                                   rxns_to_track=rxns_to_track,
                                                                   output=output,
                                                                   output_sheet=output_sheet,
                                                                   minimum_growth=minimum_growth)

    for del_analysis in analysis_to_drop:
        del analysis_configuration[del_analysis]

    return analysis_configuration


def lab_models(directory,
               results_directory,
               conditions_directory,
               icc390=None,
               icc431=None,
               icc464=None,
               icc644=None,
               ):

    _models_analysis = []
    _analysis = []

    if icc390:
        icc390 = ModelAnalysis(directory=directory,
                               model=icc390[0],
                               biomass_reaction=icc390[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        _models_analysis.append(icc390)

        icc390_analysis = build_analysis(model_analysis=icc390,
                                         model_id=icc390.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_drop=['growth_atp_tuning', 'atp_tuning',
                                                           'maintenance_atp_tuning'], )

        _analysis.append(icc390_analysis)

    if icc431:
        icc431 = ModelAnalysis(directory=directory,
                               model=icc431[0],
                               biomass_reaction=icc431[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        _models_analysis.append(icc431)

        icc431_analysis = build_analysis(model_analysis=icc431,
                                         model_id=icc431.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_drop=['growth_atp_tuning', 'atp_tuning',
                                                           'maintenance_atp_tuning'], )

        _analysis.append(icc431_analysis)

    if icc464:
        icc464 = ModelAnalysis(directory=directory,
                               model=icc464[0],
                               biomass_reaction=icc464[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        _models_analysis.append(icc464)

        icc464_analysis = build_analysis(model_analysis=icc464,
                                         model_id=icc464.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_drop=['growth_atp_tuning', 'atp_tuning',
                                                           'maintenance_atp_tuning'], )

        _analysis.append(icc464_analysis)

    if icc644:
        icc644 = ModelAnalysis(directory=directory,
                               model=icc644[0],
                               biomass_reaction=icc644[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        _models_analysis.append(icc644)

        icc644_analysis = build_analysis(model_analysis=icc644,
                                         model_id=icc644.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_drop=['growth_atp_tuning', 'atp_tuning',
                                                           'maintenance_atp_tuning'], )

        _analysis.append(icc644_analysis)

    return analysis_pipeline(models_analysis=_models_analysis, analysis=_analysis)


if __name__ == '__main__':
    _directory = 'C:/Users/ferna/OneDrive - Universidade do Minho/Pycharm_projects/lab_models/'
    _results_directory = _directory + 'results/'
    _conditions_directory = _directory + 'environmental_conditions/'

    icc390_biomass_rxn = 'e_Biomass'
    icc390_model = 'models/iCC390.xml'
    _icc390 = (icc390_model, icc390_biomass_rxn)

    icc431_biomass_rxn = 'e_Biomass'
    icc431_model = 'models/iCC431.xml'
    _icc431 = (icc431_model, icc431_biomass_rxn)

    icc464_biomass_rxn = 'e_Biomass'
    icc464_model = 'models/iCC464.xml'
    _icc464 = (icc464_model, icc464_biomass_rxn)

    icc644_biomass_rxn = 'e_Biomass'
    icc644_model = 'models/iCC644.xml'
    _icc644 = (icc644_model, icc644_biomass_rxn)

    lab_models(directory=_directory,
               results_directory=_results_directory,
               conditions_directory=_conditions_directory,
               # icc390=_icc390,
               icc431=_icc431,
               # icc464=_icc464,
               # icc644=_icc644,
               )
