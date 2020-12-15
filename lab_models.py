from collections import namedtuple
from functools import partial

from cobra.flux_analysis import pfba
from numpy import arange, linspace

from flux_analysis import ModelAnalysis

Analysis = namedtuple('Analysis', ['model_analysis', 'configuration'])


def analysis_pipeline(analysis):
    Solution = namedtuple(typename='Solution', field_names=[_analysis.model_analysis.model.id
                                                            for _analysis in analysis])

    _results = {}

    for _analysis in analysis:
        _results[_analysis.model_analysis.model.id] = run_analysis(analysis=_analysis)

    # noinspection PyArgumentList
    return Solution(**_results)


def run_analysis(analysis: Analysis):
    _results = {'analysis': analysis}

    for _analysis, method in analysis.configuration.items():

        # debug
        # print(f'Running {_analysis} for {analysis.model_analysis.model.id}')
        # _results[_analysis] = method()

        try:
            print(f'Running {_analysis} for {analysis.model_analysis.model.id}')

            _results[_analysis] = method()

        except Exception as ex:
            print(f'Failed running {_analysis} for {analysis.model_analysis.model.id}')
            print(ex)
            print(ex.args)
            with open('Logger.txt', 'a') as text_file:
                text_file.write(f'Failed running {_analysis} for {analysis.model_analysis.model.id} \n')

    ModelSolution = namedtuple(typename='ModelSolution', field_names=list(_results.keys()))

    # noinspection PyArgumentList
    return ModelSolution(**_results)


def build_analysis(model_analysis,
                   model_id,
                   conditions=None,
                   atp_tuning_conditions=None,
                   carbon_sources_conditions=None,
                   ppp_conditions=None,
                   analysis_to_drop=None,
                   analysis_to_build=None,
                   atp='atp_c',
                   h2o='h2o_c',
                   adp='adp_c',
                   h='h_c',
                   pi='pi_c',
                   growth_atp_linspace=None,
                   atp_m='ATPM',
                   atp_m_linspace=None,
                   atp_hydrolysis='ATPM',
                   main_metabolites=None,
                   carbon_sources=None,
                   amino_acids=None,
                   carbon_source='EX_lcts_e',
                   carbon_source_linspace=None,
                   reactions=None,
                   points=20,
                   dense_output=True,
                   columns_to_drop=None,
                   production_exchanges=None,
                   oxygen_exchange=None,
                   oxygen_linspace=None,
                   growth_fraction=0.9,
                   substrates=None,
                   objective=None,
                   objective_linspace=None,
                   enzyme='PFLh',
                   rxns_to_track=None,
                   tol=1E-5,
                   special_conditions=None,
                   minimum_growth=None,
                   ):

    if not conditions:
        conditions = 'wild_type_conditions.xlsx'

    if not atp_tuning_conditions:
        atp_tuning_conditions = conditions

    if not carbon_sources_conditions:
        carbon_sources_conditions = conditions

    if not ppp_conditions:
        ppp_conditions = conditions

    if not analysis_to_build:
        to_build = {key: val for key, val in model_analysis.get_analysis().items()}

    else:
        registered = model_analysis.get_analysis()
        to_build = {key: registered[key] for key in analysis_to_build}

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

    if not reactions:
        reactions = []

    if not columns_to_drop:
        columns_to_drop = []

    if not production_exchanges:
        production_exchanges = [
            'EX_co2_e',
            'EX_lac__L_e',
            'EX_ac_e',
            'EX_for_e',
            'EX_etoh_e',
            'EX_fum_e',
            'EX_succ_e',
            'EX_acald_e',
            'EX_alac__S_e',
            'EX_diact_e',
            'EX_actn__R_e',
        ]

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
        rxns_to_track = []

    if not special_conditions:
        special_conditions = {}

    if minimum_growth is None:
        minimum_growth = pfba(model_analysis.model)[model_analysis.biomass_reaction.id] * 0.1

    conditions_sheet = model_id
    output_sheet = model_id

    for analysis in analysis_to_drop:
        del to_build[analysis]

    configuration = {}

    for analysis, method in to_build.items():

        output = f'results_{analysis}.xlsx'

        if analysis == 'carbon_sources':
            _conditions = carbon_sources_conditions

        elif analysis == 'atp_tuning':

            _conditions = atp_tuning_conditions

        elif analysis == 'phenotypic_phase_plane_analysis':

            _conditions = ppp_conditions

        else:

            _conditions = conditions

        configuration[analysis] = partial(method,
                                          model_analysis,
                                          conditions=_conditions,
                                          sheet=conditions_sheet,
                                          atp=atp,
                                          h2o=h2o,
                                          adp=adp,
                                          h=h,
                                          pi=pi,
                                          growth_atp_linspace=growth_atp_linspace,
                                          atp_m=atp_m,
                                          atp_m_linspace=atp_m_linspace,
                                          atp_hydrolysis=atp_hydrolysis,
                                          main_metabolites=main_metabolites,
                                          carbon_sources=carbon_sources,
                                          amino_acids=amino_acids,
                                          carbon_source=carbon_source,
                                          carbon_source_linspace=carbon_source_linspace,
                                          reactions=reactions,
                                          points=points,
                                          dense_output=dense_output,
                                          columns_to_drop=columns_to_drop,
                                          production_exchanges=production_exchanges,
                                          oxygen_exchange=oxygen_exchange,
                                          oxygen_linspace=oxygen_linspace,
                                          growth_fraction=growth_fraction,
                                          substrates=substrates,
                                          objective=objective,
                                          objective_linspace=objective_linspace,
                                          enzyme=enzyme,
                                          rxns_to_track=rxns_to_track,
                                          tol=tol,
                                          special_conditions=special_conditions,
                                          output=output,
                                          output_sheet=output_sheet,
                                          minimum_growth=minimum_growth,
                                          )

    return Analysis(model_analysis=model_analysis, configuration=configuration)


def lab_models_atp(directory,
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

        icc390_analysis = build_analysis(model_analysis=icc390,
                                         model_id=icc390.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_build=['maintenance_atp_tuning'], )

        _analysis.append(icc390_analysis)

    if icc431:
        icc431 = ModelAnalysis(directory=directory,
                               model=icc431[0],
                               biomass_reaction=icc431[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        icc431_analysis = build_analysis(model_analysis=icc431,
                                         model_id=icc431.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_build=['growth_atp_tuning'], )

        _analysis.append(icc431_analysis)

    if icc464:
        icc464 = ModelAnalysis(directory=directory,
                               model=icc464[0],
                               biomass_reaction=icc464[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        icc464_analysis = build_analysis(model_analysis=icc464,
                                         model_id=icc464.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_build=['growth_atp_tuning'], )

        _analysis.append(icc464_analysis)

    if icc644:
        icc644 = ModelAnalysis(directory=directory,
                               model=icc644[0],
                               biomass_reaction=icc644[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        icc644_analysis = build_analysis(model_analysis=icc644,
                                         model_id=icc644.model.id,
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         analysis_to_build=['growth_atp_tuning'], )

        _analysis.append(icc644_analysis)

    return analysis_pipeline(analysis=_analysis)


def lab_models(directory,
               results_directory,
               conditions_directory,
               analysis_to_build=None,
               analysis_to_drop=None,
               icc390=None,
               icc431=None,
               icc464=None,
               icc644=None,
               ):
    _models_analysis = []
    _analysis = []

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

    amino_acids = ['EX_leu__L_e',
                   'EX_glu__L_e',
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

    columns_to_drop = ['carbon_source',
                       'flux_minimum',
                       'carbon_yield_minimum',
                       'mass_yield_minimum',
                       'carbon_yield_maximum',
                       'mass_yield_maximum']

    oxygen_exchange = 'EX_o2_e'
    oxygen_linspace = linspace(0, 5, 3)
    growth_fraction = 0.9

    if icc390:
        icc390 = ModelAnalysis(directory=directory,
                               model=icc390[0],
                               biomass_reaction=icc390[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        carbon_sources = {
            'EX_lcts_e': (-4.69, 999999),
            'EX_sucr_e': (-4.69, 999999),
            'EX_glc__aD_e': (-9.38, 999999),
            'EX_glc_D_B_e': (-9.38, 999999),
            'EX_fru_B_e': (-9.38, 999999),
            'EX_man_e': (-9.38, 999999),
            'EX_a_gal__D_e': (-9.38, 999999),
            'EX_gal_bD_e': (-9.38, 999999),
            'EX_arbt_e': (-4.69, 999999),
            'EX_malt_e': (-4.69, 999999),
            'EX_tre_e': (-4.69, 999999),
        }

        production_exchanges = [
            'EX_co2_e',
            'EX_lac__L_e',
            'EX_lac__D_e',
            'EX_ac_e',
            'EX_for_e',
            'EX_etoh_e',
            'EX_cit_e',
            'EX_succ_e',
            'EX_acald_e',
        ]

        substrates = [
            'EX_glc__aD_e',
            'EX_leu__L_e',
            'EX_glu__L_e',
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

        reactions = ['EX_glc__aD_e', 'EX_ac_e']

        carbon_source_linspace = list(linspace(1, 20, 20))

        objective_linspace = list(linspace(0.01, 0.3, 15))

        rxns_to_track = []

        icc390_analysis = build_analysis(model_analysis=icc390,
                                         model_id=icc390.model.id,
                                         conditions='wild_type_conditions.xlsx',
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         ppp_conditions='ppp_conditions.xlsx',
                                         analysis_to_build=analysis_to_build,
                                         analysis_to_drop=analysis_to_drop,
                                         main_metabolites=main_metabolites,
                                         carbon_sources=carbon_sources,
                                         amino_acids=amino_acids,
                                         carbon_source='EX_glc__aD_e',
                                         carbon_source_linspace=carbon_source_linspace,
                                         reactions=reactions,
                                         points=20,
                                         dense_output=True,
                                         columns_to_drop=columns_to_drop,
                                         production_exchanges=production_exchanges,
                                         oxygen_exchange=oxygen_exchange,
                                         oxygen_linspace=oxygen_linspace,
                                         growth_fraction=growth_fraction,
                                         substrates=substrates,
                                         objective=icc390.biomass_reaction.id,
                                         objective_linspace=objective_linspace,
                                         enzyme='PKETX',
                                         rxns_to_track=rxns_to_track,
                                         special_conditions={'EX_lac__D_e': (0.0, 0.0)},
                                         tol=1E-4,
                                         )

        _analysis.append(icc390_analysis)

    if icc431:
        icc431 = ModelAnalysis(directory=directory,
                               model=icc431[0],
                               biomass_reaction=icc431[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

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

        production_exchanges = [
            'EX_co2_e',
            'EX_lac__L_e',
            'EX_ac_e',
            'EX_for_e',
            'EX_etoh_e',
            'EX_fum_e',
            'EX_succ_e',
            'EX_acald_e',
            'EX_alac__S_e',
            'EX_diact_e',
            'EX_actn__R_e',
        ]

        substrates = [
            'EX_lcts_e',
            'EX_leu__L_e',
            'EX_glu__L_e',
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

        reactions = ['EX_lcts_e', 'EX_ac_e']

        carbon_source_linspace = list(range(2, 42, 2))

        objective_linspace = list(linspace(0.1, 1.5, 15))

        rxns_to_track = ['PDHm']

        icc431_analysis = build_analysis(model_analysis=icc431,
                                         model_id=icc431.model.id,
                                         conditions='wild_type_conditions.xlsx',
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         ppp_conditions='ppp_conditions.xlsx',
                                         analysis_to_build=analysis_to_build,
                                         analysis_to_drop=analysis_to_drop,
                                         main_metabolites=main_metabolites,
                                         carbon_sources=carbon_sources,
                                         amino_acids=amino_acids,
                                         carbon_source='EX_lcts_e',
                                         carbon_source_linspace=carbon_source_linspace,
                                         reactions=reactions,
                                         points=20,
                                         dense_output=True,
                                         columns_to_drop=columns_to_drop,
                                         production_exchanges=production_exchanges,
                                         oxygen_exchange=oxygen_exchange,
                                         oxygen_linspace=oxygen_linspace,
                                         growth_fraction=growth_fraction,
                                         substrates=substrates,
                                         objective=icc431.biomass_reaction.id,
                                         objective_linspace=objective_linspace,
                                         enzyme='PFL',
                                         rxns_to_track=rxns_to_track,
                                         tol=1E-4,
                                         )

        _analysis.append(icc431_analysis)

    if icc464:
        icc464 = ModelAnalysis(directory=directory,
                               model=icc464[0],
                               biomass_reaction=icc464[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        carbon_sources = {
            'EX_lcts_e': (-8.2, 999999),
            'EX_sucr_e': (-8.2, 999999),
            'EX_glc__aD_e': (-16.4, 999999),
            'EX_glc_D_B_e': (-16.4, 999999),
            'EX_fru_B_e': (-16.4, 999999),
            'EX_man_e': (-16.4, 999999),
            'EX_a_gal__D_e': (-16.4, 999999),
            'EX_gal_bD_e': (-16.4, 999999),
            'EX_arbt_e': (-8.2, 999999),
            'EX_malt_e': (-8.2, 999999),
            'EX_acgam_e': (-8.2, 999999),
            'EX_stys_e': (-4.1, 999999),
            'EX_tre_e': (-8.2, 999999),
            'EX_cellb_e': (-8.2, 999999),
            'EX_melib_e': (-8.2, 999999),
            'EX_raffin_e': (-5.5, 999999),
            'EX_rib__D_e': (-16.4, 999999),
            'EX_C16639_e': (-16.4, 999999),
        }

        production_exchanges = [
            'EX_co2_e',
            'EX_lac__L_e',
            'EX_lac__D_e',
            'EX_ac_e',
            'EX_for_e',
            'EX_etoh_e',
            'EX_cit_e',
            'EX_succ_e',
            'EX_acald_e',
            'EX_mal__L_e'
        ]

        substrates = [
            'EX_glc__aD_e',
            'EX_leu__L_e',
            'EX_glu__L_e',
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

        reactions = ['EX_glc__aD_e', 'EX_ac_e']

        carbon_source_linspace = list(linspace(1, 25, 20))

        objective_linspace = list(linspace(0.1, 1.5, 15))

        rxns_to_track = []

        icc464_analysis = build_analysis(model_analysis=icc464,
                                         model_id=icc464.model.id,
                                         conditions='wild_type_conditions.xlsx',
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         ppp_conditions='ppp_conditions.xlsx',
                                         analysis_to_build=analysis_to_build,
                                         analysis_to_drop=analysis_to_drop,
                                         main_metabolites=main_metabolites,
                                         carbon_sources=carbon_sources,
                                         amino_acids=amino_acids,
                                         carbon_source='EX_glc__aD_e',
                                         carbon_source_linspace=carbon_source_linspace,
                                         reactions=reactions,
                                         points=20,
                                         dense_output=True,
                                         columns_to_drop=columns_to_drop,
                                         production_exchanges=production_exchanges,
                                         oxygen_exchange=oxygen_exchange,
                                         oxygen_linspace=oxygen_linspace,
                                         growth_fraction=growth_fraction,
                                         substrates=substrates,
                                         objective=icc464.biomass_reaction.id,
                                         objective_linspace=objective_linspace,
                                         enzyme='PKETX',
                                         rxns_to_track=rxns_to_track,
                                         special_conditions={'EX_lac__D_e': (0.0, 0.0)},
                                         tol=1E-4,
                                         )

        _analysis.append(icc464_analysis)

    if icc644:
        icc644 = ModelAnalysis(directory=directory,
                               model=icc644[0],
                               biomass_reaction=icc644[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        carbon_sources = {
            'EX_lcts_e': (-4.415, 999999),
            'EX_sucr_e': (-4.415, 999999),
            'EX_glc__aD_e': (-8.83, 999999),
            'EX_glc_D_B_e': (-8.83, 999999),
            'EX_fru_B_e': (-8.83, 999999),
            'EX_man_e': (-8.83, 999999),
            'EX_a_gal__D_e': (-8.83, 999999),
            'EX_gal_bD_e': (-8.83, 999999),
            'EX_arbt_e': (-4.415, 999999),
            'EX_malt_e': (-4.415, 999999),
            'EX_acgam_e': (-4.415, 999999),
            'EX_tre_e': (-4.415, 999999),
            'EX_cellb_e': (-4.415, 999999),
            'EX_melib_e': (-4.415, 999999),
            'EX_rib__D_e': (-8.83, 999999),
            'EX_HC00832_e': (-8.83, 999999),
            'EX_tag__D_e': (-8.83, 999999),
            'EX_srb__L_e': (-8.83, 999999),
            'EX_HC00822_e': (-4.415, 999999),
            'EX_rmn_e': (-8.83, 999999),
        }

        production_exchanges = [
            'EX_co2_e',
            'EX_lac__L_e',
            'EX_lac__D_e',
            'EX_ac_e',
            'EX_for_e',
            'EX_etoh_e',
            'EX_cit_e',
            'EX_fum_e',
            'EX_succ_e',
            'EX_acald_e',
            'EX_mal__L_e',
            'EX_diact_e',
            'EX_actn__R_e',
            'EX_lald__L_e',
        ]

        substrates = [
            'EX_glc__aD_e',
            'EX_leu__L_e',
            'EX_glu__L_e',
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

        reactions = ['EX_glc__aD_e', 'EX_ac_e']

        carbon_source_linspace = list(linspace(1, 25, 20))

        objective_linspace = list(linspace(0.01, 0.6, 15))

        rxns_to_track = []

        icc644_analysis = build_analysis(model_analysis=icc644,
                                         model_id=icc644.model.id,
                                         conditions='wild_type_conditions.xlsx',
                                         carbon_sources_conditions='carbon_sources_conditions.xlsx',
                                         ppp_conditions='ppp_conditions.xlsx',
                                         analysis_to_build=analysis_to_build,
                                         analysis_to_drop=analysis_to_drop,
                                         main_metabolites=main_metabolites,
                                         carbon_sources=carbon_sources,
                                         amino_acids=amino_acids,
                                         carbon_source='EX_glc__aD_e',
                                         carbon_source_linspace=carbon_source_linspace,
                                         reactions=reactions,
                                         points=20,
                                         dense_output=True,
                                         columns_to_drop=columns_to_drop,
                                         production_exchanges=production_exchanges,
                                         oxygen_exchange=oxygen_exchange,
                                         oxygen_linspace=oxygen_linspace,
                                         growth_fraction=growth_fraction,
                                         substrates=substrates,
                                         objective=icc644.biomass_reaction.id,
                                         objective_linspace=objective_linspace,
                                         enzyme='PFL',
                                         rxns_to_track=rxns_to_track,
                                         special_conditions={'EX_lac__D_e': (0.0, 0.0)},
                                         tol=1E-4,
                                         )

        _analysis.append(icc644_analysis)

    return analysis_pipeline(analysis=_analysis)


if __name__ == '__main__':
    import os

    _directory = os.getcwd()
    _results_directory = os.path.join(_directory, 'results')
    _conditions_directory = os.path.join(_directory, 'environmental_conditions')

    biomass_rxn = 'e_Biomass'

    icc390_model = 'models/iCC390.xml'
    _icc390 = (icc390_model, biomass_rxn)

    icc431_model = 'models/iCC431.xml'
    _icc431 = (icc431_model, biomass_rxn)

    icc464_model = 'models/iCC464.xml'
    _icc464 = (icc464_model, biomass_rxn)

    icc644_model = 'models/iCC644.xml'
    _icc644 = (icc644_model, biomass_rxn)

    # _analysis_to_drop = ['growth_atp_tuning',
    #                      'atp_tuning',
    #                      'maintenance_atp_tuning',
    #                      'robustness_analysis']

    _analysis_to_build = ['enzyme_fva_analysis_no_growth', 'enzyme_fva_analysis_growth']

    _ = lab_models(directory=_directory,
                   results_directory=_results_directory,
                   conditions_directory=_conditions_directory,
                   analysis_to_build=_analysis_to_build,
                   # analysis_to_drop=_analysis_to_drop,
                   icc390=_icc390,
                   icc431=_icc431,
                   icc464=_icc464,
                   icc644=_icc644,
                   )
