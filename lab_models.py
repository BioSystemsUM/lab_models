from numpy import arange
from cobra.flux_analysis import pfba

from flux_analysis import ModelAnalysis


def s_thermophilus(directory,
                   biomass_rxn,
                   model,
                   results_directory=None,
                   model_atp=None,
                   atp_tuning=True,
                   summary=True,
                   connectivity=True,
                   topological_analysis=True,
                   carbon_sources=True,
                   amino_acids=True,
                   minimal_requirements=True,
                   robustness=True,
                   ppp=True,
                   carbon_source=True,
                   minimum_substrate=True,
                   enzyme_fva=True,
                   ):

    print("Streptococcus thermophilus")
    print()

    if atp_tuning:
        sth_tuning = ModelAnalysis(directory=directory, model=model_atp, biomass_reaction=biomass_rxn,
                                   results_directory=results_directory)

        print("ATP tuning")
        atp_tuning_conditions = 'validation/atp_tuning_conditions.xlsx'
        atp_tuning_sheet_name = 'Input'
        atp_tuning_output = 'validation/sth_atp_tuning.xlsx'
        sth_tuning.atp_tuning(conditions=atp_tuning_conditions,
                              sheet=atp_tuning_sheet_name,
                              atp_hydrolysis='ATP_Maintenance',
                              output=atp_tuning_output)

    sth = ModelAnalysis(directory=directory, model=model, biomass_reaction=biomass_rxn,
                        results_directory=results_directory)

    if summary:
        print("Summary")
        conditions = 'validation/wild_type_conditions.xlsx'
        conditions_sheet = 'Input'
        summary_output = 'sth_wild_type_simulation.xlsx'
        sth.summary(conditions=conditions, sheet=conditions_sheet, output=summary_output)

    if connectivity:
        print("Connectivity")
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
                            'pyr',]

        sth.connectivity(main_metabolites=main_metabolites)

    if topological_analysis:
        print("Topological analysis")
        sth.topological_analysis(output='sth_topological_analysis.xlsx')

    if carbon_sources:
        print("Carbon sources")
        carbon_sources_conditions = 'validation/carbon_sources_conditions.xlsx'
        carbon_sources_sheet_name = 'Input'
        carbon_sources_results = 'validation/sth_carbon_sources.xlsx'
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
            'EX_rib__D_e': (-27.624, 999999),
        }

        sth.carbon_sources(conditions=carbon_sources_conditions,
                           sheet=carbon_sources_sheet_name,
                           carbon_sources=carbon_sources,
                           output=carbon_sources_results,
                           minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

    if amino_acids:
        print("Amino acids")
        amino_acids_conditions = 'validation/amino_acids_conditions.xlsx'
        amino_acids_sheet_name = 'Input'
        amino_acids_results = 'sth_amino_acids.xlsx'
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

        sth.amino_acids(conditions=amino_acids_conditions,
                        sheet=amino_acids_sheet_name,
                        amino_acids=amino_acids,
                        output=amino_acids_results,
                        minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

    if minimal_requirements:
        print("Minimal requirements")
        minimal_requirements_conditions = 'validation/minimal_requirements_conditions.xlsx'
        minimal_requirements_sheet_name = 'Input'
        minimal_requirements_results = 'sth_minimal_requirements.xlsx'

        sth.minimal_requirements(conditions=minimal_requirements_conditions,
                                 sheet=minimal_requirements_sheet_name,
                                 output=minimal_requirements_results,
                                 minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

    if robustness:
        print("Robustness analysis")
        robustness_analysis_conditions = 'validation/robustness_analysis_conditions.xlsx'
        analysis = ['unconstrained', 'constrained', 'constrained_biomass', 'constrained_o2']
        robustness_analysis_results = 'sth_robustness_analysis.xlsx'
        columns_to_drop = ['carbon_source',
                           'flux_minimum',
                           'carbon_yield_minimum',
                           'mass_yield_minimum',
                           'carbon_yield_maximum',
                           'mass_yield_maximum']

        for sheet in analysis:

            reaction = 'EX_lcts_e'

            if sheet == 'constrained_o2':
                reaction = 'EX_o2_e'

            sth.robustness_analysis(conditions=robustness_analysis_conditions,
                                    sheet=sheet,
                                    reaction=reaction,
                                    output=robustness_analysis_results,
                                    output_sheet=sheet,
                                    columns_to_drop=columns_to_drop,
                                    minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

        print("PFL KO Robustness analysis")
        with sth.model:
            sth.get_reaction('PFLh').bounds = (0.0, 0.0)
            robustness_analysis_conditions = 'validation/robustness_analysis_conditions.xlsx'
            analysis = ['unconstrained', 'constrained', 'constrained_biomass', 'constrained_o2']
            robustness_analysis_results = 'sth_robustness_analysis_pfl_ko.xlsx'
            columns_to_drop = ['carbon_source',
                               'flux_minimum',
                               'carbon_yield_minimum',
                               'mass_yield_minimum',
                               'carbon_yield_maximum',
                               'mass_yield_maximum']

            for sheet in analysis:

                reaction = 'EX_lcts_e'

                if sheet == 'constrained_o2':
                    reaction = 'EX_o2_e'

                sth.robustness_analysis(conditions=robustness_analysis_conditions,
                                        sheet=sheet,
                                        reaction=reaction,
                                        output=robustness_analysis_results,
                                        output_sheet=sheet,
                                        columns_to_drop=columns_to_drop,
                                        minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

    if ppp:
        print("Phenotypic phase plane analysis")
        ppp_analysis_conditions = 'validation/ppp_analysis_conditions.xlsx'
        ppp_analysis_sheet = 'ppp'
        ppp_analysis_results = 'sth_ppp_analysis.xlsx'
        columns_to_drop = ['carbon_source',
                           'flux_minimum',
                           'carbon_yield_minimum',
                           'mass_yield_minimum',
                           'carbon_yield_maximum',
                           'mass_yield_maximum']

        sth.phenotypic_phase_plane_analysis(conditions=ppp_analysis_conditions,
                                            sheet=ppp_analysis_sheet,
                                            reactions=['EX_lcts_e', 'EX_o2_e'],
                                            output=ppp_analysis_results,
                                            output_sheet=ppp_analysis_sheet,
                                            columns_to_drop=columns_to_drop,
                                            minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

        print("PFL KO Phenotypic phase plane analysis")
        with sth.model:
            sth.get_reaction('PFLh').bounds = (0.0, 0.0)
            ppp_analysis_conditions = 'validation/ppp_analysis_conditions.xlsx'
            ppp_analysis_sheet = 'ppp'
            ppp_analysis_results = 'validation/sth_ppp_analysis_pfl_ko.xlsx'
            columns_to_drop = ['carbon_source',
                               'flux_minimum',
                               'carbon_yield_minimum',
                               'mass_yield_minimum',
                               'carbon_yield_maximum',
                               'mass_yield_maximum']

            sth.phenotypic_phase_plane_analysis(conditions=ppp_analysis_conditions,
                                                sheet=ppp_analysis_sheet,
                                                reactions=['EX_lcts_e', 'EX_o2_e'],
                                                output=ppp_analysis_results,
                                                output_sheet=ppp_analysis_sheet,
                                                columns_to_drop=columns_to_drop,
                                                minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

    if carbon_source:
        print("Carbon source analysis")
        carbon_source_analysis_conditions = 'validation/carbon_source_analysis_conditions.xlsx'
        carbon_source_analysis_sheet_name = 'Input'
        carbon_source_analysis_results = 'sth_carbon_source_analysis.xlsx'

        sth.carbon_source_analysis(conditions=carbon_source_analysis_conditions,
                                   sheet=carbon_source_analysis_sheet_name,
                                   carbon_source='EX_lcts_e',
                                   carbon_source_linspace=list(range(2, 42, 2)),
                                   output=carbon_source_analysis_results,
                                   minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

    if minimum_substrate:
        print("Minimum substrate analysis")
        minimum_substrate_analysis_conditions = 'validation/minimum_substrate_analysis_conditions.xlsx'
        minimum_substrate_analysis_sheet_name = 'Input'
        minimum_substrate_analysis_results = 'sth_minimum_substrate_analysis.xlsx'
        objective_linspace = arange(0.1, 1.6, 0.1)

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

        sth.minimum_substrate_analysis(conditions=minimum_substrate_analysis_conditions,
                                       sheet=minimum_substrate_analysis_sheet_name,
                                       substrates=substrates,
                                       objective=sth.biomass_reaction.id,
                                       objective_linspace=objective_linspace,
                                       output=minimum_substrate_analysis_results,
                                       minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

    if enzyme_fva:
        print("Enzyme FVA analysis no constraint")
        enzyme_fva_analysis_conditions = 'validation/enzyme_fva_analysis_conditions.xlsx'
        enzyme_fva_analysis_sheet_name = 'Input'
        enzyme_fva_analysis_results = 'sth_enzyme_fva_analysis_no_constraint.xlsx'

        rxns_to_track = [
            'PDHm',
        ]

        sth.enzyme_fva_analysis(conditions=enzyme_fva_analysis_conditions,
                                sheet=enzyme_fva_analysis_sheet_name,
                                enzyme='PFLh',
                                constraint_growth=False,
                                rxns_to_track=rxns_to_track,
                                output=enzyme_fva_analysis_results,
                                minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)

        print("Enzyme FVA analysis constraint")
        enzyme_fva_analysis_conditions = 'validation/enzyme_fva_analysis_conditions.xlsx'
        enzyme_fva_analysis_sheet_name = 'Input'
        enzyme_fva_analysis_results = 'sth_enzyme_fva_analysis_constraint.xlsx'

        rxns_to_track = [
            'PDHm',
        ]

        sth.enzyme_fva_analysis(conditions=enzyme_fva_analysis_conditions,
                                sheet=enzyme_fva_analysis_sheet_name,
                                enzyme='PFLh',
                                constraint_growth=True,
                                rxns_to_track=rxns_to_track,
                                output=enzyme_fva_analysis_results,
                                minimum_growth=pfba(sth.model)[biomass_rxn] * 0.1)


if __name__ == "__main__":
    # _directory = "C:/Users/ferna/OneDrive - Universidade do Minho/Pycharm_projects/sth/"
    _directory = "C:/Users/BiSBII/OneDrive - Universidade do Minho/Pycharm_projects/sth/"
    _results_directory = "C:/Users/BiSBII/OneDrive - Universidade do Minho/Pycharm_projects/sth/validation"
    _biomass_rxn = 'e-Biomass'
    _model = 'sth_model_26_11_2020_v4.xml'
    _model_atp = 'sth_model_26_11_2020_v2_atp_tuning.xml'

    s_thermophilus(directory=_directory, biomass_rxn=_biomass_rxn, model=_model, results_directory=_results_directory,
                   model_atp=_model_atp,
                   atp_tuning=False,
                   summary=False,
                   connectivity=True,
                   topological_analysis=False,
                   carbon_sources=False,
                   amino_acids=False,
                   minimal_requirements=False,
                   robustness=False,
                   ppp=False,
                   carbon_source=False,
                   minimum_substrate=False,
                   enzyme_fva=False)
