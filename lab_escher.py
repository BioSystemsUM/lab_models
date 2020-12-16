from os.path import join
from cobra.io import read_sbml_model, save_json_model
from pandas import DataFrame

from flux_analysis import ModelAnalysis


def to_json(directory,
            models_directory,
            icc389=None,
            icc431=None,
            icc470=None,
            icc651=None, ):
    if icc389:
        model = read_sbml_model(join(directory, models_directory, icc389))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))

    if icc431:
        model = read_sbml_model(join(directory, models_directory, icc431))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))

    if icc470:
        model = read_sbml_model(join(directory, models_directory, icc470))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))

    if icc651:
        model = read_sbml_model(join(directory, models_directory, icc651))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))


def fluxes_distribution(directory,
                        results_directory,
                        conditions_directory,
                        icc389=None,
                        icc431=None,
                        icc470=None,
                        icc651=None, ):
    conditions = 'wild_type_conditions.xlsx'

    if icc389:
        icc389 = ModelAnalysis(directory=directory,
                               model=icc389[0],
                               biomass_reaction=icc389[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        res = icc389.fluxes_distribution(conditions=conditions,
                                         sheet=icc389.model.id)

        res = DataFrame(res)

        res.to_csv(join(directory, results_directory, f'results_flux_distribution_{icc389.model.id}.csv'),
                   index_label='ID')

    if icc431:
        icc431 = ModelAnalysis(directory=directory,
                               model=icc431[0],
                               biomass_reaction=icc431[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        res = icc431.fluxes_distribution(conditions=conditions,
                                         sheet=icc431.model.id)

        res = DataFrame(res)

        res.to_csv(join(directory, results_directory, f'results_flux_distribution_{icc431.model.id}.csv'),
                   index_label='ID')

    if icc470:
        icc470 = ModelAnalysis(directory=directory,
                               model=icc470[0],
                               biomass_reaction=icc470[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        res = icc470.fluxes_distribution(conditions=conditions,
                                         sheet=icc470.model.id)

        res = DataFrame(res)

        res.to_csv(join(directory, results_directory, f'results_flux_distribution_{icc470.model.id}.csv'),
                   index_label='ID')

    if icc651:
        icc651 = ModelAnalysis(directory=directory,
                               model=icc651[0],
                               biomass_reaction=icc651[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        res = icc651.fluxes_distribution(conditions=conditions,
                                         sheet=icc651.model.id)

        res = DataFrame(res)

        res.to_csv(join(directory, results_directory, f'results_flux_distribution_{icc651.model.id}.csv'),
                   index_label='ID')


if __name__ == '__main__':
    import os

    _directory = os.getcwd()
    _results_directory = os.path.join(_directory, 'results')
    _conditions_directory = os.path.join(_directory, 'environmental_conditions')

    biomass_rxn = 'e_Biomass'

    icc389_model = 'iCC389.xml'
    icc431_model = 'iCC431.xml'
    icc470_model = 'iCC470.xml'
    icc651_model = 'iCC651.xml'

    to_json(directory=_directory,
            models_directory='models/',
            icc389=icc389_model,
            icc431=icc431_model,
            icc470=icc470_model,
            icc651=icc651_model)

    icc389_model = 'models/iCC389.xml'
    _icc389 = (icc389_model, biomass_rxn)

    icc431_model = 'models/iCC431.xml'
    _icc431 = (icc431_model, biomass_rxn)

    icc470_model = 'models/iCC470.xml'
    _icc470 = (icc470_model, biomass_rxn)

    icc651_model = 'models/iCC651.xml'
    _icc651 = (icc651_model, biomass_rxn)

    fluxes_distribution(directory=_directory, results_directory=_results_directory,
                        conditions_directory=_conditions_directory,
                        icc389=_icc389,
                        icc431=_icc431,
                        icc470=_icc470,
                        icc651=_icc651)
