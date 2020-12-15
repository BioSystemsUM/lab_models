from os.path import join
from cobra.io import read_sbml_model, save_json_model
from pandas import DataFrame

from flux_analysis import ModelAnalysis


def to_json(directory,
            models_directory,
            icc390=None,
            icc431=None,
            icc464=None,
            icc644=None, ):
    if icc390:
        model = read_sbml_model(join(directory, models_directory, icc390))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))

    if icc431:
        model = read_sbml_model(join(directory, models_directory, icc431))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))

    if icc464:
        model = read_sbml_model(join(directory, models_directory, icc464))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))

    if icc644:
        model = read_sbml_model(join(directory, models_directory, icc644))
        save_json_model(model, join(directory, models_directory, f'{model.id}.json'))


def fluxes_distribution(directory,
                        results_directory,
                        conditions_directory,
                        icc390=None,
                        icc431=None,
                        icc464=None,
                        icc644=None, ):
    conditions = 'wild_type_conditions.xlsx'

    if icc390:
        icc390 = ModelAnalysis(directory=directory,
                               model=icc390[0],
                               biomass_reaction=icc390[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        res = icc390.fluxes_distribution(conditions=conditions,
                                         sheet=icc390.model.id)

        res = DataFrame(res)

        res.to_csv(join(directory, results_directory, f'results_flux_distribution_{icc390.model.id}.csv'),
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

    if icc464:
        icc464 = ModelAnalysis(directory=directory,
                               model=icc464[0],
                               biomass_reaction=icc464[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        res = icc464.fluxes_distribution(conditions=conditions,
                                         sheet=icc464.model.id)

        res = DataFrame(res)

        res.to_csv(join(directory, results_directory, f'results_flux_distribution_{icc464.model.id}.csv'),
                   index_label='ID')

    if icc644:
        icc644 = ModelAnalysis(directory=directory,
                               model=icc644[0],
                               biomass_reaction=icc644[1],
                               results_directory=results_directory,
                               conditions_directory=conditions_directory)

        res = icc644.fluxes_distribution(conditions=conditions,
                                         sheet=icc644.model.id)

        res = DataFrame(res)

        res.to_csv(join(directory, results_directory, f'results_flux_distribution_{icc644.model.id}.csv'),
                   index_label='ID')


if __name__ == '__main__':
    import os

    _directory = os.getcwd()
    _results_directory = os.path.join(_directory, 'results')
    _conditions_directory = os.path.join(_directory, 'environmental_conditions')

    biomass_rxn = 'e_Biomass'

    icc390_model = 'iCC390.xml'
    icc431_model = 'iCC431.xml'
    icc464_model = 'iCC464.xml'
    icc644_model = 'iCC644.xml'

    to_json(directory=_directory,
            models_directory='models/',
            icc390=icc390_model,
            icc431=icc431_model,
            icc464=icc464_model,
            icc644=icc644_model)

    icc390_model = 'models/iCC390.xml'
    _icc390 = (icc390_model, biomass_rxn)

    icc431_model = 'models/iCC431.xml'
    _icc431 = (icc431_model, biomass_rxn)

    icc464_model = 'models/iCC464.xml'
    _icc464 = (icc464_model, biomass_rxn)

    icc644_model = 'models/iCC644.xml'
    _icc644 = (icc644_model, biomass_rxn)

    fluxes_distribution(directory=_directory, results_directory=_results_directory,
                        conditions_directory=_conditions_directory,
                        icc390=_icc390,
                        icc431=_icc431,
                        icc464=_icc464,
                        icc644=_icc644)
