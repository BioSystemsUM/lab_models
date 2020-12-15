from os.path import join
from cobra.io import read_sbml_model, save_json_model


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


if __name__ == '__main__':
    import os

    _directory = os.getcwd()

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