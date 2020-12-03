import os
from cobra import util, Model, Solution
from cobra.io import read_sbml_model
from cobra.exceptions import Infeasible
from cobra.flux_analysis import flux_variability_analysis, pfba, production_envelope
from numpy import linspace
from pandas import Series, read_excel, DataFrame, concat, ExcelWriter
from sympy import Add


def writer(to_write):

    def writer_decorator(fn):

        def wrapper(*args, **kwargs):

            sol = fn(*args, **kwargs)

            output = kwargs.get('output', None)

            if output is None and to_write is False:
                return sol

            self = args[0]

            output = kwargs.get('output', f'{self.model.id}_{fn.__name__}_results.xlsx')
            output_sheet = kwargs.get('output_sheet', f'{self.model.id}')

            if output is False:
                return sol

            # the writer is needed for appending sheets
            path = os.path.join(self.results_directory, output)

            if os.path.exists(path):
                with ExcelWriter(path, mode='a') as _writer:
                    sol.to_excel(_writer, sheet_name=output_sheet)

                return sol

            with ExcelWriter(path) as _writer:
                sol.to_excel(_writer, sheet_name=output_sheet)

            return sol
        return wrapper
    return writer_decorator


class ModelAnalysis:

    def __init__(self,
                 directory,
                 model,
                 biomass_reaction,
                 set_objective=True,
                 results_directory=None):

        """ Model analysis class
        When an object is created, this function automatically loads the model (xml file) using cobrapy built-in
        methods.
        It also creates the attribute biomass_reaction (a cobrapy object of the biomass reaction),
        It sets the biomass_reaction as objective function of the model.

        """

        if not results_directory:
            results_directory = directory

        os.chdir(directory)

        self.directory = directory
        self.results_directory = results_directory
        self.model_file = model

        self.model: Model = read_sbml_model(os.path.join(self.directory, self.model_file))

        self.model.solver = 'cplex'

        self.biomass_reaction = self.get_reaction(biomass_reaction)

        if set_objective:
            self.model.objective = self.biomass_reaction.id

    def get_reaction(self, reaction_id):
        return self.model.reactions.get_by_id(reaction_id)

    def get_metabolite(self, metabolite_id):
        return self.model.metabolites.get_by_id(metabolite_id)

    def maximize(self, is_pfba=False, growth=True):

        if is_pfba:

            if growth:
                try:
                    return pfba(self.model).fluxes[list(util.solver.linear_reaction_coefficients(self.model))[0].id]

                except Infeasible:
                    return 0

            try:
                return pfba(self.model)

            except Infeasible:

                return Solution(objective_value=0, status='infeasible', fluxes=Series())

        if growth:
            try:
                return self.model.optimize(objective_sense="maximize").objective_value

            except Infeasible:
                return 0

        try:
            return self.model.optimize(objective_sense="maximize")

        except Infeasible:
            return Solution(objective_value=0, status='infeasible', fluxes=Series())

    def minimize(self, biomass_value=True):

        if biomass_value:

            try:
                return self.model.optimize(objective_sense="minimize").objective_value

            except Infeasible:
                return 0

        try:
            return self.model.optimize(objective_sense="minimize")

        except Infeasible:
            return Solution(objective_value=0, status='Infeasible', fluxes=Series())

    def apply_conditions(self, data_frame):

        exchanges = data_frame['exchange'].unique()

        for exchange in exchanges:
            rxn = self.get_reaction(exchange)

            ex_mask = data_frame['exchange'] == exchange
            lb = -data_frame.loc[ex_mask, 'qs'].iloc[0]
            rxn.lower_bound = lb

    @writer(to_write=False)
    def summary(self,
                conditions=False,
                sheet=False,
                fva=None,
                filter_val=0.0,
                objectives=True,
                metabolite_index=True,
                output=None,
                output_sheet=None,
                ):

        if conditions and sheet:
            df = read_excel(os.path.join(self.directory, conditions), sheet)

            with self.model:
                self.apply_conditions(data_frame=df)

                sol = self.model.summary(fva=fva)

        else:
            sol = self.model.summary(fva=fva)

        solution_frame = sol.to_frame()

        if objectives:

            # noinspection PyProtectedMember
            objs = list(sol._objective.keys())
            objective_id = 'objective'

            if objs:
                objective_id = objs[0].id

            # noinspection PyProtectedMember
            objective_value = sol._objective_value

            if fva:
                solution_frame.loc[objective_id] = [objective_id, objective_id, 1,
                                                    objective_value, objective_value, objective_value]

            else:
                solution_frame.loc[objective_id] = [objective_id, objective_id, 1, objective_value]

        if metabolite_index:
            solution_frame.index = list(solution_frame.loc[:, 'metabolite'])

        if filter_val is not None:
            if fva:
                filter_mask = solution_frame.loc[:, ['flux', 'minimum', 'maximum']] != filter_val

            else:
                filter_mask = solution_frame.loc[:, ['flux']] != filter_val

            filter_mask = filter_mask.any(1)
            solution_frame = solution_frame.loc[filter_mask, :]

        # if output:
        #     solution_frame.to_excel(os.path.join(self.directory, output), 'Output')

        return solution_frame

    @writer(to_write=False)
    def summary_from_solution(self,
                              solution,
                              fva=False,
                              filter_val=0.0,
                              objectives=True,
                              metabolite_index=True,
                              cols_to_drop=None,
                              output=None,
                              output_sheet=None,
                              ):

        if fva:
            solution = self.model.summary(fva=solution)

        else:
            solution = self.model.summary(solution=solution)

        # reaction, metabolite, factor, flux, minimum, maximum
        # index are exchange ids
        solution_frame = solution.to_frame()

        if objectives:

            # noinspection PyProtectedMember
            objs = list(solution._objective.keys())
            objective_id = 'objective'

            if objs:
                objective_id = objs[0].id

            # noinspection PyProtectedMember
            objective_value = solution._objective_value

            if fva:
                solution_frame.loc[objective_id] = [objective_id, objective_id, 1, objective_value,
                                                    objective_value, objective_value]

            else:
                solution_frame.loc[objective_id] = [objective_id, objective_id, 1, objective_value]

        if metabolite_index:
            solution_frame.index = list(solution_frame.loc[:, 'metabolite'])

        if cols_to_drop:
            solution_frame.drop(cols_to_drop, axis=1, inplace=True)

        if filter_val is not None:
            if fva:
                filter_mask = solution_frame.loc[:, ['flux', 'minimum', 'maximum']] != filter_val

            else:
                filter_mask = solution_frame.loc[:, ['flux']] != filter_val

            filter_mask = filter_mask.any(1)
            solution_frame = solution_frame.loc[filter_mask, :]

        return solution_frame

    # old cobra version
    # def summary_from_solution(self,
    #                           solution,
    #                           fva=False,
    #                           flatten=False,
    #                           inputs=True,
    #                           outputs=True,
    #                           objectives=True,
    #                           ):
    #
    #     if flatten:
    #
    #         if fva:
    #             solution = self.model.summary(fva=solution).to_frame()
    #
    #         else:
    #             solution = self.model.summary(solution=solution).to_frame()
    #
    #         index = []
    #         flux_values = []
    #         min_flux_values = []
    #         max_flux_values = []
    #
    #         if objectives:
    #             _objectives = list(solution.loc[:, 'OBJECTIVES'].loc[:, 'ID'])
    #             index.extend(_objectives)
    #
    #             objective_fluxes = list(solution.loc[:, 'OBJECTIVES'].loc[:, 'FLUX'])
    #             flux_values.extend(objective_fluxes)
    #             min_flux_values.extend(objective_fluxes)
    #             max_flux_values.extend(objective_fluxes)
    #
    #         if inputs:
    #             _inputs = list(solution.loc[:, 'IN_FLUXES'].loc[:, 'ID'])
    #             index.extend(_inputs)
    #
    #             in_fluxes = [-val for val in solution.loc[:, 'IN_FLUXES'].loc[:, 'FLUX']]
    #             flux_values.extend(in_fluxes)
    #
    #             if fva:
    #                 min_in_fluxes = [-val for val in solution.loc[:, 'IN_FLUXES'].loc[:, 'FLUX_MIN']]
    #                 min_flux_values.extend(min_in_fluxes)
    #
    #                 max_in_fluxes = [-val for val in solution.loc[:, 'IN_FLUXES'].loc[:, 'FLUX_MAX']]
    #                 max_flux_values.extend(max_in_fluxes)
    #
    #         if outputs:
    #             _outputs = list(solution.loc[:, 'OUT_FLUXES'].loc[:, 'ID'])
    #             index.extend(_outputs)
    #
    #             out_fluxes = list(solution.loc[:, 'OUT_FLUXES'].loc[:, 'FLUX'])
    #             flux_values.extend(out_fluxes)
    #
    #             if fva:
    #                 min_out_fluxes = list(solution.loc[:, 'OUT_FLUXES'].loc[:, 'FLUX_MIN'])
    #                 min_flux_values.extend(min_out_fluxes)
    #
    #                 max_out_fluxes = list(solution.loc[:, 'OUT_FLUXES'].loc[:, 'FLUX_MAX'])
    #                 max_flux_values.extend(max_out_fluxes)
    #
    #         if fva:
    #             values = {'flux': flux_values, 'min_flux': min_flux_values, 'max_flux': max_flux_values}
    #
    #         else:
    #             values = flux_values
    #
    #         return DataFrame(data=values, index=index).dropna()
    #
    #     if fva:
    #         return self.model.summary(fva=solution).to_frame()
    #
    #     return self.model.summary(solution=solution).to_frame()

    @writer(to_write=True)
    def connectivity(self,
                     compartments=None,
                     main_metabolites=None,
                     output=None,
                     output_sheet=None,
                     ):

        if not compartments:
            compartments = self.model.compartments

        index = {met.id[:-2] for met in self.model.metabolites}

        data = DataFrame(index=index, columns=list(compartments.keys()))

        for met in self.model.metabolites:

            _idd = met.id[:-2]

            data.loc[_idd, met.compartment] = len(met.reactions)

        cols = list(compartments.values())
        cols.append('sum')
        data.loc[:, 'sum'] = data.sum(axis=1)
        data.columns = cols

        if main_metabolites:
            _new_index = set(data.index) - set(main_metabolites)
            _new_index = list(main_metabolites) + list(_new_index)
            data = data.reindex(_new_index)

        data.fillna(0, inplace=True)
        # data.to_excel(os.path.join(self.directory, output), 'output')
        return data

    @writer(to_write=True)
    def topological_analysis(self,
                             compartments=None,
                             output=None,
                             output_sheet=None,
                             ):

        if not compartments:
            compartments = self.model.compartments

        # open output file
        cols = []
        _data = []
        for compartment in compartments.keys():
            cols.append(f'{compartment} reactions')
            _data.append(0)

        for compartment in compartments.keys():
            cols.append(f'{compartment} metabolites')
            _data.append(0)

        cols.extend(['exchange reactions', 'genes', 'GPRs'])
        _data.extend([0, 0, 0])

        data = DataFrame(data=[_data], columns=cols, index=['total'])

        for rxn in self.model.reactions:

            if rxn.boundary:
                data.loc['total', 'exchange reactions'] += 1

            for compartment in rxn.compartments:

                data.loc['total', f'{compartment} reactions'] += 1

            if rxn.gene_reaction_rule.strip() != '':
                data.loc['total', 'GPRs'] += 1

        for met in self.model.metabolites:

            data.loc['total', f'{met.compartment} metabolites'] += 1

        data.loc['total', 'genes'] = len(self.model.genes)

        # df = DataFrame.from_dict(data, orient='index', columns=col)
        # df.to_excel(os.path.join(self.directory, output), 'output')

        return data

    @writer(to_write=True)
    def atp_tuning(self,
                   conditions,
                   sheet,
                   atp_hydrolysis,
                   output=True,
                   output_sheet=True,
                   ):

        """
        Simulating the qATP of each correlation between qS and experimental growth rate, for 14 growth rates that
        correspond to experimental (0.05, 0.96, ...)

        The input must contain a column named "Exchanges" followed by all points
        of growth rate.

        For the amino acids calculate their content in the biomass reaction for each growth rate
        and full fill under the respective column. Leave all the others lower bounds open
        Exchange   0.1  0.6  0.97
        lactose    2    6    10
        lactose    5    12   27
        cysteine   0,4  0,8  1,2
        alanine    0.6  0.9  1.5

        """

        df = read_excel(os.path.join(self.directory, conditions), sheet)

        df.columns = [str(col) for col in df.columns]

        col = ['qATP']

        data = {}

        with self.model:

            for column in df:

                if column != 'exchange':

                    self.apply_conditions(data_frame=df)

                    try:
                        growth_rate = float(column)

                    except ValueError:

                        growth_rate = float(column.replace('.1', ''))

                    self.biomass_reaction.bounds = (growth_rate, growth_rate)

                    self.model.objective = atp_hydrolysis
                    value = self.maximize(is_pfba=True)

                    data[column] = [value]

        df = DataFrame.from_dict(data, orient='index', columns=col)
        return df
        # df.to_excel(os.path.join(self.directory, output), 'output')

    @writer(to_write=False)
    def carbon_sources(self,
                       conditions,
                       sheet,
                       carbon_sources,
                       output,
                       minimum_growth=0.1):

        df = read_excel(os.path.join(self.directory, conditions), sheet)

        with self.model:

            self.apply_conditions(data_frame=df)

            _growth = self.maximize(is_pfba=True)
            if _growth > minimum_growth:
                raise ValueError(f'Wrong Environmental Conditions, as growth rate is {_growth}')

            data = {}

            for carbon_source, bds in carbon_sources.items():
                with self.model:
                    self.get_reaction(carbon_source).bounds = bds

                    growth = self.maximize(is_pfba=True)

                    data[carbon_source] = [growth]

        cols = ['growth_rate']
        df = DataFrame.from_dict(data, orient='index', columns=cols)
        df.to_excel(os.path.join(self.directory, output), 'output')

    @writer(to_write=False)
    def amino_acids(self,
                    conditions,
                    sheet,
                    amino_acids,
                    output,
                    minimum_growth=0.1):

        df = read_excel(os.path.join(self.directory, conditions), sheet)

        with self.model:

            self.apply_conditions(data_frame=df)

            _growth = self.maximize(is_pfba=True)
            if _growth < minimum_growth:
                raise ValueError(f'Wrong Environmental Conditions, as growth rate is {_growth}')

            data = {}

            for amino_acid in amino_acids:
                with self.model:
                    self.get_reaction(amino_acid).lower_bound = 0.0

                    growth = self.maximize(is_pfba=True)

                    data[amino_acid] = [growth]

        cols = ['growth_rate']
        df = DataFrame.from_dict(data, orient='index', columns=cols)
        df.to_excel(os.path.join(self.directory, output), 'output')

    @writer(to_write=False)
    def minimal_requirements(self,
                             conditions,
                             sheet,
                             output,
                             minimum_growth=0.1):

        df = read_excel(os.path.join(self.directory, conditions), sheet)

        with self.model:

            self.apply_conditions(data_frame=df)

            _growth = self.maximize(is_pfba=True)
            if _growth < minimum_growth:
                raise ValueError(f'Wrong Environmental Conditions, as growth rate is {_growth}')

            exchanges = df["exchange"].unique()
            data = {}

            for nutrient in exchanges:
                with self.model:
                    self.get_reaction(nutrient).lower_bound = 0.0

                    growth = self.maximize(is_pfba=True)

                    data[nutrient] = [growth]

        cols = ['growth_rate']
        df = DataFrame.from_dict(data, orient='index', columns=cols)
        df.to_excel(os.path.join(self.directory, output), 'output')

    @writer(to_write=False)
    def _robustness_ppp(self,
                        conditions,
                        sheet,
                        reactions,
                        output,
                        output_sheet,
                        points=30,
                        carbon_source=None,
                        dense_output=False,
                        columns_to_drop=None,
                        minimum_growth=0.1,
                        ):

        if not columns_to_drop:
            columns_to_drop = []

        # just a wrapper for cobrapy production_envelope function
        df = read_excel(os.path.join(self.directory, conditions), sheet)

        with self.model:

            self.apply_conditions(data_frame=df)

            _growth = self.maximize(is_pfba=True)
            if _growth < minimum_growth:
                raise ValueError(f'Wrong Environmental Conditions, as growth rate is {_growth}')

            sol = production_envelope(model=self.model,
                                      reactions=reactions,
                                      points=points,
                                      carbon_sources=carbon_source)

        sol.drop(columns_to_drop, axis=1, inplace=True)

        # transform the sparse ppp dataframe to dense ppp dataframe
        if dense_output:

            new_sol = DataFrame()

            for row in sol.values:
                _val, _col, _row = row
                new_sol.loc[_row, _col] = _val

            sol = new_sol

            sol.fillna(0, inplace=True)

        # the writer is needed for appending sheets
        path = os.path.join(self.directory, output)

        if os.path.exists(path):
            with ExcelWriter(path, mode='a') as writer:
                sol.to_excel(writer, sheet_name=output_sheet)

            return

        with ExcelWriter(path) as writer:
            sol.to_excel(writer, sheet_name=output_sheet)

    def robustness_analysis(self,
                            conditions,
                            sheet,
                            reaction,
                            output,
                            output_sheet,
                            points=30,
                            carbon_source=None,
                            columns_to_drop=None,
                            minimum_growth=0.1,
                            ):

        return self._robustness_ppp(conditions=conditions,
                                    sheet=sheet,
                                    reactions=[reaction],
                                    output=output,
                                    output_sheet=output_sheet,
                                    points=points,
                                    carbon_source=carbon_source,
                                    dense_output=False,
                                    columns_to_drop=columns_to_drop,
                                    minimum_growth=minimum_growth)

    def phenotypic_phase_plane_analysis(self,
                                        conditions,
                                        sheet,
                                        reactions,
                                        output,
                                        output_sheet,
                                        points=30,
                                        carbon_source=None,
                                        dense_output=True,
                                        columns_to_drop=None,
                                        minimum_growth=0.1,
                                        ):

        return self._robustness_ppp(conditions=conditions,
                                    sheet=sheet,
                                    reactions=reactions,
                                    output=output,
                                    output_sheet=output_sheet,
                                    points=points,
                                    carbon_source=carbon_source,
                                    dense_output=dense_output,
                                    columns_to_drop=columns_to_drop,
                                    minimum_growth=minimum_growth)

    @writer(to_write=False)
    def carbon_source_analysis(self,
                               conditions,
                               sheet,
                               carbon_source,
                               carbon_source_linspace,
                               output,
                               minimum_growth=0.1):

        df = read_excel(os.path.join(self.directory, conditions), sheet)

        with self.model:

            self.apply_conditions(data_frame=df)

            _growth = self.maximize(is_pfba=True)
            if _growth < minimum_growth:
                raise ValueError(f'Wrong Environmental Conditions, as growth rate is {_growth}')

            data = []

            for uptake in carbon_source_linspace:
                with self.model:
                    self.get_reaction(carbon_source).lower_bound = -uptake

                    pfba_sol = self.maximize(is_pfba=True, growth=False)

                    sol = self.summary_from_solution(solution=pfba_sol,
                                                     cols_to_drop=['reaction', 'metabolite', 'factor'])

                    sol.columns = [f'{carbon_source}_{uptake}']

                    data.append(sol)

        data = concat(data, axis=1, join='outer')
        data.to_excel(os.path.join(self.directory, output), 'output')

    @writer(to_write=False)
    def minimum_substrate_analysis(self,
                                   conditions,
                                   sheet,
                                   substrates,
                                   objective,
                                   objective_linspace,
                                   output,
                                   minimum_growth=0.1,
                                   ):

        if not substrates:
            substrates = []

        df = read_excel(os.path.join(self.directory, conditions), sheet)

        with self.model:

            self.apply_conditions(data_frame=df)

            _growth = self.maximize(is_pfba=True)
            if _growth < minimum_growth:
                raise ValueError(f'Wrong Environmental Conditions, as growth rate is {_growth}')

            qp_expression = [arg
                             for substrate in substrates
                             for arg in self.get_reaction(substrate).flux_expression.args]

            qp_expression = Add(*qp_expression)

            qp_objective = self.model.problem.Objective(qp_expression, direction='max')

            self.model.objective = qp_objective
            self.maximize(is_pfba=False, growth=False)

            data = []

            for objective_rate in objective_linspace:
                objective_rate = round(objective_rate, 3)

                with self.model:
                    self.get_reaction(objective).bounds = (objective_rate, objective_rate)

                    qp_sol = self.maximize(is_pfba=False, growth=False)

                    sol = self.summary_from_solution(solution=qp_sol,
                                                     objectives=False,
                                                     cols_to_drop=['reaction', 'metabolite', 'factor']
                                                     )

                    sol.columns = [f'{objective}_{objective_rate}']

                    sol.loc[self.biomass_reaction.id] = qp_sol.fluxes[self.biomass_reaction.id]
                    sol.loc['objective_function'] = qp_sol.objective_value
                    sol.loc['status'] = qp_sol.status

                    data.append(sol)

        data = concat(data, axis=1, join='outer')
        data.to_excel(os.path.join(self.directory, output), 'output')

    @writer(to_write=False)
    def enzyme_fva_analysis(self,
                            conditions,
                            sheet,
                            enzyme,
                            output,
                            constraint_growth=False,
                            rxns_to_track=None,
                            minimum_growth=0.1,
                            ):

        if not rxns_to_track:
            rxns_to_track = []

        rxns_to_track = set(rxns_to_track + [enzyme])

        df = read_excel(os.path.join(self.directory, conditions), sheet)

        with self.model:

            self.apply_conditions(data_frame=df)

            _growth = self.maximize(is_pfba=True)
            if _growth < minimum_growth:
                raise ValueError(f'Wrong Environmental Conditions, as growth rate is {_growth}')

            # enzyme knock-out
            # it is important to get a new object and not the pointer
            enzyme_bounds = tuple(self.get_reaction(enzyme).bounds)

            # KO
            self.get_reaction(enzyme).bounds = (0.0, 0.0)

            # Growth rate for PFL knock-out phenotype
            growth = self.maximize(is_pfba=True, growth=True)

            # restore the bounds
            self.get_reaction(enzyme).bounds = enzyme_bounds

            # Forcing the growth rate for PFL knock-out phenotype
            biomass_rxn_bounds = tuple(self.get_reaction(self.biomass_reaction.id).bounds)
            self.get_reaction(self.biomass_reaction.id).bounds = (growth, growth)

            # Get min and max values of pfl when growth rate limited to PFL knock-out phenotype
            fva_sol = flux_variability_analysis(model=self.model,
                                                reaction_list=[self.get_reaction(enzyme)], )

            enzyme_min_flux = fva_sol.loc[enzyme, 'minimum']
            enzyme_max_flux = fva_sol.loc[enzyme, 'maximum']

            # enzyme flux/expression 10 levels (around 10% increment)
            enzyme_range = linspace(start=enzyme_min_flux,
                                    stop=enzyme_max_flux,
                                    num=11,
                                    endpoint=True)

            # growth rate can be constrained to the PFL knock-out phenotype or not
            if not constraint_growth:
                self.get_reaction(self.biomass_reaction.id).bounds = biomass_rxn_bounds

            data = []
            columns = []

            for enzyme_rate in enzyme_range:

                enzyme_rate = round(enzyme_rate, 4)
                columns.append(f'PFL_rate_{enzyme_rate}')

                with self.model:

                    # limiting enzyme rate to the previously determined
                    self.get_reaction(enzyme).bounds = (enzyme_rate, enzyme_rate)

                    # test if there is a solution
                    pfba_sol = self.maximize(is_pfba=True, growth=False)

                    reaction_list = [self.get_reaction(rxn) for rxn in rxns_to_track]

                    # fva for the reactions to keep track
                    fva_sol = flux_variability_analysis(model=self.model,
                                                        reaction_list=reaction_list)

                    # cobrapy model summary method can only accept fva solutions performed for the exchange reactions
                    fva_sol_exchanges = flux_variability_analysis(model=self.model,
                                                                  reaction_list=self.model.exchanges,
                                                                  pfba_factor=1.1,
                                                                  processes=1)

                    # final solution for the inputs and outputs flattened
                    sol = self.summary_from_solution(solution=fva_sol_exchanges,
                                                     fva=True,
                                                     objectives=True,
                                                     cols_to_drop=['reaction', 'metabolite', 'factor']
                                                     )

                    # add additional information
                    for rxn in rxns_to_track:

                        if pfba_sol.status == 'infeasible':
                            rxn_sol = [0, 0, 0]
                            sol.loc[rxn] = rxn_sol

                        else:
                            rxn_sol = [pfba_sol.fluxes.loc[rxn],
                                       fva_sol.loc[rxn, 'minimum'],
                                       fva_sol.loc[rxn, 'maximum']]

                            sol.loc[rxn] = rxn_sol

                        sol.loc['status'] = pfba_sol.status

                    data.append(sol)

        data = concat(data, axis=1, join='outer', keys=columns)
        data.to_excel(os.path.join(self.directory, output), 'output')
