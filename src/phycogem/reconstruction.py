from __future__ import annotations
import json
import re
from pathlib import Path

from cobra import Reaction, Model


def remove_shuttle_reactions(
    model: Model, allowed_compartments: set = {"c", "e", "p"}
) -> Model:
    """
    Remove shuttle reactions between unwanted compartments.

    Args:
        model (Model): _description_
        allowed_compartments (set, optional): _description_. Defaults to {"c", "e", "p"}.

    Returns:
        Model: _description_
    """
    shuttle_rxns_in_unwanted_compartments = [
        rxn
        for rxn in model.reactions
        if (
            (len(rxn.compartments) > 1)
            and (not rxn.compartments.issubset(allowed_compartments))
        )
    ]
    model.remove_reactions(shuttle_rxns_in_unwanted_compartments, remove_orphans=True)
    return model


def move_reactions_to_cytoplasm(
    model: Model, allowed_compartments: set = {"c", "e", "p"}
) -> Model:
    """
    Update the metabolites of a reaction to include a new set of metabolites
    Args:
        reaction (Reaction): _description_
        allowed_compartments (set, optional): _description_. Defaults to {"c", "e", "p"}.

    Returns:
        Reaction: _description_
    """
    reactions_to_add = []
    reactions_to_remove = []
    for reaction in model.reactions:
        if not reaction.compartments.issubset(allowed_compartments):

            new_metabolites = {}
            for metabolite, stoich in reaction.metabolites.items():
                new_met_id = metabolite.id[:-1] + "c"
                new_metabolite = (
                    model.metabolites.get_by_id(new_met_id)
                    if new_met_id in model.metabolites
                    else metabolite.copy()
                )
                new_metabolite.compartment = "c"
                new_metabolite.id = new_met_id
                if new_met_id not in model.metabolites:
                    model.add_metabolites([new_metabolite])
                    model.remove_metabolites([metabolite])
                new_metabolites[new_metabolite] = stoich

            new_reaction = Reaction(
                id=reaction.id,
                name=reaction.name,
                lower_bound=reaction.lower_bound,
                upper_bound=reaction.upper_bound,
                subsystem=reaction.subsystem,
            )
            new_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            new_reaction.add_metabolites(new_metabolites)
            reactions_to_add.append(new_reaction)
            reactions_to_remove.append(reaction)

    model.remove_reactions(reactions_to_remove, remove_orphans=True)
    model.add_reactions(reactions_to_add)
    return model


def get_duplicated_reactions(memote_report: Path) -> list[list]:
    """Retrieve duplicated reactions from memote report

    Args:
        memote_report (Path): path to report (JSON file)

    Returns:
        list: list of lists of duplicated reaction pair IDs
    """
    with open(memote_report, "r") as f:
        report = json.load(f)
    tests = report["tests"]
    return tests["test_find_duplicated_reactions"]["data"]


def remove_duplicated_reactions(
    model: Model, duplicated_reactions: list[list], inplace: bool = True
) -> Model:
    """Remove duplicated reactions from model

    Args:
        model (Model): cobra model object
        duplicated_reactions (list[list]):list lists of duplicated pairs
        inplace (bool, optional): removes reactions in the input model or return a new one.
            Defaults to True.

    Returns:
        Model: model object without duplicated reactions
    """
    if not inplace:
        result_model = model.copy()
    else:
        result_model = model
    reactions_to_remove = [rxn_pair[0] for rxn_pair in duplicated_reactions]
    result_model.remove_reactions(reactions_to_remove, remove_orphans=True)
    return result_model


def extract_chemical_elements(formula: str) -> dict:
    """
    Extract the chemical components from a chemical formula.
    """
    components = re.findall(r"([A-Z][a-z]*)(\d*)", formula)
    component_dict = {}
    for element, count in components:
        component_dict[element] = int(count) if count else 1
    return component_dict


def is_co2(chemical_elements: dict) -> bool:
    """Check whether molecule is CO2.

    Args:
        chemical_elements (dict): dictionary of chemical elements

    Returns:
        bool: True if CO2, False otherwise
    """
    return chemical_elements == {"C": 1, "O": 2}


def is_hco3(chemical_elements: dict) -> bool:
    """Check whether molecule is HCO3-.

    Args:
        chemical_elements (dict): dictionary of chemical elements

    Returns:
        bool: True if HCO3-, False otherwise
    """
    return chemical_elements == {"C": 1, "H": 1, "O": 3}


def get_organic_exchanges(model: Model) -> list[str]:
    """
    Get IDs of all organic exchanges in a model.

    Parameters
    ----------
    model : cobra.Model
        Model to get exchanges from.

    Returns
    -------
    list
        List of exchange IDs.
    """
    organic_exchanges = []
    for rxn in model.exchanges:
        met = [met for met in rxn.metabolites][0]
        if met.formula is not None:
            chemical_elements = extract_chemical_elements(met.formula)
            if (
                (not is_co2(chemical_elements))
                or (not is_hco3(chemical_elements))
                or ("C" in chemical_elements)
            ):
                organic_exchanges.append(rxn.id)
    return organic_exchanges


def get_inorganic_exchanges(model: Model) -> list[str]:
    """
    Get IDs of all inorganic exchanges in a model.

    Parameters
    ----------
    model : cobra.Model
        Model to get exchanges from.

    Returns
    -------
    list
        List of exchange IDs.
    """
    inorganic_exchanges = []
    for rxn in model.exchanges:
        met = [met for met in rxn.metabolites][0]
        if met.formula is not None:
            chemical_elements = extract_chemical_elements(met.formula)
            if (
                (is_co2(chemical_elements))
                or (is_hco3(chemical_elements))
                or ("C" not in chemical_elements)
            ):
                inorganic_exchanges.append(rxn.id)
    return inorganic_exchanges


def open_inorganic_exchanges(
    model: Model,
    lower_bound: float = -1000,
    upper_bound: float = 1000,
    include: list = None,
    inplace: bool = True,
) -> Model:
    """
    Open all inorganic exchanges in a model and close organic exchanges,
    except for indicated ones.

    Parameters
    ----------
    model : cobra.Model
        Model to open exchanges in.
    include: list, optional
        List of exchange IDs to be opened besides inorganic ones. Default is None.
    inplace : bool, optional
        If True, open exchanges in place. If False, return a copy of the model with opened
        exchanges. Default is True.

    Returns
    -------
    cobra.Model
        Model with opened exchanges.
    """
    if inplace:
        model = model
    else:
        model = model.copy()
    for rxn in model.exchanges:
        if rxn.id in include:
            rxn.lower_bound = lower_bound
            rxn.upper_bound = upper_bound
        else:
            met = [met for met in rxn.metabolites][0]
            if met.formula is not None:
                chemical_elements = extract_chemical_elements(met.formula)
                if "C" not in chemical_elements:
                    rxn.lower_bound = lower_bound
                    rxn.upper_bound = upper_bound
                else:
                    rxn.lower_bound = 0
    return model
