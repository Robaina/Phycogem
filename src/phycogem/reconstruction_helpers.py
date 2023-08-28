import re


def remove_compartment(met_id: str) -> str:
    return re.sub(r"_[a-z]$", "", met_id)


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
