import molmass


def get_massnumber(elements: str) -> int:
    """
    Return the mass number of the given molecule's elements
    """
    return molmass.Formula(elements).isotope.massnumber
