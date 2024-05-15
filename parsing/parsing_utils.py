import argparse


def validate_arg_bounds(arg: str, lower: int, upper: int) -> int:
    """
    Validate whether argument is within specified bounds.
    :param arg: argument to be validated
    :param lower: lower bound
    :param upper: upper bound
    :return:
    """
    number = int(arg)

    if number < lower or upper < number:
        raise argparse.ArgumentTypeError(
            f"{number} is not from interval [{lower}, {upper}]")

    return number
