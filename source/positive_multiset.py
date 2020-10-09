from typing import TypeVar, Dict, Union, Generic

T = TypeVar('T')


class PositiveMultiset(Dict[T, Union[int, float]]):
    def __init__(self, base_type: Generic[T], dictionary_of_values: Dict[T, Union[int, float]], allow_infinity=False):
        if dictionary_of_values:
            for item in dictionary_of_values.keys():
                if type(item) is not base_type:
                    raise AssertionError(f"Expected item of type {base_type} and instead got {type(item)}")
                count = dictionary_of_values[item]
                if (type(count) is int and count >= 1) or (allow_infinity and count == float("inf")):
                    pass
                else:
                    raise AssertionError(
                        f"count of {type(item)} '{item}' is not a positive integer or infinity: {count}"
                    )

        super().__init__(dictionary_of_values)
