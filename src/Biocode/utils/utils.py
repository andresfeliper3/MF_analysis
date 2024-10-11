import json
import re

def list_to_str(lst):
    return json.dumps(lst)


def str_to_list(string):
    return json.loads(string)

def remove_region_part(input_string: str) -> str:
    match = re.match(r"(.+)_region_\d+$", input_string)
    if match:
        return match.group(1)
    else:
        raise ValueError("Input string does not match the required pattern.")


