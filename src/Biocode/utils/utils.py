import json


def list_to_str(lst):
    return json.dumps(lst)


def str_to_list(string):
    return json.loads(string)
