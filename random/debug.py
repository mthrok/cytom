import logging
_LG = logging.getLogger(__name__)


def print_list(list_, indent=0, brackets=True):
    if not list_:
        if brackets:
            _LG.info('{}[]'.format(' '*indent))
        return
    if brackets:
        _LG.info('{}['.format(' '*indent))
    for item in list_:
        if isinstance(item, dict):
            print_dict(item, indent+2)
        elif isinstance(item, list):
            print_list(item, indent+2)
        else:
            _LG.info('{}{}'.format(' '*(indent+2), item))
    if brackets:
        _LG.info('{}]'.format(' '*indent))


def print_dict(dict_, indent=0, brackets=True):
    if not dict_:
        if brackets:
            _LG.info('{}{{}}'.format(' '*indent))
        return
    if brackets:
        _LG.info('{}{{'.format(' '*indent))
    for key, value in dict_.items():
        if isinstance(value, dict):
            if value:
                _LG.info('{}{}: {{'.format(' '*(indent+2), key))
                print_dict(value, indent=indent+2, brackets=False)
                _LG.info('{}}}'.format(' '*(indent+2)))
            else:
                _LG.info('{}{}: {{}}'.format(' '*(indent+2), key))
        elif isinstance(value, list):
            if value:
                _LG.info('{}{}: ['.format(' '*(indent+2), key))
                print_list(value, indent=indent+2, brackets=False)
                _LG.info('{}]'.format(' '*(indent+2)))
            else:
                _LG.info('{}{}: []'.format(' '*(indent+2), key))
        else:
            _LG.info('{}{}: {}'.format(' '*(indent+2), key, value))
    if brackets:
        _LG.info('{}}}'.format(' '*indent))
