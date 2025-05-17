"""
This script contains all the functions that are defined for use
in the scripts for the project.
"""

def split_by_batch_prefix(name):
    """
    Splits the afca indices into two parts:
    - Part before 'AFCA' or 'FCA'
    - The rest starting with 'AFCA' or 'FCA'
    """
    match = re.search(r'(AFCA|FCA)', name)
    i = match.start()
    print(name[i:])

