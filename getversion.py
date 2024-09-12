import os
def get_version(changelog):
    '''automatically retrieves the version from changelog.md
    Parameters :
    ------------
    folder where the changelog.md is'''
    with open(changelog,'r') as f:
        version=''
        while not version:
            line = f.readline()
            if line.startswith('##') and 'dev' not in line:
                version = line.replace('##','').strip()
    return version