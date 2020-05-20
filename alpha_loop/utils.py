import os

#===============================================================================
# Class for string coloring
#===============================================================================
class bcolors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    WARNING = YELLOW
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    RED = '\033[91m'
    END = ENDC

def format_path(path):
    """Format the path in local format taking in entry a unix format"""
    if path[0] != '/':
        return os.path.join(*path.split('/'))
    else:
        return os.path.sep + os.path.join(*path.split('/'))

def ln(file_pos, starting_dir='.', name='', log=True, cwd=None, abspath=False):
    """a simple way to have a symbolic link without to have to change directory
    starting_point is the directory where to write the link
    file_pos is the file to link
    WARNING: not the linux convention
    """
    file_pos = format_path(file_pos)
    starting_dir = format_path(starting_dir)
    if not name:
        name = os.path.split(file_pos)[1]    

    if cwd:
        if not os.path.isabs(file_pos):
            file_pos = os.path.join(cwd, file_pos)
        if not os.path.isabs(starting_dir):
            starting_dir = os.path.join(cwd, starting_dir)        

    # Remove existing link if necessary
    path = os.path.join(starting_dir, name)
    if os.path.exists(path):
        if os.path.realpath(path) != os.path.realpath(file_pos):
            os.remove(os.path.join(starting_dir, name))
        else:
            return

    if not abspath:
        target = os.path.relpath(file_pos, starting_dir)
    else:
        target = file_pos

    try:
        os.symlink(target, os.path.join(starting_dir, name))
    except Exception as error:
        if log:
            print(error)
            print('Could not link %s at position: %s' % (file_pos, \
                                                os.path.realpath(starting_dir)))