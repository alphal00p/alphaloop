from scipy import stats
import math
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

def compute_dod(sequence):
    """ from a sequence of results in the format (scaling, complex(result)), determine the asymptotic scaling."""

    _THRESHOLD = 0.01

    xs = [math.log(abs(s[0])) for s in sequence]
    ys = [math.log(abs(s[1])) for s in sequence]

    slope, std_err = None, 0.0

    n_points_considered = 0

    if len(xs)<=5:
        n_points_considered = len(xs)
        slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
    else:
        n_entries = 5
        new_slope, new_std_err = None, 0.0
        while new_std_err<_THRESHOLD and n_entries<=len(xs):
            slope = new_slope
            std_err = new_std_err
            new_slope, intercept, r_value, p_value, new_std_err = stats.linregress(xs[-n_entries:],ys[-n_entries:])
            n_entries+=1
        n_points_considered = n_entries

        if slope is None or std_err>_THRESHOLD:
            n_entries = 5
            slope, std_err = None, 0.0
            new_slope, new_std_err = None, 0.0
            while new_std_err<_THRESHOLD and n_entries<=len(xs):
                slope = new_slope
                std_err = new_std_err
                new_slope, intercept, r_value, p_value, new_std_err = stats.linregress(xs[:n_entries],ys[:n_entries])
                n_entries+=1
                if slope is not None and new_std_err > std_err:
                    break
            n_points_considered = n_entries
        
        if slope is None or std_err>_THRESHOLD:
            n_points_considered = len(xs)
            slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)

    # Return measured_slope, standard_error, number_of_points_considered, successful_fit
    return slope, std_err, n_points_considered, (std_err < _THRESHOLD)

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