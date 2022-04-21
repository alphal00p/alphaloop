try:
    from scipy import stats
except ImportError:
    pass
import math
import os
import sys

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

class NO_colors:
    HEADER = ''
    BLUE = ''
    GREEN = ''
    YELLOW = ''
    WARNING = ''
    FAIL = ''
    ENDC = ''
    BOLD = ''
    UNDERLINE = ''
    PURPLE = ''
    CYAN = ''
    DARKCYAN = ''
    RED = ''
    END = ENDC

# An overwritable and streamable Canvas
# Example usage:
#
#    start_canvas = "A #ln0\nA #ln1\nA #ln2\nA #ln3\nA #ln4\nA #ln5\nA #ln6\nA #ln7\nA #ln8\nA #ln9"
#    canvas = Canvas(start_canvas, stream=True, overwrite=False)
#    i = 0
#    canvas_length = start_canvas.count('\n')+1
#    while True:
#        time.sleep(0.1)
#        canvas.print('\n'.join("A" for _ in range(canvas_length)),0)
#        canvas.print('\n'.join("B\nB\nB".split('\n')[-(canvas_length-i%canvas_length):]),i%canvas_length)
#        i+=1
#        if i==100:
#            break

class Canvas(object):

    cursor_up = lambda lines: '\x1b[{0}A'.format(lines)
    cursor_down = lambda lines: '\x1b[{0}B'.format(lines)

    def __init__(self, canvas_length, stream=True, overwrite=True, print_func=None):
        if isinstance(canvas_length, int):
            self.current_state = ['',]*canvas_length
        else:
            self.current_state = canvas_length.split('\n')

        print('\n'.join(self.current_state)) 
        self.canvas_length = len(self.current_state)
        self.stream = stream
        self.overwrite = overwrite

        if print_func is None:
            self.print_func = (lambda msg: print(msg))
        else:
            self.print_func = print_func

    def print(self, msg, line_number=0, offset=0):
        
        msg_n_lines = msg.count('\n')
        for i_line, line in enumerate(msg.split('\n')):
            if i_line+line_number >= self.canvas_length:
                raise Exception("Number of lines in message exceeds canvas length of %d."%self.canvas_length)
            if self.overwrite:
                prev_length = len(self.current_state[i_line+line_number])
                self.current_state[i_line+line_number] = ' '*offset+line
                self.current_state[i_line+line_number] += ' '*max(prev_length-len(self.current_state[i_line+line_number]),0)
            else:
                prefix = ''.join(list(self.current_state[line_number+i_line])[:offset])
                prefix += ' '*(len(prefix)-offset)
                self.current_state[i_line+line_number] = prefix+line+''.join(list(self.current_state[line_number+i_line])[offset+len(line):])

        if self.stream:
            print(Canvas.cursor_up(self.canvas_length-line_number),end='')
            print('\n'.join(self.current_state[line_number:line_number+msg_n_lines+1]),end='')
            print(Canvas.cursor_down(self.canvas_length-line_number - msg_n_lines ),end='\r')
        else:
            self.print_func('\n'.join(self.current_state))

    def clear(self):
        self.current_state = ['',]*self.canvas_length
        self.print('\n'.join(self.current_state),0)

    def refresh(self):
        if self.stream:
            print(Canvas.cursor_up(self.canvas_length),end='')        
            print('\n'.join(self.current_state),0)

def compute_dod(sequence, threshold=0.01):
    """ from a sequence of results in the format (scaling, complex(result)), determine the asymptotic scaling."""

    _THRESHOLD = threshold

    # First filer zeros for the sequence
    sequence = [s for s in sequence if abs(s[1])>0.]

    # We should not perform a linear regression on less than three points
    if len(sequence)<3:
        return 0.0, 0.0, 0, False

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

class suppress_output: 
    def __init__(self,suppress_stdout=False,suppress_stderr=False,active=None):
        if active is not None:
            suppress_stdout = active
            suppress_stderr = active
        self.suppress_stdout = suppress_stdout 
        self.suppress_stderr = suppress_stderr 
        self._stdout = None 
        self._stderr = None
    def __enter__(self): 
        with open(os.devnull, "w") as devnull:
            if self.suppress_stdout: 
                self._stdout = sys.stdout 
                sys.stdout = devnull        
            if self.suppress_stderr: 
                self._stderr = sys.stderr 
                sys.stderr = devnull 
    def __exit__(self, *args): 
        if self.suppress_stdout: 
            sys.stdout = self._stdout 
        if self.suppress_stderr: 
            sys.stderr = self._stderr
