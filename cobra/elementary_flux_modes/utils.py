import contextlib
import tempfile
import shutil
import subprocess
import os
import sys

from warnings import warn

@contextlib.contextmanager
def make_temp_directory(prefix=None):
    """ Create a temporary working directory for efmtool input and output files

    """
    temp_dir = tempfile.mkdtemp(prefix=prefix + '_tmp', dir='.')
    try:
        yield temp_dir
    except Exception as e:
        warn(str(e))
    finally:
        shutil.rmtree(temp_dir)


def run_process(process, verbose=True):
    """ Run a bash process in python, printing lines from STDOUT to the python
    shell as the become available

    """

    if verbose:
        stdout = subprocess.PIPE 
        stderr = subprocess.PIPE 
    else:
        stdout = open(os.devnull, 'w')
        stderr = subprocess.PIPE 

    process = subprocess.Popen(process, stdout=stdout, stderr=stderr, bufsize=1)
    if verbose:
        for line in iter(process.stdout.readline, b''):
            sys.stdout.write(line)
            sys.stdout.flush()
        process.stdout.close()
    process.wait()

def opt_gen(opts_dict):
    """ Iterate over a dictionary, returning a generator for subprocess.Popen.
    Assumes the dictionary is of the form {option : value}, where the desired
    list is [..., '-option', 'value', ...].

    Call as list(opt_gen(opts_dict)) to generate the list.

    """
    for opt, val in opts_dict.iteritems():
        yield '-' + opt
        if val is not None: yield str(val)

