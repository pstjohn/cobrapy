"""
A python wrapper around EFMTool, a java-based library for calculating
elementary flux modes. See http://www.csb.ethz.ch/tools/software/efmtool.html;
Terzer M, Stelling J (2008) Large-scale computation of elementary flux modes
with bit pattern trees. Bioinformatics 24: 2229-2235. 
"""

import os
import contextlib
import tempfile
import shutil
import subprocess

efm_lib_dir = os.path.dirname(os.path.abspath(__file__)) + '/efmtool'
efm_command =  ['java', '-jar', efm_lib_dir + '/metabolic-efm-all.jar']

import numpy as np

try:
    import pandas as pd
    pandas = True
except ImportError:
    pandas = False


from cobra import ArrayBasedModel

@contextlib.contextmanager
def make_temp_directory():
    """ Create a temporary working directory for efmtool input and output files

    """
    temp_dir = tempfile.mkdtemp(prefix='efmtool_tmp', dir='.')
    yield temp_dir
    shutil.rmtree(temp_dir)


def create_model_files(cobra_model, temp_dir, deepcopy_model=True):
    """ Write stochiometry data, metabolite, and reaction names to temporary
    files in preparation for calling efmtool

    """
    
    try: stoich_mat = cobra_model.S.toarray()
    except AttributeError:
        cobra_model = ArrayBasedModel(
            cobra_model, deepcopy_model=deepcopy_model)
        stoich_mat = cobra_model.S.toarray()

    np.savetxt(temp_dir + '/stoich.txt', stoich_mat, delimiter='\t')

    np.savetxt(temp_dir + '/revs.txt', 
               np.array([r.reversibility for r in cobra_model.reactions]),
               delimiter='\t', fmt='%d', newline='\t')

    with open(temp_dir + '/rnames.txt', 'w') as f:
        f.write('\t'.join(('"{}"'.format(r.id) for r in cobra_model.reactions)))

    with open(temp_dir + '/mnames.txt', 'w') as f:
        f.write('\t'.join(('"{}"'.format(m.id) for m in cobra_model.metabolites)))



def calculate_elementary_modes(cobra_model, opts=None):
    """ Run the java efmtool on the given cobra_model. Opts is a dictionary
    which overwrites the default options """

    if opts is None: opts = {}

    with make_temp_directory() as temp_dir:

        create_model_files(cobra_model, temp_dir)

        try: out_file = opts.pop('out_file')
        except KeyError: out_file = temp_dir + '/out.bin'

        default_opts = {
            'kind'             : 'stoichiometry',
            'stoich'           : temp_dir + '/stoich.txt',
            'rev'              : temp_dir + '/revs.txt',
            'reac'             : temp_dir + '/rnames.txt',
            'meta'             : temp_dir + '/mnames.txt',
            'arithmetic'       : 'double',
            'zero'             : 1E-10,
            'out'              : 'binary-doubles', # TODO support binary output
            'compression'      : 'default',
            'log'              : 'console',
            'level'            : 'INFO',
            'maxthreads'       : -1,
            'normalize'        : 'min',
            'adjacency-method' : 'pattern-tree-minzero',
            'rowordering'      : 'MostZerosOrAbsLexMin',
        }

        default_opts.update(opts)


        # Create a list of arguments to pass to the python subprocess module
        def opt_gen():
            for opt, val in default_opts.iteritems():
                yield '-' + opt
                yield str(val)

                # This kw takes two arguments
                if opt == 'out': yield out_file

        # Run the EFMtool, outputting STDOUT to python.
        process = subprocess.Popen(efm_command + list(opt_gen()),
                                   stdout=subprocess.PIPE)
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print output.strip()

        if   'binary-doubles' in default_opts['out']: dtype='>d'
        elif 'binary-boolean' in default_opts['out']:
            raise RuntimeError('binary-boolean not supported')
        with open(out_file, 'rb') as f:
            out_arr = np.fromstring(f.read()[13:], dtype=dtype).reshape(
                (-1, len(cobra_model.reactions))).T
            out_arr = np.array(out_arr, dtype=np.float64)

        if pandas:
            out_arr = pd.DataFrame(
                out_arr, index=(r.id for r in cobra_model.reactions), 
                columns=('EM{}'.format(i) for i in 
                         xrange(1, out_arr.shape[1] + 1)))

        return out_arr.T # I prefer to have model reactions as columns.



if __name__ == "__main__":

    # TODO: These should get moved to the test suite

    from cobra.core import Metabolite, Reaction, Model

    model = Model('simple_model')

    A = Metabolite('A')
    B = Metabolite('B')
    C = Metabolite('C')
    D = Metabolite('D')
    E = Metabolite('E')
    P = Metabolite('P')

    R1 = Reaction('R1')
    R2 = Reaction('R2')
    R3 = Reaction('R3')
    R4 = Reaction('R4')
    R5 = Reaction('R5')
    R6 = Reaction('R6')
    R7 = Reaction('R7')
    R8 = Reaction('R8')
    R9 = Reaction('R9')
    R10 = Reaction('R10')

    model.add_metabolites([A, B, C, D, E, P])
    model.add_reactions([R1, R2, R3, R4, R5, R6, R7, R8, R9, R10])

    model.reactions.R1.build_reaction_from_string('--> A')
    model.reactions.R2.build_reaction_from_string('<--> B')
    model.reactions.R3.build_reaction_from_string('P -->')
    model.reactions.R4.build_reaction_from_string('E -->')
    model.reactions.R5.build_reaction_from_string('A --> B')
    model.reactions.R6.build_reaction_from_string('A --> C')
    model.reactions.R7.build_reaction_from_string('A --> D')
    model.reactions.R8.build_reaction_from_string('B <--> C')
    model.reactions.R9.build_reaction_from_string('B --> P')
    model.reactions.R10.build_reaction_from_string('C + D --> E + P')


    calculate_elementary_modes(model, {'out' : 'binary-boolean'})




