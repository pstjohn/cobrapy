"""
A python wrapper for computing minimal cut sets, based on the berge algorithm
discussed in:
Jungreuthmayer, C., Nair, G., Klamt, S., & Zanghellini, J. (2013). Comparison
and improvement of algorithms for computing minimal cut sets. BMC
Bioinformatics, 14(1), 318. http://doi.org/10.1186/1471-2105-14-318

C code was graciously provided by C. Jungreuthmayer, but not included in this
repository. Contact the original authors for their berge algorithm
implementation. (if this works, I may also reach out for permission to create
and distribute better wrapper)
"""

import os

mcs_lib_dir = os.path.dirname(os.path.abspath(__file__)) + '/mhsCalculator'
mcs_convert_cmd = mcs_lib_dir + '/convert_text_to_bin'
mcs_berge_cmd = mcs_lib_dir + '/mhsCalculator_berge'

# We're going to need pandas here, no real way around it.
import pandas as pd

#TODO: fix me
from cobra.elementary_flux_modes.utils import (
    run_process, make_temp_directory, opt_gen)

def calculate_minimum_cut_sets(good, bad, keep=None, opts=None, verbose=True):
    """ A function to compute minimal cut sets from the input elementary flux
    modes. Will find sets of reaction knockouts which remove all of the flux
    modes in bad while keeping the flux modes in good.

    Parameters:

    good: pandas.dataframe
        A dataframe containing the reactions as columns and elementary flux
        modes as rows. Column labels should be the names of the individual
        reactions.

    bad: pandas.dataframe
        A similar dataframe to `good`, but with the elementary modes which
        should be deleted from the system. Columns must be the same as the good
        dataframe.

    keep: integer or None
        Number of good elementary modes which must be kept by the minimum cut
        set. ($n$ in the original manuscript.) If None, defaults to keeping all
        of the passed 'good' nodes.


    opts: dict or None
        Options for the mhsCalculator script. These likely shouldn't be
        changed, as all options have not necessarily been tested.

        From the mhsCalculator_berge documentation:

        -m ..... filename containing elementary flux modes (in binary form!)
        -r ..... filename containing names of reactions
        -e ..... filename containing essential reactions
        -o ..... filename (output) of file containing computed minimal cutsets
        -n ..... number of good/keeper modes
        -w ..... number of wanted good modes (number of modes that must survive
                                              if cutset is applied)
        -t ..... number of parallel threads
        -k ..... bail out of duplicate mode check if it seems to be ineffective
        -b ..... print cutsets as bitvector string
        -c ..... maximum number of cancelations
        -l ..... use linear approach to find subsets, by default
                     a tree search approach is used
        -s ..... seed value for random number generation
        -h ..... print this help message

    Note:
        I'm assuming that neutral EFMs should be left out of this class, as the
        algorithm doesn't seem to ask them explicitly.


    For details on the approach or implementation, see Jungreuthmayer, C.,
    Nair, G., Klamt, S., & Zanghellini, J. (2013). Comparison and improvement
    of algorithms for computing minimal cut sets. BMC Bioinformatics, 14(1),
    318. http://doi.org/10.1186/1471-2105-14-318

    """

    # Process Input Arguments
    efms = pd.concat((good, bad))
    if not keep: keep = len(good)
    if opts is None: opts = {}
    
    # Set some options based on input arguments
    opts['n'] = len(good)
    opts['w'] = int(keep) # Lets just make sure we're an `int1`

    with make_temp_directory('mhs') as temp_dir:

        # Create relevant temporary files and run mhsCalculator
        write_input_files(efms, temp_dir)
        assert os.path.isfile(temp_dir + '/modes.txt'), \
            "Error creating input files"

        convert_efms_to_binary(temp_dir, len(efms), verbose=verbose)
        assert os.path.isfile(temp_dir + '/modes.bin'), \
            "Error converting EFMs to binary"

        run_mhsCalculator_berge(temp_dir, efms, keep, opts, verbose=verbose)
        assert os.path.isfile(temp_dir + '/cutsets.txt'), \
            "Error calculating minimum cut sets"

        # Parse the bitvector output
        cutsets = read_bitvector_file(temp_dir + '/cutsets.txt', efms.columns)

    return cutsets



        
    
def write_input_files(efms, temp_dir):
    """ Write the elementary flux modes to the format requested by
    mhsCalculator
    
    """
    # Write EFMs to a text file
    efms.to_csv(temp_dir + '/modes.txt', sep='\t', header=False, index=False,
                float_format='%.16g')

    # write reaction names
    with open(temp_dir + '/rfile', 'w') as f:
        f.write(' '.join(('"{}"'.format(r) for r in efms.columns)))
    

def convert_efms_to_binary(temp_dir, number_of_modes, verbose):
    """ Run the 'convert_text_to_bin' program provided by mhsCalculator to
    prepare the input EFMs.

    convert_text_to_bin documentation:
    usage: convert_text_to_bin -m in_modes.txt -o out_modes.bin 
                               -n number_of_modes -r rfile

    """
    opts = {
        'm' : temp_dir + '/modes.txt',
        'o' : temp_dir + '/modes.bin',
        'n' : number_of_modes,
        'r' : temp_dir + '/rfile'
    }

    run_process([mcs_convert_cmd] + list(opt_gen(opts)), verbose=verbose)

def run_mhsCalculator_berge(temp_dir, efms, keep, opts, verbose):
    """ Run the mhsCalculator script
    
    """
    default_opts = {
        'm' : temp_dir + '/modes.bin',
        'r' : temp_dir + '/rfile',
        'o' : temp_dir + '/cutsets.txt',
        't' : 1,
        'b' : None,
    }

    default_opts.update(opts)

    run_process([mcs_berge_cmd] + list(opt_gen(default_opts)), verbose=verbose)


def read_bitvector_file(file_name, rnames):
    """ Read the bitvector output file. Still searching for a more efficient
    implementation -- cython if needed? """

    to_bool = lambda i: i == '1'
    
    outlist = []
    with open(file_name, 'r') as f:
        for line in f:
            outlist += [[to_bool(i) for i in line.strip('\n')]]

    return pd.DataFrame(
        outlist, columns=rnames, dtype=bool,
        index=('CS{}'.format(i) for i in xrange(1, len(outlist) + 1)))

            
    



if __name__ == "__main__":
    from cobra.core import Metabolite, Reaction, Model
    from cobra.elementary_flux_modes import calculate_elementary_modes

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


    out = calculate_elementary_modes(model, verbose=False)

    # Lets try to produce 'P'(R3) without producing 'E'(R4)
    good = out[(out['R3'] > 0) & (out['R4'] <= 0)]
    bad = out[out['R4'] > 0]
    
    cutsets = calculate_minimum_cut_sets(good, bad, 2, verbose=False)
