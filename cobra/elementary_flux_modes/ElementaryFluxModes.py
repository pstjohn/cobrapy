"""
A python wrapper around EFMTool, a java-based library for calculating
elementary flux modes. See http://www.csb.ethz.ch/tools/software/efmtool.html;
Terzer M, Stelling J (2008) Large-scale computation of elementary flux modes
with bit pattern trees. Bioinformatics 24: 2229-2235. 
"""

import os

efm_lib_dir = os.path.dirname(os.path.abspath(__file__)) + '/efmtool'
efm_command =  ['java', '-jar', efm_lib_dir + '/metabolic-efm-all.jar']

import numpy as np

try:
    import pandas as pd
    pandas = True
except ImportError:
    pandas = False


from ..core import ArrayBasedModel
from .utils import make_temp_directory, run_process


def create_model_files(cobra_model, temp_dir, deepcopy_model=True):
    """ Write stochiometry data, reaction reversibilities, metabolite, and
    reaction names to temporary files in preparation for calling efmtool

    """
    
    try: stoich_mat = cobra_model.S.toarray()
    except AttributeError:
        cobra_model = ArrayBasedModel(
            cobra_model, deepcopy_model=deepcopy_model)
        stoich_mat = cobra_model.S.toarray()

    # Stoichiometric Matrix
    np.savetxt(temp_dir + '/stoich.txt', stoich_mat, delimiter='\t')

    # Reaction reversibilities
    np.savetxt(temp_dir + '/revs.txt', 
               np.array([r.reversibility for r in cobra_model.reactions]),
               delimiter='\t', fmt='%d', newline='\t')

    # Reaction Names
    with open(temp_dir + '/rnames.txt', 'w') as f:
        f.write('\t'.join(('"{}"'.format(r.id) for r in cobra_model.reactions)))

    # Metabolite Names
    with open(temp_dir + '/mnames.txt', 'w') as f:
        f.write('\t'.join(('"{}"'.format(m.id) for m in cobra_model.metabolites)))



def read_double_out(cobra_model, out_file):
    """ Read the output file generated from EMFTool. Returns a numpy array or
    pandas dataframe (if pandas can be loaded)
    
    """
    with open(out_file, 'rb') as f:
        out_arr = np.fromstring(f.read()[13:], dtype='>d').reshape(
            (-1, len(cobra_model.reactions))).T
        out_arr = np.array(out_arr, dtype=np.float64)

    if pandas:
        out_arr = pd.DataFrame(
            out_arr, index=(r.id for r in cobra_model.reactions), 
            columns=('EM{}'.format(i) for i in 
                     xrange(1, out_arr.shape[1] + 1)))

    return out_arr

def calculate_elementary_modes(cobra_model, opts=None, verbose=True):
    """ A python wrapper around EFMTool, a java-based library for calculating
    elementary flux modes.
    
    Parameters:
    ===========

    cobra_model: cobra.Model or similar
        The model for which to calculate elementary flux modes. Reaction
        reversibilities are inferred from their corresponding bounds.

    opts: dict or None
        Additional dictionary options to pass to the emftool solver. See Option
        Descriptions below for additional information on arguments.

    verbose: bool
        Whether or not to redirect stdout to the python shell.


    Option Descriptions from the MATLAB Documentation:
    ==================================================

    The following options are supported (any others are silently ignored):
       'arithmetic'     Which arithmetic to use, one of
                          'double'     Double precision floating point (default)
                          'fractional' Fraction numbers, fixed or infinite
                                       precision
       'precision'      Bits to use for fraction numerator/denominator,
                        -1 for infinite precision, or 64, 128, ... bits for
                        fixed precision (only for arithmetic=fractional)
       'zero'           Values to be treated as zero
                        defaults are 1e-10 for double arithmetic and 0 for
                        fractional arithmetic
       'compression'    Which compression techniques should be applied to the
                        network before efm computation, one of
                          'default'                    Default compression
                          'default-no-combine'         No reaction merging
                          'default-no-couple-combine'  Reduced react. merging
                          'default-no-duplicate'       No duplicate gene com.
                          'default-all-duplicate'      Extended duplicate 
                                                       gene compression
                          'off'                        No compression at all
                          'unique'                     only unique flows
                                                       compression
                          'nullspace'                  Only nullspace based
                                                       compression
                          'unique-no-recursion'        Only one round of 
                                                       unique flows compr.
                          'nullspace-no-recursion'     Only one round of 
                                                       nullspace compression

       'level'          The log level, default is INFO, supported:
                          'WARNING', 'INFO', 'FINE', 'FINER', 'FINEST'

       'maxthreads'     Maximum threads to use, default is 1 for single core
                        systems and k for multi core systems with k cores.

       'normalize'      Normalization method of the output efms, one of
                          'min'     default, minimum absolute value is 1
                          'max'     maximum absolute value is 1
                          'norm2'   norm 2, that is, vector length is 1
                          'squared' like norm 2, but all values are squared, 
                                    possibly negative. To get the original 
                                    value, take the square root of the 
                                    absolute value and keep the sign, i.e.
                                    val = sign(valsq) * sqrt(abs(valsq)).
                                    This method allows error free values
                                    if fraction number arithmetic is used.
                          'none'    no normalization

       'count-only'     If true, the (binary) modes are computed and the 
                          number of efms is returned, but not the efm vectors
       'sign-only'      If true, the (binary) modes are computed, but only 
                          sign values of reaction fluxes are returned as EFMs.
                          The values +/-1 stand for forward/reverse flux, and 0
                          for no flux, respectively. Such sign-valued EFMs can
                          be converted back into double valued EFMs by using
                          the SignToDouble function.
       'parse-only'     If true, input data is processed and parsed, but no
                          elementary modes are computed. The returned 
                          structure contains no efms field
       'convert-only'   If true, input data is processed and converted into
                          the program's internal format. The efm computation
                          itself is not started, but call options (for manual
                          invocation of the Java program) are reported as
                          return string. The converted files are stored in 
                          the tmp directory.
                          Depending on the input format, these files are:
                              - stoich.txt    the stoichiometric matrix
                              - revs.txt      the reaction reversibilities
                              - mnames.txt    the metabolite names
                              - rnames.txt    the reaction names
                              - rlist.txt     the reaction list file
                              - xbml.xml      the sbml input file
                                    
       'suppress'       List of reaction names to suppress, separated by 
                          whitespace. All resulting flux modes will have a 
                          zero flux value for suppressed reactions

       'enforce'       List of reaction names to enforce, separated by 
                          whitespace. All resulting flux modes will have a 
                          non-zero flux value for enforced reactions

       'adjacency-method' Many different methods, all using pattern trees, 
                          the fastest methods are
                          'pattern-tree-minzero'   default, combinatorial
                                                   test with pattern trees
                          'pattern-tree-rank'      standard rank test with
                                                   doubles precision
                          'pattern-tree-mod-rank'  rank test with residue
                                                   arithmetic
                          'rankup-modpi-incore'    rank updating with residue 
                                                   int32 arithmetic, prime
                                                   p <= sqrt((2^31-1)/2),
                                                   pattern tree is stored in
                                                   memory (i.e. in-core)
                          'rankup-modpi-outcore'   rank updating with residue 
                                                   int32 arithmetic, prime
                                                   p <= sqrt((2^31-1)/2),
                                                   pattern tree is stored out
                                                   of the main memory (i.e. 
                                                   out-of-core)
                          'pattern-tree-rank-update-modpi' rank updating with
                                                           residue int32
                                                           arithmetic, prime
                                                           p <= sqrt((2^31-1)/2)
                          'pattern-tree-rank-update-modp'  rank updating with
                                                           residue int64
                                                           arithmetic, prime
                                                           p <= 2^31-1
                          'pattern-tree-rank-update-frac'  rank updating with
                                                           fraction numbers
                          'pattern-tree-rank-update-frac2' rank updating with
                                                           fraction numbers
                                                           and copying

       'rowordering'    Nullspace row ordering, affecting growth of
                        intermediary modes during the iteration phase.
                        The fastest row orderings are
                          'MostZerosOrAbsLexMin'    default, most zeros in a 
                                                    row, or absoulte 
                                                    lexicographical if equal
                          'MostZeros'               most zeros in a row
                          'AbsLexMin'               absoulte lexicographical
                          'LexMin'                  lexicographical
                          'FewestNegPos'            lowest product of 
                                                    negative/positive counts
                          'MostZerosOrFewestNegPos' the two orderings, the
                                                    second ordering is used 
                                                    if sorting equal with 1st
                          'MostZerosOrLexMin'       the two orderings, the
                                                    second ordering is used 
                                                    if sorting equal with 1st
                          'FewestNegPosOrMostZeros' the two orderings, the
                                                    second ordering is used 
                                                    if sorting equal with 1st
                          'FewestNegPosOrAbsLexMin' the two orderings, the
                                                    second ordering is used 
                                                    if sorting equal with 1st
                          'Random'                  random sorting

       'model'          Algorithm variant for EFM computation. Default model
                        is a nullspace-model, which uses the binary nullspace
                        approach (J. Gagneur & S. Klamt / C. Wagner).
                        Supported algorithm model variants are:
                          'nullspace'              default, binary nullspace
                                                   approach
                          'canonical'              canonical approach

       'memory'         Memory type to use, which is usually in-core (main
                        memory, RAM) or out-of-core (store intermediary 
                        results on disk).
                        Supported memory variants are:
                          'in-core'                default, keep all modes in
                                                   main memory
                          'out-core'               store intermediary results
                                                   in files on disk, in the
                                                   temp directory (see 
                                                   'tmpdir' option). For 
                                                   extreme cases, the bit 
                                                   pattern trees should also
                                                   stored on disk, which is
                                                   currently supported by the
                                                   'adjacency-method' option
                                                   'rankup-modpi-outcore'
                          'sort-out-core'          like 'out-core', but also
                                                   sorting to create the 
                                                   pattern trees is performed
                                                   out of core memory

       'impl'           Algorithm implementation for EFM computation. Default
                        is a sequential double description implementation.
                        Supported algorithm implementations:
                          'SequentialDoubleDescriptionImpl'	default

       'tmpdir'         Directory for temporary files, e.g. data files for
                        intermediary results if out core implementation is
                        used (see memory option). Fast, non-server 
                        directories with large storage capaci'models/asuc_v2_sugars.json')ty are 
                        preferable. Default is tmp subdirectory in the 
                        installation path

    For more details on the algorithm or implementation, see:
    http://www.csb.ethz.ch/tools/software/efmtool.html; Terzer M, Stelling J
    (2008) Large-scale computation of elementary flux modes with bit pattern
    trees. Bioinformatics 24: 2229-2235. 
    
    """

    if opts is None: opts = {}

    with make_temp_directory('efmtool') as temp_dir:

        create_model_files(cobra_model, temp_dir)

        try: out_file = opts.pop('out_file')
        except KeyError: out_file = temp_dir + '/out.bin'

        # Default options for the EMFtool. I don't recommend changing any of
        # these, even though the function is set up for this. This tool hasn't
        # been tested for any options other than these
        default_opts = {
            'kind'             : 'stoichiometry',
            'stoich'           : temp_dir + '/stoich.txt',
            'rev'              : temp_dir + '/revs.txt',
            'reac'             : temp_dir + '/rnames.txt',
            'meta'             : temp_dir + '/mnames.txt',
            'arithmetic'       : 'double',
            'zero'             : 1E-10,
            'out'              : 'binary-doubles', # TODO support binary output?
            'compression'      : 'default',
            'log'              : 'console',
            'level'            : 'INFO',
            'maxthreads'       : -1,
            'normalize'        : 'max',
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
        run_process(efm_command + list(opt_gen()), verbose=verbose)

        if 'binary-doubles' in default_opts['out']:
            out_arr = read_double_out(cobra_model, out_file)

        # I prefer to have model reactions as columns.
        return out_arr.T

def boundary_efms(cobra_model, efms, tol=1E-4):
    """Shortcut method to condense elementary flux modes to a unique set of
    input-output modes.

    cobra_model: a cobra.Model object

    efms: a pandas.DataFrame 
        represents the elementary flux modes of the system (output of
        calculate_elementary_modes)

    tol: float
        Tolerance for combining similar elementary flux modes

    """

    def scale_df(df):
        df_ss = np.sqrt((df**2).sum(1))
        to_drop = df_ss[df_ss == 0].index
        df_t = df.drop(to_drop)
        df_ss = df_ss.drop(to_drop)
        return df_t.divide(df_ss.values, axis='rows'), to_drop

    def remove_duplicates(df, tol=1E-4):
        df_round = df.fillna(0)
        df_round, to_drop = scale_df(df_round)
        df = df.drop(to_drop)
        df_round = np.round(df_round, decimals=int(-np.log10(1E-4)))
        return df.ix[~(df_round*(1/tol)).astype('int64').duplicated()]
    
    # Get the boundary reactions from the cobra.Model
    boundary_rxns = [r.id for r in cobra_model.reactions.query(
        'system_boundary', 'boundary') if r.id in efms.columns]

    ioefms = efms[boundary_rxns]

    return remove_duplicates(ioefms)








if __name__ == "__main__":

    # TODO: These should get moved to the test suite.
    # This just builds a test model (from
    # http://www.csb.ethz.ch/tools/software/efmtool/documentation.html) and
    # calculates the elementary flux modes.

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


    out = calculate_elementary_modes(model, verbose=False)

    print out




