from six.moves import zip_longest
from six import iterkeys, print_, text_type

import pandas as pd

from .variability import flux_variability_analysis


def format_long_string(string, max_length):
    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string


def metabolite_summary(met, threshold=0.01, fva=False):
    """Print a summary of the reactions which produce and consume this
    metabolite

    threshold: float
    a value below which to ignore reaction fluxes

    fva: float (0->1), or None
    Whether or not to include flux variability analysis in the output.
    If given, fva should be a float between 0 and 1, representing the
    fraction of the optimum objective to be searched.

    """

    def rxn_summary(r):
        return {
            'id': r.id,
            'flux': r.x * r.metabolites[met],
            'reaction': r.reaction,
        }

    flux_summary = pd.DataFrame((rxn_summary(r) for r in met.reactions))
    assert flux_summary.flux.sum() < 1E-6, "Error in flux balance"
    producing = flux_summary[flux_summary.flux > 0].copy()
    consuming = flux_summary[flux_summary.flux < 0].copy()

    for df in [producing, consuming]:
        df['percent'] = df.flux / df.flux.sum()
        df.drop(df[df['percent'] < threshold].index, axis=0,
                inplace=True)

    producing.sort_values('percent', ascending=False, inplace=True)
    consuming.sort_values('percent', ascending=False, inplace=True)

    if not fva:

        producing.flux = producing.flux.apply(
            lambda x: '{:6.2g}'.format(x))
        consuming.flux = consuming.flux.apply(
            lambda x: '{:6.2g}'.format(x))

        flux_len = 6

    else:
        fva_results = pd.DataFrame(
            flux_variability_analysis(met.model, met.reactions,
                                      fraction_of_optimum=fva)).T
        half_span = (fva_results.maximum - fva_results.minimum) / 2
        median = fva_results.minimum + half_span

        producing.flux = producing.id.apply(
            lambda x, median=median, err=half_span:
            u'{0:0.2f} \u00B1 {1:0.2f}'.format(median[x], err[x]))
        consuming.flux = consuming.id.apply(
            lambda x, median=median, err=half_span:
            u'{0:0.2f} \u00B1 {1:0.2f}'.format(median[x], err[x]))

        flux_len = max(producing.flux.apply(len).max(),
                       consuming.flux.apply(len).max()) + 1

    for df in [producing, consuming]:

        df['reaction'] = df['reaction'].map(
            lambda x: format_long_string(x, 52))
        df['id'] = df['id'].map(
            lambda x: format_long_string(x, 8))

    head = "PRODUCING REACTIONS -- " + format_long_string(met.name, 55)
    print_(head)
    print_("-" * len(head))
    print_(("{0:^6} {1:>" + str(flux_len) + "} {2:>8} {3:^54}").format(
        '%', 'FLUX', 'RXN ID', 'REACTION'))

    for row in producing.iterrows():
        print_((u"{0.percent:6.1%} {0.flux:>" + str(flux_len) +
                "} {0.id:>8} {0.reaction:>54}").format(row[1]))

    print_()
    print_("CONSUMING REACTIONS -- " + format_long_string(met.name, 55))
    print_("-" * len(head))
    print_(("{0:^6} {1:>" + str(flux_len) + "} {2:>8} {3:^54}").format(
        '%', 'FLUX', 'RXN ID', 'REACTION'))

    for row in consuming.iterrows():
        print_((u"{0.percent:6.1%} {0.flux:>" + str(flux_len) +
                "} {0.id:>8} {0.reaction:>54}").format(row[1]))


def model_summary(model, threshold=1E-8, fva=None, decimals=2):
    """Print a summary of the input and output fluxes of the model.

    threshold: float
        tolerance for determining if a flux is zero (not printed)
    fva: int or None
        Whether or not to calculate and report flux variability in the
        output summary
    decimals: int
        number of digits after the decimal place to print

    """
    obj_fluxes = pd.Series({'{:<15}'.format(r.id): '{:.3f}'.format(r.x)
                            for r in iterkeys(model.objective)})

    if not fva:

        out_rxns = model.reactions.query(
            lambda rxn: rxn.x > threshold, None
        ).query(lambda x: x, 'boundary')

        in_rxns = model.reactions.query(
            lambda rxn: rxn.x < -threshold, None
        ).query(lambda x: x, 'boundary')

        out_fluxes = pd.Series({r.reactants[0]: r.x for r in out_rxns})
        in_fluxes = pd.Series({r.reactants[0]: r.x for r in in_rxns})

        out_fluxes = out_fluxes.sort_values(ascending=False).round(decimals)
        in_fluxes = in_fluxes.sort_values().round(decimals)

        table = pd.np.array(
            [((a if a else ''), (b if b else ''), (c if c else ''))
             for a, b, c in zip_longest(
                ['IN FLUXES'] + in_fluxes.to_string().split('\n'),
                ['OUT FLUXES'] + out_fluxes.to_string().split('\n'),
                ['OBJECTIVES'] + obj_fluxes.to_string().split('\n'))])

    else:
        fva_results = pd.DataFrame(
            flux_variability_analysis(model, fraction_of_optimum=fva)).T

        half_span = (fva_results.maximum - fva_results.minimum) / 2
        median = fva_results.minimum + half_span

        out_rxns = model.reactions.query(
            lambda rxn: median.loc[rxn.id] > threshold, None
        ).query(lambda x: x, 'boundary')

        in_rxns = model.reactions.query(
            lambda rxn: median.loc[rxn.id] < -threshold, None
        ).query(lambda x: x, 'boundary')

        out_fluxes = pd.DataFrame(
            {r.reactants[0]: {'x': median.loc[r.id],
                              'err': half_span.loc[r.id]}
             for r in out_rxns}).T

        in_fluxes = pd.DataFrame(
            {r.reactants[0]: {'x': median.loc[r.id],
                              'err': half_span.loc[r.id]}
             for r in in_rxns}).T

        out_fluxes = out_fluxes.sort_values(
            by='x', ascending=False).round(decimals)
        in_fluxes = in_fluxes.sort_values(by='x').round(decimals)

        in_fluxes_s = in_fluxes.apply(
            lambda x: u'{0:0.2f} \u00B1 {1:0.2f}'.format(x.x, x.err),
            axis=1)
        out_fluxes_s = out_fluxes.apply(
            lambda x: u'{0:0.2f} \u00B1 {1:0.2f}'.format(x.x, x.err),
            axis=1)
        out_fluxes_s = out_fluxes.apply(lambda x: text_type(x.x) +
                                        u" \u00B1 " + text_type(x.err), axis=1)

        table = pd.np.array(
            [((a if a else ''), (b if b else ''), (c if c else ''))
             for a, b, c in zip_longest(
                ['IN FLUXES'] + in_fluxes_s.to_string().split('\n'),
                ['OUT FLUXES'] + out_fluxes_s.to_string().split('\n'),
                ['OBJECTIVES'] + obj_fluxes.to_string().split('\n'))])

    print_(u'\n'.join([u"{a:<30}{b:<30}{c:<20}".format(a=a, b=b, c=c) for
                       a, b, c in table]))
