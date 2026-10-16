"""Classify genes as essential / ambiguous / non-essential from sites2genes output
(Python re-port of tradis_essentiality_wscore.R). Writes QC / changepoint plots as PDF.

Ported verbatim from the CAPTn-seq Snakemake pipeline; only the entrypoint changed, plus an
explicit Agg backend so it runs headless in containers / HPC jobs.
"""

# Faithful re-port of tradis_essentiality_wscore.R
#
# Fixes vs. the previous python port (tradis_essentiality_wscore.py):
#   1. lower/upper threshold indexing was off by one bin (0-based numpy index
#      divided directly by 10000, instead of mapped through the same 1-indexed
#      x-grid R uses) -> thresholds were systematically ~0.0001 too low.
#   2. loess() (local weighted quadratic regression) was replaced by a global
#      degree-3 np.polyfit, which finds a different minimum -> different m.
#   3. scipy's expon.fit / gamma.fit float the `loc` parameter by default;
#      R's fitdistr assumes loc=0 for both. Now fit with floc=0 to match.
#   4. First histogram (finding the second maxima) used bins='fd' where R
#      uses breaks=200 equal-width bins.
#   5. Second histogram (finding the inter-mode minimum m) used a fixed
#      99-bin linspace where R uses r bins of fixed width 1/2000 (r depends
#      on the data via maxval).
#
# ValueError guards (abort_and_output early exits) and the epsilon in the
# log-odds ratio are kept from the previous python port; R has no equivalent
# and would simply error on this pathological input.

import argparse

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def loess_predict(x, y, span=0.75, degree=2):
    """Local weighted polynomial regression (R's loess default: span=0.75,
    degree=2, tricube weights), evaluated at the original x points -
    equivalent to R's predict(loess(y ~ x)) with no newdata."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    k = min(max(int(np.floor(span * n)), degree + 1), n)
    fitted = np.empty(n)
    for i in range(n):
        d = np.abs(x - x[i])
        order = np.argsort(d, kind="mergesort")[:k]
        dmax = d[order].max()
        w = np.ones(k) if dmax == 0 else np.clip(1 - (d[order] / dmax) ** 3, 0, None) ** 3
        dx = x[order] - x[i]
        X = np.vstack([dx ** p for p in range(degree + 1)]).T
        sw = np.sqrt(w)
        beta, *_ = np.linalg.lstsq(X * sw[:, None], y[order] * sw, rcond=None)
        fitted[i] = beta[0]
    return fitted


def run(input_file, output_all, output_essen, output_ambig, output_nonessen):
    STM_baseline = pd.read_csv(input_file, sep="\t", header=0)
    input_base = input_file.split(".")[0]
    ii = STM_baseline['ins_index'].values
    STM_baseline['read_index'] = STM_baseline['read_count'] / STM_baseline['gene_length']

    def abort_and_output(score_E=0, lower=1, upper=99):
        STM_baseline['score_E'] = score_E
        STM_baseline['F_threshold'] = lower
        STM_baseline['N_threshold'] = upper
        STM_baseline.to_csv(output_all, index=False)
        STM_baseline[STM_baseline['ins_index'] < lower].to_csv(output_essen, index=False)
        STM_baseline[STM_baseline['ins_index'] >= upper].to_csv(output_nonessen, index=False)
        if lower < upper:
            STM_baseline[(STM_baseline['ins_index'] >= lower) & (STM_baseline['ins_index'] < upper)].to_csv(output_ambig, index=False)
        else:
            STM_baseline[STM_baseline['ins_index'] >= upper].to_csv(output_ambig, index=False)
        return

    # Identify second maxima (R: hist(ii, breaks=200))
    h, bin_edges = np.histogram(ii, bins=200, density=True)
    if len(h) < 10:
        abort_and_output()
        return
    maxindex = np.argmax(h[10:])
    midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
    maxval = midpoints[maxindex + 3]

    nG = len(STM_baseline['read_count'])
    r = int(np.floor(maxval * 2000))
    if r <= 0:
        abort_and_output()
        return

    fig1 = plt.figure()

    # R: hist(ii[I], breaks=(0:r/2000)) -> r bins of fixed width 1/2000
    I = ii < r / 2000
    h1, edges1 = np.histogram(ii[I], bins=np.linspace(0, r / 2000, r + 1))
    loess_curve = loess_predict(np.arange(1, len(h1) + 1), h1)

    plt.plot(h1, label="Density")
    plt.plot(loess_curve, color='red', lw=2, label="Loess Curve")
    plt.legend()
    plt.title("Density")

    mid1 = (edges1[:-1] + edges1[1:]) / 2
    m = mid1[np.argmin(loess_curve)]
    I1 = (ii < m) & (ii >= 0)

    h_full, edges_full = np.histogram(ii, bins='fd')
    I2 = (ii >= m) & (ii < edges_full[np.max(np.where(h_full > 5))])

    f1 = (np.sum(I1) + np.sum(ii == 0)) / nG
    f2 = np.sum(I2) / nG

    # Fit exponential and gamma distributions, anchored at loc=0 like R's fitdistr
    try:
        d1 = stats.expon.fit(ii[I1], floc=0)
        d2 = stats.gamma.fit(ii[I2], floc=0)
    except ValueError:
        abort_and_output()
        return

    # Generate and save the plot
    x_vals = np.linspace(0, max(ii), 500)
    fig2 = plt.figure()
    plt.hist(ii, bins='fd', density=True, label="Data")
    plt.plot(x_vals, f1 * stats.expon.pdf(x_vals, *d1), label="Exponential fit")
    plt.plot(x_vals, f2 * stats.gamma.pdf(x_vals, *d2), label="Gamma fit")

    # Calculate log-odds ratios (x grid matches R's 1:1000/10000 exactly)
    epsilon = 1e-10  # prevents log(0)/log(inf) on pathological input; no R equivalent
    x_grid = np.arange(1, 1001) / 10000.0
    log_odds = np.log2(
        ((stats.gamma.cdf(x_grid, *d2) *
          (1 - stats.expon.sf(x_grid, *d1))) + epsilon) /
        ((stats.expon.sf(x_grid, *d1) *
          (1 - stats.gamma.cdf(x_grid, *d2))) + epsilon)
    )

    where_lower = np.where(log_odds < -2)[0]
    where_upper = np.where(log_odds > 2)[0]
    if len(where_lower) == 0 or len(where_upper) == 0:
        abort_and_output()
        return
    lower = x_grid[np.max(where_lower)]
    upper = x_grid[np.min(where_upper)]

    # Adding lines to the plot for thresholds
    plt.legend()
    plt.axvline(lower, color='red', label="Lower Threshold")
    plt.axvline(upper, color='red', label="Upper Threshold")
    plt.text(lower, 30, f"Essential changepoint: {lower:.4f}", color='red')
    plt.text(upper, 10, f"Ambiguous changepoint: {upper:.4f}", color='red')
    plt.title("Gamma fits")

    with plt.rc_context():
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(f"{input_base}_QC_and_changepoint_plots.pdf") as pdf:
            pdf.savefig(fig1)
            pdf.savefig(fig2)
    plt.close(fig1)
    plt.close(fig2)

    # Calculate 'read_index' and 'score_E'
    score_E = []
    for iindex in STM_baseline['ins_index']:
        if iindex <= lower:
            score_E.append(iindex / lower)
        elif lower < iindex <= upper:
            score_E.append((iindex - lower) / lower + 1)
        else:
            score_E.append((iindex - upper) / upper + 2)

    abort_and_output(score_E, lower, upper)


def add_arguments(parser):
    """Register this subcommand's options on *parser*."""
    parser.add_argument("--input", required=True, help="sites2genes output (genes_insertions.tsv)")
    parser.add_argument("--output-all", required=True)
    parser.add_argument("--output-essen", required=True)
    parser.add_argument("--output-ambig", required=True)
    parser.add_argument("--output-nonessen", required=True)


def run_args(args):
    """Run the subcommand from an ``argparse.Namespace``."""
    run(
        input_file=args.input,
        output_all=args.output_all,
        output_essen=args.output_essen,
        output_ambig=args.output_ambig,
        output_nonessen=args.output_nonessen,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="captnseq-toolkit essentiality",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_arguments(parser)
    run_args(parser.parse_args(argv))


if __name__ == "__main__":
    main()
