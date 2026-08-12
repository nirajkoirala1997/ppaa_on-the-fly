#!/usr/bin/env python3
"""Merge several independent (seed) runs into one report.

    tools/merge_runs.py ../summary/run_00/LO_*.dat -o ../summary/run_00/LO.dat
    tools/merge_runs.py ../summary/run_00/LO_*.dat          # -> LO.dat beside them

THE SAMPLE UNIT IS A VEGAS ITERATION, NOT AN EVENT.  This is the one thing that
does not carry over from the KALEU tree's merge, and getting it wrong would
produce a plausible number with a wrong error bar.

That code draws independent points, so it merges per event: sum_w over sum of
npoints.  VEGAS cannot be merged that way.  With MDS = 1 the points inside one
iteration are STRATIFIED -- deliberately anti-correlated -- so their spread is
not a variance of anything.  What is independent is one iteration from the next.
So each report stores, per bin,

    hw  = sum_i u_i        over n = niter iterations, u_i the bin's
    hww = sum_i u_i^2      contribution from iteration i

where u_i already carries the 1/ITMX from the vwgt the integrand was called
with, so hw is itself the mean.  The independent per-iteration estimates of the
bin integral are therefore v_i = n * u_i, and merging over runs r is a pooled
mean over ALL the iterations of ALL the runs:

    N    = sum_r n_r
    S1   = sum_r n_r   * hw_r        =  sum over every v_i
    S2   = sum_r n_r^2 * hww_r       =  sum over every v_i^2
    mean = S1 / N
    err  = sqrt( (S2 - S1^2/N) / (N*(N-1)) )        [per bin: / binwidth]

Applied to a single run this reproduces that run's own numbers exactly, which is
the identity test worth running after any change here.  Unlike averaging per-run
means it stays correct when the runs contribute unequal numbers of iterations.

WHAT THE MERGE REFUSES.  Pooling the v_i treats them as identically
distributed, and a v_i's variance depends on how many points went into its
iteration.  So runs with different npoints_per_iter are NOT mergeable, however
tempting it looks -- that is a hard error, not a warning.  Same for a different
ecm, PDF, cut, scale, or physics card: those are different quantities.  A
crashed or truncated seed, by contrast, is merely skipped, because one dead job
must not cost the campaign.

WHAT MERGING CANNOT FIX.  Summing accumulators is exact but blind.  A seed whose
VEGAS grid adapted badly still contributes its full weight, drags the mean and
SHRINKS the pooled error -- a confidently quoted wrong number.  The per-seed
table and chi2/dof printed at the end are what expose that, and they are the
reason to run a campaign rather than one long job.
"""

import argparse
import math
import os
import re
import sys
from collections import OrderedDict

# Blocks that carry data rows. RUN and INPUT CARD are comment-only and are
# handled separately: the first supplies the metadata, the second is echoed.
DATA_KINDS = ("TOTAL", "HISTOGRAM2D", "HISTOGRAM")

# Header keys that must agree across every input, with the value assumed for a
# report written before that key existed. Two runs differing in any of these
# describe different quantities and must never be pooled.
CONSISTENCY = OrderedDict([
    ("mode",             None),
    ("generator",        ""),
    ("ecm",              None),
    ("pdf",              None),
    ("amplitude",        ""),
    ("cuts",             "none"),
    ("scales",           None),
    ("nvar",             "1"),
    ("npoints_per_iter", None),
])


class Block(object):
    """A TOTAL or HISTOGRAM section of one report."""

    def __init__(self, kind, contrib, name):
        self.kind = kind
        self.contrib = contrib          # 'LO'
        self.name = name                # 'mee1', or '' for a TOTAL
        self.header = []                # comment lines, kept from the first file
        self.cols = []                  # from '# columns:', units stripped
        self.colline = ""               # ... and the same line verbatim
        self.meta = {}                  # '# key = value' inside the block
        self.rows = []

    @property
    def niter(self):
        return int(self.meta["niter"])

    def col(self, name, path):
        try:
            return self.cols.index(name)
        except ValueError:
            raise ValueError("%s: block '%s %s' has no %r column; it carries: %s"
                             % (path, self.kind, self.name, name,
                                " ".join(self.cols) or "(no '# columns:' line)"))


def parse_report(path):
    """Return (meta, blocks, card_lines).

    Deliberately mirrors tools/extract_data.py: block starts are recognised only
    for the data kinds, so the RUN block's 'key = value' lines land in the
    top-level meta, and the INPUT CARD block's '# [run] seed = ...' lines match
    nothing and are collected verbatim for the echo.
    """
    meta, blocks, card = {}, OrderedDict(), []
    cur, in_card = None, False

    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue

            if re.match(r"#\s*====\s*END", s):
                cur, in_card = None, False
                continue

            if re.match(r"#\s*====\s*INPUT CARD\s*====", s):
                in_card = True
                continue

            m = re.match(r"#\s*====\s*(%s)\s*(.*?)\s*====" % "|".join(DATA_KINDS), s)
            if m:
                kind, rest = m.group(1), m.group(2).strip()
                parts = rest.split()
                # 'HISTOGRAM LO mee1' -> contrib 'LO', name 'mee1'
                # 'TOTAL LO'          -> contrib 'LO', name ''
                name = parts[-1] if kind != "TOTAL" and len(parts) > 1 else ""
                contrib = " ".join(parts[:-1]) if name else rest
                cur = "%s %s" % (kind, name)
                blocks[cur] = Block(kind, contrib, name)
                continue

            if in_card:
                card.append(line.rstrip("\n"))
                continue

            if s.startswith("#"):
                mc = re.match(r"#\s*columns:\s*(.*)$", s)
                if mc and cur:
                    blocks[cur].colline = mc.group(1)
                    blocks[cur].cols = [re.sub(r"\[.*?\]", "", c)
                                        for c in mc.group(1).split()]
                    continue
                mm = re.match(r"#\s*(\w+)\s*=\s*(.*)$", s)
                if mm:
                    (blocks[cur].meta if cur else meta)[mm.group(1)] = mm.group(2).strip()
                    continue
                if cur:
                    blocks[cur].header.append(line.rstrip("\n"))
                continue

            if cur is None:
                raise ValueError("%s: data row outside any block: %r" % (path, s))
            blocks[cur].rows.append([float(x) for x in s.split()])

    return meta, blocks, card


def n_fixed_columns(kind):
    """Leading columns of a histogram row that describe the bin, not the value."""
    return 6 if kind == "HISTOGRAM2D" else 3


def accumulator_index(blk, iv, nvar, path):
    """Column of (sum_w, sum_w2) for variation `iv` in a histogram row.

    For the single-variation report this build writes today the answer is read
    by NAME, so a column added ahead of the accumulators cannot silently shift
    what gets summed.  Beyond that it follows the layout the report format
    documents -- fixed bin columns, then nvar (value, error) pairs, then nvar
    (sum_w, sum_w2) pairs -- which is what tools/extract_data.py also assumes.
    """
    if nvar == 1:
        i = blk.col("sum_w", path)
        if blk.col("sum_w2", path) != i + 1:
            raise ValueError("%s: sum_w2 does not follow sum_w in '%s %s'"
                             % (path, blk.kind, blk.name))
        return i
    return n_fixed_columns(blk.kind) + 2 * nvar + 2 * iv


def bin_width(kind, row):
    if kind == "HISTOGRAM2D":
        return (row[1] - row[0]) * (row[4] - row[3])
    return row[1] - row[0]


def pooled(S1, S2, N):
    """(mean, error) of the pooled per-iteration estimates. See the module doc."""
    mean = S1 / N
    if N < 2:
        return mean, 0.0
    return mean, math.sqrt(max(0.0, S2 - S1 * S1 / N) / (N * (N - 1.0)))


def repack(S1, S2, N):
    """(hw, hww) that a reader of the merged file recovers (mean, err) from.

    The file format's contract is that a row's sum_w/sum_w2 are per-iteration
    accumulators over the block's own niter, and that mean and error follow from
    them by histo1d_mean/histo1d_error.  A merged row must satisfy the same
    contract with niter = N, or the merged file would be readable only by this
    script and could not itself be merged again.  Inverting those two formulas
    for (hw, hww) at n = N gives what is below; for a single input it returns
    that input's own hw and hww unchanged.
    """
    hw = S1 / N
    hww = hw * hw / N + (S2 - S1 * S1 / N) / (N * N)
    return hw, hww


def weighted_mean(vals, errs):
    """Inverse-variance combination, and the chi2/dof of the inputs about it.

    Used for the VEGAS column, which cannot be pooled from accumulators: VEGAS
    quotes a variance-WEIGHTED average over its own iterations, and each seed
    adapted its own grid, so the seeds are independent measurements of the same
    number with individually meaningful errors -- exactly the case this is for.
    """
    pairs = [(v, e) for v, e in zip(vals, errs) if e > 0.0]
    if not pairs:
        return (vals[0] if vals else 0.0), 0.0, 0.0
    sw = sum(1.0 / (e * e) for _, e in pairs)
    mean = sum(v / (e * e) for v, e in pairs) / sw
    err = math.sqrt(1.0 / sw)
    ndof = len(pairs) - 1
    chi2 = (sum(((v - mean) / e) ** 2 for v, e in pairs) / ndof) if ndof >= 1 else 0.0
    return mean, err, chi2


def card_signature(card_lines):
    """The physics half of the echoed card.

    [run] and [output] are dropped on purpose: the seed differs between inputs
    by design, n_adapt/it_adapt only shape a grid that is thrown away, and where
    a job wrote its report says nothing about what it computed.  Everything else
    -- masses, ew, cuts, the histogram definitions -- must match, and the report
    header alone does not cover it: two runs at different mZ have identical
    headers.
    """
    sig = []
    for line in card_lines:
        m = re.match(r"#\s*(\[\w+\])\s*(.*)$", line.strip())
        if not m:
            continue
        if m.group(1) in ("[run]", "[output]"):
            continue
        sig.append("%s %s" % (m.group(1), m.group(2).strip()))
    return sig


def read_inputs(paths):
    """(usable, skipped): parse every input and drop the ones that cannot be used.

    Two failure classes, on purpose:
      unusable  (missing, truncated, unparseable) -> warn, skip, carry on.
      wrong     (different physics, or an already-merged input) -> hard error.
                Combining those would be wrong, and quietly dropping one would
                be just as bad.
    """
    usable, skipped = [], []

    for p in paths:
        try:
            meta, blocks, card = parse_report(p)
        except (IOError, OSError, ValueError) as exc:
            skipped.append((p, "unreadable: %s" % exc))
            continue

        # An already-merged report has its iterations counted once already;
        # folding it back in would inflate N and shrink the error bar.
        if "nruns" in meta or str(meta.get("seed", "")).startswith("merged"):
            sys.exit("ERROR: %s is itself a merged report -- merging it again would\n"
                     "       double-count its iterations. Glob the seed-tagged files\n"
                     "       (LO_*.dat); the merged LO.dat has no _<seed> suffix and\n"
                     "       is deliberately not matched by them." % p)

        data = [b for b in blocks.values() if b.rows]
        if not data:
            skipped.append((p, "no data blocks (the run may have died early)"))
            continue
        if not meta.get("npoints_per_iter"):
            skipped.append((p, "no npoints_per_iter in header (truncated?)"))
            continue
        missing = [k for k, b in blocks.items() if "niter" not in b.meta]
        if missing:
            skipped.append((p, "block %r has no '# niter' line" % missing[0]))
            continue

        usable.append((p, meta, blocks, card))

    return usable, skipped


def check_consistent(usable):
    """Hard-error on any input that describes a different quantity."""
    ref_path, ref_meta, ref_blocks, ref_card = usable[0]

    for key, dflt in CONSISTENCY.items():
        ref = ref_meta.get(key, dflt)
        for p, m, _, _ in usable[1:]:
            got = m.get(key, dflt)
            if got != ref:
                extra = ""
                if key == "npoints_per_iter":
                    extra = ("\n       Points per iteration IS part of the sample unit's"
                             "\n       variance, so these cannot be pooled at all.")
                sys.exit("ERROR: %s differs between %s (%r) and %s (%r)%s"
                         % (key, ref_path, ref, p, got, extra))

    ref_sig = card_signature(ref_card)
    for p, _, _, card in usable[1:]:
        sig = card_signature(card)
        if sig != ref_sig:
            diff = ([x for x in sig if x not in ref_sig] +
                    [x for x in ref_sig if x not in sig])[:4]
            sys.exit("ERROR: the physics section of the run card differs between\n"
                     "       %s and %s, e.g. %s"
                     % (ref_path, p, "; ".join(diff) or "(ordering)"))

    seeds = [m.get("seed") for _, m, _, _ in usable]
    dupes = sorted(set(s for s in seeds if seeds.count(s) > 1))
    if dupes:
        sys.exit("ERROR: duplicate seed(s) among the inputs: %s\n"
                 "       Merging the same run twice would understate the error."
                 % ", ".join(str(d) for d in dupes))


def check_shapes(usable):
    """Drop inputs whose blocks or binning do not line up with the first."""
    ref_path, _, ref_blocks, _ = usable[0]
    keep, skipped = [usable[0]], []

    for entry in usable[1:]:
        p, _, blocks, _ = entry
        if list(blocks.keys()) != list(ref_blocks.keys()):
            skipped.append((p, "different sections than %s" % ref_path))
            continue

        bad = None
        for k, ref in ref_blocks.items():
            if len(blocks[k].rows) != len(ref.rows):
                bad = ("%s has %d rows, expected %d (truncated?)"
                       % (k, len(blocks[k].rows), len(ref.rows)))
                break
            if ref.kind == "TOTAL":
                continue
            nfix = n_fixed_columns(ref.kind)
            # Same bin count at different edges is the one shape mismatch a row
            # count cannot see, and it would merge silently into nonsense.
            for r_row, ref_row in zip(blocks[k].rows, ref.rows):
                if any(abs(r_row[i] - ref_row[i]) > 1e-9 * max(1.0, abs(ref_row[i]))
                       for i in range(nfix)):
                    bad = "%s is binned differently than in %s" % (k, ref_path)
                    break
            if bad:
                break

        if bad:
            skipped.append((p, bad))
            continue
        keep.append(entry)

    return keep, skipped


def merge(paths, outfile):
    usable, skipped = read_inputs(paths)
    if not usable:
        for p, why in skipped:
            sys.stderr.write("  skipped %s -- %s\n" % (p, why))
        sys.exit("ERROR: no usable reports among %d input file(s)" % len(paths))

    check_consistent(usable)
    usable, more = check_shapes(usable)
    skipped += more

    ref_path, ref_meta, ref_blocks, ref_card = usable[0]
    nvar = int(ref_meta.get("nvar", "1"))
    seeds = [m.get("seed", "?") for _, m, _, _ in usable]

    # --- accumulate ---------------------------------------------------------
    # acc[key] = (S1, S2, N) laid out like the block: one entry per TOTAL row,
    # or one entry per (bin, variation) for a histogram.
    acc = OrderedDict()
    for key, ref in ref_blocks.items():
        nrow = len(ref.rows)
        if ref.kind == "TOTAL":
            S1, S2 = [0.0] * nrow, [0.0] * nrow
        else:
            S1 = [[0.0] * nvar for _ in range(nrow)]
            S2 = [[0.0] * nvar for _ in range(nrow)]
        N = 0

        for p, _, blocks, _ in usable:
            blk = blocks[key]
            n = blk.niter
            N += n
            if ref.kind == "TOTAL":
                isw = blk.col("sum_w", p)
                for ir, row in enumerate(blk.rows):
                    S1[ir] += n * row[isw]
                    S2[ir] += n * n * row[isw + 1]
            else:
                idx = [accumulator_index(blk, iv, nvar, p) for iv in range(nvar)]
                for ir, row in enumerate(blk.rows):
                    for iv, i in enumerate(idx):
                        S1[ir][iv] += n * row[i]
                        S2[ir][iv] += n * n * row[i + 1]

        acc[key] = (S1, S2, N)

    # --- write --------------------------------------------------------------
    # Every block of a run closes on the same iteration boundary, so the pooled
    # N is a property of the campaign rather than of a block. Checked rather
    # than assumed: a file where they disagree is one where some block was
    # written from a different pass, and its accumulators would not be pooled
    # against the same denominator as the rest.
    counts = set(N for _, _, N in acc.values())
    if len(counts) != 1:
        sys.exit("ERROR: the blocks disagree on their iteration counts (%s).\n"
                 "       These reports were not all written by one production pass."
                 % ", ".join(str(c) for c in sorted(counts)))
    total_iter = counts.pop()
    npi = int(ref_meta.get("npoints_per_iter", "0"))

    with open(outfile, "w") as out:
        out.write("# ==== RUN ====\n")
        out.write("# mode = %s\n" % ref_meta.get("mode", ""))
        out.write("# generator = %s\n" % ref_meta.get("generator", ""))
        out.write("# seed = merged(%s)\n" % ",".join(str(s) for s in seeds))
        out.write("# nruns = %d\n" % len(usable))
        out.write("# npoints_per_iter = %d\n" % npi)
        out.write("# niter = %d\n" % total_iter)
        out.write("# npoints_total = %d\n" % (npi * total_iter))
        out.write("# ecm = %s\n" % ref_meta.get("ecm", ""))
        out.write("# pdf = %s\n" % ref_meta.get("pdf", ""))
        for key in ("amplitude", "cuts", "scales", "nvar", "nscale",
                    "variation_order"):
            if key in ref_meta:
                out.write("# %s = %s\n" % (key, ref_meta[key]))
        out.write("# ==== END RUN ====\n\n")

        if ref_card:
            out.write("# ==== INPUT CARD ====\n")
            # The shared card, echoed from the first input. [run] seed and
            # [output] file below describe THAT job, not this file: the run
            # header above is what describes this one.
            for line in ref_card:
                out.write(line + "\n")
            out.write("# ==== END INPUT CARD ====\n\n")

        for key, ref in ref_blocks.items():
            S1, S2, N = acc[key]
            title = ("%s %s" % (ref.contrib, ref.name)).strip()
            out.write("# ==== %s %s ====\n" % (ref.kind, title))
            out.write("# niter = %d\n" % N)
            out.write("# nruns = %d\n" % len(usable))
            for k in ("nbins", "npoints_per_iter", "xlabel"):
                if k in ref.meta:
                    v = str(npi) if k == "npoints_per_iter" else ref.meta[k]
                    out.write("# %s = %s\n" % (k, v))
            out.write("# columns: %s\n" % ref.colline)

            if ref.kind == "TOTAL":
                out.write("# chi2_per_it is the chi2/dof of the %d SEEDS about the\n"
                          "# weighted mean here, not of one run's iterations.\n"
                          % len(usable))
                _write_total(out, key, ref, usable, S1, S2, N, ref_path)
            else:
                for line in ref.header:
                    out.write(line + "\n")
                _write_hist(out, ref, S1, S2, N, nvar)

            out.write("# ==== END %s %s ====\n\n" % (ref.kind, title))

    _report(usable, skipped, paths, acc, ref_blocks, outfile, total_iter, npi)


def _write_total(out, key, ref, usable, S1, S2, N, ref_path):
    imur = ref.col("mur", ref_path)
    imuf = ref.col("muf", ref_path)
    iveg = ref.col("vegas", ref_path)

    for ir, ref_row in enumerate(ref.rows):
        vals = [blocks[key].rows[ir][iveg] for _, _, blocks, _ in usable]
        errs = [blocks[key].rows[ir][iveg + 1] for _, _, blocks, _ in usable]
        vegas, vegas_err, chi2 = weighted_mean(vals, errs)
        # htot is a one-bin histogram of unit width, so its accumulators need
        # no width division -- the same convention the Fortran writer uses.
        histo, histo_err = pooled(S1[ir], S2[ir], N)
        hw, hww = repack(S1[ir], S2[ir], N)
        out.write("%10.4f%10.4f%24.16E%24.16E%24.16E%24.16E%24.16E%24.16E%24.16E\n"
                  % (ref_row[imur], ref_row[imuf], vegas, vegas_err, chi2,
                     histo, histo_err, hw, hww))


def _write_hist(out, ref, S1, S2, N, nvar):
    nfix = n_fixed_columns(ref.kind)
    for ir, row in enumerate(ref.rows):
        width = bin_width(ref.kind, row)
        out.write("".join("%24.16E" % row[i] for i in range(nfix)))
        for iv in range(nvar):
            mean, err = pooled(S1[ir][iv], S2[ir][iv], N)
            out.write("%24.16E%24.16E" % (mean / width, err / width))
        for iv in range(nvar):
            hw, hww = repack(S1[ir][iv], S2[ir][iv], N)
            out.write("%24.16E%24.16E" % (hw, hww))
        out.write("\n")


def _report(usable, skipped, paths, acc, ref_blocks, outfile, total_iter, npi):
    print("merged %d of %d file(s): %d iterations x %d points = %d points -> %s"
          % (len(usable), len(paths), total_iter, npi, total_iter * npi, outfile))
    for p, why in skipped:
        sys.stderr.write("  WARNING skipped %s -- %s\n" % (p, why))

    key = next((k for k, b in ref_blocks.items() if b.kind == "TOTAL"), None)
    if key is None:
        return
    ref = ref_blocks[key]
    ref_path = usable[0][0]
    iveg = ref.col("vegas", ref_path)
    S1, S2, N = acc[key]

    # Row 0 is the central variation while nvar is 1; with a scale grid it is
    # the (1,1) factor pair, which is where the writer puts the central point.
    ir = 0
    vals = [blocks[key].rows[ir][iveg] for _, _, blocks, _ in usable]
    errs = [blocks[key].rows[ir][iveg + 1] for _, _, blocks, _ in usable]
    vegas, vegas_err, chi2 = weighted_mean(vals, errs)
    histo, histo_err = pooled(S1[ir], S2[ir], N)

    print("  sigma [VEGAS, weighted] = %.8E +- %.4E fb" % (vegas, vegas_err))
    print("  sigma [histograms, pooled] = %.8E +- %.4E fb" % (histo, histo_err))

    report_seed_spread(usable, key, iveg, vegas, chi2)


def report_seed_spread(usable, key, iveg, mean, chi2):
    """Compare the seeds against each other, not just against nothing.

    Summing accumulators is exact but blind: if one seed's VEGAS grid adapted
    badly and covered only part of the region, it still contributes its full
    weight, which drags the mean and SHRINKS the pooled error.  The result is a
    confidently quoted wrong number.  Comparing each seed's own value against
    its own error is what exposes that, and it is the reason to run a campaign
    at all.
    """
    rows = []
    for p, meta, blocks, _ in usable:
        v = blocks[key].rows[0][iveg]
        e = blocks[key].rows[0][iveg + 1]
        if e <= 0.0:
            continue        # skip this seed only; a degenerate error is
        rows.append((meta.get("seed", "?"), v, e,   # exactly when the rest
                     (v - mean) / e, p))            # of the check matters

    if len(rows) < 2:
        return

    print()
    print("  %-12s %22s %12s %8s" % ("seed", "sigma [fb]", "error", "pull"))
    for seed, v, e, pull, _ in rows:
        print("  %-12s %22.8E %12.4E %8.2f" % (seed, v, e, pull))
    print("  seed consistency: chi2/dof = %.2f  (%d seeds)" % (chi2, len(rows)))

    if chi2 > 3.0:
        worst = max(rows, key=lambda r: abs(r[3]))
        print()
        print("  WARNING: the seeds disagree by far more than their quoted errors.")
        print("           Worst is %s: %.6E +- %.4E fb, %.0f sigma from the mean."
              % (worst[4], worst[1], worst[2], abs(worst[3])))
        print("           The merged error does NOT cover this -- treat it as")
        print("           unconverged VEGAS adaptation, not as statistics.")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="seed-tagged reports (LO_*.dat)")
    ap.add_argument("-o", "--output", default=None,
                    help="merged report (default: the untagged name beside the inputs)")
    a = ap.parse_args()

    out = a.output
    if out is None:
        # Strip the _<seed> the writer added, which lands exactly on the name
        # the campaign reserves for the merged file.
        base = os.path.basename(a.files[0])
        out = os.path.join(os.path.dirname(a.files[0]) or ".",
                           re.sub(r"_\d+(\.[^.]*)$", r"\1", base))
        if os.path.abspath(out) == os.path.abspath(a.files[0]):
            sys.exit("ERROR: cannot infer an output name from %r -- pass -o" % base)

    merge(a.files, out)


if __name__ == "__main__":
    main()
