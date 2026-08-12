#!/usr/bin/env python3
"""Extract columns from a run report.

    tools/extract_data.py ../summary/run_00/LO.dat inv
    tools/extract_data.py ../summary/run_00/LO.dat ptll
    tools/extract_data.py ../summary/run_00/LO.dat tot     ('total' also works)
    tools/extract_data.py ../summary/run_00/LO.dat list

The two arguments may also be given the other way round -- the kaleu tree's
extract_data.py takes <what> first, and that habit should not cost you a rerun.
Whichever argument names a file that exists is the report.

The distribution can be given as the histogram name (mee, mee1, ptll, ...), a
case-insensitive prefix, or an alias against the recorded xlabel:
inv/invmass/mass/m -> the invariant-mass distribution, pt -> the first pt
histogram, y/rap -> the rapidity one, dphi/phi -> the azimuthal one.

    #        Q                  sigma                  error
        60.500    3.8236862948825e+03    1.1296605642249e+01
        61.500    3.7070284193348e+03    1.5477519169375e+01
        ...

Q is the bin centre; sigma is d(sigma)/d(observable) in fb per unit of that
observable, so multiply by the bin width to get the cross section in a bin.

Options:
    --scales      the whole scale grid instead of just the central point
    --with-pdf    also the PDF member variations

WHY THOSE TWO FLAGS EXIST TODAY.  This build integrates one scale point:
write_report_header stamps '# nvar = 1' and '# variation_order = central', so
both flags currently print a note and give you the central column anyway.  They
are here because the run is going to be extended to the 7-point scale variation
plus PDF members in one pass, and this reader is already written for that file.

The extension contract, so the Fortran side and this script stay in step:

  * '# nvar = N' counts the variations, '# nscale = M' says how many of them
    (the leading ones) are the scale grid; the rest are PDF members.  Without
    an explicit nscale every variation is taken to be a scale point.
  * the TOTAL block carries one row per variation, its mur/muf columns
    identifying each -- as factors, once there is a grid for a factor to index.
    A 'mem' column there labels the PDF members; without one they are labelled
    by their position in the file, which works but is less readable.
  * a histogram row carries bin_low bin_high bin_center, then N (sigma, error)
    pairs in the same order, then whatever trailing accumulator columns the
    writer keeps.  Column positions are read from the block's '# columns:'
    line, not assumed, so a new trailing column costs nothing here.

With --scales the output is grouped by bin: every scale choice of one bin, a
blank line, then the next bin.  A single variation prints one row per bin, so
there is nothing to group and no blank line is written.
"""

import argparse
import os
import re
import signal
import sys

# Let `... | head` close the pipe quietly instead of raising BrokenPipeError.
try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass


def parse(path):
    """Return (meta, blocks) where blocks maps (kind, name) -> dict."""
    meta, blocks, cur = {}, {}, None
    try:
        fh = open(path)
    except IOError as e:
        sys.exit("ERROR: %s" % e)
    with fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if re.match(r"#\s*====\s*END", s):
                cur = None
                continue
            m = re.match(r"#\s*====\s*(TOTAL|HISTOGRAM2D|HISTOGRAM)\s*(.*?)\s*====", s)
            if m:
                kind, rest = m.group(1), m.group(2).strip()
                # rest is "<contribution> <name>" for histograms,
                # "<contribution>" for totals
                parts = rest.split()
                name = parts[-1] if kind != "TOTAL" and len(parts) > 1 else ""
                contrib = " ".join(parts[:-1]) if name else rest
                cur = (kind, name)
                blocks[cur] = {"contrib": contrib, "rows": [], "meta": {}, "cols": []}
                continue
            mc = re.match(r"#\s*columns:\s*(.*)$", s)
            if mc and cur:
                # Names are stripped of their unit suffix, so 'vegas[fb]' is
                # looked up as 'vegas' and the writer can relabel units freely.
                blocks[cur]["cols"] = [re.sub(r"\[.*?\]", "", c)
                                       for c in mc.group(1).split()]
                continue
            # Only '# key = value' is metadata.  The INPUT CARD block's lines
            # read '# [run] seed = 1002' and are deliberately NOT matched here:
            # the card is echoed for provenance, and its keys collide with the
            # run header's (seed, file, ...).  Do not "fix" this regex to
            # accept them.
            mm = re.match(r"#\s*(\w+)\s*=\s*(.*)$", s)
            if mm:
                (blocks[cur]["meta"] if cur else meta)[mm.group(1)] = mm.group(2).strip()
                continue
            if s.startswith("#"):
                continue
            if cur:
                blocks[cur]["rows"].append([float(x) for x in s.split()])
    return meta, blocks


def column(blk, name, what):
    """Index of a named column in a block, or a hard error.

    Reading by name rather than by position is what lets the report grow a
    column without silently shifting the number this script prints.
    """
    try:
        return blk["cols"].index(name)
    except ValueError:
        sys.exit("ERROR: %s block has no %r column. It carries: %s"
                 % (what, name, " ".join(blk["cols"]) or "(no '# columns:' line)"))


def total_block(blocks):
    key = [k for k in blocks if k[0] == "TOTAL"]
    if not key:
        sys.exit("ERROR: no TOTAL block in the report.")
    return blocks[key[0]]


def variations(meta, blocks):
    """(nvar, nscale, icentral) for this report.

    icentral is the row of the central scale: the (1,1) factor pair once there
    is a grid, and simply the only row while nvar is 1.  A fixed-scale report
    writes mur/muf as absolute GeV, so no (1,1) is found there and the first
    row -- the only one -- is the central one by construction.
    """
    tot = total_block(blocks)
    nvar = int(meta.get("nvar", len(tot["rows"])))
    if nvar != len(tot["rows"]):
        sys.exit("ERROR: header says nvar = %d but the TOTAL block has %d rows."
                 % (nvar, len(tot["rows"])))
    nscale = int(meta.get("nscale", nvar))
    imur = column(tot, "mur", "TOTAL")
    imuf = column(tot, "muf", "TOTAL")
    icentral = next((i for i, r in enumerate(tot["rows"][:nscale])
                     if r[imur] == 1.0 and r[imuf] == 1.0), 0)
    return nvar, nscale, icentral


def selected(nvar, nscale, icentral, a):
    """Variation indices to print, in the report's own order.

    The central scale alone is the default: --scales widens that to the scale
    grid and --with-pdf appends the members, so both together are the whole
    file and --with-pdf on its own is the central scale plus members, which is
    exactly the comparison the member rows are there for.
    """
    which = [icentral]
    if a.scales:
        if nscale < 2:
            sys.stderr.write("  note: this report has a single scale point.\n")
        else:
            which = list(range(nscale))
    if a.with_pdf:
        if nvar <= nscale:
            sys.stderr.write("  note: this report has no PDF members.\n")
        which = which + list(range(nscale, nvar))
    return which


def variation_tag(tot, which, nscale):
    """(header, format function) for the column that tells variations apart.

    Scale rows identify themselves by their mur/muf factors, but PDF members
    all sit at mur = muf = 1 and would otherwise print as several identical
    rows.  The writer's own 'mem' column is used when it exists; failing that
    the variation's index in the report's order is already unambiguous and
    needs no cooperation from the Fortran side.  Returns (None, None) when no
    member rows are being printed, which is the whole of today's output.
    """
    if not any(i >= nscale for i in which):
        return None, None
    if "mem" in tot["cols"]:
        imem = tot["cols"].index("mem")
        return "mem", lambda i: "%7d" % int(round(tot["rows"][i][imem]))
    return "var", lambda i: "%7d" % i


def pick(blocks, want):
    """Resolve a user string to a histogram block key."""
    hists = [k for k in blocks if k[0] in ("HISTOGRAM", "HISTOGRAM2D")]
    w = want.lower()
    aliases = {"inv": "mass", "invmass": "mass", "mass": "mass", "m": "mass",
               "pt": "pt", "y": "y", "rap": "y", "dphi": "dphi", "phi": "dphi"}
    # exact name, then prefix, then alias against the recorded xlabel
    for k in hists:
        if k[1].lower() == w:
            return k
    for k in hists:
        if k[1].lower().startswith(w):
            return k
    if w in aliases:
        for k in hists:
            if blocks[k]["meta"].get("xlabel", "").lower() == aliases[w]:
                return k
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", help="the run report to read")
    ap.add_argument("what", help="distribution name/alias, or 'tot'/'total', or 'list'")
    ap.add_argument("--scales", action="store_true",
                    help="the whole scale grid, not just the central point")
    ap.add_argument("--with-pdf", action="store_true", dest="with_pdf",
                    help="also the PDF member variations")
    a = ap.parse_args()

    # Accept the kaleu tool's <what> <file> order as well. No observable is
    # ever the name of an existing file, so the test cannot pick wrong; if
    # neither argument names a file, the order is left alone and parse()
    # reports the missing one the user actually meant.
    if not os.path.exists(a.file) and os.path.exists(a.what):
        a.file, a.what = a.what, a.file

    meta, blocks = parse(a.file)
    nvar, nscale, icentral = variations(meta, blocks)

    want = a.what.lower()

    if want == "list":
        print("# mode = %s   seed = %s   npoints = %s x %s iterations"
              % (meta.get("mode", "?"), meta.get("seed", "?"),
                 meta.get("npoints_per_iter", "?"), meta.get("niter", "?")))
        print("# scales = %s" % meta.get("scales", "?"))
        print("# pdf = %s   nvar = %d" % (meta.get("pdf", "?"), nvar))
        print("# available:")
        print("   %-10s %-10s" % ("total", ""))
        for k in blocks:
            if k[0] != "TOTAL":
                b = blocks[k]
                print("   %-10s %-10s %d bins" %
                      (k[1], b["meta"].get("xlabel", ""), len(b["rows"])))
        return

    if want in ("total", "tot"):
        tot = total_block(blocks)
        # Both estimators, on purpose: VEGAS quotes a variance-weighted average
        # over iterations, the histograms an unweighted one, and their
        # agreement is the cheapest signal that the run has converged. Printing
        # only one of them would throw that check away.
        cols = ["mur", "muf", "vegas", "vegas_err", "chi2_per_it",
                "histo", "histo_err"]
        idx = [column(tot, c, "TOTAL") for c in cols]
        which = selected(nvar, nscale, icentral, a)
        tag_hdr, tag = variation_tag(tot, which, nscale)
        print("#" + "%9s%9s" % ("mur", "muf") + ("%7s" % tag_hdr if tag else "") +
              "".join("%22s" % c for c in cols[2:]))
        for i in which:
            r = tot["rows"][i]
            print(" %9.4f%9.4f" % (r[idx[0]], r[idx[1]]) + (tag(i) if tag else "") +
                  "".join("%22.13e" % r[j] for j in idx[2:]))
        return

    key = pick(blocks, a.what)
    if key is None:
        names = ", ".join(sorted(k[1] for k in blocks if k[0] != "TOTAL"))
        sys.exit("ERROR: no distribution matching %r. Available: %s\n"
                 "       (or use 'tot' / 'list')" % (a.what, names))
    if key[0] == "HISTOGRAM2D":
        sys.exit("ERROR: %s is 2D; this extractor handles 1D distributions." % key[1])

    blk = blocks[key]
    isig = column(blk, "dsigma_dx", key[1])
    ierr = column(blk, "error", key[1])
    icen = column(blk, "bin_center", key[1])
    if ierr != isig + 1:
        sys.exit("ERROR: %s has error at column %d, not next to dsigma_dx at %d; "
                 "the (sigma, error) pairing this reader assumes is gone."
                 % (key[1], ierr, isig))

    which = selected(nvar, nscale, icentral, a)
    tot = total_block(blocks)
    imur, imuf = column(tot, "mur", "TOTAL"), column(tot, "muf", "TOTAL")
    tag_hdr, tag = variation_tag(tot, which, nscale)
    grid = len(which) > 1

    print("# %s: %s   [%s]" % (key[1], blk["meta"].get("xlabel", ""), blk["contrib"]))
    print("# mode = %s   seed = %s   npoints = %s x %s iterations"
          % (meta.get("mode", "?"), meta.get("seed", "?"),
             meta.get("npoints_per_iter", "?"), meta.get("niter", "?")))
    # The scale columns only say something once more than one variation is
    # printed; with a single one they would repeat the header's fixed scales on
    # every line and push sigma out of the column the plotting scripts read.
    print("#%9s" % "Q" + ("%9s%9s" % ("mur", "muf") if grid else "") +
          ("%7s" % tag_hdr if tag else "") + "%22s%22s" % ("sigma", "error"))

    for r in blk["rows"]:
        for i in which:
            scales = ("%9.4f%9.4f" % (tot["rows"][i][imur], tot["rows"][i][imuf])
                      if grid else "")
            print(" %9.3f" % r[icen] + scales + (tag(i) if tag else "") +
                  "%22.13e%22.13e" % (r[isig + 2*i], r[isig + 2*i + 1]))
        # Separates one bin's variations from the next bin's. With a single row
        # per bin there is nothing to separate and it only doubles the output.
        if grid:
            print()


if __name__ == "__main__":
    main()
