#!/usr/bin/env python3
"""Extract columns from a run report.

    tools/extract_data.py summary/run_00/LO.dat inv
    tools/extract_data.py summary/run_00/LO.dat ptll
    tools/extract_data.py summary/run_00/LO.dat tot     ('total' also works)
    tools/extract_data.py summary/run_00/LO.dat list

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

With 'tot', each flag also prints the uncertainty its axis exists to produce,
so neither number needs a second tool:

    --scales    the 7-point table in the 'Q sc1 sc2 sigma error' form, then the
                7-point / muF-only / muR-only envelopes around the central scale.
    --with-pdf  the PDF band (the symmetric-Hessian sum over the error members)
                and, for a set carrying a trailing alpha_s pair, the alpha_s
                band and the two in quadrature.

A report with a single variation has neither axis; both flags then print a note
and give you the central column anyway.

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
import math
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
    rows.  The writer's own 'pdf_member' column (common/histogram_mod.f90's
    write_total_xsec) is used when it exists; failing that the variation's
    index in the report's order is already unambiguous and needs no
    cooperation from the Fortran side.  Returns (None, None) when no member
    rows are being printed, which is the whole of today's output.
    """
    if not any(i >= nscale for i in which):
        return None, None
    if "pdf_member" in tot["cols"]:
        imem = tot["cols"].index("pdf_member")
        return "mem", lambda i: "%7d" % int(round(tot["rows"][i][imem]))
    return "var", lambda i: "%7d" % i


def format_sigma(val):
    """Format sigma_central as readable fixed-point for moderate magnitudes,
    falling back to scientific notation for very small or very large values."""
    aval = abs(val)
    if aval == 0:
        return "0.0"
    if 1e-2 <= aval < 1e6:
        return f"{val:.2f}"
    else:
        return f"{val:.4e}"


def compute_scale_uncertainties(Q, scale_data, nscale):
    """Compute 7-point scale uncertainties (sigma_central, 7pt/muF/muR bands).

    scale_data: dict mapping (sc1, sc2) -> sigma for this Q value.
    nscale: number of scale variations in the grid.

    Returns (sigma_central, unc_7pt, unc_muF, unc_muR) or (None, None, None, None)
    if the required 7-point grid is not present.
    """
    required_scales = {(0.5, 0.5), (0.5, 1.0), (1.0, 0.5), (1.0, 1.0),
                       (1.0, 2.0), (2.0, 1.0), (2.0, 2.0)}
    required_scales1 = {(0, 0), (0, 1), (1, 0), (1, 1),
                        (1, 2), (2, 1), (2, 2)}

    if nscale < 7:
        return None, None, None, None

    scales_present = set(scale_data.keys())
    active_scales = None
    if required_scales.issubset(scales_present):
        active_scales = required_scales
    elif required_scales1.issubset(scales_present):
        active_scales = required_scales1
    else:
        return None, None, None, None

    sigma_central = scale_data.get((1.0, 1.0))
    if sigma_central is None:
        return None, None, None, None

    sigmas_7pt = [scale_data[scale] for scale in active_scales]
    max_sigma_7pt = max(sigmas_7pt)
    min_sigma_7pt = min(sigmas_7pt)
    max_unc_7pt = ((max_sigma_7pt - sigma_central) / sigma_central) * 100.0
    min_unc_7pt = ((min_sigma_7pt - sigma_central) / sigma_central) * 100.0

    muF_scales = [s for s in active_scales if s[0] == 1.0]
    sigmas_muF = [scale_data[s] for s in muF_scales]
    max_unc_muF = ((max(sigmas_muF) - sigma_central) / sigma_central) * 100.0
    min_unc_muF = ((min(sigmas_muF) - sigma_central) / sigma_central) * 100.0

    muR_scales = [s for s in active_scales if s[1] == 1.0]
    sigmas_muR = [scale_data[s] for s in muR_scales]
    max_unc_muR = ((max(sigmas_muR) - sigma_central) / sigma_central) * 100.0
    min_unc_muR = ((min(sigmas_muR) - sigma_central) / sigma_central) * 100.0

    return sigma_central, (max_unc_7pt, min_unc_7pt), (max_unc_muF, min_unc_muF), (max_unc_muR, min_unc_muR)


def compute_pdf_uncertainties(pdf_rows, ihist, imem, pdf_set):
    """PDF and alpha_s bands from the member rows of a TOTAL block.

    pdf_rows: the rows after the scale grid, member 0 first.
    Returns (central, pdf_band, as_band, combined, as_pair, n_hessian), with
    as_band and as_pair None for a set that carries no alpha_s members, or
    None when there are no member rows at all.

    PDF4LHC21_40 is a symmetric-Hessian set, so LHAPDF's GetPDFuncertainty
    reduces to the quadrature sum below and does not need to be called.
    """
    if not pdf_rows:
        return None

    central = pdf_rows[0][ihist]

    # The trailing alpha_s pair is named by the set, not by a fixed member
    # count: the '_pdfas' suffix is LHAPDF's own convention, and bare
    # PDF4LHC21_40 also ships members 41/42 (see pdfvar_mod.f90's npdf + 2).
    name = pdf_set.strip().lower()
    has_as = name.endswith("_pdfas") or name.startswith("pdf4lhc21_40")

    hessian = pdf_rows[1:-2] if has_as else pdf_rows[1:]
    pdf_band = math.sqrt(sum((r[ihist] - central) ** 2 for r in hessian))

    if has_as:
        down, up = pdf_rows[-2][ihist], pdf_rows[-1][ihist]
        # Members 41/42 sit at alpha_s = 0.117/0.119, so half their spread is
        # the response to d(alpha_s) = 0.001. PDF4LHC21 quotes the alpha_s
        # uncertainty at 0.0015, hence the 1.5x: |up-down|/2 * 3/2.
        as_band = abs(up - down) * 0.75
        combined = math.sqrt(pdf_band ** 2 + as_band ** 2)
        as_pair = (int(pdf_rows[-2][imem]), int(pdf_rows[-1][imem]))
    else:
        as_band = None
        combined = pdf_band
        as_pair = None

    return central, pdf_band, as_band, combined, as_pair, len(hessian)


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
        cols = ["mur", "muf", "vegas", "vegas_err", "chi2_per_it",
                "histo", "histo_err"]
        idx = [column(tot, c, "TOTAL") for c in cols]
        which = selected(nvar, nscale, icentral, a)
        tag_hdr, tag = variation_tag(tot, which, nscale)

        # Every PDF member sits at mur = muf = 1, so feeding the member rows
        # into the scale table would overwrite its central point with the last
        # member's cross section. The two axes are read separately.
        scale_idx = [i for i in which if i < nscale]
        pdf_idx = [i for i in which if i >= nscale]

        if a.scales and nscale > 1:
            # Print 7-point scale table in Q sc1 sc2 sigma error format (for 7scale_extract)
            Q = float(meta.get("ecm", "0"))
            if Q == 0:
                sys.exit("ERROR: report has no 'ecm' metadata; cannot create scale table.")

            imur, imuf = idx[0], idx[1]
            ihist, ihist_err = idx[5], idx[6]

            print("# 5-column scale table for 7scale_extract")
            print("# columns: Q sc1[factor] sc2[factor] sigma[fb] error[fb]")

            scale_data = {}
            for i in scale_idx:
                r = tot["rows"][i]
                sc1, sc2 = r[imur], r[imuf]
                sigma, error = r[ihist], r[ihist_err]
                scale_data[(sc1, sc2)] = sigma
                print(" %10d %8.4f %8.4f %24.16E %24.16E" % (int(Q), sc1, sc2, sigma, error))

            result = compute_scale_uncertainties(Q, scale_data, nscale)
            if result[0] is not None:
                sigma_central, unc_7pt, unc_muF, unc_muR = result
                print()
                print("# Scale uncertainties")
                print(f"{'Q_value':<12} {'Sigma_Central':<18} {'7pt Unc(%)':<18} {'muF Unc(%)':<18} {'muR Unc(%)':<18}")
                print("-" * 84)
                sigma_str = format_sigma(sigma_central)
                unc_7pt_str = f"{unc_7pt[0]:>+5.1f} {unc_7pt[1]:>5.1f}"
                unc_muF_str = f"{unc_muF[0]:>+5.1f} {unc_muF[1]:>5.1f}"
                unc_muR_str = f"{unc_muR[0]:>+5.1f} {unc_muR[1]:>5.1f}"
                print(f"{int(Q):<12} {sigma_str:<18} {unc_7pt_str:<18} {unc_muF_str:<18} {unc_muR_str:<18}")
            else:
                print()
                print(f"{int(Q):<12} [Skipped: Missing 1 or more of the 7 scale variations]")
        else:
            # Original output format for single variation or non-scales mode
            print("#" + "%9s%9s" % ("mur", "muf") + ("%7s" % tag_hdr if tag else "") +
                  "".join("%22s" % c for c in cols[2:]))
            for i in which:
                r = tot["rows"][i]
                print(" %9.4f%9.4f" % (r[idx[0]], r[idx[1]]) + (tag(i) if tag else "") +
                      "".join("%22.13e" % r[j] for j in idx[2:]))

        if a.with_pdf and pdf_idx:
            imem = column(tot, "pdf_member", "TOTAL")
            pdf_set = meta.get("pdf", "?")
            result = compute_pdf_uncertainties([tot["rows"][i] for i in pdf_idx],
                                               idx[5], imem, pdf_set)
            if result is not None:
                central, pdf_band, as_band, combined, as_pair, nhess = result
                pct = lambda v: 100.0 * v / central if central else 0.0
                print()
                print("# PDF / alpha_s uncertainties   (%d Hessian members, %s)"
                      % (nhess,
                         "alpha_s pair = members %d,%d" % as_pair if as_pair
                         else "no alpha_s pair"))
                print(f"{'PDF_set':<22} {'Central':<14} {'PDF(%)':<12} "
                      f"{'alphas(%)':<12} {'Combined(%)':<12}")
                print("-" * 74)
                as_str = f"{pct(as_band):.2f}" if as_band is not None else "n/a"
                print(f"{pdf_set:<22} {format_sigma(central):<14} "
                      f"{pct(pdf_band):<12.2f} {as_str:<12} {pct(combined):<12.2f}")
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
