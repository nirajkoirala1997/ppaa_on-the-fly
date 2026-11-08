#!/usr/bin/env python3
"""Rebuild a per-seed report (<contrib>_<seed>.dat) from a VEGAS checkpoint.

The Fortran programs write their report only when a run finishes, but
checkpoint_<seed>.dat is rewritten after every production iteration and
holds every accumulator the report is built from. This script applies the
same formulas as the Fortran writer (histogram_mod.f90, analysis_mod.f90,
vsup.f) so an interrupted or still-running campaign can be inspected and
merged with tools/merge_runs.py.

Usage:
    checkpoint_to_report.py CHECKPOINT_OR_DIR [...] [-o OUT | --outdir DIR]

Default output: partial/<contrib>_<seed>.dat next to each checkpoint.
The whole report, header included, is derived from the checkpoint alone;
only the amplitude string comes from a per-contribution table (amplitude_label).
"""

import argparse
import math
import re
import sys
from pathlib import Path

DEFAULT_SCALE_ORDER = "(0.5,0.5) (0.5,1) (1,0.5) (1,1) (1,2) (2,1) (2,2)"


def fe(x, w=24, d=16):
    """Fortran Ew.d, e.g. 0.5994265130202753E+02."""
    if x == 0.0 or not math.isfinite(x):
        mant, e = "0." + "0" * d, 0
        if not math.isfinite(x):
            return ("%" + str(w) + "s") % repr(x)
    else:
        s = "%.*e" % (d - 1, abs(x))
        digits, ex = s.split("e")
        mant = "0." + digits.replace(".", "")
        e = int(ex) + 1
    exp = "E%+03d" % e if abs(e) <= 99 else "%+04d" % e
    s = ("-" if x < 0 else "") + mant + exp
    return s.rjust(w)


# ---------------------------------------------------------------- parsing

class Histo:
    def __init__(self, name, nbins, nvar, nit, hw, hww):
        self.name, self.nbins, self.nvar, self.nit = name, nbins, nvar, nit
        self.hw, self.hww = hw, hww  # indexed [iv][ib]


def parse_checkpoint(path):
    lines = Path(path).read_text().splitlines()
    ck = {"path": Path(path)}
    i = 0

    def find(marker, start):
        for j in range(start, len(lines)):
            if lines[j].strip() == marker:
                return j
        raise ValueError("%s: missing '%s' -- truncated checkpoint?" % (path, marker))

    i = find("# ==== CHECKPOINT ====", 0)
    for key in ("seed", "ndim", "n_prod", "it_prod", "iterations_done"):
        i += 1
        m = re.match(r"#\s*%s\s*=\s*(-?\d+)" % key, lines[i])
        if not m:
            raise ValueError("%s: expected '# %s = ' at line %d" % (path, key, i + 1))
        ck[key] = int(m.group(1))
    find("# ==== END CHECKPOINT ====", i)

    c0 = find("# ==== INPUT CARD ====", i)
    c1 = find("# ==== END INPUT CARD ====", c0)
    card_lines = lines[c0:c1 + 1]
    if c1 + 1 < len(lines) and lines[c1 + 1].strip() == "":
        card_lines.append("")
    ck["card_lines"] = card_lines
    card = {}
    for ln in card_lines:
        m = re.match(r"#\s*\[(\w+)\]\s*(.*?)\s*=\s*(.*?)\s*$", ln)
        if m:
            card[(m.group(1), m.group(2))] = m.group(3)
    ck["card"] = card

    v = find("# ==== VEGAS ====", c1)
    ints = [int(t) for t in lines[v + 1].split()]
    reals = [float(lines[v + 2 + k].replace("D", "E")) for k in range(9)]
    ck["vegas_it"] = ints[6]
    (_, _, _, _, _, ck["si"], ck["si2"], ck["swgt"], ck["schi"]) = reals

    h0 = find("# ==== HISTOGRAMS ====", v)
    tok = iter(lines[h0 + 1:])
    nhists = int(next(tok))
    hists = []
    for _ in range(nhists + 1):  # booked histograms, then 'total'
        name = next(tok).strip()
        nbins, nvar, nit = (int(t) for t in next(tok).split())
        n = nbins * nvar
        flat_w = [float(next(tok)) for _ in range(n)]
        flat_ww = [float(next(tok)) for _ in range(n)]
        # Fortran column order: bin index runs fastest.
        hw = [flat_w[iv * nbins:(iv + 1) * nbins] for iv in range(nvar)]
        hww = [flat_ww[iv * nbins:(iv + 1) * nbins] for iv in range(nvar)]
        hists.append(Histo(name, nbins, nvar, nit, hw, hww))
    ck["hists"], ck["htot"] = hists[:-1], hists[-1]
    return ck


def header_value(run_lines, key):
    for ln in run_lines or ():
        m = re.match(r"#\s*%s\s*=\s*(.*?)\s*$" % re.escape(key), ln)
        if m:
            return m.group(1)
    return None


# ------------------------------------------------------------ computation

def rescale_ceiling(h, it_prod):
    """histo1d_rescale_ceiling: undo the WGT/it_prod fold of an early stop."""
    if h.nit <= 0 or h.nit >= it_prod:
        return
    c = it_prod / h.nit
    h.hw = [[x * c for x in col] for col in h.hw]
    h.hww = [[x * c * c for x in col] for col in h.hww]


def h_mean(h, ib, iv, bw):
    return h.hw[iv][ib] / bw


def h_error(h, ib, iv, bw):
    if h.nit < 2:
        return 0.0
    n = float(h.nit)
    return math.sqrt(max(0.0, n * h.hww[iv][ib] - h.hw[iv][ib] ** 2) / (n - 1.0)) / bw


def vegas_result(ck):
    """vsup.f end-of-iteration: AVGI, SD, CHI2A from the running sums."""
    it, si, si2, swgt, schi = (ck["vegas_it"], ck["si"], ck["si2"],
                               ck["swgt"], ck["schi"])
    avgi = si / swgt
    t = swgt * it / si2
    chi2 = t * (schi / swgt - avgi * avgi) / (it - 1) if it > 1 else 0.0
    return avgi, math.sqrt(1.0 / t), chi2


# ---------------------------------------------------------------- writing

def f0(x, d):
    """Fortran F0.d: gfortran drops the leading zero (0.5 -> .5)."""
    s = "%.*f" % (d, x)
    if s.startswith("0."):          # gfortran writes 0 as '.0' too
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return s


def card_float(card, sec, key, default):
    v = card.get((sec, key))
    return float(v.replace("d", "e").replace("D", "e")) if v else default


def card_bool(card, sec, key, default):
    v = card.get((sec, key))
    if v is None:
        return default
    return v.strip().strip(".").lower() in ("yes", "true", "t", "1", "on")


def cuts_label(card, ecm):
    """cuts_mod.f90 cuts_label: the m(ee) window."""
    low = card_float(card, "cuts", "mee_min", 0.0)
    high = card_float(card, "cuts", "mee_max", 0.0)
    if high <= 0.0:
        high = ecm
    return "mee=%s:%s" % (f0(low, 3), f0(high, 3))


def scales_label(card):
    """common/run_card_mod.f90 scales_label."""
    cf = card_float(card, "scales", "central_factor", 1.0)
    cfs = ("%12.4f" % cf).strip() if abs(cf - 1.0) > 1e-12 else ""
    if card_bool(card, "scales", "dynamic", False):
        mid = cfs + " x " if cfs else ""
        return "dynamic: muR = muF = %sm(final state, partons excluded), 7-point grid" % mid
    mur = card_float(card, "scales", "mur_fixed", 0.0)
    muf = card_float(card, "scales", "muf_fixed", 0.0)
    return "fixed: muR = %s, muF = %s, 7-point grid" % (f0(cf * mur, 4), f0(cf * muf, 4))


# The '# amplitude =' strings, copied verbatim from each program's main.f90
# (lo_sa/kp_sa: provider + proc_str; virt_sa: a literal; rs_sa: proc_str).
# merge_runs.py checks this header field across seeds, so a report rebuilt
# here must carry exactly what the Fortran writer would -- update this table
# whenever one of those strings changes.
AMPLITUDE = {
    "born": "recola: q q~ -> e+ e-",
    "kp":   "recola: KP on q q~ -> e+ e-",
    "virt": "recola: q q~ -> e+ e- one-loop QCD + I",
    "rs":   "recola: q q~, q g, g q -> e+ e- + parton",
}


def amplitude_label(prefix, card):
    return AMPLITUDE.get(prefix, "unknown")


def build_run_header(ck, prefix, contrib, nvar):
    n_prod, niter = ck["n_prod"], ck["htot"].nit
    card = ck["card"]
    ecm = card_float(card, "beam", "ecm", 0.0)
    source = "# source = checkpoint (iterations_done = %d of it_prod = %d)" % (
        ck["iterations_done"], ck["it_prod"])
    nscale = len(DEFAULT_SCALE_ORDER.split())
    npdf = nvar - nscale if nvar > nscale else 0
    var_order = "scale:1-%d" % nscale + (" pdfmem:0-%d" % (npdf - 1) if npdf else "")
    return [
        "# ==== RUN ====",
        "# mode = " + contrib,
        "# generator = PSGen + VEGAS (vsup)",
        "# seed = %d" % ck["seed"],
        "# npoints_per_iter = %d" % n_prod,
        "# niter = %d" % niter,
        "# npoints_total = %d" % (n_prod * niter),
        "# ecm = " + fe(ecm),
        "# pdf = " + card.get(("pdf", "set"), ""),
        "# amplitude = " + amplitude_label(prefix, card),
        "# cuts = " + cuts_label(card, ecm),
        "# scales = " + scales_label(card),
        "# nscale = %d" % nscale,
        "# scale_order = " + DEFAULT_SCALE_ORDER,
        "# npdfvar = %d" % npdf,
        "# pdf_central_member = " + card.get(("pdf", "member"), "0"),
        "# nvar = %d" % nvar,
        "# variation_order = " + var_order,
        source,
        "# ==== END RUN ====",
    ]


def write_report(ck, prefix, contrib, outfile):
    card = ck["card"]
    htot = ck["htot"]
    nvar = htot.nvar
    for h in ck["hists"] + [htot]:
        rescale_ceiling(h, ck["it_prod"])

    run = build_run_header(ck, prefix, contrib, nvar)
    scale_order = header_value(run, "scale_order") or DEFAULT_SCALE_ORDER
    var_order = header_value(run, "variation_order")
    central_mem = int(header_value(run, "pdf_central_member")
                      or card.get(("pdf", "member"), "0"))
    scales = [tuple(float(x) for x in p.strip("()").split(","))
              for p in scale_order.split()]
    nscale = len(scales)
    icentral = scales.index((1.0, 1.0))

    ans, sd, chi2 = vegas_result(ck)
    o = []
    o.extend(run)
    o.append("")
    o.extend(ck["card_lines"])

    o.append("# ==== TOTAL %s ====" % contrib)
    o.append("# niter = %d" % htot.nit)
    o.append("# columns: mur[factor] muf[factor] pdf_member vegas[fb] "
             "vegas_err[fb] chi2_per_it histo[fb] histo_err[fb] sum_w sum_w2")
    for iv in range(nvar):
        if iv < nscale:
            fmur, fmuf = scales[iv]
            imem = central_mem
        else:
            fmur = fmuf = 1.0
            imem = iv - nscale
        hm, he = h_mean(htot, 0, iv, 1.0), h_error(htot, 0, iv, 1.0)
        v, verr, c2 = (ans, sd, chi2) if iv == icentral else (hm, he, 0.0)
        o.append("%10.4f%10.4f%6d" % (fmur, fmuf, imem) + "".join(
            fe(x) for x in (v, verr, c2, hm, he, htot.hw[iv][0], htot.hww[iv][0])))
    o.append("# ==== END TOTAL %s ====" % contrib)
    o.append("")

    for h in ck["hists"]:
        if h.nbins == 0:
            continue
        key = lambda k: card.get(("histograms", "%s %s" % (h.name, k)))
        xmin, xmax = float(key("xmin")), float(key("xmax"))
        bw = (xmax - xmin) / h.nbins
        obs = key("obs") or h.name
        legs = key("legs") or key("leg")
        xlabel = "%s(%s)" % (obs, legs.replace(" ", "")) if legs else obs

        o.append("# ==== HISTOGRAM %s %s ====" % (contrib, h.name))
        o.append("# nbins = %d" % h.nbins)
        o.append("# niter = %d" % h.nit)
        o.append("# npoints_per_iter = %d" % ck["n_prod"])
        o.append("# xlabel = " + xlabel)
        o.append("# scale_order = " + scale_order)
        if var_order:
            o.append("# variation_order = " + var_order)
        o.append("# columns 1-3: bin_low bin_high bin_center ; "
                 "columns 4-%d: (dsigma_dx error) per variation" % (3 + 2 * nvar))
        o.append("# columns %d-%d: (sum_w sum_w2) per variation -- raw "
                 "per-iteration accumulators (NOT divided by the bin width)"
                 % (4 + 2 * nvar, 3 + 4 * nvar))
        o.append("# columns: bin_low bin_high bin_center"
                 + " dsigma_dx[fb] error[fb]" * nvar + " sum_w sum_w2" * nvar)
        for ib in range(h.nbins):
            lo = xmin + ib * bw
            hi = lo + bw
            row = [fe(lo), fe(hi), fe(0.5 * (lo + hi))]
            for iv in range(nvar):
                row.append(fe(h_mean(h, ib, iv, bw)) + fe(h_error(h, ib, iv, bw)))
            for iv in range(nvar):
                row.append(fe(h.hw[iv][ib]) + fe(h.hww[iv][ib]))
            o.append("".join(row))
        o.append("# ==== END HISTOGRAM %s %s ====" % (contrib, h.name))
        o.append("")

    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_text("\n".join(o) + "\n")


# ------------------------------------------------------------------- main

def file_prefix(ck_path):
    """Report file prefix: the program's directory."""
    d = ck_path.resolve().parent
    return d.name


def mode_label(prefix, card):
    return prefix


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inputs", nargs="+", help="checkpoint files or directories")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("-o", "--output", help="output file (single checkpoint only)")
    g.add_argument("--outdir", help="write every report into this directory")
    args = ap.parse_args()

    paths = []
    for p in map(Path, args.inputs):
        paths.extend(sorted(p.glob("checkpoint_*.dat")) if p.is_dir() else [p])
    if not paths:
        sys.exit("no checkpoint_*.dat found")
    if args.output and len(paths) > 1:
        sys.exit("-o takes a single checkpoint; use --outdir for several")

    for p in paths:
        ck = parse_checkpoint(p)
        prefix = file_prefix(p)
        contrib = mode_label(prefix, ck["card"])
        name = "%s_%d.dat" % (prefix, ck["seed"])
        if args.output:
            out = Path(args.output)
        elif args.outdir:
            out = Path(args.outdir) / name
        else:
            out = p.parent / "partial" / name
        write_report(ck, prefix, contrib, out)
        print("%s -> %s  (%d/%d iterations)" % (p, out, ck["htot"].nit, ck["it_prod"]))


if __name__ == "__main__":
    main()
