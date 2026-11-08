#!/usr/bin/env python3
"""Live monitor for a run_nlo_seeds.sh campaign.

Reads flushed VEGAS logs and the coordinator's marker files.  Each row is one
seed and has a progress/status line followed by a value/error line.  Adaptation
values are intentionally shown (in red) and participate in the live partial
NLO sum, as requested; production values are green and completed jobs cyan.
"""
import argparse
import curses
import math
import re
import time
from pathlib import Path

REFRESH = 2.0
TARGETS = ("lo", "virt", "kp", "rs")
LABELS = {"lo":"LO", "virt":"VIRT", "kp":"KP", "rs":"RS"}
PATHS = {"lo":("born","born"), "virt":("virt","virt"), "kp":("kp","kp"),
         "rs":("rs","rs")}
PASS_RE = re.compile(r"\[ VEGAS \]\s+(adaptation|production) pass", re.I)
ITER_RE = re.compile(r"^\s*(\d+)\s+([+-]?[\d.]+(?:[EDed][+-]?\d+)?)\s+([+-]?[\d.]+(?:[EDed][+-]?\d+)?)\s+([+-]?[\d.]+(?:[EDed][+-]?\d+)?)\s+([+-]?[\d.]+(?:[EDed][+-]?\d+)?)")

# input.input has one [stats_*] block per program -- see its own header --
# giving each pass's iteration count (it_adapt/it_prod) directly, which is
# never echoed to the log. rsgg/rsqg/rsqq all share [stats_rs], since run_rs
# takes its channel on the command line, not from a separate stats block.
STATS_SECTION = {"lo":"stats_lo", "kp":"stats_kp", "virt":"stats_virt",
                  "rs":"stats_rs"}

def info(path):
    out = {}
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1); out[k.strip()] = v.strip()
    return out

def parse_card_stats(path):
    """Return {target: (it_adapt, it_prod)} read from input.input's [stats_*] blocks."""
    sections = {}
    section = None
    for line in path.read_text(errors="replace").splitlines():
        line = line.split("#", 1)[0].split("!", 1)[0].strip()
        if not line: continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip(); sections[section] = {}
            continue
        if section is not None and "=" in line:
            k, v = line.split("=", 1); sections[section][k.strip()] = v.strip()
    totals = {}
    for target, name in STATS_SECTION.items():
        block = sections.get(name, {})
        try:
            totals[target] = (int(block["it_adapt"]), int(block["it_prod"]))
        except (KeyError, ValueError):
            totals[target] = (None, None)
    return totals

def number(s): return float(s.replace("D", "E").replace("d", "e"))

def parse_log(path):
    """Return latest VEGAS cumulative estimate and its pass/iteration."""
    if not path.exists(): return None
    phase, latest = None, None
    try:
        lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return None
    for line in lines:
        m = PASS_RE.search(line)
        if m: phase = m.group(1).lower(); continue
        m = ITER_RE.match(line)
        if m and phase:
            try: latest = (phase, int(m.group(1)), number(m.group(4)), number(m.group(5)))
            except ValueError: pass
    return latest

def state(status, seed, target):
    d = status / str(seed)
    for name in ("failed", "done", "running", "queued"):
        if (d / f"{target}.{name}").exists(): return name
    return "queued"

def fmt(x, width=15):
    return f"{x:.5g}".rjust(width)

def estimate(value, error):
    return f"{value:.5g}+/-{error:.3g}"

def put(scr, y, x, text, width, attr=0):
    """Draw a clipped string; a narrow terminal must not stop monitoring."""
    if x >= width or y < 0:
        return
    try:
        scr.addnstr(y, x, text, max(0, width - x - 1), attr)
    except curses.error:
        pass

def draw(scr, cfg):
    curses.curs_set(0); scr.nodelay(True); scr.timeout(100)
    curses.start_color(); curses.use_default_colors()
    curses.init_pair(1, curses.COLOR_RED, -1); curses.init_pair(2, curses.COLOR_GREEN, -1)
    curses.init_pair(3, curses.COLOR_CYAN, -1); curses.init_pair(4, curses.COLOR_YELLOW, -1)
    base = Path(cfg["base_out"]); status = base / ".nlo_campaign/status"
    totals = parse_card_stats(Path(cfg["card"]))
    seeds = cfg["seeds"].split(); colw = 17; scroll = 0
    while True:
        rows = [] ; counts = {"done":0, "running":0, "queued":0, "failed":0}
        by_target = {t:0 for t in TARGETS}
        for seed in seeds:
            cells, usable = {}, []
            for target in TARGETS:
                subdir, stem = PATHS[target]; st = state(status, seed, target); counts[st] += 1
                if st == "done": by_target[target] += 1
                val = parse_log(base / subdir / "logs" / f"{stem}_{seed}.log")
                cells[target] = (st, val)
                if val: usable.append(val)
            if usable:
                total = sum(x[2] for x in usable); err = math.sqrt(sum(x[3]**2 for x in usable))
                adap = any(x[0] == "adaptation" for x in usable)
                nlo = (len(usable), total, err, adap)
            else: nlo = (0, None, None, False)
            rows.append((seed, cells, nlo))
        scr.erase(); h, w = scr.getmaxyx()
        put(scr, 0, 0, f"On-shell NLO campaign  {base}", w, curses.A_BOLD)
        put(scr, 1, 0, "q quit   r refresh   arrows scroll", w, curses.A_DIM)
        seedw = 7
        colw = 15
        def cell_text(text):
            return "|" + text.center(colw)
        header = "seed".ljust(seedw) + "".join(cell_text(LABELS.get(t, t)) for t in TARGETS + ("NLO",))
        put(scr, 3, 0, header, w, curses.A_UNDERLINE)
        visible = max(1, (h - 6) // 2)
        scroll = min(scroll, max(0, len(rows)-visible))
        for i, (seed, cells, nlo) in enumerate(rows[scroll:scroll+visible]):
            y = 4 + 2*i; put(scr, y, 0, str(seed).ljust(seedw), w, curses.A_BOLD)
            for j, target in enumerate(TARGETS):
                st, val = cells[target]; x = seedw + j*(colw+1)
                if val:
                    phase, it, _, _ = val
                    it_adapt, it_prod = totals.get(target, (None, None))
                    denom = it_adapt if phase == "adaptation" else it_prod
                    text = f"{('A' if phase=='adaptation' else 'P')} {it}/{denom or '?'}"
                    color = 1 if phase == "adaptation" else 2
                else:
                    text = "queued" if st == "queued" else st; color = 4 if st == "running" else (3 if st == "done" else 1 if st == "failed" else 0)
                put(scr, y, x, text.center(colw), w, curses.color_pair(color))
                if val: put(scr, y+1, x, estimate(val[2], val[3]).center(colw), w, curses.color_pair(color))
            n, total, err, adap = nlo; x = seedw + len(TARGETS)*(colw+1)
            color = 1 if adap else (2 if n == 6 else 4)
            put(scr, y, x, f"{n}/6".center(colw), w, curses.color_pair(color))
            if n: put(scr, y+1, x, estimate(total, err).center(colw), w, curses.color_pair(color))
        footer = f"Jobs: {len(seeds)*6} total | {counts['done']} finished | {counts['running']} running | {counts['queued']} queued | {counts['failed']} failed"
        put(scr, h-2, 0, footer, w, curses.A_BOLD)
        detail = "Completed: " + " · ".join(f"{LABELS[t]} {by_target[t]}/{len(seeds)}" for t in TARGETS)
        put(scr, h-1, 0, detail, w)
        scr.refresh()
        end = time.monotonic() + REFRESH
        while time.monotonic() < end:
            key = scr.getch()
            if key in (ord('q'), 27): return
            if key == curses.KEY_DOWN: scroll = min(scroll+1, max(0, len(rows)-visible))
            if key == curses.KEY_UP: scroll = max(0, scroll-1)
            if key == ord('r'): break
            time.sleep(.05)

def main():
    p = argparse.ArgumentParser(); p.add_argument("manifest", nargs="?", default="campaign.info")
    args = p.parse_args(); path = Path(args.manifest).resolve()
    if not path.is_file(): p.error(f"manifest not found: {path}")
    cfg = info(path)
    if "base_out" not in cfg or "seeds" not in cfg: p.error("invalid campaign manifest")
    curses.wrapper(draw, cfg)

if __name__ == "__main__": main()
