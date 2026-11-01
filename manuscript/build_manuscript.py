#!/usr/bin/env python3
"""
build_manuscript.py -- tag-based citation management.

Pass 1 (--convert): rewrite existing superscript numerals as {{TAG}} markers.
Pass 2 (--build)  : assign numbers in order of first appearance, emit final
                    markdown with superscripts + an ordered reference list,
                    and report uncited / undefined tags.

Renumbering is automatic, so inserting a citation anywhere cannot corrupt
the ordering -- the failure mode this script exists to prevent.
"""
import json, re, sys, argparse

SUP = {'0':'⁰','1':'¹','2':'²','3':'³','4':'⁴','5':'⁵',
       '6':'⁶','7':'⁷','8':'⁸','9':'⁹'}
def sup(n): return ''.join(SUP[c] for c in str(n))

# original numbering in Manuscript_SLC1A2_rs4354668.md
LEGACY = {1:"RAJKOWSKA13", 2:"CHOUDARY05", 3:"MEDINA13", 4:"CHANDLEY13",
          5:"BECHTHOLT10", 6:"JOHN12", 7:"ALS23", 8:"MALLOLAS06",
          9:"RODEK22", 10:"BORDER19", 11:"PGCMDD25", 12:"WACHOLDER04",
          13:"IOANNIDIS07", 14:"THAWEETHEE19", 15:"MERK19", 16:"ZHANG15",
          17:"NAGAI09", 18:"DENG04", 19:"GELMAN06", 20:"HARDINGHAM10"}


def convert(src, dst):
    t = open(src).read()
    t = t.split("## References")[0] + "{{REFERENCES}}\n"
    # longest numerals first so ¹⁰ is not read as ¹ then ⁰
    for n in sorted(LEGACY, reverse=True):
        t = t.replace(sup(n), "{{%s}}" % LEGACY[n])
    t = t.replace("˒", "")           # citation separator, regenerated on build
    open(dst, "w").write(t)
    print(f"converted -> {dst}")


def build(src, bib_path, dst):
    t = open(src).read()
    bib = json.load(open(bib_path))

    body = t.split("{{REFERENCES}}")[0]
    order, seen = [], set()
    for m in re.finditer(r"\{\{([A-Z0-9_]+)\}\}", body):
        tag = m.group(1)
        if tag not in seen:
            seen.add(tag); order.append(tag)

    undefined = [x for x in order if x not in bib]
    uncited   = [k for k in bib if k not in seen]
    if undefined:
        sys.exit(f"ERROR undefined tags: {undefined}")

    num = {tag: i + 1 for i, tag in enumerate(order)}

    # collapse adjacent markers into a single superscript run
    def collapse(m):
        tags = re.findall(r"\{\{([A-Z0-9_]+)\}\}", m.group(0))
        ns = sorted(num[x] for x in tags)
        out, i = [], 0
        while i < len(ns):                      # compress runs into a–b
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                j += 1
            out.append(sup(ns[i]) if j == i else f"{sup(ns[i])}–{sup(ns[j])}")
            i = j + 1
        return "˒".join(out)

    body = re.sub(r"(?:\{\{[A-Z0-9_]+\}\})+", collapse, body)

    refs = "## References\n\n" + "\n".join(
        f"{num[tag]}. {bib[tag]}" for tag in order) + "\n"
    tail = t.split("{{REFERENCES}}")[1]
    open(dst, "w").write(body + refs + tail)

    print(f"built -> {dst}")
    print(f"  {len(order)} references, numbered in order of first citation")
    print(f"  uncited entries in refs.json: {uncited or 'none'}")
    print("\n  final order:")
    for tag in order:
        print(f"   {num[tag]:>2}. {tag}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--convert", nargs=2, metavar=("SRC", "DST"))
    ap.add_argument("--build", nargs=3, metavar=("SRC", "BIB", "DST"))
    a = ap.parse_args()
    if a.convert: convert(*a.convert)
    if a.build:   build(*a.build)
