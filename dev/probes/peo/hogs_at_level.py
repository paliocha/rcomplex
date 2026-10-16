"""Cut the FastOMA hog hierarchy at the MRCA of each experiment's species.

Usage: python3 hogs_at_level.py <orthoxml> <tree.nwk> <sets.tsv> <root>
Writes <root>/complex_<exp>/inputs/orthologs.tsv (species, gene, hog),
hog_age.tsv (hog, root_hog, root_taxrange: the TaxRange of the hog's
root orthologGroup), age_order.txt (tree nodes from the root down to
the level, then the nodes below it) and level.txt, for every
experiment in sets.tsv, in one streaming pass. An existing
orthologs.tsv is kept, not overwritten.

A gene's hog at level L is the innermost enclosing orthologGroup whose
TaxRange is L or an ancestor of L, with ids <L>:<n> in file order. A
gene with no such group sits in a root hog born below L; that root hog
is one young hog, <L>:y<n>. Genes in no group get no row. Gene ids are
the protId up to the first whitespace.
"""
import collections
import os
import sys
import xml.etree.ElementTree as ET

NS = "{http://orthoXML.org/2011/}"


def parse_newick(text):
    """Return (labels, parent) lists; parent[root] is None."""
    labels, parent, stack = [], [], []
    last, i, n = None, 0, len(text)
    while i < n:
        c = text[i]
        if c == "(":
            labels.append(None)
            parent.append(stack[-1] if stack else None)
            stack.append(len(labels) - 1)
            last = None
            i += 1
        elif c == ",":
            last = None
            i += 1
        elif c == ")":
            last = stack.pop()
            i += 1
        elif c == ";":
            break
        elif c.isspace():
            i += 1
        else:
            j = i
            while j < n and text[j] not in "(),;":
                j += 1
            name = text[i:j].strip().split(":")[0]
            if last is not None:
                labels[last] = name
            else:
                labels.append(name)
                parent.append(stack[-1] if stack else None)
            last = None
            i = j
    return labels, parent


def path_up(k, parent):
    out = []
    while k is not None:
        out.append(k)
        k = parent[k]
    return out


def level_of(species, labels, parent):
    """MRCA label of species, and the labels of it and its ancestors."""
    idx = {lab: k for k, lab in enumerate(labels)}
    miss = [s for s in species if s not in idx]
    if miss:
        sys.exit("species not in tree: " + ", ".join(miss))
    paths = [path_up(idx[s], parent) for s in species]
    common = set(paths[0]).intersection(*map(set, paths[1:]))
    mrca = next(k for k in paths[0] if k in common)
    return labels[mrca], [labels[k] for k in path_up(mrca, parent)][::-1]


def below(k, labels, parent):
    """Labels of the nodes under node k, in preorder."""
    kids = collections.defaultdict(list)
    for c, p in enumerate(parent):
        kids[p].append(c)
    out, todo = [], list(reversed(kids[k]))
    while todo:
        c = todo.pop()
        out.append(labels[c])
        todo.extend(reversed(kids[c]))
    return out


def read_sets(path):
    sets = collections.OrderedDict()
    with open(path) as f:
        next(f)
        for line in f:
            exp, sp, _ = line.rstrip("\n").split("\t")
            sets.setdefault(exp, []).append(sp)
    return sets


def main(xml_path, tree_path, sets_path, root):
    with open(tree_path) as f:
        labels, parent = parse_newick(f.read())
    sets = read_sets(sets_path)
    level, ok = {}, {}
    for exp, sp in sets.items():
        level[exp], ok[exp] = level_of(sp, labels, parent)
        print(f"{exp}: {len(sp)} species, level {level[exp]}", flush=True)
    by_sp = collections.defaultdict(list)
    for exp, sp in sets.items():
        for s in sp:
            by_sp[s].append(exp)

    genes = {}  # orthoxml gene id -> (species, gene)
    rows = {exp: {} for exp in sets}  # gene -> (species, hog, root hog)
    hog_ids = {exp: {} for exp in sets}
    young_ids = {exp: {} for exp in sets}
    stack, cur_sp, n_notax = [], None, 0
    for ev, el in ET.iterparse(xml_path, events=("start", "end")):
        tag = el.tag[len(NS):] if el.tag.startswith(NS) else el.tag
        if ev == "start":
            if tag == "species":
                cur_sp = el.get("name")
            elif tag in ("orthologGroup", "paralogGroup"):
                stack.append([tag, el.get("id"), None])
            elif tag == "geneRef":
                hit = genes.get(el.get("id"))
                if hit is None or hit[0] not in by_sp:
                    continue
                for exp in by_sp[hit[0]]:
                    if hit[1] in rows[exp]:
                        sys.exit(f"gene twice in groups: {hit[1]}")
                    for g in reversed(stack):
                        if g[0] != "orthologGroup":
                            continue
                        if g[2] is None:
                            n_notax += 1
                        if g[2] in ok[exp]:
                            ids = hog_ids[exp]
                            h = ids.setdefault(
                                g[1], f"{level[exp]}:{len(ids) + 1}")
                            break
                    else:
                        # root hog born below the level: one young hog
                        ids = young_ids[exp]
                        h = ids.setdefault(
                            stack[0][1], f"{level[exp]}:y{len(ids) + 1}")
                    rows[exp][hit[1]] = (
                        hit[0], h, stack[0][1].split("_")[0], stack[0][2])
            continue
        if tag == "gene" and cur_sp in by_sp:
            genes[el.get("id")] = (cur_sp, el.get("protId").split()[0])
            el.clear()
        elif tag == "gene":
            el.clear()
        elif tag == "property" and el.get("name") == "TaxRange":
            stack[-1][2] = el.get("value")
        elif tag in ("orthologGroup", "paralogGroup"):
            stack.pop()
            if not stack:
                el.clear()
    if n_notax:
        print(f"warning: {n_notax} lookups met a group without TaxRange")

    for exp, sp in sets.items():
        d = os.path.join(root, f"complex_{exp}", "inputs")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "level.txt"), "w") as f:
            f.write(level[exp] + "\n")
        k = labels.index(level[exp])
        ages = ok[exp] + below(k, labels, parent)
        with open(os.path.join(d, "age_order.txt"), "w") as f:
            f.write("\n".join(ages) + "\n")
        r = rows[exp]
        out = os.path.join(d, "orthologs.tsv")
        if os.path.exists(out):
            print(f"{out} exists; kept")
        else:
            with open(out, "w") as f:
                f.write("species\tgene\thog\n")
                for gene, (s, h, _, _) in r.items():
                    f.write(f"{s}\t{gene}\t{h}\n")
        age = {h: (rh, tr) for (_, h, rh, tr) in r.values()}
        with open(os.path.join(d, "hog_age.tsv"), "w") as f:
            f.write("hog\troot_hog\troot_taxrange\n")
            for h, (rh, tr) in age.items():
                f.write(f"{h}\t{rh}\t{tr}\n")
        print(f"\n{exp} at {level[exp]}: {len(r)} genes, "
              f"{len(hog_ids[exp]) + len(young_ids[exp])} hogs in "
              f"{len({v[2] for v in r.values()})} root hogs; by age:")
        n_age = collections.Counter(v[1] for v in age.values())
        for a in ages:
            if n_age.get(a, 0):
                print(f"  {a}\t{n_age[a]}")
        print("species\tgenes_with_hog\thogs\tshare_multi_copy")
        for s in sp:
            cnt = collections.Counter(h for (t, h, _, _) in r.values()
                                      if t == s)
            multi = sum(v > 1 for v in cnt.values())
            share = multi / len(cnt) if cnt else float("nan")
            print(f"{s}\t{sum(cnt.values())}\t{len(cnt)}\t{share:.3f}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
