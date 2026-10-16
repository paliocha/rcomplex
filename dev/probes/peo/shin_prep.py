"""Build the FastOMA input of the Shin 2021 benchmark.

Usage: python3 shin_prep.py <fastoma dir> <mendeley dir>
Reads <fastoma>/pep/*.pep.all.fa.gz (Ensembl Plants release 38), keeps
the longest protein per gene, writes it as
<fastoma>/in_folder/proteome/<Species>.fa with the bare gene id (the
gene: tag) as header, and writes in_folder/species_tree.nwk. Prints
the share of each compendium's gene ids found in its proteome, and the
share of proteome genes found in the compendium. Exits 1 when both are
below 95 %. The Arabidopsis compendium also holds TAIR10 non-coding
and transposon genes, which have no protein (about 80 % one way, but
nearly every protein-coding gene is in the compendium).
"""
import glob
import gzip
import os
import sys

SPECIES = {
    "ath": "Arabidopsis_thaliana", "mtr": "Medicago_truncatula",
    "osa": "Oryza_sativa", "ppa": "Physcomitrella_patens",
    "stu": "Solanum_tuberosum", "zma": "Zea_mays",
}
TREE = ("((((Medicago_truncatula,Arabidopsis_thaliana)Anc4,"
        "Solanum_tuberosum)Anc3,(Oryza_sativa,Zea_mays)Anc5)Anc2,"
        "Physcomitrella_patens)Anc1;\n")


def read_pep(path):
    """Longest sequence per gene id from an Ensembl pep fasta."""
    best, gene, seq = {}, None, []

    def keep():
        if gene is not None:
            s = "".join(seq)
            if len(s) > len(best.get(gene, "")):
                best[gene] = s
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                keep()
                tags = [t for t in line.split() if t.startswith("gene:")]
                if not tags:
                    sys.exit(f"no gene: tag in {path}: {line.strip()}")
                gene, seq = tags[0][5:], []
            else:
                seq.append(line.strip())
        keep()
    return best


def main(fdir, mdir):
    out = os.path.join(fdir, "in_folder", "proteome")
    os.makedirs(out, exist_ok=True)
    # FastOMA stages hogmap_in (may be empty) and splice (must exist and
    # hold a file of 2+ lines). One protein per gene already, so a dummy
    # splice file maps two genes to themselves, as the PEO run did.
    os.makedirs(os.path.join(fdir, "in_folder", "hogmap_in"), exist_ok=True)
    sd = os.path.join(fdir, "in_folder", "splice")
    os.makedirs(sd, exist_ok=True)
    with open(os.path.join(sd, "Arabidopsis_thaliana.splice"), "w") as f:
        f.write("AT1G01010;AT1G01010\nAT1G01020;AT1G01020\n")
    low = []
    for code, sp in SPECIES.items():
        dirname = sp.lower()
        pep = glob.glob(os.path.join(fdir, "pep",
                                     sp + ".*.pep.all.fa.gz"))
        if len(pep) != 1:
            sys.exit(f"{dirname}: expected one pep file, found {pep}")
        best = read_pep(pep[0])
        with open(os.path.join(out, sp + ".fa"), "w") as f:
            for g, s in best.items():
                f.write(f">{g}\n")
                for i in range(0, len(s), 60):
                    f.write(s[i:i + 60] + "\n")
        with open(os.path.join(mdir, f"{code}_compendium_matrix.txt")) as f:
            next(f)
            ids = [line.split("\t", 1)[0] for line in f]
        share = sum(i in best for i in ids) / len(ids)
        back = len(best.keys() & set(ids)) / len(best)
        print(f"{sp}\tgenes {len(best)}\tcompendium {len(ids)}\t"
              f"in proteome {share:.4f}\tproteome in compendium "
              f"{back:.4f}", flush=True)
        if share < 0.95 and back < 0.95:
            low.append(sp)
    with open(os.path.join(fdir, "in_folder", "species_tree.nwk"), "w") as f:
        f.write(TREE)
    if low:
        sys.exit("below 95 %: " + ", ".join(low))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
