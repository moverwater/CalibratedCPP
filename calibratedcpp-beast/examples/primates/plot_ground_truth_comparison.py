#!/usr/bin/env python3
"""
Compare posterior node-age estimates from BEAST runs against the independent dos Reis et
al. calibration ground truth (data/dosReis.csv), for each of the primates_remove8,
primates_10pct_temp_node5, and primates_10pct_temp_severalCali experiments.

Works entirely from each run's BEAST output (.trees + .log) — no XML or LPhy script is
read. Each calibrated node is identified directly from the log's mrca.age(TaxonSetN)
columns; its identity (taxa membership, and a representative tipA/tipB pair for
labelling) comes from locating that same node, by height, in the first sampled tree.
Calibrated clades are monophyletic for the whole chain (a hard MRCAPrior constraint), so
identifying membership from one tree is valid for every sample.

Each dos Reis clade is matched to a calibrated node by exact taxa-set equality (dos
Reis's tipsA | tipsB must equal the node's full leaf set). A looser "smallest common
ancestor + bipartition" check was tried and rejected: the first sampled tree is an
early/unconverged state, so its *uncalibrated* splits are close to arbitrary, and that
looser check produced false matches (e.g. Catarrhini spuriously matching Primates).
Truncated ("+N more") dos Reis rows can't be verified this way and are skipped.

Usage:
    python plot_ground_truth_comparison.py

Expects (relative to this script's directory):
    data/dosReis.csv
    data/<run>.trees or data/<run>.trees  (BEAST posterior trees)
    data/<run>.txt                                   (BEAST .log, tab-separated)
"""
import csv
import os
import re

import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(BASE, "data", "dosReis.csv")

RUNS = [
    {
        "name": "primates_remove8",
        "trees": os.path.join(BASE, "data", "primates_remove8.trees"),
        "log": os.path.join(BASE, "data", "primates_remove8.txt"),
        "out": os.path.join(BASE, "primates_remove8_ground_truth_comparison.pdf"),
    },
    {
        "name": "primates_10pct_temp_node5",
        "trees": os.path.join(BASE, "data", "primates_10pct_temp_node5.trees"),
        "log": os.path.join(BASE, "data", "primates_10pct_temp_node5.txt"),
        "out": os.path.join(BASE, "primates_10pct_temp_node5_ground_truth_comparison.pdf"),
    },
    {
        "name": "primates_10pct_temp_severalCali",
        "trees": os.path.join(BASE, "data", "primates_10pct_temp_severalCali.trees"),
        "log": os.path.join(BASE, "data", "primates_10pct_temp_severalCali.txt"),
        "out": os.path.join(BASE, "primates_10pct_temp_severalCali_ground_truth_comparison.pdf"),
    },
]

BURNIN_FRACTION = 0.1
GROUND_TRUTH_COLOR = "#e15759"

# Node numbers as in the de Vries & Beck fossil calibration table
# (~/Downloads/primates/primates.caibrationPriors.pdf, Table 1), keyed by each node's
# exact taxon membership in this 29-taxon dataset. Used only for labelling plot panels
# with the paper's own node number instead of an arbitrary sequential count — matching
# by taxa (not by dos Reis's own clade name) matters because dos Reis's "Cercopithecidae"
# row is actually the 8-taxon Cercopithecinae grouping, not the true 12-taxon
# Cercopithecidae (see match_node_to_dosreis's docstring for the exact-match rationale).
_EUARCHONTOGLIRES = frozenset({
    "Tupaia_chinensis", "Galeopterus_variegatus", "Otolemur_garnettii", "Microcebus_murinus",
    "Propithecus_coquereli", "Carlito_syrichta", "Cebus_capucinus_imitator", "Saimiri_boliviensis",
    "Aotus_nancymaae", "Callithrix_jacchus", "Nomascus_leucogenys", "Pongo_abelii", "Gorilla_gorilla",
    "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes", "Colobus_angolensis_palliatus",
    "Piliocolobus_tephrosceles", "Rhinopithecus_bieti", "Rhinopithecus_roxellana",
    "Chlorocebus_sabaeus", "Macaca_nemestrina", "Macaca_fascicularis", "Macaca_mulatta",
    "Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada",
    "Mus_musculus",
})
_EUARCHONTA = _EUARCHONTOGLIRES - {"Mus_musculus"}
_PRIMATES = _EUARCHONTA - {"Tupaia_chinensis", "Galeopterus_variegatus"}
_STREPSIRRHINI = frozenset({"Otolemur_garnettii", "Microcebus_murinus", "Propithecus_coquereli"})
_HAPLORRHINI = _PRIMATES - _STREPSIRRHINI
_ANTHROPOIDEA = _HAPLORRHINI - {"Carlito_syrichta"}
_CERCOPITHECIDAE = frozenset({
    "Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
    "Rhinopithecus_roxellana", "Chlorocebus_sabaeus", "Macaca_nemestrina", "Macaca_fascicularis",
    "Macaca_mulatta", "Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis",
    "Theropithecus_gelada",
})
_COLOBINAE = frozenset({
    "Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
    "Rhinopithecus_roxellana",
})
_CERCOPITHECINAE = _CERCOPITHECIDAE - _COLOBINAE
_PAPIONINI = frozenset({
    "Macaca_nemestrina", "Macaca_fascicularis", "Macaca_mulatta", "Cercocebus_atys",
    "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada",
})
_HOMINOIDEA = frozenset({
    "Nomascus_leucogenys", "Pongo_abelii", "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus",
    "Pan_troglodytes",
})
_HOMINIDAE = _HOMINOIDEA - {"Nomascus_leucogenys"}
_HOMO_PAN = frozenset({"Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"})
_CATARRHINI = _CERCOPITHECIDAE | _HOMINOIDEA
_CALLITRICHIDAE_CEBIDAE = frozenset({
    "Cebus_capucinus_imitator", "Saimiri_boliviensis", "Aotus_nancymaae", "Callithrix_jacchus",
})
_CEBIDAE = frozenset({"Cebus_capucinus_imitator", "Saimiri_boliviensis"})

NODE_NUMBER_BY_TAXA = {
    _EUARCHONTOGLIRES: 1,
    _EUARCHONTA: 3,
    _PRIMATES: 5,
    _STREPSIRRHINI: 6,
    _HAPLORRHINI: 10,
    _ANTHROPOIDEA: 11,
    _CATARRHINI: 12,
    _CERCOPITHECIDAE: 13,
    _COLOBINAE: 14,
    _CERCOPITHECINAE: 15,
    _PAPIONINI: 16,
    _HOMINOIDEA: 18,
    _HOMINIDAE: 19,
    _HOMO_PAN: 20,
    _CALLITRICHIDAE_CEBIDAE: 23,
    _CEBIDAE: 24,
}


class Node:
    __slots__ = ("name", "children", "branch_length", "height")

    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.branch_length = 0.0
        self.height = None

    def is_leaf(self):
        return not self.children

    def leaf_names(self):
        if self.is_leaf():
            return frozenset({self.name})
        out = frozenset()
        for c in self.children:
            out |= c.leaf_names()
        return out

    def representative_leaf(self):
        node = self
        while not node.is_leaf():
            node = node.children[0]
        return node.name


def parse_translate_table(trees_path):
    with open(trees_path) as f:
        content = f.read()
    m = re.search(r"Translate\s*(.*?);", content, re.DOTALL)
    table = {}
    for entry in m.group(1).split(","):
        entry = entry.strip()
        if not entry:
            continue
        num, name = entry.split(None, 1)
        table[num.strip()] = name.strip()
    return table


def parse_newick_string(newick, translate):
    newick = re.sub(r"\[&[^\]]*\]", "", newick)  # strip [&rate=...] annotations

    pos = 0

    def parse_node():
        nonlocal pos
        node = Node()
        if newick[pos] == "(":
            pos += 1
            node.children.append(parse_node())
            while newick[pos] == ",":
                pos += 1
                node.children.append(parse_node())
            assert newick[pos] == ")"
            pos += 1
        else:
            start = pos
            while newick[pos] not in ",(): ;":
                pos += 1
            label = newick[start:pos]
            node.name = translate.get(label, label)
        if pos < len(newick) and newick[pos] == ":":
            pos += 1
            start = pos
            while newick[pos] not in ",();":
                pos += 1
            node.branch_length = float(newick[start:pos])
        return node

    root = parse_node()
    compute_heights(root)
    return root


def parse_first_tree(trees_path, translate):
    with open(trees_path) as f:
        content = f.read()
    m = re.search(r"tree\s+STATE_\d+\s*=\s*(?:\[&R\]\s*)?(\(.*?;)", content, re.DOTALL)
    return parse_newick_string(m.group(1), translate)


def iter_all_trees(trees_path, translate):
    """Yield the root Node of every posterior tree in the file, in file order."""
    with open(trees_path) as f:
        content = f.read()
    for m in re.finditer(r"tree\s+STATE_\d+\s*=\s*(?:\[&R\]\s*)?(\(.*?;)", content, re.DOTALL):
        yield parse_newick_string(m.group(1), translate)


def find_exact_clade(root, target_leaves):
    """Return the node whose leaf set exactly equals target_leaves, or None if
    target_leaves isn't monophyletic in this particular tree."""
    stack = [root]
    while stack:
        node = stack.pop()
        if not node.is_leaf():
            if node.leaf_names() == target_leaves:
                return node
            stack.extend(node.children)
    return None


def collect_ages_across_trees(post_burnin_roots, target_leaves):
    """For a clade with no logged mrca.age column, get its age posterior by checking, in
    every post-burnin tree, whether target_leaves forms a monophyletic clade — if so,
    recording that node's height. Returns (ages, n_trees_checked); ages may be shorter
    than n_trees_checked if the clade isn't monophyletic in every sample (since it's not
    under a hard MRCAPrior constraint here)."""
    ages = [node.height for root in post_burnin_roots
            if (node := find_exact_clade(root, target_leaves)) is not None]
    return np.array(ages), len(post_burnin_roots)


def compute_heights(node):
    """Post-order: leaf height = 0, internal height = child's height + child's branch length
    (tree is ultrametric — all children give the same answer, so just use the first)."""
    if node.is_leaf():
        node.height = 0.0
        return 0.0
    child_height = compute_heights(node.children[0]) + node.children[0].branch_length
    for c in node.children[1:]:
        compute_heights(c)
    node.height = child_height
    return node.height


def find_node_by_height(root, target_height, tol=1e-4):
    best = None
    best_diff = None
    stack = [root]
    while stack:
        node = stack.pop()
        if not node.is_leaf():
            diff = abs(node.height - target_height)
            if diff < tol and (best_diff is None or diff < best_diff):
                best, best_diff = node, diff
            stack.extend(node.children)
    return best


def parse_dos_reis(csv_path):
    rows = []
    with open(csv_path, encoding="latin-1", newline="") as f:
        reader = csv.reader(f)
        header = [h.strip().lower() for h in next(reader)]
        idx = {name: header.index(name) for name in ("clade", "lower", "upper", "tipsa", "tipsb")}
        for r in reader:
            if len(r) < len(header):
                continue
            clade = r[idx["clade"]]
            lower, upper = r[idx["lower"]], r[idx["upper"]]
            tips_a, tips_b = r[idx["tipsa"]], r[idx["tipsb"]]
            if "not represented" in tips_a.lower():
                continue
            if upper.strip().lower() in ("n/s", "", "na"):
                continue

            truncated = "(+" in tips_a or "(+" in tips_b

            def clean(s):
                s = re.sub(r",?\s*\.\.\.\s*\(\+\d+ more\)", "", s)
                return frozenset(t.strip() for t in s.split(",") if t.strip() and t.strip() != "...")

            a, b = clean(tips_a), clean(tips_b)
            rows.append({"clade": clade, "lower": float(lower), "upper": float(upper),
                         "tipsA": a, "tipsB": b, "truncated": truncated})
    return rows


def match_node_to_dosreis(calibrated_nodes, dos_reis_rows):
    """Return {node_id(python id of Node): dos_reis_row}.

    Only the nodes BEAST actually calibrates (monophyletic for the entire chain — passed
    in via calibrated_nodes) are considered as match candidates. A row matches if its full
    taxa (tipsA | tipsB) exactly equal that node's leaf set. Exact equality — rather than
    a looser "these two tips end up on opposite sides of some node" check — matters here:
    the STATE_0 tree used for node identification is an early/unconverged sample, so
    UNCALIBRATED splits in it are close to arbitrary. A looser check can spuriously match
    an uncalibrated dos Reis clade (e.g. Catarrhini) to an unrelated calibrated node whose
    two children happen, by chance in this one tree, to separate the same two tips.
    Requiring exact equality against a node whose monophyly is a hard constraint for the
    whole chain avoids that failure mode; truncated ("+N more") dos Reis rows can't be
    verified this way and are skipped rather than guessed at.
    """
    matches = {}
    for row in dos_reis_rows:
        if row["truncated"]:
            continue
        target = row["tipsA"] | row["tipsB"]
        for node in calibrated_nodes:
            if node.leaf_names() == target:
                matches[id(node)] = row
                break
    return matches


def read_log_columns(log_path, columns, burnin_fraction):
    with open(log_path) as f:
        header = f.readline().strip().split("\t")
    col_idx = {c: header.index(c) for c in columns if c in header}
    missing = [c for c in columns if c not in col_idx]
    if missing:
        raise ValueError(f"Columns not found in log file: {missing}")

    data = {c: [] for c in columns}
    with open(log_path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            for c in columns:
                data[c].append(float(parts[col_idx[c]]))

    n = len(next(iter(data.values())))
    burnin = int(n * burnin_fraction)
    return {c: np.array(v[burnin:]) for c, v in data.items()}


def run_analysis(name, trees_path, log_path, out_pdf, dos_reis_rows):
    print(f"=== {name} ===")
    translate = parse_translate_table(trees_path)
    root = parse_first_tree(trees_path, translate)

    with open(log_path) as f:
        header = f.readline().strip().split("\t")
        row0 = dict(zip(header, f.readline().strip().split("\t")))

    mrca_cols = [c for c in header if c.startswith("mrca.age(")]

    calibrated_nodes = {}
    for col in mrca_cols:
        target_height = float(row0[col])
        node = find_node_by_height(root, target_height)
        if node is None:
            print(f"WARNING: could not locate tree node for {col} (height={target_height})")
            continue
        calibrated_nodes[col] = node

    matches_by_node_id = match_node_to_dosreis(list(calibrated_nodes.values()), dos_reis_rows)

    calibrated_cols_needed = [col for col in mrca_cols
                              if calibrated_nodes.get(col) is not None
                              and id(calibrated_nodes[col]) in matches_by_node_id]
    log_data = read_log_columns(log_path, calibrated_cols_needed, BURNIN_FRACTION) if calibrated_cols_needed else {}

    plot_entries = []  # (node_number, label, row, samples, source_note)
    for col in mrca_cols:
        node = calibrated_nodes.get(col)
        if node is None:
            continue
        row = matches_by_node_id.get(id(node))
        if row is None:
            continue
        tip_a, tip_b = node.children[0].representative_leaf(), node.children[1].representative_leaf()
        label = f"{tip_a} / {tip_b}"
        node_number = NODE_NUMBER_BY_TAXA.get(node.leaf_names())
        plot_entries.append((node_number, label, row, log_data[col], "calibrated"))

    # For dos Reis clades that aren't one of this model's calibrated (and hence always
    # monophyletic) nodes, there's no mrca.age log column to read — but we can still
    # recover a posterior by checking every sampled tree for whether that exact taxon
    # set happens to be monophyletic there, and using its height where it is.
    matched_rows = {id(row) for row in matches_by_node_id.values()}
    unmatched_rows = [r for r in dos_reis_rows if not r["truncated"] and id(r) not in matched_rows]
    tree_scan_entries = []
    if unmatched_rows:
        all_roots = list(iter_all_trees(trees_path, translate))
        burnin = int(len(all_roots) * BURNIN_FRACTION)
        post_burnin_roots = all_roots[burnin:]
        for row in unmatched_rows:
            target = row["tipsA"] | row["tipsB"]
            ages, n_trees = collect_ages_across_trees(post_burnin_roots, target)
            if len(ages) == 0:
                continue
            tip_a = sorted(row["tipsA"])[0]
            tip_b = sorted(row["tipsB"])[0]
            label = f"{tip_a} / {tip_b}"
            node_number = NODE_NUMBER_BY_TAXA.get(target)
            note = f"tree scan, monophyletic in {len(ages)}/{n_trees} trees"
            tree_scan_entries.append((node_number, label, row, ages, note))

    def node_label(node_number):
        return f"Node {node_number}" if node_number is not None else "Node ?"

    print(f"{len(plot_entries)} calibrated nodes matched a dos Reis clade "
          f"(via mrca.age log columns):")
    for node_number, label, row, _, _ in plot_entries:
        print(f"  {node_label(node_number)}: {label} -> dos Reis '{row['clade']}' [{row['lower']}, {row['upper']}]")

    print(f"{len(tree_scan_entries)} more dos Reis clades recovered by scanning trees "
          f"directly (not an explicit calibration):")
    for node_number, label, row, ages, note in tree_scan_entries:
        print(f"  {node_label(node_number)}: {label} -> dos Reis '{row['clade']}' [{row['lower']}, {row['upper']}] ({note})")

    plot_entries += tree_scan_entries
    plot_entries.sort(key=lambda e: (e[0] is None, e[0]))

    if not plot_entries:
        print("No matched nodes to plot; skipping PDF.\n")
        return

    n = len(plot_entries)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, (node_number, label, row, samples, source) in zip(axes, plot_entries):
        ax.hist(samples, bins=50, density=True, color="#4e79a7", alpha=0.7, label="Posterior estimate")
        # Thin marker bar at the bottom of the panel instead of a full-height overlay,
        # so it doesn't visually dominate the posterior histogram.
        ax.axvspan(row["lower"], row["upper"], ymin=0, ymax=0.05,
                   color=GROUND_TRUTH_COLOR, alpha=0.9, label="dos Reis ground truth")
        ax.axvline(row["lower"], color=GROUND_TRUTH_COLOR, linestyle=":", linewidth=1, alpha=0.6)
        ax.axvline(row["upper"], color=GROUND_TRUTH_COLOR, linestyle=":", linewidth=1, alpha=0.6)
        prefix = f"{node_label(node_number)}" if source == "calibrated" else f"{node_label(node_number)} (uncalibrated)"
        ax.set_title(f"{prefix}: {label}", fontsize=9)
        ax.set_xlabel("Age (Ma)")
        ax.set_yticks([])
        ax.legend(fontsize=7, loc="upper right")

    for ax in axes[len(plot_entries):]:
        ax.axis("off")

    fig.suptitle(f"{name}: posterior node ages vs. dos Reis calibration ground truth", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Wrote {out_pdf}\n")


def main():
    dos_reis_rows = parse_dos_reis(CSV_PATH)
    for run in RUNS:
        run_analysis(run["name"], run["trees"], run["log"], run["out"], dos_reis_rows)


if __name__ == "__main__":
    main()
