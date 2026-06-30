# The Constraint Tree Format

**Reference implementation:** [ConstraintTree.java](ConstraintTree.java)
**Consumers:** [CalibrationPrior.java](../calibrationprior/CalibrationPrior.java) (`constraintTree` input), [CalibratedCoalescentPointProcess.java](../calibratedcpp/CalibratedCoalescentPointProcess.java) (`constraintTree` input)

## 1. Purpose

A constraint tree is a compact, single-string alternative to declaring calibration
clades one at a time as explicit `<taxonset>`/`<calibration>` XML blocks. It is an
**extended Newick** string in which internal nodes may carry a metadata annotation
block giving the node a name and/or a soft age interval. The tree's topology is used
only to derive **taxon sets** (which leaves are grouped under which internal node);
branch lengths and the overall tree shape carry no other biological or statistical
meaning to the parser.

A constraint tree is not itself a probability model — it is a data-entry convenience.
Given a `newick` string, `ConstraintTree` produces:

- a `TaxonSet` for every internal node that is named and/or bears bounds,
- a `CalibrationCladePrior` for every internal node that bears **both** a lower and
  an upper bound,
- the flat set of all leaf taxon names appearing in the tree.

The distribution family used to model each calibrated age (log-normal, Beta, etc.)
is **not** part of this format — that is a downstream modelling choice made by
`CalibrationPrior`.

## 2. Grammar

```
constraint-tree  ::= WS* [ subtree ] WS* [ ';' ] WS*
                    | WS*                                   ; empty ⇒ no calibrations

subtree          ::= internal | leaf

internal         ::= '(' WS* branch-list WS* ')' label? annotation? branch-length?
branch-list      ::= subtree ( WS* ',' WS* subtree )*

leaf             ::= taxon-name annotation? branch-length?

label            ::= any run of characters not in { ':' '[' ',' ')' ';' }
                      ; equivalent to standard Newick internal-node labels;
                      ; overridden by an explicit name=... key inside annotation, if present

taxon-name       ::= any run of characters not in { ',' ')' ':' '[' }

annotation       ::= '[&' key-value-list ']'
key-value-list   ::= key-value ( ',' key-value )*
key-value        ::= key '=' value
key              ::= /[^=,]+/            ; matched case-insensitively, trimmed
value            ::= /[^,\]]+/           ; trimmed

branch-length    ::= ':' /[^,)\[;]*/     ; parsed and discarded (no semantic effect)

WS               ::= whitespace (space, tab, newline)
```

Comment lines are stripped before parsing: any input line whose trimmed content
starts with `#` is removed in its entirety, and all remaining lines are
concatenated (see `stripComments`). `#` is therefore only a line-comment marker
when it is the first non-whitespace character on a line — it is not recognized
mid-line.

### 2.1 Recognized annotation keys

Keys are matched case-insensitively after trimming. Unrecognized keys are
silently ignored (this includes any distribution-type hint — see §1).

| Key | Aliases | Type | Effect |
|---|---|---|---|
| `name` | — | string | Sets/overrides the node's ID. Used as the `TaxonSet` ID and, if bounds are also present, the `CalibrationCladePrior` ID. |
| `lower` | `lowerAge` | double | Soft lower age bound for the clade's MRCA. |
| `upper` | `upperAge` | double | Soft upper age bound for the clade's MRCA. |
| `virtualRoot` | — | boolean (`true`/`false` via `Boolean.parseBoolean`) | Marks the node as a non-calibrating organizational root (see §4). |

Values for `lower`/`upper` that fail to parse as a `double` are silently ignored
(the bound is left unset, `NaN`), rather than raising a parse error.

### 2.2 Leaf vs. internal-node labels

- A **leaf**'s label (the text before `:` / `[` / `,` / `)`) is always the taxon
  name and is added to the tree's flat taxon set.
- An **internal node**'s label — the standard bare Newick label written directly
  after the closing `)` — is captured as a candidate name, but is overridden by
  an explicit `name=...` key inside `[&...]` if one is present.
- Internal nodes with **no** label and **no** `name=` annotation and **no**
  bounds are structural only: they group taxa but produce neither a `TaxonSet`
  nor a `CalibrationCladePrior`.

## 3. Semantics — what each node produces

For every internal node `N` (excluding the virtual root, see §4), let
`taxa(N)` be the set of leaf names in `N`'s subtree.

| Node has... | Produces |
|---|---|
| neither `name` nor bounds | nothing (pure grouping node) |
| `name` only | a `TaxonSet` with ID = `name`, taxa = `taxa(N)` |
| `lower` **and** `upper` (name optional) | a `TaxonSet` (auto-named if `name` absent — BEAST assigns a default ID) **and** a `CalibrationCladePrior` with the *same* `TaxonSet` instance, `lowerAge = lower`, `upperAge = upper`, and ID = `name` if given |
| only one of `lower`/`upper` | bounds are incomplete (`hasCalibrationBounds()` is false since `NaN` fails the check) → treated as if no bounds were given; only a `TaxonSet` is produced if `name` is present, otherwise nothing |

Leaves never produce a `TaxonSet` or `CalibrationCladePrior` regardless of
annotation — only internal (clade) nodes do.

The `CalibrationCladePrior` for a node reuses the exact same `TaxonSet` object
returned in `getTaxonSets()` for that node (verified by
`priorSharesTaxonSetWithTaxonSetList` in the test suite) — they are not
independently constructed.

Branch lengths (`:1.0`, `:0.5`, ...) are parsed only to be skipped; they have
no effect on output.

### 3.1 Downstream nesting/disjointness constraint

`ConstraintTree` itself does not enforce any relationship between the taxon
sets it emits. However, the flat list of `CalibrationCladePrior`s it produces
is subsequently consumed by `CalibrationForest.fromPriors(...)`
([CalibrationForest.java](CalibrationForest.java)), which **does** validate
that every pair of calibration clades is either fully nested (one's taxa are
a subset of the other's) or fully disjoint — a partial overlap between two
clades' taxa, or two clades with identical taxa, both throw
`IllegalArgumentException`. Because a constraint tree's topology guarantees
this by construction (any two subtrees of a tree are automatically nested or
disjoint), this failure mode only matters if calibrations from a
`constraintTree` are mixed with — or a `ConstraintTree`'s output is
combined with — an independently supplied list of clades.

## 4. The virtual root

`[&virtualRoot=true]` on a node (conventionally the outermost node of the
string) marks that node as **excluded** from both the `TaxonSet` list and the
`CalibrationCladePrior` list, even if it also carries `name`/`lower`/`upper`.
This lets the outermost parenthesis simply bundle several independent,
non-overlapping calibrated clades into one Newick string without being
mistaken for a calibration on the whole tree:

```
((A,B)[&name=Clade1,lower=10.0,upper=20.0],
 (C,D)[&name=Clade2,lower=5.0,upper=15.0])[&virtualRoot=true];
```

Here `Clade1` and `Clade2` each yield a `TaxonSet` + `CalibrationCladePrior`,
the outer node yields neither, and `getAllTaxa()` still returns `{A, B, C, D}`.

Absent an explicit `virtualRoot=true`, the outermost node is treated like any
other internal node — if it has `name` and/or bounds, it produces output
exactly as described in §3.

## 5. Empty input

If `newickInput` is the empty string, or trims to `""` or `";"`, the tree is
treated as having **no calibrations at all**: `getTaxonSets()`,
`getCalibrationCladePriors()`, and `getAllTaxa()` all return empty
collections. This is the default value of the `newick` input, so omitting a
constraint tree entirely is equivalent to supplying one with no calibrations.

## 6. Worked examples

Single calibrated clade:
```
(A,B)[&name=Root,lower=10.0,upper=20.0];
```
→ one `TaxonSet` `Root = {A, B}`, one `CalibrationCladePrior` `Root` with bounds `[10, 20]`.

Nested calibrations (inner clade's taxa are a subset of the outer's):
```
((A,B)[&name=Inner,lower=5.0,upper=10.0],C)[&name=Outer,lower=15.0,upper=30.0];
```
→ `Inner = {A, B}` bounded `[5, 10]`; `Outer = {A, B, C}` bounded `[15, 30]`.
`CalibrationPrior` treats `Inner` as nested inside `Outer` and fits the joint
age distribution accordingly (a separate concern from this format — see
`CalibrationPrior.initAndValidate`).

Named clade with no bounds (topology-only constraint, e.g. for use as a plain
`TaxonSet` elsewhere, or for future MRCA-only conditioning):
```
(A,B)[&name=TopClade];
```
→ one `TaxonSet` `TopClade = {A, B}`; zero `CalibrationCladePrior`s.

Multiple independent calibrations bundled under a virtual root:
```
((A,B)[&name=Clade1,lower=10.0,upper=20.0],
 (C,D)[&name=Clade2,lower=5.0,upper=15.0])[&virtualRoot=true];
```
→ two `TaxonSet`/`CalibrationCladePrior` pairs, no calibration on the root itself.

Alternate key spelling and ignored branch lengths/comments:
```
# comment lines starting with '#' are stripped
((A:1.0,B:2.0):0.5)[&name=Cl,lowerAge=8.0,upperAge=16.0];
```
→ `Cl = {A, B}` bounded `[8, 16]`; all `:length` tokens discarded.

## 7. Usage in BEAST2 XML

`ConstraintTree` is a `BEASTObject` with a single string input, `newick`, and is
wired into the two prior classes that accept it as an alternative to an explicit
list of clades/taxon sets:

```xml
<constraintTree spec="calibration.ConstraintTree" id="myConstraints"
    newick="((A,B)[&amp;name=Clade1,lower=10.0,upper=20.0],
             (C,D)[&amp;name=Clade2,lower=5.0,upper=15.0])[&amp;virtualRoot=true];"/>

<distribution id="calibrationPrior" spec="calibrationprior.CalibrationPrior"
    tree="@tree" constraintTree="@myConstraints"/>
```

or, for the tree prior itself:

```xml
<distribution id="treePrior" spec="calibratedcpp.CalibratedBirthDeathSkylineModel"
    tree="@tree" constraintTree="@myConstraints" .../>
```

`&` must be XML-escaped as `&amp;` inside the `newick` attribute since `[&...]`
is standard NEXUS/BEAST metadata-annotation syntax reused here. In both
consumer classes it is invalid to supply both `constraintTree` and the
corresponding explicit list (`calibration` / `calibrations`) simultaneously —
exactly one of the two should be provided.

The BEAUti editor ([CalibratedCPPInputEditor.java](../calibratedcpp/beauti/CalibratedCPPInputEditor.java))
provides an Import/Export "Constraint Tree (Newick)" dialog that round-trips
exactly this grammar to/from its calibration-clade table, with the in-app hint:
*"Annotations: `[&name=Label,lower=X,upper=Y]` — name and bounds are optional."*

## 8. A related but non-canonical extended grammar

The 93 example files under
`calibratedcpp-beast/validation/calibratedcpp/phylodata_calibration_forests/`,
together with the standalone parser in
[SimpleConstraintBenchmark.java](../../../../../../test/java/calibratedcpp/SimpleConstraintBenchmark.java)
(a performance-benchmark harness, not part of the production model), use a
**superset** of this annotation vocabulary, additionally recognizing `dist`,
`M`, `S`, `offset`, `monophyletic`, and `meanInRealSpace` — metadata preserved
from the original BEAST2 `MRCAPrior` XML these files were converted from, used
to describe a fitted log-normal distribution in real space. This extended
vocabulary is **not** understood by the production `calibration.ConstraintTree`
parser: any such keys are simply ignored as unrecognized. This is intentional
— per §1, distribution shape is a modelling choice made by `CalibrationPrior`,
not data the constraint tree format carries — so these files should be treated
as legacy/benchmark inputs, not as valid `constraintTree=` values for
production XML.
