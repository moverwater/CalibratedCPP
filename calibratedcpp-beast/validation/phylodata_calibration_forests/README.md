# PhyloData Calibration Forest Constraint Trees

Calibration constraint trees extracted from BEAST2 XML files in the PhyloData repository.
These represent real-world calibration structures used in published phylogenetic analyses.

## Format

Annotated Newick format with BEAST-style metadata:
```
((A,B)[&name=Clade1,dist=LogNormal,M=5.2,S=0.5,offset=33.9],(C,D)[&name=Clade2])[&name=Root];
```

## Files

| File | Taxa | Calibrations | Source |
|------|------|--------------|--------|
| burridge-2020-*.newick | 725 | 17 | Galaxioid fish diversification |
| harrington-2024-*.newick | 5940 | 22 | Large dispersal analysis |
| jorge-2018-*.newick | 114 | 7 | - |
| marburger-2018-*.newick | 1105 | 7 | - |
| miri-2024-*.newick | 202 | 5 | Zygoptera |
| presslee-2019-*.newick | 108 | 16 | Xenarthra palaeoproteomics |
| saladin-2017-*.newick | 230 | 12-15 | Fossil calibrations |
| serrano-2024-*.newick | 314 | 7 | - |
| trotta-2018-*.newick | 1080 | 50 | Pine Rockland community |

## Usage

Run benchmark:
```bash
cd calibratedcpp-beast
java -cp "lib/*:target/classes" calibratedcpp.SimpleConstraintBenchmark [iterations] [warmup]
```

## Notes

Some files have overlapping (non-nested) calibrations that are incompatible with
CalibratedCoalescentPointProcess, which requires calibrations to be nested or disjoint.
These files will produce errors during benchmarking.

## Source

Generated from `calibrated-time-tree-dist` analysis of PhyloData BEAST2 XML files.
