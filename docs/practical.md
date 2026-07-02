# Practical guide

This document gives practical workflows for running Willi26.

## Input files

Willi26 reads Hennig86/TNT-style character data using `xread`.

Example:

```text
xread
'SALA -- BOLICK ADV. CLAD. 1:115-125 SALMEA (COMPOSITAE: HELIANTHEAE)'
24 11 
ANCT 000000000000000000000000
SCAN 000000100001210000000000
ORTH 000000100001101000000000
OLIG 100000100002010111000000
PALM 000000100002000110100000
PETR 001100100000000000011000
CALE 100021011110000000000000
GLAB 000200111100000000000001
MONT 321010111100000000000000
INSI 220000011110000000000110
PAUC 210010011110000000000101
 ;
proc / ;
```


## Small datasets

For small datasets, use implicit enumeration:

```text
procedure example.tnt;
ie;
tsave mpts.tre;
nelsen;
tchoose /;
tsave con.tre;
log consensus_plot.txt;
tplot;
log /;
yama;
```

This is the recommended workflow when the dataset is small enough for exact search.

Use:

```text
ie;
```

to find all shortest trees by implicit enumeration.

Use:

```text
ie -;
```

to stop after finding one shortest tree.

Use:

```text
ie *;
```

to use all available tree space.

## Large datasets

For larger datasets, use a heuristic search first:

```text
procedure example.tnt;
mhennig;
bb;
tsave mpts.tre;
nelsen;
tchoose /;
tsave con.tre;
log consensus_plot.txt;
tplot;
log /;
yama;
```

This workflow runs multiple-tree search with `mhennig`, improves the result with branch breaking using `bb`, calculates a Nelson consensus tree with `nelsen`, and saves the tree plot to `consensus_plot.txt`.

## More intensive large-dataset search

If the number of retained trees is still small, for example fewer than 100 trees, and a more extensive search is desired, continue with:

```text
bb *;
```

A more intensive workflow is:

```text
procedure example.tnt;
mhennig;
bb *;
tsave mpts.tre;
nelsen;
tchoose /;
tsave con.tre;
log consensus_plot.txt;
tplot;
log /;
yama;
```

`bb *` uses all available tree space and may take substantially longer than `bb`.

## Saving a consensus tree plot

To save the consensus tree plot:

```text
nelsen;
tchoose /;
log consensus_plot.txt;
tplot;
log /;
```

The output will be written to:

```text
consensus_plot.txt
```

## Listing trees

To list trees in parenthetical notation:

```text
tlist;
```

## Saving trees

To save the current tree set to a file:

```text
tsave  n;
```

where `n` is the file name.

## Diagnosing trees

To diagnose the current tree file:

```text
xsteps;
```

Useful variants include:

```text
xsteps l;   list tree lengths
xsteps w;   set character weights according to fits
```

## Common commands

```text
assist;          list available commands
xread;           read character data
procedure n;     open procedure file n
hennig;          calculate a single tree
mhennig;         calculate multiple trees
ie;              implicit enumeration
bb;              branch breaking
nelsen;          Nelson consensus tree
tlist;           list trees in parenthetical notation
tplot;           plot trees as text diagrams
tsave  n;        save current tree file on dos file n
xsteps;          diagnose trees
log n;           write output to log file n
log /;           close log file
yama;            exit Willi26
```

