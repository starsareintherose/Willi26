```
       Willi26 Version Copyright (c) Guoyi Zhang 2026
                     All rights reserved.
         This copy produced for the exclusive use of
                       Every  Cladist.
```

```
assist;         list available commands
assist *;       show full help for all commands
assist xx;      show help for commands starting with 'xx'

batch;          turn on batch switch (DOS legacy; not implemented)
batch -;        turn it off (DOS legacy, not implemented)

bb;             produce multiple trees by branch breaking
bb *;           use all available tree space

bytes;          display bytes of free ram

ccode [opt];    control character coding
                [opt] / set weight [ activate ] deactivate
                [opt] + additive - nonadditive
                [opt] * connect different opts
ccode ;         display codings

cget x;         set coding from code file x

ckeep x;        save current coding in code file x

display;        short listings to display
display -;      no listing to display
display *;      all listing to display

erase  s;       delete tree files in scope s

files;          display directory of treefiles

get x;          make tree file x current

hennig;         calculate single tree
hennig *;       use branch breaking

ie;             find trees by implicit enumeration
ie -;           find just one tree
ie *;           use all available tree space

keep x;         save current tree file as tree file x

log  n;         open dos file n as new log file
log  -;         deactivate log file
log  *;         activate it
log  /;         close it

mhennig;        calculate multiple trees
mhennig *;      use branch breaking

nelsen;         calculate nelson consensus tree

outgroup xx;    control outgroup
outgroup = x;   set outgroup to list
outgroup ;      list outgroup alphabetically

procedure n;    open dos file n as procedure file
procedure -;    deactivate procedure file
procedure *;    activate it
procedure /;    close it

quote;          copy message

reroot;         current treefile according to current outgroup

steps;          display max/min steps per character

tchoose  s;     select trees in scope s from current tree file

tlist;          display trees in parenthetical notation

tplot;          produce tree diagrams

tread;          read trees

tsave  n;       save current tree file on dos file n

txascii;        use extended ascii characters in tree plots
txascii -;      don't

view  n;        inspect dos file n
view  *;        close and inspect current log file

watch;          turn on stopwatch
watch -;        turn it off

xread;          read character data

xsteps;         diagnose trees in current tree file
xsteps h;       list possible states for hypothetical ancestors
xsteps c;       list character fits
xsteps m;       list best/worst fits
xsteps l;       list tree lengths
xsteps u;       produce file of distinct trees
xsteps w;       set character weights according to fits

xx              display and modify diagnosed tree
                number switchs to specific character
                /[]+- follow ccode grammar
                \(branch1) (branch2) move branch
                \\(branch) delete a branch
                = save current settings and exit
                ; exit without saving

yama;           return to dos
```
