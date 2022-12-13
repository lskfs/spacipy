# stomics tools for custom users
collection of modules for some basic STOmic data processing

## Description of the Gem Format

Gem files are usually four columns, tab-delimited, plain text files generated from STOmics (stereo-seq: https://www.stomics.tech/) analysis workflow. The file contents include the coordination of DNB (the basic unit of stereo-seq), as well as gene name found in the spatial location and the expression number of this gene.

    geneID  x       y       MIDCounts
    sox2    6232    8293    3
    Cr2     15949   9839    1
    Cr1l    11011   2334    1
    Cd46    7149    9671    1
    Cd46    22900   9266    2
    Plxna2  13842   18100   1
    Kcnh1   19118   1036    1

## gem.py
#### import Gem class and read in the gem file
```python
from .gem import Gem
gem = Gem.readin('./data/example.gem.gz')
```
