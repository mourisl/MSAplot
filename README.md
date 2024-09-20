MSAplot
======

### What is it?
A Python package to visualize multiple sequence alignment (MSA). MASplot is inspired by the package [pyMSAviz](https://github.com/moshi4/pyMSAviz) and [ggmsa](https://yulab-smu.top/ggmsa/). The major feature of MSAplot is customizable panels for plotting MSA, seqlogo, annotation, and consensus histogram. 

### Installation

1. Clone the [GitHub repo](https://github.com/mourisl/msaplot), e.g. with `git clone https://github.com/mourisl/msaplot`.
2. Copy "msaplot" folder to your project folder.

I will try to add msaplot to PyPi in future.

### Usage
Here is a minimal example:

```python
import matplotlib.pyplot as plt
from msaplot import msaplot

msaplot.DrawComplexMSA(msa=["AC-GAT", "A-CGT-"],
               panels=[msaplot.DrawSeqLogo, msaplot.DrawMSA, msaplot.DrawConsensusHisto])
```

There are four drawinging functions that can be used in the panels option:

- DrawMSA: visualize multiple sequence alignment
- DrawSeqLogo: represent the seqlogo of the alignment
- DrawConsensusHisto: represent the fraction of the consensus character for each position
- DrawAnnotation: need to take a list to specify the name for each annotated region and region coordinate.

See more examples in [example.ipynb](https://github.com/mourisl/msaplot/blob/main/example.ipynb).

### Requirements
+ Python >= 3.6
+ seaborn >= 0.9
+ matplotlib >= 2.2.2
