# SARS-CoV-2

Integrate multi-omics data for SARS-Cov-2.

Data obtained from:

- [SARS-CoV-2 infected host cell proteomics reveal potential therapy targets](DOI:10.21203/rs.3.rs-17218/v1)
  Denisa Bojkova, Kevin Klann, Benjamin Koch, Marek Widera, David Krause, Sandra Ciesek, Jindrich Cinatl, Christian Münch
- [Supp table 1](https://assets.researchsquare.com/files/rs-17218/v1/Supplementary%20Table%2001.xlsx)
- [Supp table 2](https://assets.researchsquare.com/files/rs-17218/v1/Supplementary%20Table%2002.xlsx)

Method reference:

- [DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays](https://doi.org/10.1093/bioinformatics/bty1054) Amrit Singh, Casey P Shannon, Benoît Gautier, Florian Rohart, Michaël Vacher, Scott J Tebbutt, Kim-Anh Lê Cao

Experimental design (n=3 for all classes):

| Infection state | Timepoint (hour) |
|-----------------|------------------|
| Control         | 2                |
| Control         | 6                |
| Control         | 10               |
| Control         | 24               |
| Infected        | 2                |
| Infected        | 6                |
| Infected        | 10               |
| Infected        | 24               |

Software libraries involved are:

```
python 3.8.2
  jupyter-notebook 1.0.0
  pandas 1.0.2

R 3.6.2
  argparser 0.6
  igraph_1.2.5
  mixOmics 6.13.3
```

Acknowledgements:

[We thank Kim-Anh Lê Cao](https://lecao-lab.science.unimelb.edu.au/) for contributions to the code and analysis.
