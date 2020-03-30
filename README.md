# SARS-CoV-2

Integrate multi-omics data for SARS-Cov-2.

Data obtained from:

- [SARS-CoV-2 infected host cell proteomics reveal potential therapy targets](DOI:10.21203/rs.3.rs-17218/v1)
  Denisa Bojkova, Kevin Klann, Benjamin Koch, Marek Widera, David Krause, Sandra Ciesek, Jindrich Cinatl, Christian Münch
- [Supp table 1](https://assets.researchsquare.com/files/rs-17218/v1/Supplementary%20Table%2001.xlsx)
- [Supp table 2](https://assets.researchsquare.com/files/rs-17218/v1/Supplementary%20Table%2002.xlsx)

Method reference:

- [DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays](https://doi.org/10.1093/bioinformatics/bty1054) Amrit Singh, Casey P Shannon, Benoît Gautier, Florian Rohart, Michaël Vacher, Scott J Tebbutt, Kim-Anh Lê Cao

Experimental design:

| Design  |               |              |             |             |          |          |           |           |
|---------|---------------|--------------|-------------|-------------|----------|----------|-----------|-----------|
| Analysis 1 | Control       | Virus        |             |             |          |          |           |           |
|         | 12            | 12           |             |             |          |          |           |           |
| Analysis 2 | Control early | Control late | Virus early | Virus late  |          |          |           |           |
|         | 6             | 6            | 6           | 6           |          |          |           |           |
| Analysis 3 | Control_2h    | Control_6h   | Control_10h | Control_24h | Virus_2h | Virus_6h | Virus_10h | Virus_24h |
|         | 3             | 3            | 3           | 3           | 3        | 3        | 3         | 3         |
|         |               |              |             |             |          |          |           |           |

Software libraries involved are:

- python 3.8.2
  - jupyter-notebook 1.0.0
  - pandas 1.0.2

- R 3.6.2
  - argparser 0.6
  - BiocManager 1.30.10
  - mixOmics 6.11.11
