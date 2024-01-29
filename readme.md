# MelanoDB
This repository is a tool to charge MelanoDB, a database collection of patients with metastatic melanoma treated with MAPK inhibitors. We formatted data from 8 different studies for a total of 420 cases to gather common clinical and molecular features. Whole or partial exome sequencing is available for 191 cases and gene expression for 132 cases. Web application to explore the database is available is this repository: . We whare this dataset to the scientific community according to FAIR principles <br />

MelanoDB created using python 3.9

### Dependencies:

```peewee``` version 3.13.3 <br />
```starlette``` version 0.13.4 <br />
```playhouse``` <br />
```numpy``` <br />
```ujson``` <br />

### Installation

#### Mac and Linux
Create a virtual environment in melano-py root and load required package:
```bash env.sh```

#### Windows

Create a virtual environment in melano-py root:
```python -m virtualenv .venv```

Activate the virtual environment:
```.\.venv\Scripts\activate```

Install required packages:
```pip -m install -r requirements.txt```

### Class diagram
![Alt text](melanoDB_relational_scheme.png)


### Sources
This database has been created using present article:

\cite{blateauTERTPromoterMutation2020}
\cite{vanallenGeneticLandscapeClinical2014}
\cite{louveauBaselineGenomicFeatures2019}
\cite{yanGenomicFeaturesExceptional2019}
\cite{catalanotti2017pten}