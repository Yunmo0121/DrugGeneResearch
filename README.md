# DrugGeneResearch

**A Python proof-of-concept for literature-augmented drug–gene relationship analysis in cancer research.**

The program accepts a gene and a drug name, retrieves relevant PubMed article records and abstracts through NCBI E-utilities, and asks a GLM language model to produce an evidence-based summary. It is an educational implementation of the **core workflow**, not a full reproduction or validation of the research paper *Inferring Drug–Gene Relationships in Cancer Using Literature-Augmented Large Language Models*.

## Features

- Search PubMed using gene and drug terms (`ESearch`).
- Retrieve titles, PMIDs, and available abstracts (`EFetch`).
- Analyze retrieved abstracts using the Zhipu AI / GLM API.
- Report relationship type, possible mechanism, supporting PMIDs, and limitations.
- Save timestamped literature JSON and analysis text files.
- Resume analysis from a local JSON file if the model API is temporarily rate-limited.

## Requirements

- Python 3.10 or later (a virtual environment is recommended)
- Python packages: `requests`, `python-dotenv`, and `zai-sdk`
- A Zhipu BigModel **standard API key** with access to the configured model
- Internet access to NCBI E-utilities and the model API

## Setup

Create and activate a virtual environment (macOS/Linux):

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Create a local `.env` file in the project directory:

```dotenv
ZHIPU_API_KEY=your_own_standard_api_key
```

The basic PubMed retrieval in `main_v2.py` does not require an NCBI API key. **Never upload `.env` or real API keys to GitHub or include them in a submission.**

## Usage

Run a new analysis:

```bash
python main_v2.py
```

Example input:

```text
Enter gene name: EGFR
Enter drug name: Osimertinib
```

Check previously saved literature **without making API calls**:

```bash
python main_v2.py --resume articles_CDK4_Osimertinib.json --check
```

Resume model analysis from a saved literature file (calls only the model API):

```bash
python main_v2.py --resume articles_CDK4_Osimertinib.json
```

## Outputs

- `articles_<gene>_<drug>_<timestamp>.json`: retrieved literature records and abstracts
- `analysis_<gene>_<drug>_<timestamp>.txt`: model-generated relationship analysis, written **only after a successful model response**

The model may encounter rate limits or account usage limits. A failed model request should not be treated as a completed relationship analysis.

## Preliminary tests

| Gene | Drug | Outcome |
| --- | --- | --- |
| EGFR | Osimertinib | End-to-end analysis ran; the model identified a direct-target relationship. |
| CDK4 | Palbociclib | End-to-end analysis ran; the model identified CDK4/6 inhibition. |
| CDK4 | Osimertinib | PubMed retrieval succeeded; model analysis was not completed because the API reached a usage/rate limit. |

These are **exploratory functional tests**, not a validated accuracy benchmark. The third pair is a *non-direct-target control*, not proof of biological irrelevance.

## Limitations

- Only a small number of retrieved abstracts (currently up to five) are supplied to the model; selection bias and incomplete coverage are possible.
- An abstract may not contain enough information to establish a direct molecular mechanism.
- The model may generate incorrect or overly strong claims. All PMIDs and claim-to-evidence links should be manually verified.
- The program does not reproduce the full paper's methods, experiments, or performance evaluation.
- The output is for research/education only, not clinical decision-making.

## Data and attribution

PubMed records are accessed through NCBI E-utilities. Check the original article licenses before redistributing full abstracts or other copyrighted material. Do not publish private research material without the supervisor's permission.

### Related publication

*Inferring Drug–Gene Relationships in Cancer Using Literature-Augmented Large Language Models* — https://pmc.ncbi.nlm.nih.gov/articles/PMC12036822/
