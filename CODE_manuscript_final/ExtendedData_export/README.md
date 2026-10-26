# Extended Data publication export

Pipeline intermediates live under `DATA_manuscript_final/ExtendedData/`.
Publication sheets are written to `RESULT_manuscript_final/ExtendedData/`
(one file per sheet name) plus a combined `ExtendedData.xlsx`.

## Run

```bash
python3 /home2/sjim/openKoGES/CODE_manuscript_final/ExtendedData_export/export_ExtendedData_publication.py
```

## Inputs → outputs

| Publication sheet (RESULT) | Intermediate source (DATA) |
|---|---|
| Extended Data1.xlsx | `ExtendedData1_{continuous,categorical}.csv` + `ED1_phenotype_metadata.tsv` |
| Extended Data2A.xlsx | `ExtendedData2_continuous.xlsx` → `Total_summary_Novelty_Repl` |
| Extended Data2B.xlsx | `ExtendedData2_categorical.xlsx` → `Total_summary_Novelty_Repl` |
| Extended Data3.xlsx | `ExtendedData_conditionally_independent_172.csv` |
| Extended Data4A/4B.xlsx | `ExtendedData2_continuous.xlsx` → `Replication_{EAS,EUR}_meta` |
| Extended Data5A/5B.xlsx | `ExtendedData2_categorical.xlsx` → `Replication_{EAS,EUR}_meta` |
| Extended Data6A.xlsx | `ExtendedData3.csv` (winner’s-curse continuous dump) |
| Extended Data6B.xlsx | WC recomputed from `ExtendedData2_categorical.xlsx` |
| Extended Data8A.xlsx / Extended Data 8B.xlsx | Curated multi-header layout from combined workbook |

`ExtendedData2_*.xlsx` and related intermediates are also inputs to Table1,
SuppleTable7, Figure4A, SuppleFigure8, and winner’s-curse scripts.
