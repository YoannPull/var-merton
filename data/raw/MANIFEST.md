# MANIFEST — data/raw/

Inventaire des fichiers d'entrée attendus dans `data/raw/`. Renseignez la colonne
"Récupéré le" lors du dépôt final. Aucun de ces fichiers n'est produit par le code :
ce sont tous des **entrées**.

## Séries FRED (https://fred.stlouisfed.org) — CSV, colonne `observation_date`

| Fichier attendu | Mnémonique FRED | Variable | Récupéré le |
|---|---|---|---|
| `VIXCLS.csv`    | VIXCLS    | VIX (volatilité implicite) | |
| `CPIAUCSL.csv`  | CPIAUCSL  | IPC (déflateur)            | |
| `CNP16OV.csv`   | CNP16OV   | Population 16+             | |
| `USPRIV.csv`    | USPRIV    | Emploi privé               | |
| `WTISPLC.csv`   | WTISPLC   | Pétrole WTI                | |
| `GS2.csv`       | GS2       | Taux 2 ans                 | |
| `T10Y2Y.csv`    | T10Y2Y    | Spread 10A–2A              | |
| `T10Y3M.csv`    | T10Y3M    | Spread 10A–3M              | |
| `NFCI.csv`      | NFCI      | Conditions financières     | |
| `GDPC1.csv`     | GDPC1     | PIB réel                   | |
| `GPDIC1.csv`    | GPDIC1    | Investissement réel        | |
| `UNRATE.csv`    | UNRATE    | Taux de chômage            | |
| `PAYEMS.csv`    | PAYEMS    | Emploi non agricole        | |

> Le papier mentionne aussi **HOANBS** (heures non agricoles) et **DRALACBN**
> (delinquency, Annexe A.8). Ajoutez les CSV correspondants si vous reproduisez
> ces parties.

## Indice GPR — Caldara & Iacoviello (2022)

| Fichier attendu | Source | Notes |
|---|---|---|
| `data_gpr_daily.csv` | https://www.matteoiacoviello.com/gpr.htm | quotidien ; séparateur décimal "," (lu via `read.csv2`/`fread(dec=",")`). Colonnes GPRD, GPRD_ACT, GPRD_THREAT, DAY (AAAAMMJJ). |

## Incertitude de politique économique — Baker, Bloom & Davis (2016)

| Fichier attendu | Source | Notes |
|---|---|---|
| `EPU_US.xlsx` | https://www.policyuncertainty.com | feuille 1, colonnes `Year`, `Month`, `News_Based_Policy_Uncert_Index` |

## Marché actions

| Fichier attendu | Source | Notes |
|---|---|---|
| `sp500_GSPC_snapshot.csv` | Yahoo Finance (^GSPC) via `tidyquant` | **généré au 1er run** par `02_build_dataset.R` ; à committer pour la réplication exacte. |

## Données de défaut (entrée du satellite)

| Fichier attendu | Source | Notes |
|---|---|---|
| `default/DRALACBN.csv` | FRED — Delinquency Rate on All Loans and Leases, All Commercial Banks | colonnes `observation_date`, `DRALACBN` ; trimestriel, en \%. |

## Comparaison / validation

| Fichier attendu | Source | Notes |
|---|---|---|
| `vix_caldara.csv` | jeu de comparaison interne | colonnes `quarter`, `SPVXO` (complète le VIX avant l'ère VIXCLS) |
| `comparaison_data_caldara.csv` | réplication Caldara | séparateur décimal "," ; préfixes `C_*` ; utilisé par `91_validate_caldara.R` |
