import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import argparse
import json

VARIANT_ID_SEP = ":"
FILL_NA = 0.0
MISSING_VALUE_INDICATOR = -999.0

VCF_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def _parse_gt(gt_str):
    sep = "|" if "|" in gt_str else "/"
    parts = gt_str.split(sep)
    try:
        a1 = int(parts[0]) if parts[0] != "." else 0
        a2 = int(parts[1]) if len(parts) > 1 and parts[1] != "." else 0
    except (ValueError, IndexError):
        a1, a2 = 0, 0
    return a1, a2


def _parse_mtd_count(info_str):
    for field in info_str.split(";"):
        if field.startswith("MTD="):
            methods = field[4:].split(",")
            return len([m for m in methods if m])
    return 0


def _extract_features(row, sample_col):
    ref = str(row["REF"])
    alt = str(row["ALT"]).split(",")[0]
    flt = str(row["FILTER"])
    info = str(row["INFO"])

    ref_len = len(ref)
    alt_len = len(alt)
    is_snp = int(ref_len == 1 and alt_len == 1)
    is_indel = int(ref_len != alt_len)

    filter_pass       = int(flt == "PASS")
    filter_suspicious = int("Suspicious" in flt)
    filter_overlap    = int("OverlapConflict" in flt)

    mtd_count = _parse_mtd_count(info)

    gt_raw = str(row.get(sample_col, "0/0"))
    a1, a2 = _parse_gt(gt_raw)
    gt_is_het     = int(a1 != a2)
    gt_is_hom_alt = int(a1 == 1 and a2 == 1)
    gt_is_hom_ref = int(a1 == 0 and a2 == 0)

    return {
        "ref_len":           float(ref_len),
        "alt_len":           float(alt_len),
        "is_snp":            float(is_snp),
        "is_indel":          float(is_indel),
        "filter_pass":       float(filter_pass),
        "filter_suspicious": float(filter_suspicious),
        "filter_overlap":    float(filter_overlap),
        "mtd_count":         float(mtd_count),
        "gt_allele1":        float(a1),
        "gt_allele2":        float(a2),
        "gt_is_het":         float(gt_is_het),
        "gt_is_hom_alt":     float(gt_is_hom_alt),
        "gt_is_hom_ref":     float(gt_is_hom_ref),
    }


def load_vcf(path):
    path = Path(path)
    records = []
    col_names = None
    sample_col = None

    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM") or col_names == None:
                if line.startswith("#CHROM"):
                    col_names = line.lstrip("#").split("\t")
                else:
                    col_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SRR5131688']
                format_idx = col_names.index("FORMAT")
                sample_cols = col_names[format_idx + 1:]
                if not sample_cols:
                    raise ValueError(f"VCF '{path}' has no sample column after FORMAT.")
                sample_col = sample_cols[0]
                continue
            parts = line.split("\t")
            row = dict(zip(col_names, parts))
            records.append(row)

    if not records:
        raise ValueError(f"VCF '{path}' contains no data records.")

    raw_df = pd.DataFrame(records)

    raw_df["variant_id"] = raw_df["CHROM"].str.strip() + VARIANT_ID_SEP + raw_df["POS"].str.strip()

    alleles_dict = {}
    feature_rows = []
    for _, row in raw_df.iterrows():
        ref = str(row["REF"])
        alt = str(row["ALT"]).split(",")[0]
        gt_raw = str(row.get(sample_col, "0/0"))
        a1, a2 = _parse_gt(gt_raw)
        
        if a1 == 0 and a2 == 0:
            allele = ref
        elif a1 == 1 or a2 == 1:
            allele = alt
        else:
            allele = ref
        
        alleles_dict[row["variant_id"]] = (ref, alt, allele)
        
        feats = _extract_features(row, sample_col)
        feats["variant_id"] = row["variant_id"]
        feature_rows.append(feats)

    feat_df = pd.DataFrame(feature_rows).set_index("variant_id")
    feat_df = feat_df.fillna(FILL_NA)

    return feat_df, alleles_dict


def validate_vcf_tables(
    t1: pd.DataFrame,
    t2: pd.DataFrame,
    t3: pd.DataFrame,
) -> None:
    if not (t1.columns.tolist() == t2.columns.tolist() == t3.columns.tolist()):
        raise ValueError("All three VCF tables must produce identical feature columns.")

    ids1, ids2, ids3 = set(t1.index), set(t2.index), set(t3.index)
    if not (ids1 == ids2 == ids3):
        only_in_1 = ids1 - ids2 - ids3
        only_in_2 = ids2 - ids1 - ids3
        only_in_3 = ids3 - ids1 - ids2
        raise ValueError(
            f"Variant ID mismatch across VCF files.\n"
            f"  Only in vcf1: {list(only_in_1)[:5]} ...\n"
            f"  Only in vcf2: {list(only_in_2)[:5]} ...\n"
            f"  Only in vcf3: {list(only_in_3)[:5]} ..."
        )


def preprocess_tables(
    path1,
    path2,
    path3,
    scale = True,
    intersect_ids = False,
    fill_missing = True,
):
    t1, alleles1 = load_vcf(path1)
    t2, alleles2 = load_vcf(path2)
    t3, alleles3 = load_vcf(path3)

    if not (t1.columns.tolist() == t2.columns.tolist() == t3.columns.tolist()):
        raise ValueError("All three VCF tables must produce identical feature columns.")

    feature_names = t1.columns.tolist()

    if intersect_ids:
        common_ids = sorted(set(t1.index) & set(t2.index) & set(t3.index))
        if not common_ids:
            raise ValueError("No common variant IDs found across the three VCF files.")
        t1 = t1.loc[common_ids]
        t2 = t2.loc[common_ids]
        t3 = t3.loc[common_ids]
    elif fill_missing:
        all_ids = sorted(set(t1.index) | set(t2.index) | set(t3.index))
        if not all_ids:
            raise ValueError("No variant IDs found in any of the VCF files.")
        
        t1 = t1.reindex(all_ids)
        t2 = t2.reindex(all_ids)
        t3 = t3.reindex(all_ids)
        
        t1 = t1.fillna(MISSING_VALUE_INDICATOR)
        t2 = t2.fillna(MISSING_VALUE_INDICATOR)
        t3 = t3.fillna(MISSING_VALUE_INDICATOR)
        
        missing_in_1 = set(all_ids) - set(t1.index[t1.iloc[:, 0] != MISSING_VALUE_INDICATOR])
        missing_in_2 = set(all_ids) - set(t2.index[t2.iloc[:, 0] != MISSING_VALUE_INDICATOR])
        missing_in_3 = set(all_ids) - set(t3.index[t3.iloc[:, 0] != MISSING_VALUE_INDICATOR])
        
        print(f"Total unique variant IDs: {len(all_ids)}")
        print(f"Missing in VCF1: {len(missing_in_1)} variants")
        print(f"Missing in VCF2: {len(missing_in_2)} variants")
        print(f"Missing in VCF3: {len(missing_in_3)} variants")
        print(f"Missing values filled with indicator: {MISSING_VALUE_INDICATOR}")
        
        common_ids = all_ids
    else:
        validate_vcf_tables(t1, t2, t3)
        common_ids = sorted(t1.index.tolist())
        t1 = t1.loc[common_ids]
        t2 = t2.loc[common_ids]
        t3 = t3.loc[common_ids]

    allele_match = []
    for vid in common_ids:
        alleles_list = []
        if vid in alleles1:
            alleles_list.append(alleles1[vid][2])
        if vid in alleles2:
            alleles_list.append(alleles2[vid][2])
        if vid in alleles3:
            alleles_list.append(alleles3[vid][2])
        
        if len(alleles_list) < 2:
            allele_match.append(0.0)
        else:
            from collections import Counter
            counts = Counter(alleles_list)
            max_count = max(counts.values())
            allele_match.append(1.0 if max_count >= 2 else 0.0)
    
    allele_match_array = np.array(allele_match).reshape(-1, 1)

    scaler = None
    if scale:
        all_data = np.vstack([t1.values, t2.values, t3.values])
        
        missing_mask = (all_data == MISSING_VALUE_INDICATOR).any(axis=1)
        valid_mask = ~missing_mask
        
        if valid_mask.sum() > 0:
            scaler = StandardScaler()
            scaler.fit(all_data[valid_mask])
            
            def transform_with_missing(data, scaler, missing_val):
                """Transform data while preserving missing value indicators."""
                result = data.copy()
                non_missing = (data != missing_val).all(axis=1)
                if non_missing.any():
                    result[non_missing] = scaler.transform(data[non_missing])
                return result
            
            X1 = transform_with_missing(t1.values, scaler, MISSING_VALUE_INDICATOR)
            X2 = transform_with_missing(t2.values, scaler, MISSING_VALUE_INDICATOR)
            X3 = transform_with_missing(t3.values, scaler, MISSING_VALUE_INDICATOR)
        else:
            raise ValueError("No valid (non-missing) data points found for scaling.")
    else:
        X1 = t1.values.astype(float)
        X2 = t2.values.astype(float)
        X3 = t3.values.astype(float)

    X = np.hstack([X1, X2, X3, allele_match_array])
    print(f"Allele matches: {len(allele_match_array == 1.0)}")

    return X, common_ids, feature_names, scaler


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Preprocess 3 VCF files for IsolationForest."
    )
    parser.add_argument("vcf1", help="Path to first VCF file")
    parser.add_argument("vcf2", help="Path to second VCF file")
    parser.add_argument("vcf3", help="Path to third VCF file")
    parser.add_argument("--no-scale", action="store_true", help="Skip StandardScaler")
    parser.add_argument(
        "--intersect",
        action="store_true",
        help="Use intersection of variant IDs instead of requiring exact match",
    )
    parser.add_argument(
        "--fill-missing",
        action="store_true",
        default=True,
        help="Use union of variant IDs and fill missing data with special indicator (-999.0) (default: True)",
    )
    parser.add_argument("--output", default="preprocessed.npy", help="Output .npy file for X")
    parser.add_argument("--ids-output", default="ids.json", help="Output JSON file for variant IDs")
    args = parser.parse_args()

    X, ids, feature_names, scaler = preprocess_tables(
        args.vcf1,
        args.vcf2,
        args.vcf3,
        scale=not args.no_scale,
        intersect_ids=args.intersect,
        fill_missing=args.fill_missing,
    )

    np.save(args.output, X)
    with open(args.ids_output, "w") as f:
        json.dump(ids, f, indent=2)

    print(f"Preprocessed matrix shape  : {X.shape}")
    print(f"Number of variant IDs      : {len(ids)}")
    print(f"Features per view ({len(feature_names)}): {feature_names}")
    print(f"Saved X to '{args.output}' and IDs to '{args.ids_output}'")
