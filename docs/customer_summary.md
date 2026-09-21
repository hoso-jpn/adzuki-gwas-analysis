# 顧客summary statistics入力 v1

既存のDryad専用6データセット契約を変更せず、別adapterとして受け入れる。
一つのTSVとmetadata TOML、一つの検証済みreference bundleを使う。
[GWAS-SSF](https://github.com/gwas-catalog/summary-statistics-standard)のdata/metadata分離を参考にした
独自の植物向け契約であり、GWAS-SSF完全準拠とは称しない。

```toml
schema_version = 1
dataset_id = "customer_trait"
cohort_id = "cohort-1"
analysis_id = "analysis-1"
reference = "VERIFIED_REFERENCE"
assembly_id = "VERIFIED_ASSEMBLY_VERSION"
species = "Vigna angularis"
trait = "seed_weight"
trait_unit = "g"
trait_coding = "mass of 100 seeds"
test = "wald"
effect_type = "beta"
effect_scale = "g"
standard_error_scale = "g"
effect_strand = "forward"
target_effect_allele = "ALT"
sample_size_definition = "analyzed_individuals"
invalid_row_policy = "error"

[columns]
chr = "CHROM"
pos = "BP"
ref = "REF"
alt = "ALT"
effect_allele = "EA"
other_allele = "OA"
effect = "BETA"
standard_error = "SE"
pvalue = "P"
sample_size = "N"
effect_allele_frequency = "EAF" # optional
```

```bash
uv run python -m adzuki_gwas_analysis.customer_summary \
  --summary /case/summary.tsv --metadata-path /case/summary.toml \
  --bundle-path /case/reference/reference_bundle.toml \
  --clustering-distance 1000 --output-dir /case/summary-report
```

REF/ALTは同一bundleの正鎖・1-based座標で検証する。このv1は明示的なREF/ALTを要求し、
effect alleleをALTと推測しない。`effect_strand=forward`が確認済みの場合だけ、ALTへの
統一でbetaの符号とEAFを変換する。`unknown`ではpalindromic SNPを含めstrandを推測せず、
原方向のeffectとunresolved状態を保持する。ref/altと同じ組合せでもforward宣言は必要。
有利アレルは出力から推定しない。betaとSEの尺度は一致が必要。
OR入力ではeffect_typeをodds_ratio、effect_scaleとstandard_error_scaleをlog_oddsにする。
入力ORをlogへ変換し、SEは既にlog-odds尺度の値を要求する。

各検定familyはdataset×cohort×trait×analysis×test。補正のmは保持された有効行数。
不正値は既定でエラー。`invalid_row_policy=exclude`を明示した場合だけ不正行を除外し、
source row番号と理由を記録する。重複variant identityは除外設定でもエラー。
サンプル数は正整数、SE>0、EAFは[0,1]、nonfiniteは拒否する。空欄・NAを勝手に補完しない。

pはDecimalで読み、`1e-400`等を保持する。0は原計算の-log10(p)がある場合だけ許可し、
columnsに`neg_log10_pvalue`を追加する。正のpとlog-pが両方ある場合は整合性を検査する。
補正・候補順位・図はlog空間で処理する。`source_pvalue`は入力値を保持し、0に丸められた
値を復元値と偽らない。`neg_log10_pvalue`と`pvalue_source`を併せて読む。

χ²への変換根拠がある場合だけmetadataに`chi_square_df`と`chi_square_basis`を指定する。
指定がない場合、検定不明の場合、中央値が対応可能な数値範囲を超える場合はλGCを未評価にする。
旧DryadのLRT df=1は流用しない。λGCでp値を再補正しない。

成果物はnormalized_summary/candidate_snps/excluded_rows TSV、Manhattan/QQ PNG、
statistical_diagnostics.json、analysis_report.md、生成時来歴。candidate_snpsは周辺配列抽出へ
直接渡せる。物理距離signalはLD/QTLの証拠ではない。初期行数上限は100万行で、超過は
部分結果を納品せず停止する。変更には`--max-rows`を明示し、実案件規模の性能は別途検証する。
個体別genotypeのない入力からGWAS・kinship・LD・GSを再実行しない。
