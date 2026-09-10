# 個体別データによる量的形質 GWAS v1

Issue #15。公開 Dryad summary の再解析とは別の入口で、利用許可のある個体データを
入力する。CI の入力は一時ディレクトリに生成する合成データだけである。
実コホート、Miyagi/Shumari の実配列、Seedcore-01 の実機性能は未検証。

## 入力と実行

```bash
uv sync --locked
MPLBACKEND=Agg uv run python -m adzuki_gwas_analysis.individual_gwas \
  --genotypes cases/example/genotypes.tsv --phenotypes cases/example/phenotypes.tsv \
  --samples cases/example/samples.tsv --variants cases/example/variants.tsv \
  --config-path cases/example/run.toml --bundle-path cases/example/reference/reference_bundle.toml \
  --output-dir deliveries/example-gwas
```

TSV はUTF-8、IDは文字列。`001` と `1` は別ID、`NA` はsample IDとして有効。
行の位置ではなくIDで結合する。3つのsample ID集合は完全一致が必要で、黙って内結合しない。

| 入力 | 必須列・条件 |
|---|---|
| genotypes.tsv | `sample_id`、その後に各 `marker_id` の列。セルはリテラル `0`/`1`/`2`/`NA`（ALT dosage）。VCF直読、dosage確率、倍数性変更は未対応 |
| phenotypes.tsv | `sample_id,trait,value,unit`。個体ごとに1観測、有限数、同一trait/unit。欠測や重複、2値形質は拒否 |
| samples.tsv | `sample_id,cohort_id` と指定した数値共変量。コホートは1つ。欠測を拒否 |
| variants.tsv | `marker_id,chr,pos,ref,alt,assembly_id`。1-based SNP、明示REF/ALTをFASTAで確認。重複変異・マーカー・別assemblyは拒否 |
| reference bundle | #12 の検証済みバンドル。dataset/reference/assembly/speciesが一致すること |

共変量は解析担当者が試験設計から指定する。カテゴリ変数は基準水準を除くダミー列を
事前に作る。切片は自動追加するため、全てのカテゴリ水準を追加して共線性を作らない。
反復測定、多環境のランダム効果、二値・順序形質、G×E、空間補正はこのモデルの対象外。
異なる試験を同一個体の独立観測として投入しない。試験計画からの自動モデル選択は行わない。

TOML設定例（`synthetic` は合成データの例であり、実案件は確認済み情報を入力する）：

```toml
schema_version = 1
dataset_id = "synthetic_trait"
cohort_id = "synthetic-cohort"
analysis_id = "analysis-1"
reference = "Synthetic"
assembly_id = "synthetic-v1"
species = "synthetic plant"
trait = "mass"
trait_unit = "g"
trait_coding = "untransformed mass"
trait_type = "quantitative"
genotype_encoding = "ALT_0_1_2"
design_note = "one observation per sample; numeric trial block covariate"
data_scope = "synthetic"
data_use_confirmed = true
permission_reference = "synthetic fixture"
transfer_method = "local generation"
retention_policy = "test duration"
deletion_policy = "temporary directory cleanup"
output_ownership = "synthetic fixture"
covariates = ["block_numeric"]
n_pcs = 0
min_maf = 0.05
max_marker_missing = 0.1
max_sample_missing = 0.1
max_samples = 1000
max_markers = 100000
max_genotype_cells = 5000000
clustering_distance = 10000
alpha = 0.05
fdr_level = 0.05
```

利用許可、転送、保管、削除、成果物帰属の確認記録を必須とする。文字列を入力したことは
その内容をシステムが法的に確認したことを意味しない。実際の保管・削除は案件運用で行う。
設定には秘密情報そのものを記載しない。

## QC・モデル・検定

sample欠測率超過は実行を中止する。markerは観測コールからMAFと欠測率を計算し、
全欠測・単型・閾値違反を除外会計に残す。残存markerの欠測はそのコホートの平均ALT dosage
`2p`で補完する。2 marker未満、定数形質、共変量rank不足、残差自由度不足は中止する。

エンジンは `numpy-scipy-null-reml-lmm-v1`。依存パッケージは既存の `uv.lock` を使用する。
モデルは `y = Cb + xβ + u + e`、`u ~ N(0, σg² K)`、`e ~ N(0, σe² I)`。
`K = ZZ' / (2 Σ p(1-p))`、`Z = G - 2p` で全QC通過markerから作る。
指定数の上位kinship PCと、平均・標準偏差で標準化した指定共変量をCへ加える。
上限は密行列実装の受付制限であり、性能保証ではない。メモリは概ねO(nm+n²)、
計算はO(n³+n²m)。設定上限を上げる前に実機で測定する。

帰無モデルのREML profile objectiveを `log(σe²/σg²) ∈ [-12,12]` で数値最適化し、
共分散比を全markerで固定して一般化最小二乗を行う。p値は二側t検定、自由度は
`n - rank(C) - 1`。共変量と共線なmarkerは検定せず理由を残す。
境界解を `null_reml_boundary` に記録する。ランダムな処理・seedはない。

この検定は**共分散推定を固定する近似**であり、GEMMAのmarker別REML/LRTでも、EMMAX本体でもない。
全markerのKには検定markerも入り、proximal contaminationにより検出力を落とす可能性がある。
LOCO、外部エンジンとの実コホート照合、較正の保証は未実装。正規性、等分散、試験設計、
境界解の妥当性は案件ごとに確認する。PC固有値が重複する場合は回転が一意でなく、
異なる数値ライブラリ・ハードウェア間のbitwise一致を保証しない。

背景となる分散成分法は [Kang et al. 2010](https://pubmed.ncbi.nlm.nih.gov/20208533/)、
relationship matrixの構成は [VanRaden 2008](https://pubmed.ncbi.nlm.nih.gov/18946147/) を参照。
これらの論文の性能は本実装の性能検証を代替しない。

## 出力・接続・検証

- `association_results.tsv`: エンジンのbeta、SE、pと元の-log10(p)、ALT方向、n。
- `normalized_summary.tsv`, `candidate_snps.tsv`: #20と同じ列契約。候補は#13のcontextへ接続可能。
- `sample_qc.tsv`, `marker_qc.tsv`, `excluded_rows.tsv`: 分母と除外理由。nは補完後の解析個体数、観測コール数はmarker QCに別記。
- `kinship.npy`, `sample_order.tsv`, `pc_scores.tsv`, `model.json`: 個体順、行列、PC、変換、エンジン条件。
- `manhattan.png`, `qq.png`, `statistical_diagnostics.json`, `analysis_report.md`: 独立した1 familyのBonferroni/BHと物理距離候補。t検定に根拠のないchi-square dfを付けず、λGCは未評価。
- `analysis_run.json`, `bundle_contract.json`: 入力hash・実装/環境・成果物hash。原データはコピーしない。

実行失敗時は完成品を出さず、既存の非空納品先は上書きしない。顧客ID・kinship・測定値を含む
出力は非公開の案件領域で扱い、公開Gitに追加しない。既存6 datasetへ混合しない。
正のALT betaは形質値の増加であり、有利アレルや原因変異を意味しない。

テストはK=Iでの独立OLS比較、直接行列計算のGLS比較、誤入力拒否、64個体の植込み効果、
QC除外、補正・図・候補・来歴、同条件の結果再現を確認する。Seedcore-01測定と実験検証は含まない。
