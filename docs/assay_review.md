# ARMS実験結果・検証状態の管理 v1

Issue #22。#14の計算設計を実験担当者に引き渡し、返却された測定・照合情報から
限定された条件での検証状態を記録する。CLIはオフラインであり、実験や発注を実行しない。
今回の実装・CIで使用するのは合成データだけで、実測の性能や検証済みパネルは提供しない。

## 引き渡しと返却

来歴が整合する `arms_design` バンドルを指定する。実験前は結果を付けずに実行できる。

```bash
uv run python -m adzuki_gwas_analysis.assay_review \
  --design-dir deliveries/arms-design --output-dir deliveries/assay-handoff
```

`handoff.tsv`は候補・設計・primer version・assembly・アレル・5′→3′配列・計算条件と
設計run IDを保持する。REF-specific + common、ALT-specific + commonの2反応を
区別し、wet-lab側で測定条件・対照・判定ルールを定める。
`assay_return_template.tsv` は列名だけの返却様式。空の様式は測定済みと扱わない。
結果がない全設計は `computational_candidate` のままである。

実験後は結果TSVと測定・レビューのmetadataを一緒に指定する。

```bash
uv run python -m adzuki_gwas_analysis.assay_review \
  --design-dir deliveries/arms-design --results cases/example/assay-results.tsv \
  --metadata-path cases/example/assay-review.toml --output-dir deliveries/assay-review
```

TSVの全必須列（追加列・未知の列を拒否する）：

| 列 | 契約 |
|---|---|
| candidate_id, design_id, primer_version, assembly_id, ref, alt | 引き渡し表と完全一致。REF/ALTはゲノムforward方向。反転・旧版を暗黙変換しない |
| allele_encoding | 常に `ALT_0_1_2_forward`。primer strandに関わらずゲノムALT dosageで返す |
| conditions_id | metadataの測定条件IDに一致 |
| plate_id, batch_id, sample_id, replicate_id | 空でない文字列。先頭ゼロやID `NA`を保持。design/batch/plate/sample/replicateの重複を拒否 |
| control_type | `sample`, `positive`, `negative` |
| call | `0`, `1`, `2`, `NA`。生のバンド像・弱いバンドの解釈は実験側で実施。判定不能はNA |
| no_call_reason | sample/positiveのNAでは必須。call済みの場合は空 |
| expected_genotype | positive対照の期待0/1/2。negative/sampleではNA |
| reference_genotype | sampleの独立した照合genotype 0/1/2、未知ならNA。対照ではNA |

同じcandidate/sampleの対照種別・期待値・照合genotypeが反復間や設計代替間で変わる場合は
入力を拒否する。正しい形式の不良対照やno-callは除去せず、採否に使う。
各design/batch/plateに0/1/2のpositive対照とnegative対照を求め、欠けたplateを数える。
違う測定条件は別reviewとして扱う。

metadata例。以下の件数・率は**様式例**であり、汎用の検証基準ではない。
実案件の基準は用途と試験計画に応じて事前に決め、review記録を残す。

```toml
schema_version = 1
validation_dataset_id = "assay-validation-1"
data_scope = "synthetic"
conditions_id = "conditions-1"
measurement_conditions = "synthetic example; record chemistry, cycling, DNA input and call protocol"
target_population = "synthetic breeding population"
reference_method = "synthetic known genotype; real cases identify independent reference assay"
analyst_review_reference = "synthetic-review-record"
reviewer_id = "synthetic-reviewer"
analyst_approved = true
requested_state = "assay_validated"
population_scope_confirmed = false
population_evidence_reference = "not_applicable"
minimum_samples = 24
minimum_comparable_samples = 24
minimum_replicated_samples = 3
minimum_call_rate = 0.95
minimum_concordance = 0.99
```

`data_scope`はsynthetic/customer/public。データ件数、測定条件、照合方法、レビューを
プログラムが外部から認証するものではない。入力hashと確認記録を保持し、人が内容を確認する。

## 分母と採否

`assay_metrics.tsv`はdesignごとに以下を別々に記録する。

- call rate：callされたsample測定数 / 全sample測定数。技術反復を含む。
- complete-sample call rate：全反復でcallできたsample数 / 異なるsample数。
- concordance：照合genotypeに一致したcall数 / callと照合genotypeが両方ある測定数。
- sample数、照合可能な異なるsample数、比較可能な反復sample数、全反復sample数。
- no-call数、反復間の異なるcallがあるsample群数、no-callを含む反復群数。
- 対照件数、不良対照数、plate数、必要対照が欠けたplate数。

対照はcall rate/concordanceの分母に含めない。no-call・未知の照合genotypeはconcordanceに
含めず、専用列と元の結果で追跡する。分母0の率は空/未評価とし、100%を代入しない。
反復は両方callされたsample群を最低件数に数える。反復が片方NAの場合は不一致ではないが、
比較可能な反復としても数えない。小標本の100%一致や技術反復を母集団性能保証に使わない。

## 状態と履歴

| state | 要件 |
|---|---|
| computational_candidate | 計算候補。測定なし、基準未達、未承認、または昇格を要求していない |
| assay_validated | 測定条件・照合方法・reviewが明示され承認済み。設定した件数と率を満たし、対照不足/不良・反復不一致がない |
| population_validated | 同一design/primer/assembly/alleleの過去assay reviewを検証し、別の検証datasetと対象集団の証拠・明示承認があり、今回も基準を満たす |

合成データでは運用用`state`を**常にcomputational_candidate**とし、計算上の遷移を
`simulated_state`にのみ記録する。以前の状態を新しいprimer版へ自動継承しない。
population昇格の条件だけが不足する場合、今回の結果によるassay状態までを保存し、
要求した状態を拒否した理由を記録する。検証状態は記録された集団と測定条件に限定される。

population reviewは `--previous-review deliveries/assay-review` を追加する。
前回の実験データと同じdataset IDを育種集団検証の別証拠として再利用できない。
別IDは独立試料の証明ではなく、対象集団と測定規模の確認責任はreviewerにある。

`alternative_designs.tsv`に同じ候補の利用可能な計算代替を示す。代替の実験成功は保証しない。
再設計は `--previous-review` と `--redesign-history` で入力する。
様式は `redesign_history_template.tsv`。candidate_id、prior_design_id、replacement_design_id、
reason、evidence_reference、analyst_idが必須で、前回の設計と今回の別設計を同じ候補・
assembly・アレルとして結ぶ。履歴と前回reviewのhash/run IDを保存する。

## 納品・責任範囲

`marker_states.tsv`、`assay_review.json`、`assay_review_report.md`から採否・測定条件・
対象集団・試験規模・確認証拠を追跡する。`population_validated_panel.tsv`には運用用の
population_validated設計だけを出す。合成データと未検証候補では空である。
`design_failures.tsv`に計算設計の失敗理由も残す。
`analysis_run.json`/`bundle_contract.json`が入力・前回review・成果物のhashを検証する。

結果・handoff・kinship等は顧客情報を含み得る非公開成果物であり、公開Gitに追加しない。
特に`assay_results.tsv`は返却されたsample IDと測定結果を保持する。原データの保管・削除条件は
案件運用で管理する。処理失敗時に部分納品を出さず、既存の非空出力先を上書きしない。
返却行・履歴行の初期上限は各100,000件（`--max-rows`）。

解析側は計算設計・版・指標・監査記録を提供する。実験側は測定条件、対照、バンドの解釈と
callを担い、担当者が用途に照らした採否を確認する。この実装に実験発注、試薬購入、
外部事業者への送信、実際のwet-lab作業は含まれない。
