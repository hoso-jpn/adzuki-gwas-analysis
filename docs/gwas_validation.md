# 個体別GWASの検証経路（#31）

商用提供の承認とは別に、事前指定のsimulationと外部結果の比較を再現可能に実行する。
この変更で実コホート・実参照・外部エンジンを実行済みとはしない。合成試験だけでは
`commercial_status=not_validated`を変更しない。

## LOCO

個体別GWASのrun TOMLに `kinship_mode = "loco"` を指定すると、検定する染色体上の
全QC通過マーカーを除いてKを作り、染色体別に帰無REML共分散を推定する。元のglobal動作は
既定値として維持する。除外後に多型マーカーが2つ未満なら停止し、globalへfallbackしない。
LOCO v1は `n_pcs=0` を必要とする。集団構造を無視してよいという意味ではなく、
必要な独立共変量と適用可否を解析計画で決定する。

`model.json.loco_models`に染色体、背景マーカー数、共分散比、境界解を記録する。
`kinship.npy`は従来互換のglobal Kで、LOCO検定に使ったKそのものではない。LOCO Kは
固定された入力/QC/染色体除外規則から再構成する。単一の`covariance_ratio`はLOCOではnull。
`ploidy=2`、`observation_design="one_observation_per_sample"`以外の宣言は拒否する。
反復区画/環境調整は別の表現型前処理契約を経由し、行を独立個体として追加しない。

## 事前固定simulation

```bash
MPLBACKEND=Agg uv run python -m adzuki_gwas_analysis.gwas_validation simulate \
  --plan-path config/calibration.example.json --output-dir /case/calibration
```

設定例の数値は商用基準ではない。用途に応じて評価前にseed、反復数、標本/マーカー数、
有意水準、帰無棄却率の95% Wilson上限と検出力下限の基準を固定する。
帰無、集団構造、局所効果、欠測、低頻度、遺伝分散0の境界シナリオを全件実施する。
Gaussian random effectを背景Kから生成し、生の表現型をshuffleしない。

1反復につき事前指定した1マーカーを検定する。独立simulation反復を分母にして棄却率、
95% Wilson区間、効果bias、95%効果区間のcoverage、境界解の件数を保存する。
これはgenome-wide FWER、任意の試験設計、実集団の精度の検証ではない。
標的が単型などで検定できなかった反復は除外を報告し、そのシナリオを合格にしない。
CIの少数反復fixtureは処理/判定の回帰検査のみ。

## 外部結果比較

`gwas_validation compare --left normalized.tsv --right external-normalized.tsv
--plan-path comparison.json --output-dir /case/comparison`を使用する。
外部エンジン結果は既存`customer_summary`で明示的に正規化してから渡す。
必須列はchr/pos/ref/alt/effect_allele/other_allele/beta/se/neg_log10_pvalue。
variant集合とeffect方向の不一致は拒否し、推測flipしない。

比較計画JSONはschema_version=1、両TSVのleft_sha256/right_sha256、
beta_atol/se_atol/logp_atolとleft_contract/right_contractを要求する。
各contractはdataset_id/assembly_id/input_genotypes_sha256/input_phenotypes_sha256/
covariates_sha256/model/covariance_strategy/test/effect_encoding/engine/version/license/data_scope
を非空文字列で宣言する。data_scopeはsynthetic/public/customer。
エンジン名/version/license以外の条件が異なる場合は`not_comparable`となり、数値差だけを保存する。
同一engine/versionのコピーは独立比較としてpassにならない。

比較artifactの一致は外部エンジンの実行そのものの証明ではない。原runの入力hash・実行版・
モデル・アレル契約・ログを別途保存し、担当者が確認する。独立エンジン照合、実コホート/参照、
予定規模の資源測定、担当者reviewが揃うまでは#31の商用受入は未完了。
