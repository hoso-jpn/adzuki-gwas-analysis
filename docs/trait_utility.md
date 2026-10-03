# 遺伝型判定と形質選抜の証拠を分ける（#32）

assayの`population_validated`は、宣言した集団での**遺伝型判定**の検証を意味する。
形質との関連の再現性・マーカー選抜の有効性を意味しない。旧assay schemaを変更せず、
`trait_utility` consumerで独立した軸へmappingする。新assay出力にも`evidence_axes.tsv`を
追加し、assayだけのtrait utilityは常にnot_assessedとする。

```bash
uv run python -m adzuki_gwas_analysis.trait_utility \
  --assay-dir /case/assay-review --output-dir /case/trait-review
```

この呼び出しだけでは形質証拠が無いため、trait_supported_panel.tsvは空。
形質証拠がある場合は`--evidence-path evidence.json --discovery-path discovery.tsv
--validation-path validation.tsv --analysis-run-path external-analysis.json`を全て指定する。

## 検証結果の受入契約

これは外部の試験解析結果を照合するconsumerであり、圃場試験の混合モデルを自動実行したり、
数値結果を再計算/認証するものではない。元解析の入力、モデル・版、共変量、区間推定法、
事前指定した選抜規則と基準法を担当者が確認し、その実行記録をhashで紐付ける。

`discovery.tsv`はsample_id/family_id、`validation.tsv`はそれにpopulation/environment/
trait/unit/valueを加えた厳格な列契約。各sampleは独立した生物単位で1回だけ。
v1は独立家系の検証に限定し、発見とのsample/family重複、検証内のfamily反復は非承認とする。
relatednessをモデル化した複雑な試験はこのv1で承認しない。

evidence.jsonはschema_version=1、candidate_id/design_id/primer_version/assembly_id/ref/alt、
trait/trait_unit/target_population/environment、discovery_dataset_id/validation_dataset_id、
analysis_plan_reference/analysis_engine/analysis_version/analysis_model/analysis_run_sha256、
validation_data_sha256/discovery_data_sha256、reviewer_id/review_reference/baseline_strategy/
relatedness_review_referenceを非空文字列で要求する。さらに以下を指定する。

- data_scope: synthetic/public/customer。
- desired_direction: increase/decrease、effect_allele: ゲノム正鎖ALT。
- analyst_approved: boolean、independence_basis: independent_families。
- minimum_samples: 4以上の事前指定基準。4が十分な性能保証という意味ではない。
- minimum_effect/minimum_utility: 0以上の事前指定された実用上の最小差。
- association_effect/selection_contrast: それぞれestimate/lower95/upper95/unit。
  後者は事前選抜規則による群とbaseline_strategyとの差で、単なるGWAS betaの流用は禁止。

元run記録・発見/検証TSVのSHA-256が一致し、assayの対象集団とdesign/primer版が一致すること。
形質・単位・環境の違い、未review、sample不足、再利用があれば承認しない。
望ましい方向の95%区間下端が、効果と選抜utilityの各最低基準を**超える**場合のみ、
その範囲でのsupportを記録する。因果性や他環境での収量向上を主張しない。

## 出力と限界

evidence_axes.tsvはgenotyping_validation/association_replication/trait_utilityを分離する。
trait_review.jsonに原metadata、数値証拠、分母、採否理由を保存する。
trait_supported_panel.tsvは運用用trait utilityが支持された設計だけを含む。
合成データはsimulated_trait_utilityにだけ結果を保存し、運用用はnot_assessedのまま。
実験・圃場試験・外部送信を実行せず、他のprimer版へ状態を継承しない。
