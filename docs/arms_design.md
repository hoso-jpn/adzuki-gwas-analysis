# ARMSプライマーの計算候補

```bash
uv run python -m adzuki_gwas_analysis.arms \
  --context-dir /case/context --config config/primer_design.toml \
  --output-dir /case/arms
```

検証済みcontext bundleだけを入力とする。各設計はアレル特異的REF primer、アレル特異的
ALT primer、共通primerの3本で、REF+共通／ALT+共通の**2反応**を想定する。
tetra-primer ARMSではない。識別塩基を3'末端に置き、追加の意図的mismatchは導入しない。
参照の正負鎖、両primerの座標、アレル、長さ、GC、Tm、想定断片長を保存する。

内部engine `adzuki-arms-screen-v1`（本プロジェクトMIT、標準ライブラリのみ）は、
指定した長さ・断片長を列挙し、設定条件に合う候補をTm差などで順位付けする。
設定例はアズキで校正済みの基準ではなく、実験条件に合わせて見直す必要がある。
各候補の上位alternatives件だけを保持し、採用・不採用理由を記録する。

TmはSantaLucia 1998のDNA/DNA最近接塩基対モデル（DNA_NN3 convention）に基づく。
塩はNa濃度、鎖濃度は各鎖のnM。等濃度の完全一致二本鎖を想定し、entropyへの塩補正を使う。
[Biopython公式の公開計算例](https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html)
との一致をテストする。Mg、dNTP、DMSO、mismatch、PCR反応速度はモデル化しない。
自己相補性とdimerの値は最長連続相補長によるscreenで、hairpin/dimerの自由エネルギーではない。
詳細な熱力学・特異性検証へ進む際は[Primer3の設定と前提](https://primer3.org/manual.html)も確認する。

- `primer_candidates.tsv`: 設計ID/version、順位、3配列、座標、各設計指標、検証状況。
- `arms_marker_candidates.tsv`: 候補ごとの設計可否、条件で落ちた理由と件数、repeat情報。
- `marker_design_report.md`: 設計モデル、制約、実験前の確認と推奨検証手順。
- 来歴には設定値・engine・source hash・context run IDを保存する。

cohort VCFで確認された追加変異が結合領域に重なる設計は除く。VCFなしはunavailable、
contextの検索窓が結合領域を覆い切らない場合はpartial_window_unverifiedと表示する。
近傍変異がないことを集団内の全変異がない証拠にはしない。
全ゲノム特異性はこのengineでは常にunverified。計算候補を実験済みマーカーと呼ばない。
実験前に特異性・primer間相互作用を詳しく評価し、対照・反復を含むPCR条件検討と
対象集団でのcall rate/concordance測定を行う。試薬購入や実験発注は行わない。
