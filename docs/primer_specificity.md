# プライマーの計算確認（#33）

既存ARMS計算候補に、オフラインの参照全体検索と任意のPrimer3熱力学計算を追加する。
配列と条件を外部サービスへ送信しない。合成mini-referenceと公開された計算例を試験し、
実アズキ全ゲノムの性能・wet-lab成功・allele discriminationを実証したとはしない。

```bash
uv sync --locked --extra thermodynamics
uv run --extra thermodynamics python -m adzuki_gwas_analysis.primer_specificity \
  --design-dir /case/arms --context-dir /case/context \
  --bundle-path /case/reference/reference_bundle.toml \
  --settings-path config/specificity.example.json --output-dir /case/primer-confirmation
```

設定例のmismatch/3′ seed/断片長/Tm基準は汎用の実験合格基準ではない。
元contextのrun IDとreference FASTA checksumが設計・現bundleの両方に一致することを要求する。
案件のPCR条件をもとに担当者が決める。thermodynamics.enabled=falseではPrimer3を必要とせず、
熱力学確認はnot_checkedとなる。trueなのに依存が無い/固定版と違う場合は停止する。

## 特異性の範囲

各primerの正負両鎖について、3′側exact seedと許容mismatch数を使う**ungapped**検索を行う。
予測ampliconは、内向きでprimer結合区間が重ならないpairを対象とする。
重なるprimer区間やprimer-dimer由来の産物をこのamplicon検索で網羅するとはしない。
2 primerの向き・位置・距離からREF-specific+common、ALT-specific+commonの2反応を別々に評価する。
各反応内で同一primerが2か所に結合するspecific/specific、common/commonの予測産物も検出する。
ALT反応の意図した標的は、候補位置のみALTに替えた明示的なhaplotypeを使う。
REF参照しかないためALT primerの3′末端が一致しない問題を、標的不存在と取り違えない。
現行ARMS設計エンジンは二次的な意図的mismatchを導入しない。外部で変更したprimerを
使う場合も3′ seed内のmismatchを検索規則から除外している点を確認し、設計の来歴を更新する。

全contigを検索し、primer_hits.tsv/amplicons.tsvにhit・向き・mismatch・予測断片を残す。
意図した断片が1つで他の予測断片がなければpassed_within_declared_reference_and_rules。
余分な断片・標的不在はfailed、参照にN等の未解決塩基があればambiguous。
設定した参照サイズ/contig/設計数/hit数上限を超えれば停止し、部分的なpassを納品しない。
seed候補window数と予測amplicon数にも予算を設ける。予算超過は未解決として停止する。
`reference_search_scope=supplied_reference_contigs_only`、
`reference_completeness=not_assessed`を常に記録する。短い参照や部分的なassemblyでの
passを「全ゲノム特異性確認済み」に読み替えない。seed候補window上限はprimer・contig・
反応ごと、hit/amplicon上限は反応内の全contig累計にも適用する。

この検索ではindelを含むprimer結合、3′ seed内のmismatchを許容した増幅、未提供の
haplotype/構造変異を網羅しない。これらを含めた生物学的な特異性保証ではない。
cohort_binding_checkのunavailable/partialはそのまま残す。実ゲノム規模の資源測定は別途必要。

## 熱力学

[Primer3-py API](https://libnano.github.io/primer3-py/api/bindings.html)のcalc_hairpin、
calc_homodimer、calc_heterodimerを使用する。ref/alt/common各primerと、同一反応内の
ref+common/alt+commonを評価する（REFとALT primerは別反応なのでその相互dimerを合否に使わない）。
Na相当一価イオン/Mg相当二価イオン/dNTPはmM、oligoはnM、温度は℃。
Tm、ΔG/ΔH/ΔS、structure_found、化学条件・package版を保存する。
ΔG計算温度と反応のannealing温度を同一視しない。DMSO等はこのadapterでは扱わない。
公開docsのhairpin既知例との照合をoptional CIでも実行する。

任意依存`primer3-py==2.3.1`はGPL-2.0-or-later。プロジェクトのMIT表明をこの依存へ
適用しない。[公式license説明](https://libnano.github.io/primer3-py/quickstart.html#contributing)
と配布物のライセンスを、商用配布形態のリリース確認に含める。

## 引き渡し

specificity_review.tsvとthermodynamics.jsonを元のdesign/primer versionに結び付けて渡す。
`assay_review --specificity-dir /case/primer-confirmation`で元design run IDを照合し、
計算確認を`primer_confirmation/`へ添付できる。assayの採否には計算結果を代入しない。
元ARMS/assay bundleは変更しない。計算結果だけでassay/population/trait utilityを昇格させない。
発注・試薬購入・実験・外部送信は行わない。
