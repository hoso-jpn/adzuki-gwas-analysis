# 参照ゲノム間の候補比較 v1

Issue #21。#12の参照バンドルと、#20/#15の正規化済み候補を使用する。
同一traitの候補だけを、出典・版・checksum・確認記録のある外部UCSC chainを介して比較する。
Miyagi/Shumariの実対応資産を取得・検証したことを意味しない。

## 入力・実行

```bash
uv run python -m adzuki_gwas_analysis.cross_reference \
  --source-table deliveries/source/candidate_snps.tsv \
  --source-bundle cases/example/source/reference_bundle.toml \
  --target-table deliveries/target/candidate_snps.tsv \
  --target-bundle cases/example/target/reference_bundle.toml \
  --chain-metadata cases/example/mapping/chain.toml \
  --trait mass --output-dir deliveries/comparison
```

`source`はchainの**q/query**、`target`は**t/target**。一般的なファイル名の`AtoB`から
方向を推測しない。今回のv1はtarget strand `+`、query strand `+/-`の非圧縮chainに対応する。
whole-genome alignmentからchainへの変換、chain生成、orthology-only形式の取り込みは含まない。
対応する実資産を準備していない場合、`--chain-metadata`を省略すれば参照内注釈だけを行う。
target table/bundleを両方省略した単一参照の注釈も可能。

候補表は `cross_reference.CANDIDATE_FIELDS` の列を持つTSV。dataset/cohort/analysis/trait/testを
保持した1 familyに限定し、dataset/reference/assemblyをバンドルに照合する。
候補ID・変異の重複、REF不一致、無効なbeta/SE/log-p、別trait、異なるspeciesを拒否する。
空の候補表は有効で、0件の会計を保存する。旧Dryadの候補表はそのまま投入せず、明示的な
REF/ALTと形質・effect契約を持つ形式へ変換する必要がある。

`chain.toml`の例：

```toml
schema_version = 1
source_assembly = "source-assembly-accession-version"
target_assembly = "target-assembly-accession-version"
source_uri = "alignment publication or internal asset reference"
version = "alignment-v1"
license = "documented terms"
validation_reference = "review record identifying assemblies and alignment QC"
reviewer_id = "review-record-owner"
reviewed = true
direction = "query_to_target"
chain_path = "alignment.chain"
chain_sha256 = "REPLACE_WITH_64_HEX_SHA256"
```

source/target assemblyが一致し、ファイルhash、contig名と長さ、block構造とheaderの終点が
整合することを検証する。対応資産への確認記録を必須にするが、確認記録の真偽や
全ゲノムalignment品質をプログラムが認証するものではない。

## 対応とアレル

[UCSC chain仕様](https://genome.ucsc.edu/goldenPath/help/chain.html)に従い、block内を
0-based half-openで解釈し、入出力のSNP座標は1-basedとする。
負鎖queryはcontig長から逆変換する。query gapに落ちるSNPを周囲のoffsetで埋めない。
target gapを超えた場合は次blockの正しいoffsetへ移る。chainはストリーム処理し、候補点だけを保存する。

- `mapped`: 入力chain内で1つの対応があり、target FASTAの塩基がstrand変換後のREF/ALTに整合。
- `unmapped`: 対応blockなし。未対応と記録し、関連なしと解釈しない。
- `multimap`: 複数chain hit。scoreで自動的に1つを採用せず、全hitを記録する。
- `ambiguous`: 座標は1つでもtarget塩基をアレルとして整合できない。
- `not_assessed`: 対応資産なし。参照内注釈は実行できる。

一致表は、単一対応と両側FASTA検証に加えて、target候補のREF/ALT集合が一致するものに限定する。
反転とREF/ALT交換は別の列に記録する。effect alleleを対応させ、同じ単位・coding・scale・testで
方向が判明している場合だけ符号を比較する。元のbetaと検定値を書き換えない。
未解決strandや不明な検定契約では符号比較も未評価とする。
`mapped`でも対応先候補が表にないことは、未検定・候補選抜外などを区別できず「関連なし」ではない。
uniqueは**供給されたchain内**の一意性であり、相互一意性・全ゲノムorthologyの証明ではない。

## 遺伝子注釈と出力

[GFF3仕様](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)のgene feature
のID、またはGTF gene featureのgene_idを使用する。overlapは距離0、非重複時はgene intervalの
近い端までの塩基距離を記録し、同距離の全geneを残す。strandとgene intervalも保持する。
gene featureを持たない注釈では `no_gene_features_on_contig`、資産なしでは `not_assessed`。
transcriptをgeneとして推測・統合せず、機能注釈・因果性・調節領域を断定しない。

| 出力 | 内容 |
|---|---|
| candidate_annotations.tsv | 両参照の候補とgene overlap/nearest。近傍であること自体は候補geneの確証ではない |
| coordinate_mappings.tsv | 元/先assembly・座標・アレル、strand、chain ID/score/block長、対応状態と理由、全hit |
| cross_reference_comparison.tsv | 同一traitで対応アレルを確認できた候補の検定・effect・geneを横並びで保持 |
| candidate_accounting.tsv | source/target全候補を1候補1行で会計。未対応理由や対応先未選抜を保持 |
| comparison_summary.json / comparison_report.md | 件数・比較の範囲・解釈上の条件 |
| analysis_run.json / bundle_contract.json | FASTA/GFF/chain/入力のhash、対応資産の出典とreview記録、成果物の完全性 |

上限は各候補表10,000件（`--max-candidates`で変更可能）、chain 1,000,000 blocks、
候補への総hit 100,000件。上限超過は部分結果を出さず失敗する。既存の非空出力先は上書きしない。

元の参照別補正familyは維持し、p値の再補正・集約・meta-analysisを行わない。
同じcohortを2参照で解析したものを独立replicationとしない。cohort IDが違うだけでも独立性は
確認できない。raw p値の差を直接、生物学的差と解釈しない。

合成mini-assemblyで反転、両側gap、contig端、1対多、未対応、第三アレル、REF/ALT交換、
遺伝子距離、形質/assembly不一致、全候補会計を検証する。実配列・顧客データはGitに含めない。
