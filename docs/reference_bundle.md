# 参照ゲノム資産の契約 v1

FASTA本体は案件のローカル領域に置き、Gitへ登録しません。Miyagi、Shumari、
Longxiaodou 4の実資産はこの変更では取得・検証していません。合成配列で契約を検証します。

```toml
schema_version = 1
reference = "Miyagi"
assembly_id = "REPLACE_WITH_VERIFIED_ACCESSION_AND_VERSION"
species = "Vigna angularis"
source = "REPLACE_WITH_VERIFIED_SOURCE"
license = "REPLACE_WITH_CONFIRMED_TERMS"
data_scope = "public" # public / customer / synthetic
dataset_ids = ["miyagi_water_permeability"]

[fasta]
path = "reference.fasta"
assembly_id = "REPLACE_WITH_VERIFIED_ACCESSION_AND_VERSION"
sha256 = "REPLACE_WITH_ACTUAL_SHA256"

[fai]
path = "reference.fasta.fai"
assembly_id = "REPLACE_WITH_VERIFIED_ACCESSION_AND_VERSION"
sha256 = "REPLACE_WITH_ACTUAL_SHA256"
```

この例はテンプレートであり、未置換では検証を通過しません。`annotation`（GFF3/GTF）と
`mask`（BED、0-based half-open）は任意。同じpath/assembly_id/sha256の3項目が必要です。
配列は非圧縮FASTA、FAIは5列。実ファイルを走査してcontigの名前・順序・長さ・offset・
改行幅を再計算しFAIと照合します。FASTA全体をメモリへ保持しません。
GFFの埋込みFASTAは未対応。annotationのcontigと座標範囲も検証します。

```bash
uv run python -m adzuki_gwas_analysis.reference \
  --bundle /case/reference/reference_bundle.toml \
  --dataset-id miyagi_water_permeability --reference Miyagi \
  --assembly-id VERIFIED_ACCESSION_AND_VERSION
```

`dataset_ids`は利用者が確認したデータとアセンブリの明示的な対応です。ファイル名から
アセンブリを推測せず、利用者の申告の科学的妥当性をchecksumだけで証明したとも扱いません。
Python APIの`check_dataset`と`check_snp`を後続処理の入口で呼び出します。
座標は1-based inclusive。正負鎖を指定でき、REF/ALTは必ず参照の正鎖表記です。
相対pathのみ許可し、symlink、ディレクトリ外への参照、checksum違いは拒否します。
読取中の資産変更を避けるため、案件ごとの不変スナップショットを使用してください。

参照間のliftoverやorthology対応はこの機能に含みません。異なるアセンブリの
同じ染色体名・数値座標を同一変異と扱うことはできません。
