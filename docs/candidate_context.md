# 候補SNPの周辺配列と近傍変異

```bash
uv run python -m adzuki_gwas_analysis.context \
  --candidate-table /case/analysis/candidate_snps.tsv \
  --bundle-path /case/reference/reference_bundle.toml \
  --dataset-id miyagi_water_permeability \
  --flank-bp 200 --neighbor-window-bp 1000 --output-dir /case/context
```

窓幅は例でありLDの推定値ではない。入力は一datasetだけを含むTSVで、`chr,pos,trait`と
`ref,alt`または既存GWASの`allele1,allele0`が必要。既存candidate_snps.tsvを直接読める。
`reference,assembly_id,dataset_id`があればbundleの対応と照合する。省略された情報は
利用者が指定するbundle/datasetの対応に基づき明示的に補う。別アセンブリを自動推測しない。

GWASのallele1/allele0から配列用REF/ALTへ変換するときはFASTAとの一致だけを用いる。
**betaの方向や有利アレルは推測しない。** `candidate_id`省略時は座標・アレル・assembly・
datasetから安定したIDを生成する。重複ID/重複座位、参照・アレルの不一致は拒否する。
`strand`は省略時`+`、`-`ならFASTA出力を逆相補配列とし、SNPのoffsetも変換する。
REF/ALTは常に参照正鎖。端で切り詰められた配列は`flanks_complete=False`になる。

- `candidate_context.tsv`: 同一性、配列区間、向き、SNP offset、配列hash、入力hash、近傍件数。
- `flanking_sequences.fasta`: ID（percent encoding）とdataset/reference/座標/向き/配列hash。
- `neighboring_variants.tsv`: 候補と近傍変異、summaryまたはcohort_vcfの由来と元ファイルhash。
- `analysis_run.json` / `bundle_contract.json`: 入出力と生成環境の来歴。

`--summary-table`は明示された同じdatasetのsummary TSV（chr/posと2アレル）を読む。
候補自身を近傍件数へ含めず、候補を中心とするwindow内で情報源ごとに記録する。
`--cohort-vcf`は非圧縮VCFまたはgzipを受け付け、`##reference=`がbundleのassembly_idと
完全一致する必要がある。REF正規塩基のSNP/indelに対応し、symbolic ALTは対象window内では
拒否する。VCFのFILTER値を保持し、PASS以外を検証済みと読み替えない。

summary由来の変異は集団内の全変異ではなく、VCFについても未観測変異がないとは保証しない。
GT、LD、因果性は推定しない。maskがあれば配列区間との重なりを示し、注釈未取得は未取得と
表示する。遺伝子との距離の解釈は参照内注釈機能で別途行う。実配列はGitへ登録しない。
