# 解析生成時の来歴

旧batch/candidatesの25/43/5ファイル契約は維持する。新規の納品向け解析では
`--record-provenance`を付け、レポート生成には`--require-provenance`を付ける。

```bash
uv run adzuki-gwas-analyze batch --manifest manifest.toml --data-dir data/raw \
  --clustering-distance 1000 --record-provenance --output-dir /case/analysis
uv run adzuki-gwas-analyze report --analysis-dir /case/analysis \
  --require-provenance --output-dir /case/delivery
```

1000 bpはコマンド例であり推奨窓幅ではない。窓は案件ごとに明示して選ぶ。
`candidates`にも同じ記録オプションがある。全成果物をstaging内で作り、入力が開始時から
変更されていないことを検査してから公開する。候補出力もこの経路では一括公開になる。

- `analysis_run.json` v1: 入力のroleとSHA-256、Python/依存版、ソース内容hash、
  取得可能なGit commitとlock hash、統計family、解析引数、成果物のhashと一覧。
- `bundle_contract.json` v1: 来歴必須の宣言とanalysis_runのhash。
- 同じソース・環境・入力・設定から同じrun_idを得る。時刻、hostname、username、
  環境変数、入力の絶対pathは採取しない。wheel配布等でcommit/lockがない場合はnull。
  未コミット変更も反映するソース内容hashがコード同一性の基準になる。
- reportは来歴と成果物を照合し、コピー後も再検証する。解析環境とレポート環境を
  `software_versions.json` v2の別項目として記録し、来歴2ファイルも納品へ引き継ぐ。

旧成果物は記録なしとして読み込み、解析環境を`unavailable_from_source_artifacts`とする。
来歴ファイルが片方しかない場合、hash違い、余分/欠落成果物は拒否する。
両方の来歴ファイルが除去された旧形式との区別には、呼出側の`--require-provenance`が必要。
これは署名や悪意のある全面改ざんの認証機構ではない。

6データセットを同時にDataFrameへ保持しない。既存の数値・検定family・候補抽出規則は
変更せず、入力hashの再読取だけを追加する。RO-Crate exportは将来の任意拡張。
