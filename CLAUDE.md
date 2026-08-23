# CLAUDE.md

adzuki-gwas-analysis: 公開GWAS summary statistics（Chien et al. 2025 / Dryad）の再解析・
可視化リポジトリ。背景・データ出典・全コマンド一覧は README.md と
docs/gwas_input_contract.md が正。ここにはコードから自明でない制約と正規コマンドのみ置く。

## 環境（uv のみ）

- Python 3.11 固定。`requires-python = ">=3.11,<3.12"` なので3.12以降は依存解決自体が失敗する。
- 依存は uv.lock が正。環境構築は `uv sync --locked` のみ。
- `pip install` / `python -m venv` / `pytest` は使わない。テストは標準 unittest。
- 依存を変えるときは pyproject.toml と `uv lock` をセットで更新し、lockfile差分も必ずcommitする。

## 変更後に必ず実行する（CIと同一）

```bash
uv run ruff check src tests scripts
uv run ruff format --check src tests scripts
uv run mypy src
MPLBACKEND=Agg uv run python -m unittest discover -s tests -v
git diff --check
```

- ruff の対象は `src tests scripts` の3つ。mypy は `src` のみ（strict）。`ruff format --check` を省略しない。
- matplotlib はヘッドレス実行のため `MPLBACKEND=Agg` を付ける。`src/adzuki_gwas_analysis/analysis/plotting.py`
  は意図的に `matplotlib.use()` を呼ばないので、この環境変数の指定を省略しない。
- テストは tests/fixtures/ の合成データのみを読む。data/raw/ が無くても全テストが通る設計。
  テストが raw データを要求し始めたら設計が壊れている。

## パッケージ構造の責務分離

- `src/adzuki_gwas_analysis/`（manifest / schema / validate / loader）: Dryad入力契約の
  検証層。`loader.py` はストリーミング専用で DataFrame を返さない。
- `src/adzuki_gwas_analysis/analysis/`: 検証済みデータに対する解析・プロット層。
  `analysis/loader.py` は検証成功後にのみ呼ばれる別のDataFrame loaderで、この2層分離は
  「validation後に1データセットずつ解析する」設計方針の基盤。呼び出し側も複数DataFrameを
  同時に保持してはならない（loaderの呼び出し方自体は複数回呼べるため、保証ではなく方針）。統合しない。
- `scripts/01`-`04`: `analysis/` の薄いbackward-compatibleラッパー。tracked な
  `plots/`/`results/water_permeability/` の既存ファイル名・配置を再現する経路。
- `adzuki-gwas-analyze`（unified CLI）の `all`/`regions` は単一 `--output-dir` へまとめて
  書き出す bundle コマンドであり、tracked成果物の再現コマンドではない（README参照）。
  実データで動作確認するときは常に一時 output directory（`mktemp -d` 等）へ出力し、
  `plots/`/`results/` を直接指定しない。
- `analysis/batch.py`の`batch` subcommandは6データセットをmanifest記載順に逐次処理する。
  1データセットにつき validation 1回・DataFrame load 1回のみで、その同じDataFrameから
  Manhattan/QQ/diagnosticsを全て作る（`analysis/pipeline.py`の`compute_diagnostics_result`
  が「既にload済みのDataFrameから計算する」内部APIを提供する）。6件のDataFrameや
  adjusted p-value配列を同時保持しない。出力はstaging directoryへ全25ファイルを作った後
  `--output-dir`へ一括publishし、途中失敗時はstagingを削除して`--output-dir`を変更しない。
- `analysis/candidates.py`（Issue #10）は`significant_variants.tsv`と同じ、既にload済みの
  有意variant集合をインメモリで消費するpure関数のみ。raw fileの再読込・再検証は行わない。
  `clustering_distance`（bp）はこのリポジトリ内にLDやQTL区間の推定根拠が無いため必須引数で
  default値を持たない（`candidates` subcommandは required、`batch --clustering-distance`は
  未指定なら候補抽出自体をskipするoptional拡張で、指定しない限りbatchの既存25ファイル出力・
  `batch_summary.tsv`の`schema_version=1`は完全に不変）。物理距離クラスタは"signal"と呼び、
  LD blockやQTL区間と表現しない。lead variantをcausal variantと表現しない。
  `priority_tier`/`priority_reasons`はBonferroni/BH significance flagのみから決まる
  downstream validation priorityであり、生物学的重要度ランキングや検証済みマーカーではない。
  6データセット・Miyagi/Shumariを跨いだ候補の統合・共通rankingは行わない。
- `analysis/report.py`（Issue #11）は`batch`成果物のconsumerであり、GWAS・candidate rankingの
  再計算は一切行わない。`--analysis-dir`は`batch_summary.tsv`の`schema_version==2`
  （candidate-enabled batch）のみを受け付け、schema_version=1は
  `ReportRequiresCandidateEnabledBatchError`で明示的にfailする（暗黙のclustering distanceを
  設定しない）。`analysis/report_validation.py`が唯一のゲートで、全artifactを読み取り専用で
  cross-checkしてからでないと`report_content.py`（pure関数、ファイルI/Oなし）が呼ばれない。
  delivery packageへコピーするのは明示allowlist上のderived artifactsのみで、raw
  `.assoc.txt`やディレクトリ丸ごとcopyは行わない。`software_versions.json`の
  `analysis_generation_environment`は常に`"unavailable_from_source_artifacts"`固定文字列
  （batch/candidates成果物自体がgeneration-time software versionを記録していないため、
  report生成環境をanalysis生成環境と偽らない）。`platform.node()`（hostname）・
  `os.environ`のdump・絶対パスをoutputへ含めない。生成物に生成時刻を含めない
  （同じ入力から常に同一内容を再現できる設計）。output-dirの安全性チェックは
  `analysis/output_safety.py`で`batch`と共有。

## データ

- raw Dryad データ（data/raw/**, *.assoc.txt）は絶対にcommitしない。`git add -f` も禁止。
- .gitignore に `!tests/fixtures/*.assoc.txt` の例外があるため、raw ファイルを tests/fixtures/ に
  置くと除外が効かず追跡されてしまう。fixture は小さい合成データのみ。
- 6データセットを同時にメモリへ保持しない。解析は常に1データセットずつ処理する。
  複数データセット対応を求められても、ループ内で DataFrame を溜め込まない。
- tracked な plots/*.png と results/water_permeability/top_variants_by_region.tsv を
  不用意に上書きしない。再生成する場合は legacy wrapper を使い、実行前後で `git status`/
  `git diff` を確認してから commitする。
- results/water_permeability/top_snps.tsv の生成規則は未検証。これを生成したコードは
  git履歴上に存在しない。推測でこのファイルを変更・再生成しない。

## 統計・用語（表現の精度に注意する）

- `pval` は likelihood ratio test (LRT) p-value。`p_wald` / `p_score` も検証対象だが primary ではない。
- `1e-5`（`--threshold`）は legacy visualization threshold。Bonferroni補正や genome-wide
  significance と表現しない。`--threshold`は`adzuki-gwas-analyze diagnostics`の
  `--alpha`/`--fdr-level`とは無関係で、`diagnostics`が既存subcommandの挙動を変えることはない。
- multiple-testing familyは「1 dataset_id × manifest.pvalue_columns.primary × その
  ファイルでschema v1 validationを通過した全variant」に固定。6ファイル・Miyagi/Shumari・
  3形質・post-hoc regionを混ぜない。`m`はmanifestの`row_count`を無条件に使わず、実際に
  loadしたp-value数と一致することを確認する。
- BH-adjusted p-value(`pval_bh`)をStorey q-valueと呼ばない。BHの保証は独立検定または
  PRDS(positive regression dependency on a subset)下でのみ成立し、LDで相関するSNPに
  対して任意の依存構造での保証があるとは表現しない。
- λGCは診断値のみ。test statisticやp-valueをλGCで再補正しない。λGC単独で population
  stratification/kinship/batch effect/polygenicityの原因を断定しない。df=1の根拠は
  Dryad metadata・GEMMA univariate LMMのモデル記述であり、論文本文Methodsは未確認
  （README「Statistical Diagnostics」節参照）。
- `batch`は6データセットそれぞれを独立したfamilyとして扱う。合計8,187,994行という数字は
  処理件数の合計であり、6データセットを1つの補正対象（共通family）として扱った結果では
  ない。MiyagiとShumariの座標を比較・結合しない。
- config/water_permeability_regions.toml の窓は Manhattan plot を目視して事後に選んだ
  post-hoc visualization window。独立に定義されたQTL区間やLD blockであると断定しない。
- 本リポジトリは公開 summary statistics を可視化・再解析するだけで GWAS を再実行しない。
  "re-ran the GWAS" ではなく "re-analyzed / visualized published summary statistics" と書く。
- Miyagi / Shumari / Longxiaodou 4（姉妹repo adzuki-snp-pipeline）は互換性のない3つの座標系。
  liftover なしで座標を比較・アノテーションしない。
- liftover・LD解析・kinship補正・genomic-selection学習に必要なデータ・処理を、
  このリポジトリが現在持つsummary statisticsだけで行えるものと混同しない。

## Git / GitHub 操作

- 1 Issue = 1 branch = 1 PR。作業前後に `git status` を確認する。
- stage は明示したpathだけを使う（例: `git add -- README.md`）。
  `git add .` / `git add -A` / `git add --all` は禁止。
- push、Issue/PR の作成・更新・close、merge、force系操作（force push, reset --hard,
  branch削除等）は、明示的なユーザー指示・承認がある場合のみ行う。

## Claude Code設定の限界

- `.claude/settings.json` の `permissions` は完全なセキュリティ境界ではない。
  Read deny はビルトインのファイル読み取りツールと、認識可能な一部のBashコマンド形にのみ
  best-effortで適用され、Python/Node等の任意サブプロセス経由の読み取りは防げない。
