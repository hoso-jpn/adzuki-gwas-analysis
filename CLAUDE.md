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
  `analysis/loader.py` は検証成功後にのみ呼ばれる別のDataFrame loaderで、この2層分離が
  「6データセットを同時にメモリへ保持しない」保証の実体。統合しない。
- `scripts/01`-`04`: `analysis/` の薄いbackward-compatibleラッパー。tracked な
  `plots/`/`results/water_permeability/` の既存ファイル名・配置を再現する経路。
- `adzuki-gwas-analyze`（unified CLI）の `all`/`regions` は単一 `--output-dir` へまとめて
  書き出す bundle コマンドであり、tracked成果物の再現コマンドではない（README参照）。
  実データで動作確認するときは常に一時 output directory（`mktemp -d` 等）へ出力し、
  `plots/`/`results/` を直接指定しない。

## データ

- raw Dryad データ（data/raw/**, *.assoc.txt）は絶対にcommitしない。`git add -f` も禁止。
- .gitignore に `!tests/fixtures/*.assoc.txt` の例外があるため、raw ファイルを tests/fixtures/ に
  置くと除外が効かず追跡されてしまう。fixture は小さい合成データのみ。
- 6データセットを同時にメモリへ保持しない。解析は常に1データセットずつ処理する。
  複数データセット対応を求められても、ループ内で DataFrame を溜め込まない。
- tracked な plots/*.png と results/water_permeability/top_variants_by_region.tsv を
  不用意に上書きしない。再生成する場合は legacy wrapper を使い、実行前後で `git status`/
  `git diff` を確認してから committする。

## 統計・用語（表現の精度に注意する）

- `pval` は likelihood ratio test (LRT) p-value。`p_wald` / `p_score` も検証対象だが primary ではない。
- `1e-5`（`--threshold`）は legacy visualization threshold。Bonferroni補正や genome-wide
  significance と表現しない。本リポジトリは多重検定補正を一切計算しない。
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
