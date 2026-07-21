# HGVS notation to VCF coordinates API server

[English](README.md)

オフラインで HGVS 表記を VCF 座標（`CHROM POS REF ALT`）に変換する HTTP API サーバです。クエリ時に VEP・Mutalyzer・NCBI へはアクセスせず、事前に用意した転写産物モデル（cdot）と参照ゲノム（GRCh38）をローカルインデックスから読み出して変換します。

## できること

| 入力形式 | 例 | 転写産物の解決 |
|----------|-----|----------------|
| HGVSp（遺伝子シンボル） | `ALDH2:p.Glu504Lys` | シンボル → **MANE Select** |
| HGVSc（遺伝子シンボル） | `ALDH2:c.1510G>A` | シンボル → **MANE Select** |
| HGVSc（RefSeq） | `NM_000690:c.1510G>A` | バージョン省略時は **最新** |
| HGVSc（RefSeq、版指定） | `NM_000603.4:c.894T>G` | 指定どおり |
| HGVSc（Ensembl） | `ENST00000297494.3:c.894T>G` | 指定どおり / 最新 |

`p.` / `c.` / `n.` の置換は実装済みです。`del` / `dup` / `ins` / `delins` やイントロンオフセットはパース・マッピングまで対応し、インデルは `bcftools norm` 相当の左寄せ正規化を行います。

レスポンスには `vcf`（座標配列）、`ambiguous`（アミノ酸変化が複数の塩基変化に対応する場合）、`warnings`（参照塩基不一致など）が含まれます。

## クイックスタート

### 1. データの準備

バージョン付きファイル名はリリースごとに変わるため、ダウンロードは手動で行います。下記の取得元から最新ファイルを `data/` に置き、固定名へシンボリックリンクを張ります。

```bash
mkdir -p data && cd data

# 1) cdot 転写産物モデル — 最新版は https://github.com/SACGF/cdot/releases
#    (例: data_v0.2.33 → cdot-0.2.33.refseq.GRCh38.json.gz / .ensembl.GRCh38.json.gz)
ln -s cdot-0.2.33.refseq.GRCh38.json.gz   cdot.refseq.GRCh38.json.gz
ln -s cdot-0.2.33.ensembl.GRCh38.json.gz  cdot.ensembl.GRCh38.json.gz

# 2) MANE サマリー（監査・出典確認用。実行時には読みません）
#    最新版は https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
#    (例: MANE.GRCh38.v1.5.summary.txt.gz)
ln -s MANE.GRCh38.v1.5.summary.txt.gz     MANE.GRCh38.summary.txt.gz

# 3) 参照ゲノム — NCBI GCF_000001405.40（cdot の contig と一致する RefSeq アクセッション）
#    GCF_000001405.40_GRCh38.p14_genomic.fna.gz の取得先:
#    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
ln -s GCF_000001405.40_GRCh38.p14_genomic.fna.gz  GRCh38.fna.gz

cd ..
```

各ファイルの役割は [データソース](#データソース) を参照してください。

### 2. コンテナ（Docker / Podman）で起動する場合

```bash
docker build -t hgvs2vcf .

# LMDB インデックスを構築（データ更新のたびに再実行）
docker run -it --rm --user $(id -u):$(id -g) -v $(pwd)/data:/data hgvs2vcf index

# サーバー起動（初回は GRCh38.fna.gz の展開と samtools faidx を一度だけ実行、~3GB）
docker run -it --rm --user $(id -u):$(id -g) -v $(pwd)/data:/data -p 4567:4567 hgvs2vcf
```

Podman の場合は `podman` に置き換え、マウントに `:z` を付けます。

```bash
podman run -it --rm --userns=keep-id -v $(pwd)/data:/data:z -p 4567:4567 hgvs2vcf
```

`data/` にはシンボリックリンク、構築済み `index/`（LMDB）、展開済み `GRCh38.fna`（と `.fai`）が入ります。コンテナ内では `0.0.0.0` にバインドします（`BIND` で変更可）。

### 3. コンテナを使わずにインストールする場合

データの準備の手順は同じですが、コンテナを使わずローカル環境に Ruby と Bundler を入れて直接動かしたい場合は、こちらの手順に従ってください。

`bundle install` は Ruby gem だけを入れます。`lmdb` gem は C 拡張のため、**OS 側に LMDB ライブラリ（開発用ヘッダとリンク用ライブラリ）と C コンパイラが別途必要**です。参照 FASTA のインデックス作成には `samtools` も必要です。コンテナ（セクション 2）ではこれらがイメージに含まれています。

```bash
# 例: Debian/Ubuntu
sudo apt install build-essential liblmdb-dev samtools

# 例: macOS (Homebrew)
brew install lmdb samtools
```

```bash
bundle install
ruby bin/build_index.rb --out data/index --build GRCh38 \
     data/cdot.refseq.GRCh38.json.gz data/cdot.ensembl.GRCh38.json.gz
gunzip -kc data/GRCh38.fna.gz > data/GRCh38.fna && samtools faidx data/GRCh38.fna
INDEX_DIR=data/index GENOME_FASTA=data/GRCh38.fna bundle exec ruby app.rb
```

## API

| メソッド | パス | 説明 |
|----------|------|------|
| GET | `/health` | 生存確認（データ不要） |
| GET | `/decode?hgvs=...` | 単一 HGVS をデコード |
| POST | `/decode` | JSON `{"hgvs":["...", ...]}` でバッチデコード |

**生存確認**

```bash
curl localhost:4567/health
# -> ok
```

**単一デコード** — `>` は URL でエンコードが必要です。

```bash
curl -G localhost:4567/decode --data-urlencode 'hgvs=NM_000603.4:c.894T>G'
```

```json
{"input":"NM_000603.4:c.894T>G","transcript":"NM_000603.4","gene":"NOS3",
 "vcf":[{"chrom":"NC_000007.14","pos":150999023,"ref":"T","alt":"G"}],
 "ambiguous":false,"warnings":[]}
```

**バッチデコード** — `Content-Type: application/json` を指定してください。

```bash
curl -X POST localhost:4567/decode \
     -H 'Content-Type: application/json' \
     -d '{"hgvs":["ALDH2:p.Glu504Lys","ALDH2:c.1510G>A","NM_000690:c.1510G>A"]}'
```

項目ごとに独立してデコードし、エラーは `{"input":...,"error":...}` として返します。

## 背景

### アーキテクチャ

```
                data/ の準備（手動ダウンロード + シンボリックリンク、クイックスタート参照）
                        │
     ┌──────────────────┴───────────────────┐
     ▼                                        ▼
 cdot JSON (RefSeq+Ensembl, GRCh38)     GRCh38 FASTA + .fai
 エクソンブロック、CDS、鎖、             (contig id は RefSeq
 MANE/canonical タグ、                   アクセッション: NC_0000..)
 転写産物↔ゲノム gap CIGAR
     │                                        │
     └──────────────┬─────────────────────────┘
                    ▼
              Decoder (lib/hgvs2vcf)
   Parser → 転写産物解決 → Mapper(c.→g.) → Sequence(参照塩基 + 正規化)
                    ▼
             Sinatra API  (GET/POST /decode)
```

### 転写産物バックボーンに cdot を使う理由

`biocommons/hgvs` + UTA は参照実装ですが Postgres が必要で、配列も自前で行います。**[cdot](https://github.com/SACGF/cdot)**（`SACGF/cdot`、MIT）は RefSeq GFF と Ensembl GTF を、プロジェクタに必要な情報だけのフラット JSON に事前変換します。特に **転写産物↔ゲノム gap**（`"M196 I1 M61 …"`）を持つため、参照ゲノムと完全一致しない RefSeq 転写産物も正しくマッピングできます。素朴な GFF エクソンプロジェクタで最も起きやすい off-by-N エラーの主因を避けられます。執筆時点の最新データリリース: **`data_v0.2.33`**（2026-06-26; Ensembl 116, RefSeq RS_2025_08 / GCF_000001405.40 GRCh38.p14）。MANE タグも同梱されています。

cdot のエクソンタプル: `[g_start(0), g_end(excl), exon_number, cdna_start(1), cdna_end(1), gap|null]`。エクソンは **ゲノム座標順** に並ぶため、マッパーは配列位置ではなく cDNA 区間でインデックスします（マイナス鎖で重要）。

### データソース

いずれも定期的に更新されますが、デコード時に **ライブ照会はしません**。

| 役割 | 取得元 | 使用タイミング |
|------|--------|----------------|
| 転写産物モデル（RefSeq + Ensembl、GRCh38） | [cdot releases](https://github.com/SACGF/cdot/releases) | インデックス構築 |
| 参照ゲノム | [NCBI GCF_000001405.40_GRCh38.p14](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/) — `GCF_000001405.40_GRCh38.p14_genomic.fna.gz`（cdot の `contig` と同じ RefSeq アクセッション） | サーバー実行 |
| MANE サマリー（監査・独立したシンボル対応表） | [NCBI MANE_human/current](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/) | 監査のみ |
| HGNC complete set（任意） | シンボルエイリアス・廃止シンボル | デフォルトでは未使用 |

実行時に読むのは LMDB インデックスとゲノム FASTA のみです。cdot JSON はインデックス構築時のみ必要で、MANE ファイルは出典確認用に保持します。

### 知っておくべき 3 つの判断

1. **HGVSp は本質的に曖昧です。** アミノ酸置換は複数の塩基変化に対応することがあります。参照コドン（転写産物から読み取り）は *どの参照コドンか* を特定しますが、*どの変異コドンか* までは常に特定できません:
   - `ALDH2:p.Glu504Lys`、参照コドン `GAG` → 一意に `c.1510G>A`（GAG→AAG）。
   - `NOS3:p.Asp298Glu`、参照コドン `GAT` → 最小 SNV が **2 つ**: `c.894T>A` と `c.894T>G`。（既知の rs1799983 は T>G ですが、HGVSp だけでは断定できません。）
   デコーダーは最小候補を **すべて** 返し `ambiguous: true` にします。2 塩基以上の変化が必要な候補は `requires_mnv` とし MNV/delins として出力します。運用方針（すべて返す・1塩基を優先・拒否）は利用側で決めてください。

2. **シンボルやアクセッションのみの入力では MANE Select をデフォルトにします。** cdot タグから解決し、MANE のない遺伝子（非コードング、一部パッチ遺伝子）は Ensembl canonical にフォールバックします。`ALDH2` → `NM_000690.4`、`NOS3` → `NM_000603.x`。`NM_000690` のみ → 読み込み済みの最高バージョン。

3. **ENST も同様に扱え、むしろ単純です。** Ensembl/GENCODE 転写産物はゲノム由来のため `gap` は常に null — マイナス鎖とスプライス算術が本質です。cdot の *ensembl* ファイルを refseq と併せて読み込めば、ENST も同じ手順で解決します。

### 実装の注意点（正解集合で検証すること）

- **CDS 基点 → cDNA オフセット。** cdot の `start_codon`/`stop_codon` は転写産物上のオフセットです。コードでは `c.1 = cdna(start_codon + 1)` としています。±1 の扱い、UTR（`c.-N`、`c.*N`）、イントロン（`c.N+M`）は off-by-one が起きやすい箇所です。
- **鎖方向。** マイナス鎖では HGVS の ref/alt は転写産物座標系のため、参照塩基チェックと正規化の前にゲノム座標系へ補完します。
- **参照塩基チェック。** すべての SNV で、暗黙の参照塩基をゲノム FASTA と照合し、不一致は `warnings` に記録します（転写産物バージョン違い、ビルド混在、実際の転写産物/ゲノム不一致など）。
- **正規化。** `Normalize.left_align` は `bcftools norm`（トリム + 左シフト）相当です。HGVS は 3′ シフト、VCF は左寄せのため、インデルは座標がずれます。

### オンディスクインデックス（LMDB）— 起動時に JSON を RAM に載せない

cdot JSON は **一度だけオフラインで** LMDB 環境（`bin/build_index.rb`）に変換します。サーバーはクエリごとに小さなパック済みレコードを mmap から読み出すだけで、巨大 JSON はヒープに載りません。

1 つの環境に 2 つの名前付き DBI:

- `kv` — 解決とメタデータ: `mane:<SYMBOL>`、`canon:<SYMBOL>`、`latest:<BASE>`、`meta:contigs`（id→名前）、`meta:build`
- `tx` — `tx_id` → パック済み座標レコード

**パック済みレコードレイアウト**（リトルエンディアン。`lib/hgvs2vcf/packing.rb` 参照）:

```
magic u8 | strand u8 | coding u8 | contig_id u16 |
cds_start_i u32 | cds_end_i u32 | gene_len u16 | gene[gene_len] |
n_exons u16 |
  n_exons × { g_start u32, g_end u32, cdna_start u32, cdna_end u32, gap_len u16, gap[gap_len] }
```

クエリ = mmap からの 1 レコード読み出し + エクソン cDNA 区間への二分探索 + 算術（gap CIGAR は gap を持つ転写産物のみ）。構築は LMDB の単一書き込みトランザクション（原子的コミット）。サーバーは `nolock: true` でインデックスを開き（`lock.mdb` への書き込みも含め書き込みなし）、静的インデックスを複数プロセスで共有可能にしています。

## テスト

```bash
ruby spec/engine_test.rb   # パーサー、HGVSp 逆翻訳、マッパー、正規化
ruby spec/lmdb_test.rb     # インデックス構築 → LMDB → デコード
```

## データ更新時の確認

正しさはデータリリースに依存します。cdot や参照ゲノムを更新したら、代表例で既知の正解と突き合わせてください。

- ClinVar `variant_summary.txt` や cdot 同梱の `clinvar_hgvs_*.tsv`
- VariantValidator / Mutalyzer によるスポットチェック（マイナス鎖、gap 付き転写産物、エクソン境界、`c.-N` / `c.*N` など）

不一致率が急に上がった場合は、コードより転写産物バージョンや MANE の再割り当てを疑ってください。

## 環境変数

| 変数 | デフォルト | 説明 |
|------|------------|------|
| `INDEX_DIR` | `data/index`（ローカル）/ `/data/index`（コンテナ） | LMDB インデックスのパス |
| `GENOME_FASTA` | `data/GRCh38.fna` / `/data/GRCh38.fna` | 参照 FASTA（`.fai` 必須） |
| `BUILD` | `GRCh38` | ゲノムビルド名 |
| `BIND` | `0.0.0.0` | バインドアドレス |
| `PORT` | `4567` | 待ち受けポート |

## Author

Toshiaki Katayama が最初 Claude (Opus 4.8) で開発し、Cursor エージェントでデバッグして動くように調整した。
