---
marp: true
theme: gaia
_class: lead
paginate: true
style: |
    h1, h2, h3, h4, h5, header, footer {
        color: black;
    }
    section {
        # background-color: white;
        font-family: 'Noto Sans CJK JP Regular';
        color: black;
        font-size: 18px;
    }
# backgroundColor: #fff
# backgroundImage: url('https://marp.app/assets/hero-background.jpg')
---

![bg left:40% 80%](https://marp.app/assets/marp.svg)

# **Marp**

Markdown Presentation Ecosystem

https://marp.app/

---

# How to write slides

Split pages by horizontal ruler (`---`). It's very simple! :satisfied:

```markdown
# Slide 1

foobar

---

# Slide 2

foobar
```
---
<style>
h1 {font-size: 40px;}
{font-size: 20px;}
</style>
# 自己紹介

* yu340102
* とある高専卒
* T大M2
* 生活道路での自動運転実現に関する研究
* 研究では機械学習とかあまりしていない
* 大学院の授業でデータサイエンスに関する授業を取る

---

# きっかけ

* 就活で athletics に登録 -> Slack でコンペのことを知る
* 第4回Brain(s)に参加 -> 序盤に一瞬LB載る程度
* Award会に参加 -> めちゃくちゃ化学に詳しい人ではなく，色々試した人が入賞していた  
-> 第5回（今回）はリベンジのつもりで入賞目指して，参加した．

---

# 参加するにあたって

* テーブルデータ分析の基本を一通りやるつもりでいた
* ドメイン知識は分からないから考えないことにした
* 全力出し切って，特徴量生成～モデルのチューニングまで一通りできたら感無量というスタンス

---

# 結論

**情報収集×試行回数=>1st**
* 基礎情報
    * Kaggleで勝つデータ分析の技術
    * 各ライブラリのDocument
* ドメイン情報
    * Slackの過去メッセージ
    * 第4回解法ブログ
    * RDKit, mordred, FingerPrintに関するQiita等のまとめ記事
* 試行回数
    * とにかくsub
    * 何かやれることはないか考える
**=1st**

---

# ベストスコア

画像を挿入
* 特徴量は291
* モデルはMLPRegressor(from sklearn)
* 5-fold CV の平均値を提出

---

# 概要

* 特徴量生成
    * 記述子（RDKit, mordred）
    * フィンガープリント（Morgan FP, MACCS Keys）
    * Count（各原子や記号（=,-,+等），SMILESの文字列長）
* 特徴量選択
    * 0のカラムを除去
    * LightGBMのfeature importance(gainで高いもの)
    * RDKitの記述子 + MACCS Keys + Count = 291変数
* モデル選択
    * 基本LGBM
    * sklearnにMLPあること知って後から採用
    * パラメータはどっかのタイミングでColab + Optunaでやったもの
* アンサンブル
    * 5-fold CV の平均
    * LGBM+MLP+SVR->LinearReg のstackingとかも試したけどベストスコアではなかった

---

# 大変だったこと

* 一回も中間ランキングに載らなかった
* trackのスコアばかり良くなる  
-> 0.15台
* 採点ルール分かっていなかった  
-> 画像挿入

---

# まとめと感想

* 色々試したらたまたまいいスコアが出た
* 今後はもっと**分析**して結果を出したい
* リベンジできてめちゃくちゃ嬉しい
* 化学分からなくても何とかなった
* ルールはよく確認したい
* (データ分析コンペは無料のソシャゲ...？ 本業に支障がないようにしたい)
* 解法について，ブログと書いてみたいなあ

---

# memo

＜発表資料＞
　①自己紹介：実名は伏せてください。●●大 HN　で大丈夫です。
　　　　　　　研究テーマなど
　②このコンテストを知ったきっかけ／期待していたこと
　③解法説明
　　・使用モデル／工夫／苦労／息抜き方法など
　④まとめ／感想