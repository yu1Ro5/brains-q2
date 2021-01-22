---
marp: true
theme: future
header: "**Eye Tracker 4C 引継ぎ**"
footer: "by **Yuichiro Shimizu**"
paginate: True

---
<style>
h1 {font-size: 40px;}
{font-size: 20px;}
</style>
# 前提

* PC、Eye Tracker 4C が手元にある。
* ソフト類のインストールができる
* やっていたことをある程度知っている。

上記を満たす方向けの資料です。

---
<!--

-->
# 大まかな流れ

1. 必要なソフトのインストール
2. 実行環境の整備
3. コードの実行

---
# 確認

"GazeViz4C" フォルダを自身のPCにコピーし、以下の構成になっていることを確認する。

GAZEVIZ4C
│  base.md
│  README.md
│  sms.yml
│  
├─code
│     eyetracker4c.py
│     gazeheatplot.py
│     run.py
│     step_by_step.py
│          
├─data
│      
├─image
│      
└─reference
        

---
# 必要なソフトのインストール

* Tobii Eye Tracking Core Software
    * Eye Tracking のために必要なコア機能を提供するソフト
    * 4C からデータを取得するために必要
    * キャリブレーションもこれで行う
* Anaconda
    * Python 本体と実行のために必要なパッケージを簡単に管理できるソフト
    * コードを実行するために必要

上記二点のソフトを公式HPもしくは、"installer"フォルダの.exeファイル実行してインストールする。

---
# 実行環境の整備

## Python 環境の構築
1. Windowsキーを押し、Anaconda Prompt (anaconda3) と入力し、実行する。
2. 下記を順番に入力しエンター  
    1. cd GazeViz4C/ （指定のフォルダに移動  ）
    2. conda env create -f=sms.yml （必要な環境を再構築）  
    3. conda activate sms （実行環境をアクティベート）

## Tobii Eye Tracker のキャリブレーション
* reference\tobii_Eye_Tracking.pdf の8番までを参考に設定を完了させる．
* コードを実行する際は，有効になっていることを確認する．

---
# コードの実行1

* 場所：code/step_by_step.py
* 概要：4Cを認識し，5秒間ディスプレイ上の左右の目のXY座標をコマンドライン上に出力する.
* コマンド
cd code/
python step_by_step.py
* イメージ
![w:600](image/step_by_step.png)

---
# コードの実行2-1

* 場所：code/run.py
* 概要
    1. 4Cを認識
    2. 画面のスクリーンショットを撮影
    3. 30秒間の視線位置座標をcsvファイルに記録
    4. 2.で撮影した画像に3.で保存したデータをもとにHeatMapを描画
    （data/ にスクリーンショット画像，csvファイル，HeatMap画像を保存）
* コマンド
cd code/
python run.py
---
# コードの実行2-2

* イメージ
![w:800](image/sample_output.png)
---
# なんかうまくいかなかったら

* tobii eye tracking が実行されているか確認．
* 視行動が認識されているか確認．
タスクバーのtobii eye tracking のアイコンが××になっていると測定できる角度や位置に座れていない可能性が高い．
* pythonの仮想環境がアクティベートされているか確認．
(sms)が頭に表示されいていないと，実行できる環境がアクティベートされていない．
'conda activate sms' を実行しよう．
* カレントディレクトリを確認．
/GazeViz4C/code に移動した状態での実行を想定している．現在位置を確認しよう．