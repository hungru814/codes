#!/usr/bin/env python3

###  2019/03/22 藤木研人
## 【背景】
## R6L3陳さん向けプログラム 陳さん作成のPOSCAR内の空間にLiを敷き詰めるプログラムとの併用を前提とする。
## LiのみのPOSCARの原子座標をクラスタリング。周期境界を跨いで原子が繋がっているか判断する。
## これにより、Liのイオン伝導パスの有無を簡易的に判断する。
## 【動作内容】
## 原子のクラスタリング
## クラスタ毎に周期境界をを跨ぐ経路があるかを判断
## 【動作確認環境】
## Anaconda 仮想環境下
## Python 3.6.4 |Anaconda, Inc.| (default, Mar 13 2018, 01:15:57) [GCC 7.2.0] on linux
## pymatgen '2018.5.22'
###

import sys
import numpy as np
from pymatgen.core.structure import IStructure

fileName=sys.argv[1]

# print("reading data")

## POSCAR（Liの原子のみのデータ）読み込み
## 原子数,原子のインデックス,原子間距離の最小値を格納
structure = IStructure.from_file(fileName)
points=structure.cart_coords
atom=len(points)*["Li"]
num=np.arange(0,len(atom))
num=num.reshape(1,-1)
minDist=structure.distance_matrix[structure.distance_matrix.nonzero()].min()

## クラスタリング
# print("clustering start")

## minDist内の原子リストの行列,変数の初期化
trueList=np.array(np.where(structure.distance_matrix<minDist*4)).T
searchList=num[0]
clusterArray=[]
count=0
unionFlag=0
unionList=[]
unionIdx=[]

## searchListが空になるまで回す。(原子インデックスの順にループ)
while len(searchList) > 0:
    unionList=[]
    unionIdx=[]
    unionFlag=0
    count+=1
    nearAtomList=np.unique(trueList[np.any(trueList==searchList[0],axis=1)].flatten())

    ## クラスター行列を全探索して、同じクラスターの原子がないかチェック
    ## ある場合、unionListに同クラスターである原子を追加
    for i in range(len(clusterArray)):
        if set(nearAtomList).isdisjoint(clusterArray[i]) == False:
            unionIdx.append(i)
            unionList=np.append(unionList,clusterArray[i])
            unionFlag=1

    ## 同一クラスターのフラグがある場合   操作の手順上、インデックスリストunionIdxを逆順に
    if unionFlag==1:
        unionIdx=unionIdx[::-1]

    ## unionIdxに値がある場合    clusterArrayのunionIdx番目を削除
    while len(unionIdx) > 1:
        del clusterArray[unionIdx[0]]
        unionIdx=np.delete(unionIdx,0)

    ## unionIdxがひとつになった場合
    ## 探索中の原子からminDist内に存在する原子インデックスリストnearAtomListと
    ## nearAtomList内の原子からminDist内に存在する原子リストunionListを結合し
    ## clusterArrayに格納
    if len(unionIdx) == 1:

        nearAtomList=np.unique(np.append(nearAtomList,unionList))
        clusterArray[unionIdx[0]]=nearAtomList

    ## clusterArray内に同クラスターの原子が存在しない場合
    ## clusterArrayに格納
    else:
        clusterArray.append(nearAtomList)

    ## 探索リストsearchListを更新。クラスタリングされた原子を排除して効率化。フラグの初期化
    searchList=np.setdiff1d(searchList,nearAtomList)


# print("clustering finish")

## 経路の有無判定

# print()
# print("judge start")
count=0 #経路があると判断されたクラスタの数

## clusterArrayの要素数ループ
## クラスタ毎にabc軸に対して、Liが分割した点数分存在するかをチェック
frac=structure.frac_coords
for n in range(len(clusterArray)):
    ## frac_coordsを軸ごとに分割
    cl0=np.array([frac.T[:,int(j)] for j in clusterArray[i]])
    cl_a=cl0.T[0]
    cl_b=cl0.T[1]
    cl_c=cl0.T[2]
    ## 軸ごとのLiの最大点数を算出(切り捨て面倒なので-1にしている)
    num_a=round(structure.lattice.a/0.2)-1
    num_b=round(structure.lattice.b/0.2)-1
    num_c=round(structure.lattice.c/0.2)-1

    ## 軸ごとにLiの点数が最大数あるかチェック
    if num_a <= len(set(cl_a)):
        print(cl_a)
        count+=1
    if num_b <= len(set(cl_b)):
        print(cl_b)
        count+=1
    if num_c <= len(set(cl_c)):
        print(cl_c)
        count+=1


# print("judge finish")
# print()
if count > 0:
    print(1) # 経路が存在する
else:
    print(0)
