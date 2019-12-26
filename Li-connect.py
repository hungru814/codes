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
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
import itertools
import pickle
import os
import copy


# print("reading data")

## POSCAR（Liの原子のみのデータ）読み込み
## 原子数,原子のインデックス,原子間距離の最小値を格納


minDist = 3.0



#################################################################################
## クラスタリング
# print("clustering start")

## minDist内の原子リストの行列,変数の初期化
def clustering(structure):
    trueList=np.array(np.where(structure.distance_matrix<minDist)).T
    clusterArray=[]
    count=0
    unionFlag=0
    unionList=[]
    unionIdx=[]


    points=structure.cart_coords
    atom=len(points)*["Li"]
    num=np.arange(0,len(atom))
    num=num.reshape(1,-1)

    searchList=num[0]



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
#        searchList=np.setdiff1d(searchList,nearAtomList)
        searchList=searchList[1:]
    return(len(clusterArray))
# print("clustering finish")

#############################################################################




with open('/home/chen/MP_data_2019_Mar/mp_structure.pkl.pd_','rb') as f :
    mp_structure = pickle.load(f)


with open('/home/chen/MP_data_2019_Mar/mp_inorganic.pkl.pd_','rb') as g :
    mp_inorganic = pickle.load(g)




for i in range(0,len(mp_inorganic)):
    if 'Li' in mp_inorganic.iloc[i]['elements'] and 'O' in mp_inorganic.iloc[i]['elements']  \
       and mp_inorganic.iloc[i]['e_above_hull'] < 0.025 :

        e_above_hull = mp_inorganic.iloc[i]['e_above_hull']
        formula = mp_inorganic.iloc[i]['pretty_formula']
        m_id = mp_inorganic.iloc[i].name
        ss = SpacegroupAnalyzer(mp_structure.loc[m_id]['structure']).get_conventional_standard_structure()

        remove_list=[]

        for j in range(0,len(ss)):
            if ss.sites[j].specie.symbol != 'Li' :
                remove_list.append(j)
        ss.remove_sites(remove_list)
        super_x = copy.deepcopy(ss)
        super_y = copy.deepcopy(ss)
        super_z = copy.deepcopy(ss)

        super_x.make_supercell([2,1,1])
        super_y.make_supercell([1,2,1])
        super_z.make_supercell([1,1,2])

        if clustering(ss) == clustering(super_x) and clustering(ss) == clustering(super_y)  \
           and clustering(ss) == clustering(super_z):

            print(clustering(ss))

            mp_structure.loc[m_id]['structure'].to("poscar",'POSCAR-'+m_id)
            with open('entries.dat','at') as fff:
                fff.write(m_id+'   '+formula + '    ' + str(e_above_hull) + '\n')

