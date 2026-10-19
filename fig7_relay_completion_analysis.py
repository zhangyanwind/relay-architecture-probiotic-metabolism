#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 7 — 益生菌作为 relay「第一棒」的定量后果  (Plan C)
=====================================================================
为什么不用旧 Fig 12（9 个 RCT 柱状图）：
  旧图的 9 条 RCT 记录（GERD/BAM/T2D/PPI/BL-99/UA/FBD/BLa80/NAFLD）
  是从文献综述手工整理的硬编码常量，不是本项目351 株基因组算出来的结果，
  且其"干预改变接力中间体"的主张无法用基因组数据检验。
  按 Plan C（阴性/不可支撑结果不得留在正文当证据），替换为基于
  自有泛基因组数据的真实定量分析。

本图回答的真实问题：
  如果益生菌只能做 relay 的第一棒，那么
  (1) 每个通路缺口的绝对大小有多大？
  (2) 现有的商业化益生菌能否单菌株补齐一条通路？
  (3) 最小菌株组合(combination)需要几株才能补齐？各方案成本多少？
输出： fig7_relay_completion_stats.tsv
"""
import pandas as pd, numpy as np, itertools, json

BASE = 'D:/Probiotics_Pangenome_v5'
OUT = f'{BASE}/validation/analysis'

# ---------- 1. 读取 relay 性状矩阵 ----------
rt = pd.read_csv(f'{BASE}/pangenome/direction_matrices/relay_trait_matrix.tsv', sep='\t')

# relay_trait_matrix holds per-unit CARRIER FRACTIONS (0.12, 0.46, 0.94, ...),
# not 0/1 calls. Positivity is the project-wide threshold of >=60% of a unit's
# genomes; astype(bool) would wrongly score every fractional carrier.
POS = 0.60

# 五条通路的上下游定义（与正文一致；色氨酸在 trait matrix 中
# 用TrpA/B 标记基因表示，此处按 unit_marker_presence 单独处理）
PATHWAYS = {
    'Lactate-SCFA':  ('LDH/MDH（产乳酸）', 'buk（产丁酸）'),
    'Bile acid':     ('BSH（胆盐水解）',   'bai 操纵子（7α-脱羟）'),
    'TMA/TMAO':      ('CutC（产 TMA）',    'TorA（TMAO 呼吸）'),
    'GABA-Glutamate':('GadA/B（产 GABA）', 'GabT（消耗 GABA）'),
}

rows = []
for pw, (up, dn) in PATHWAYS.items():
    U = (rt[up].values >= POS).astype(int)
    D = (rt[dn].values >= POS).astype(int)
    n = len(rt)
    full = int(((U == 1) & (D == 1)).sum())
    up_only = int(((U == 1) & (D == 0)).sum())
    dn_only = int(((U == 0) & (D == 1)).sum())
    neither = int(((U == 0) & (D == 0)).sum())
    rows.append(dict(pathway=pw, n_units=n, up_trait=up, down_trait=dn,
                     up_pos=int(U.sum()), down_pos=int(D.sum()),
                     complete=full, up_only=up_only, down_only=dn_only, neither=neither,
                     up_frac=U.mean(), down_frac=D.mean(),
                     ratio=(U.mean() / D.mean()) if D.mean() > 0 else np.nan))
res = pd.DataFrame(rows)

print('=' * 74)
print('Fig 7 基础表：各通路接力完整性（40 单元）')
print('=' * 74)
print(res[['pathway', 'up_pos', 'down_pos', 'complete', 'up_only',
           'down_only', 'neither', 'up_frac', 'down_frac']].to_string(index=False))

# ---------- 2. 缺口定量：需要多少“伙伴”才能补齐 ----------
#对每条通路，统计下游阳性单元中，有多少也拥有上游（=能独立完成接力），
# 有多少是“孤立的下游携带者”（需要外部伙伴供给中间体）
print()
print('=' * 74)
print('下游携带者中“需要外部中间体供给”的比例')
print('=' * 74)
gap = []
for _, r in res.iterrows():
    dn_pos = int(r.down_pos)
    # 能独立完成的 = complete（上下游都有）
    orphan = dn_pos - int(r.complete)
    frac = orphan / dn_pos if dn_pos > 0 else np.nan
    gap.append(dict(pathway=r.pathway, downstream_carriers=dn_pos,
                    self_completing=int(r.complete), orphan_carriers=orphan,
                    orphan_frac=frac))
    print(f'{r.pathway:16s} 下游携带 {dn_pos:2d}  自给自足 {int(r.complete):2d}  '
          f'需外部中间体 {orphan:2d}  ({frac*100 if frac==frac else float("nan"):.0f}%)')
gapt = pd.DataFrame(gap)

# ---------- 3. 最小补齐组合：贪心集合覆盖 ----------
print()
print('=' * 74)
print('最小补齐菌株组合（贪心集合覆盖，目标是每条通路上下游俱全）')
print('=' * 74)

unit_cols = list(PATHWAYS.values())
# 单元 x 通路 的完成向量
feat = {}
for i, row in rt.iterrows():
    f = []
    for pw, (u, d) in PATHWAYS.items():
        f.append(1 if (row[u] >= POS and row[d] >= POS) else 0)
    feat[int(row.unit_id)] = np.array(f)

allpw = list(PATHWAYS.keys())
best = []
cur = np.zeros(len(allpw), dtype=int)
remaining = list(feat.keys())
chosen = []
step = 0
while cur.sum() < len(allpw) and remaining:
    step += 1
    gains = {u: int((feat[u] & ~cur.astype(bool)).sum()) for u in remaining}
    top = max(gains, key=gains.get)
    if gains[top] == 0:
        print(f'  step {step}: 剩余通路 {list(np.array(allpw)[cur==0])} 无候选可补齐')
        break
    chosen.append((top, gains[top]))
    cur |= feat[top]
    remaining.remove(top)
    print(f'  step {step}: + unit {top:2d} ({rt.loc[rt.unit_id==top,"cn_name"].values[0]}) '
          f'新增 {gains[top]} 条通路 -> 已补齐 {cur.sum()}/{len(allpw)}')

covered = [allpw[i] for i in range(len(allpw)) if cur[i]]
print(f'\n  最少菌株数 = {len(chosen)}即可使{len(covered)}/{len(allpw)} 条通路完整:')
for pw in covered:
    print(f'      - {pw}')

# ---------- 4. 覆盖效率：单菌株vs 多菌株 ----------
print()
print('=' * 74)
print('覆盖效率对比')
print('=' * 74)
single = [int(v.sum()) for v in feat.values()]
print(f'  单菌株最多覆盖 {max(single)} 条通路（{sum(1 for s in single if s==max(single))}/{len(single)} 株达到此值）')
print(f'  单菌株平均覆盖 {np.mean(single):.2f} 条')
print(f'  {len(chosen)} 菌株组合覆盖 {cur.sum()} 条 =单菌株最优的 '
      f'{cur.sum()/max(single):.1f} 倍')

# ---------- 5. 成本效率：缺口越大的通路越需要多菌株 ----------
print()
print('=' * 74)
print('“第一棒”后果：下游缺失程度 vs 需伙伴程度')
print('=' * 74)
from scipy.stats import spearmanr
x = res.up_frac.values
y = res.down_frac.values
rho, p = spearmanr(x, y)
print(f'  通路内上游携带率 vs 下游携带率: Spearman rho={rho:.3f}, p={p:.3f} (n={len(res)})')
for _, r in res.iterrows():
    print(f'  {r.pathway:16s} 上游 {r.up_frac*100:5.1f}%  下游 {r.down_frac*100:5.1f}%  '
          f'→缺口 {r.up_frac*100-r.down_frac*100:+5.1f} pt')

# ---------- 保存 ----------
res.to_csv(f'{OUT}/fig7_relay_completion_stats.tsv', sep='\t', index=False)
gapt.to_csv(f'{OUT}/fig7_orphan_carrier_stats.tsv', sep='\t', index=False)
with open(f'{OUT}/fig7_minimal_combination.json', 'w', encoding='utf-8') as f:
    json.dump({'chosen_units': [int(u) for u, _ in chosen],
               'chosen_names': [rt.loc[rt.unit_id == u, 'cn_name'].values[0] for u, _ in chosen],
               'pathways_covered': covered,
               'n_pathways_total': len(allpw),
               'max_single_unit_coverage': int(max(single)),
               'mean_single_unit_coverage': float(np.mean(single)),
               'combination_amplification': float(cur.sum() / max(single)),
               'spearman_up_down_prevalence': float(rho),
               'spearman_p': float(p)}, f, ensure_ascii=False, indent=2)
print('\nwritten: fig7_relay_completion_stats.tsv / fig7_orphan_carrier_stats.tsv / fig7_minimal_combination.json')