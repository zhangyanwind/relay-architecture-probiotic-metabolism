#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""疾病断裂方向一致性的分层统计（Plan C：真实计算，替代硬编码矩阵）"""
import pandas as pd, math, numpy as np
from scipy import stats
from scipy.stats import binomtest

A = 'D:/Probiotics_Pangenome_v5/validation/analysis'
m = pd.read_csv(f'{A}/relay_metabolite_disease_stats_ibdmdb.tsv', sep='\t')
m['lfc_cd'] = m.apply(lambda r: math.log2((r.CD_mean + 1) / (r.nonIBD_mean + 1)), axis=1)
m['lfc_uc'] = m.apply(lambda r: math.log2((r.UC_mean + 1) / (r.nonIBD_mean + 1)), axis=1)
m['pmin'] = m[['p_CD_vs_nonIBD', 'p_UC_vs_nonIBD']].min(axis=1)
sig = m[m.pmin < 0.05].copy()

print(f'=== 分层: 显著代谢物 (n={len(sig)}/{len(m)}) ===')
for role in ['substrate', 'intermediate', 'end_product']:
    s = sig[sig.role == role]
    print(f'{role:13s} sig {len(s):2d}/{len(sig[sig.role==role]):2d}  '
          f'medianCD={s.lfc_cd.median():+.2f}  medianUC={s.lfc_uc.median():+.2f}')

si = sig[sig.role.isin(['substrate', 'intermediate'])]
ep = sig[sig.role == 'end_product']
u, p = stats.mannwhitneyu(si.lfc_cd.abs(), ep.lfc_cd.abs())
print(f'\nCD |log2FC| substrate+intermed = {si.lfc_cd.abs().median():.2f} (n={len(si)}) '
      f'vs end_product = {ep.lfc_cd.abs().median():.2f} (n={len(ep)})')
print(f'  Mann-Whitney p = {p:.4f}')

for tag, col in [('CD', 'lfc_cd'), ('UC', 'lfc_uc')]:
    s2 = sig.copy()
    up_i = ((s2.role.isin(['substrate', 'intermediate'])) & (s2[col] > 0)).sum()
    n_i = s2.role.isin(['substrate', 'intermediate']).sum()
    dn_e = ((s2.role == 'end_product') & (s2[col] < 0)).sum()
    n_e = (s2.role == 'end_product').sum()
    b = up_i + dn_e
    bt = binomtest(b, len(s2), 0.5, alternative='greater').pvalue
    print(f'{tag}: intermediate/substrate UP {up_i}/{n_i} ; end_product DOWN {dn_e}/{n_e} '
          f'-> concordant {b}/{len(s2)}, binom p = {bt:.4g}')

# 按通路看方向
print('\n=== 按通路 (CD) ===')
for pw in sig.pathway.unique():
    s = sig[sig.pathway == pw]
    si2 = s[s.role.isin(['substrate', 'intermediate'])]
    ep2 = s[s.role == 'end_product']
    print(f'{pw:14s} int/sub中位={si2.lfc_cd.median() if len(si2) else float("nan"):+.2f} (n={len(si2)})  '
          f'end中位={ep2.lfc_cd.median() if len(ep2) else float("nan"):+.2f} (n={len(ep2)})')

# 胆汁酸：初级(中间体) vs 次级(终产物) 专门检验
ba = sig[sig.pathway == 'bile_acid']
prim = ba[ba.role == 'intermediate'].lfc_cd
sec = ba[ba.role == 'end_product'].lfc_cd
print(f'\n胆汁酸 CD: primary(intermediate) median={prim.median():+.2f} n={len(prim)} ; '
      f'secondary(end) median={sec.median():+.2f} n={len(sec)} ; MWU p={stats.mannwhitneyu(prim,sec).pvalue:.4g}')

sig.to_csv(f'{A}/fig5_disease_breakdown_stats.tsv', sep='\t', index=False)
print('\nwritten: fig5_disease_breakdown_stats.tsv')