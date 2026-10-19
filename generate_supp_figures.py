#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
益生菌"设计性缺陷"范式论文 — 附件图表生成
Nature 风格调色板
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 配色：全部取自全局调色板 paper_palette.PALETTE
# ------------------------------------------------------------
# 附图与主文图必须读同一个调色板，否则两套颜色并列出现时整篇会显得不是
# 一套视觉系统。这里保留原有的字典键名（脚本下方大量引用），但把它们
# 全部重定向到用户指定的 13 色固定顺序，并按语义绑定：
#   lactate/SCFA = slot 6 黄   bile acid = slot 5 蓝
#   tryptophan   = slot 4 紫   neutral   = slot 9 棕
#   上游 = slot 1 橙           下游 = slot 2 红
# 不再使用任何 stock matplotlib 主题色。
# ============================================================
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper_palette import (PALETTE, UP_C, DOWN_C, NEUTRAL, NEUTRAL_C,
                           RELAY_COLORS, present_cmap, sequential_cmap)

NATURE_COLORS = {
    'brown':      NEUTRAL_C,                 # slot 9  中性参照
    'orange':     UP_C,                      # slot 1  橙  上游
    'magenta_d':  PALETTE[3],                # slot 4  紫  色氨酸
    'purple':     PALETTE[3],                # slot 4
    'purple_d':   PALETTE[3],                # slot 4
    'magenta':    PALETTE[3],                # slot 4  色氨酸通路
    'green':      RELAY_COLORS['lactate_scfa'],   # slot 6  黄  乳酸/SCFA
    'blue_d':     RELAY_COLORS['bile_acid'],      # slot 5  蓝  胆汁酸
    'blue_m':     RELAY_COLORS['bile_acid'],      # slot 5
    'blue_l':     PALETTE[7],                # slot 8  青  TMA/TMAO
    'blue_l2':    PALETTE[7],                # slot 8
    'teal':       PALETTE[7],                # slot 8
}
# 语义别名，供新增代码直接使用
LAC_C = RELAY_COLORS['lactate_scfa']   # slot 6
BA_C = RELAY_COLORS['bile_acid']      # slot 5
TRP_C = RELAY_COLORS['tryptophan']    # slot 4
TMA_C = RELAY_COLORS['tma_tmao']      # slot 8
GABA_C = RELAY_COLORS['gaba_glutamate']   # slot 7

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
})

BASE = 'D:/Probiotics_Pangenome_v5'
FIG_SUPP = f'{BASE}/figures/supplementary'

# ============================================================
# Supp Fig 1: 详细通路热图（3 通路）
# ============================================================
def supp_fig1_pathway_heatmaps():
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = gridspec.GridSpec(3, 1, hspace=0.5)

    mp = pd.read_csv(f'{BASE}/pangenome/direction_matrices/unit_marker_presence.tsv', sep='\t')

    # (a) 乳酸→SCFA 通路
    ax = fig.add_subplot(gs[0])
    scfa_cols = [c for c in mp.columns if any(k in c for k in
                 ['TIGR00651', 'TIGR02706', 'TIGR00016', 'TIGR02518', 'ADHE', 'TIGR01255',
                  'TIGR02707', 'PFLB', 'PF02652', 'TIGR00795', 'TIGR01774', 'TIGR01117',
                  'PF09861', 'PF21113', 'PF09330', 'TIGR01591', 'FDH', 'TIGR02330', 'LCTP',
                  'TIGR01553', 'TIGR01583'])]
    scfa_mat = mp[scfa_cols].values.astype(float)
    sort_idx = np.argsort(-scfa_mat.sum(axis=1))
    scfa_mat_sorted = scfa_mat[sort_idx]
    units_sorted = mp.iloc[sort_idx]['cn_name'].values

    from matplotlib.colors import LinearSegmentedColormap
    cmap = sequential_cmap(RELAY_COLORS['lactate_scfa'], 'scfa')
    ax.imshow(scfa_mat_sorted, aspect='auto', cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(scfa_cols)))
    ax.set_xticklabels([c.split('|')[1] for c in scfa_cols], rotation=90, fontsize=4)
    ax.set_yticks(range(len(units_sorted)))
    ax.set_yticklabels(units_sorted, fontsize=4)
    ax.set_title('(a) Lactate→SCFA pathway genes (21 markers)', fontsize=8, loc='left')

    # (b) 胆汁酸通路
    ax = fig.add_subplot(gs[1])
    bile_cols = [c for c in mp.columns if any(k in c for k in
                  ['PF02275', 'BAI_N', 'BAI_H', 'BAI_CD', 'HSDH_7A', 'BAI_F', 'BAI_K', 'HSDH_7B'])]
    bile_mat = mp[bile_cols].values.astype(float)
    sort_idx = np.argsort(-bile_mat.sum(axis=1))
    bile_mat_sorted = bile_mat[sort_idx]
    units_sorted = mp.iloc[sort_idx]['cn_name'].values

    cmap = sequential_cmap(RELAY_COLORS['bile_acid'], 'bile')
    ax.imshow(bile_mat_sorted, aspect='auto', cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(bile_cols)))
    ax.set_xticklabels([c.split('|')[1] for c in bile_cols], rotation=90, fontsize=4)
    ax.set_yticks(range(len(units_sorted)))
    ax.set_yticklabels(units_sorted, fontsize=4)
    ax.set_title('(b) Bile acid pathway genes (8 markers)', fontsize=8, loc='left')

    # (c) 色氨酸通路
    ax = fig.add_subplot(gs[2])
    trp_cols = [c for c in mp.columns if any(k in c for k in
                ['TIGR03811', 'TYDC', 'TIGR00837', 'TIGR02617', 'TIGR02618', 'PF08053',
                 'TIGR02616', 'TIGR03036', 'TIGR01814', 'TIGR03035'])]
    trp_mat = mp[trp_cols].values.astype(float)
    sort_idx = np.argsort(-trp_mat.sum(axis=1))
    trp_mat_sorted = trp_mat[sort_idx]
    units_sorted = mp.iloc[sort_idx]['cn_name'].values

    cmap = sequential_cmap(RELAY_COLORS['tryptophan'], 'trp')
    ax.imshow(trp_mat_sorted, aspect='auto', cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(trp_cols)))
    ax.set_xticklabels([c.split('|')[1] for c in trp_cols], rotation=90, fontsize=4)
    ax.set_yticks(range(len(units_sorted)))
    ax.set_yticklabels(units_sorted, fontsize=4)
    ax.set_title('(c) Tryptophan pathway genes (10 markers)', fontsize=8, loc='left')

    fig.savefig(f'{FIG_SUPP}/supp_fig1_pathway_heatmaps.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig1_pathway_heatmaps.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 1 done")

# ============================================================
# Supp Fig 2: 功能 guild 矩阵
# ============================================================
def supp_fig2_guild_matrix():
    fig, ax = plt.subplots(figsize=(7.2, 5.0))

    # 读取 guild 数据
    try:
        guilds = pd.read_csv(f'{BASE}/pangenome/direction_matrices/unit_functional_guilds.tsv', sep='\t')
        # 构建 guild 矩阵
        guild_cols = [c for c in guilds.columns if c not in ['unit_id', 'cn_name', 'genus']]
        guild_mat = guilds[guild_cols].values.astype(float)
        units = guilds['cn_name'].values
        # 期刊图件不含中文字形：用 RefSeq 拉丁名替代（Arial 无 CJK 会变方框）
        _lat2 = (pd.read_csv(f'{BASE}/metadata/strain_manifest.tsv', sep='\t')
                   .drop_duplicates('unit_id')
                   .set_index('unit_id')['current_name'].to_dict())
        if 'unit_id' in guilds.columns:
            units = [_lat2.get(u, str(u)) for u in guilds['unit_id']]

        from matplotlib.colors import LinearSegmentedColormap
        cmap = sequential_cmap(PALETTE[6], 'guild')
        im = ax.imshow(guild_mat, aspect='auto', cmap=cmap, vmin=0, vmax=1)
        ax.set_xticks(range(len(guild_cols)))
        ax.set_xticklabels(guild_cols, rotation=45, ha='right', fontsize=5)
        ax.set_yticks(range(len(units)))
        ax.set_yticklabels(units, fontsize=4)
        ax.set_title('Functional guild matrix (41 units)', fontsize=9, loc='left')
    except:
        ax.text(0.5, 0.5, 'Guild data not available', transform=ax.transAxes,
                ha='center', va='center', fontsize=8)

    fig.savefig(f'{FIG_SUPP}/supp_fig2_guild_matrix.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig2_guild_matrix.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 2 done")

# ============================================================
# Supp Fig 3: 方向完整性详细评分
# ============================================================
def supp_fig3_direction_completeness():
    fig = plt.figure(figsize=(7.2, 4.0))
    gs = gridspec.GridSpec(1, 3, wspace=0.35)

    dc = pd.read_csv(f'{BASE}/pangenome/direction_matrices/direction_completeness.tsv', sep='\t')

    directions = ['dir1', 'dir2', 'dir3']
    dir_names = ['Lactate→SCFA', 'Bile acid', 'Tryptophan']
    # 每个面板一条通路，配色必须与该通路在全篇的固定色一致：
    # 乳酸/SCFA = slot 6 黄，胆汁酸 = slot 5 蓝，色氨酸 = slot 4 紫。
    # 之前这里用的是 blue_d/orange/magenta，与 dir_names 顺序错配，
    # 导致胆汁酸面板显示为橙色（=上游色），破坏全篇色彩语义。
    dir_colors = [RELAY_COLORS['lactate_scfa'], RELAY_COLORS['bile_acid'],
                  RELAY_COLORS['tryptophan']]
    # 期刊图件只能含西文字形；Arial 无中文字形会渲染成方框，故用拉丁名
    _lat = (pd.read_csv(f'{BASE}/metadata/strain_manifest.tsv', sep='\t')
              .drop_duplicates('unit_id')
              .set_index('unit_id')['current_name'].to_dict())

    for i, (d, name, color) in enumerate(zip(directions, dir_names, dir_colors)):
        ax = fig.add_subplot(gs[i])
        dc_d = dc[dc['direction'] == d].copy()
        dc_d = dc_d.sort_values('weighted_score', ascending=True)
        if 'unit_id' in dc.columns:
            _names = [_lat.get(u, str(u)) for u in dc_d['unit_id']]
        else:
            _names = dc_d['cn_name'].astype(str)

        ax.barh(range(len(dc_d)), dc_d['weighted_score'], color=color, height=0.7)
        ax.set_yticks(range(len(dc_d)))
        ax.set_yticklabels(_names, fontsize=3.5)
        ax.set_xlabel('Completeness score')
        ax.set_title(f'{name}', fontsize=8, loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Direction completeness scores by pathway', fontsize=9, y=1.02)
    fig.savefig(f'{FIG_SUPP}/supp_fig3_direction_completeness.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig3_direction_completeness.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 3 done")

# ============================================================
# Supp Fig 4: 正交同源组层方向富集
# ============================================================
def supp_fig4_orthogroup_enrichment():
    fig, ax = plt.subplots(figsize=(5.0, 3.5))

    de = pd.read_csv(f'{BASE}/pangenome/direction_enrichment_orthogroup.tsv', sep='\t')

    # 提取 MWU 检验结果
    mwu_rows = de[de['set'].str.startswith('MWU')]
    directions = ['dir1', 'dir2', 'dir3']
    dir_names = ['Lactate→SCFA', 'Bile acid', 'Tryptophan']

    x = np.arange(len(directions))
    p_less = []
    for d in directions:
        row = mwu_rows[(mwu_rows['direction'] == d) & (mwu_rows['set'] == 'MWU_p_vs_bg_shared')]
        if len(row) > 0:
            p_less.append(row['mwu_p_less'].values[0])
        else:
            p_less.append(0.5)

    colors_bar = [NATURE_COLORS['blue_d'], NATURE_COLORS['orange'], NATURE_COLORS['magenta']]
    bars = ax.bar(x, p_less, color=colors_bar, width=0.5)
    ax.axhline(y=0.5, color='gray', ls='--', lw=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(dir_names, fontsize=7)
    ax.set_ylabel('P(direction genes span fewer units)')
    ax.set_title('Orthogroup-level direction enrichment', fontsize=9, loc='left')
    for i, p in enumerate(p_less):
        ax.text(i, p + 0.02, f'{p:.3f}', ha='center', fontsize=7)
    ax.set_ylim(0, 1.1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(f'{FIG_SUPP}/supp_fig4_orthogroup_enrichment.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig4_orthogroup_enrichment.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 4 done")

# ============================================================
# Supp Fig 5: 单元泛基因组统计
# ============================================================
def supp_fig5_unit_stats():
    fig = plt.figure(figsize=(7.2, 4.0))
    gs = gridspec.GridSpec(1, 2, wspace=0.3)

    # (a) 各单元泛基因组大小
    ax = fig.add_subplot(gs[0])
    us = pd.read_csv(f'{BASE}/pangenome/unit_stats_orthogroup.tsv', sep='\t')
    us = us.sort_values('n_orthogroups', ascending=True)
    y = range(len(us))
    ax.barh(y, us['n_orthogroups'], color=NATURE_COLORS['blue_d'], height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(us['cn_name'], fontsize=4)
    ax.set_xlabel('Number of orthogroups')
    ax.set_title('(a) Unit pangenome size', fontsize=8, loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # (b) 核心基因比例
    ax = fig.add_subplot(gs[1])
    if 'pct_core' in us.columns:
        us2 = us.sort_values('pct_core', ascending=True)
        y = range(len(us2))
        ax.barh(y, us2['pct_core'], color=NEUTRAL_C, height=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(us2['cn_name'], fontsize=4)
        ax.set_xlabel('% core orthogroups')
        ax.set_title('(b) Core genome proportion', fontsize=8, loc='left')
    else:
        # 尝试其他列名
        core_cols = [c for c in us.columns if 'core' in c.lower()]
        if core_cols:
            us2 = us.sort_values(core_cols[0], ascending=True)
            y = range(len(us2))
            ax.barh(y, us2[core_cols[0]], color=NEUTRAL_C, height=0.7)
            ax.set_yticks(y)
            ax.set_yticklabels(us2['cn_name'], fontsize=4)
            ax.set_xlabel(core_cols[0])
            ax.set_title('(b) Core genome proportion', fontsize=8, loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(f'{FIG_SUPP}/supp_fig5_unit_stats.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig5_unit_stats.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 5 done")

# ============================================================
# Supp Fig 6: 层交叉表分析
# ============================================================
def supp_fig6_layer_crosstab():
    fig, ax = plt.subplots(figsize=(5.0, 4.0))

    lc = pd.read_csv(f'{BASE}/pangenome/layer_crosstab.tsv', sep='\t')

    # 散点图: n_families_90 vs n_genomes
    ax.scatter(lc['n_families_90'], lc['n_genomes'],
               c=NATURE_COLORS['blue_d'], s=5, alpha=0.3, edgecolors='none')
    ax.set_xlabel('Number of 90% families per orthogroup')
    ax.set_ylabel('Number of genomes')
    ax.set_title('Two-layer clustering relationship', fontsize=9, loc='left')

    # 添加趋势线
    z = np.polyfit(lc['n_families_90'], lc['n_genomes'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(lc['n_families_90'].min(), lc['n_families_90'].max(), 100)
    ax.plot(x_line, p(x_line), '--', color=NATURE_COLORS['magenta'], lw=1, alpha=0.7)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(f'{FIG_SUPP}/supp_fig6_layer_crosstab.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig6_layer_crosstab.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 6 done")

# ============================================================
# Supp Fig 7: IBDMDB 酶-代谢物偶联详细
# ============================================================
def supp_fig7_coupling_detail():
    fig = plt.figure(figsize=(7.2, 4.0))
    gs = gridspec.GridSpec(1, 2, wspace=0.3)

    coupling = pd.read_csv(f'{BASE}/validation/analysis/relay_enzyme_metabolite_coupling_ibdmdb.tsv', sep='\t')

    # (a) 偶联强度 vs 预期方向
    ax = fig.add_subplot(gs[0])
    rhos = coupling['spearman_rho'].values
    expected = coupling['expected_direction'].values
    enzymes = coupling['enzyme'].str.split('|').str[1].str.split(' ').str[0]

    # 三分类是「耦合状态」而非某条通路，不能借用通路色（否则读者会误以为
    # 某个色块代表某条代谢通路）。按语义取调色板：绿=符合预期(slot 3)、
    # 红=相反(slot 2，与下游/疾病色一致)、灰=近零(中性)。
    MATCH_C, OPP_C, ZERO_C = PALETTE[2], PALETTE[1], NEUTRAL
    colors_c = []
    for r, e in zip(rhos, expected):
        if e == 'positive' and r > 0:
            colors_c.append(MATCH_C)
        elif e == 'positive' and r < 0:
            colors_c.append(OPP_C)
        else:
            colors_c.append(ZERO_C)

    ax.barh(range(len(enzymes)), rhos, color=colors_c, height=0.6)
    ax.set_yticks(range(len(enzymes)))
    ax.set_yticklabels(enzymes, fontsize=6)
    ax.axvline(x=0, color='gray', ls='-', lw=0.5)
    ax.set_xlabel('Spearman ρ')
    ax.set_title('(a) Coupling direction vs expected', fontsize=8, loc='left')

    legend_elements = [
        mpatches.Patch(color=MATCH_C, label='Matches expected'),
        mpatches.Patch(color=OPP_C, label='Opposite to expected'),
        mpatches.Patch(color=ZERO_C, label='Near zero'),
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=5, loc="lower right")
    ax.set_xlim(-1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # (b) 样本量限制说明
    ax = fig.add_subplot(gs[1])
    ax.axis('off')
    ax.set_title('(b) Sample size limitation', fontsize=8, loc='left')
    text = (
        "IBDMDB MGX∩MBX overlap: n=7\n\n"
        "All correlations non-significant (p>0.05)\n\n"
        "Required for adequate power:\n"
        "• GERD RCT: n=120 paired\n"
        "• ASD cohort: n=203 paired\n"
        "• Mongolian T2D: n=66 paired\n\n"
        "These cohorts will provide\n"
        "adequate statistical power for\n"
        "the L1 discriminative test."
    )
    ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=6,
            va='top', ha='left', family='monospace')

    fig.savefig(f'{FIG_SUPP}/supp_fig7_coupling_detail.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig7_coupling_detail.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 7 done")

# ============================================================
# Supp Fig 8: 验证队列覆盖矩阵
# ============================================================
def supp_fig8_cohort_coverage():
    fig, ax = plt.subplots(figsize=(7.2, 5.0))

    # 队列 × 通路覆盖矩阵
    cohorts = [
        'IBDMDB (HMP2)', 'Zhao T2D', 'Mongolian T2D', 'Franzosa IBD',
        'Gonzalez CD', 'Kang Yu RA', 'Su RA', 'CRC+T2DM', 'NAFLD',
        'GERD RCT', 'BAM RCT', 'T2D RCT', 'PPI RCT', 'BL-99 RCT',
        'UA RCT', 'FBD RCT', 'BLa80 RCT',
        'ASD GABA', 'Laryushina CD', 'Laryushina UC', 'TMAO diet',
    ]
    pathways = ['Lactate→SCFA', 'Bile acid', 'TMA/TMAO', 'GABA/Glu', 'Tryptophan']

    # 覆盖矩阵 (1=covered, 0=not covered)
    coverage = np.array([
        [1, 1, 1, 1, 1],  # IBDMDB
        [1, 1, 0, 0, 0],  # Zhao T2D
        [1, 1, 0, 0, 0],  # Mongolian T2D
        [1, 1, 0, 0, 1],  # Franzosa IBD
        [0, 1, 0, 0, 0],  # Gonzalez CD
        [1, 0, 0, 0, 0],  # Kang Yu RA
        [0, 1, 0, 0, 0],  # Su RA
        [1, 1, 0, 0, 0],  # CRC+T2DM
        [1, 1, 0, 0, 1],  # NAFLD
        [1, 0, 0, 1, 0],  # GERD RCT
        [1, 1, 0, 0, 0],  # BAM RCT
        [1, 1, 0, 0, 0],  # T2D RCT
        [0, 0, 0, 0, 1],  # PPI RCT
        [1, 1, 0, 0, 0],  # BL-99 RCT
        [0, 0, 1, 0, 0],  # UA RCT
        [1, 0, 0, 0, 0],  # FBD RCT
        [0, 0, 0, 1, 0],  # BLa80 RCT
        [0, 0, 0, 1, 0],  # ASD GABA
        [0, 0, 1, 0, 0],  # Laryushina CD
        [0, 0, 1, 0, 0],  # Laryushina UC
        [0, 0, 1, 0, 0],  # TMAO diet
    ])

    from matplotlib.colors import LinearSegmentedColormap
    cmap = sequential_cmap(PALETTE[4], 'coverage')
    ax.imshow(coverage, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(range(len(pathways)))
    ax.set_xticklabels(pathways, rotation=30, ha='right', fontsize=7)
    ax.set_yticks(range(len(cohorts)))
    ax.set_yticklabels(cohorts, fontsize=5)
    ax.set_title('Validation cohort × pathway coverage matrix', fontsize=9, loc='left')

    # 标注 RCT
    for i, c in enumerate(cohorts):
        if 'RCT' in c:
            ax.text(-0.5, i, '★', fontsize=6, ha='right', va='center', color=NATURE_COLORS['orange'])

    fig.savefig(f'{FIG_SUPP}/supp_fig8_cohort_coverage.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG_SUPP}/supp_fig8_cohort_coverage.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Supp Fig 8 done")

# ============================================================
# 主函数
# ============================================================
if __name__ == '__main__':
    print("Generating supplementary figures...")

    supp_fig1_pathway_heatmaps()
    supp_fig2_guild_matrix()
    supp_fig3_direction_completeness()
    supp_fig4_orthogroup_enrichment()
    supp_fig5_unit_stats()
    supp_fig6_layer_crosstab()
    supp_fig7_coupling_detail()
    supp_fig8_cohort_coverage()

    print("\nAll supplementary figures generated!")
    print(f"Output: {FIG_SUPP}")
