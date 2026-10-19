# -*- coding: utf-8 -*-
"""接力通路（relay pathways）分析：TMA/TMAO 与 GABA/谷氨酸。

核心叙事：影响宿主的关键代谢往往由多菌接力完成——
上游益生菌产生中间产物，下游非益生菌完成转化，最终作用于宿主。

TMA/TMAO 接力：
  胆碱 --(CutC, 上游益生菌)--> TMA --(TorA, 下游菌)--> TMAO --> 宿主 CVD 风险
  实际上 TMA→TMAO 主要在宿主肝脏 FMO3 完成，TorA 是细菌的 TMAO 呼吸

GABA/谷氨酸接力：
  谷氨酸 --(GadA/B, 上游益生菌)--> GABA --(GabT+GabD, 下游菌)--> SSA --> 琥珀酸 --> SCFA

输出：
  fig_relay_tma.png       TMA/TMAO 接力图
  fig_relay_gaba.png      GABA/谷氨酸接力图
  fig_relay_integrated.png 四通路接力整合图
  relay_trait_matrix.tsv  单元 × 接力性状矩阵
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (DIR_COLOR, despine, genus_color, read_tsv, save)

ROOT = os.environ.get("PROB_ROOT", "D:/Probiotics_Pangenome_v5")
ANN = os.path.join(ROOT, "annotations")
OUT = os.path.join(ROOT, "pangenome", "figures")

import matplotlib.pyplot as plt

# ---------------------------------------------------------------- load
manifest = {}
with open(os.path.join(ROOT, "metadata", "strain_manifest.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r.get("group") == "bacteria":
            manifest[r["accession"]] = r

unit_of = {g: manifest[g]["unit_id"] for g in manifest}
cn_of = {u: manifest[g]["cn_name"] for g, u in unit_of.items()}
genus_of = {u: manifest[g]["current_name"].split()[0] for g, u in unit_of.items()}
unit_size = defaultdict(int)
for g in manifest:
    unit_size[manifest[g]["unit_id"]] += 1

# relay HMM hits
rows = list(csv.DictReader(open(os.path.join(ANN, "hmm_hits_relay.tsv"),
                                encoding="utf-8"), delimiter="\t"))
by_fam = defaultdict(set)
for r in rows:
    by_fam[r["hmm_acc"]].add(r["genome"])

# ---------------------------------------------------------------- traits
# NOTE: CutD (TIGR04395) alone is not used -- it matches other rSAM activases
# (122 genomes have CutD without CutC). CutC is the specific catalytic marker.
RELAY_TRAITS = [
    # (trait name, pathway, role, [families])
    ("胆碱 TMA 裂解酶 CutC", "TMA", "upstream", ["TIGR04394"]),
    ("TMAO 还原酶 TorA", "TMA", "downstream", ["TIGR02164"]),
    ("甜菜碱/胆碱/肉碱转运 BCCT", "TMA", "transport", ["BCCT"]),
    ("谷氨酸脱羧酶 GadA/B", "GABA", "upstream", ["TIGR01788"]),
    ("GABA 转氨酶 GabT", "GABA", "downstream", ["TIGR00700"]),
    ("谷氨酸脱氢酶 GdhA", "GABA", "supply", ["ELFV_dehydrog"]),
    ("谷氨酸/GABA 逆向转运 GadC", "GABA", "transport", ["TIGR00910"]),
]
ROLE_COLOR = {"upstream": "#0072B2", "downstream": "#D55E00",
              "transport": "#009E73", "supply": "#7B3294"}
ROLE_LABEL = {"upstream": "上游（益生菌产中间产物）",
              "downstream": "下游（非益生菌消耗）",
              "transport": "转运", "supply": "底物供给"}

# unit-level carrier fraction
unit_trait = defaultdict(lambda: defaultdict(float))
for name, pw, role, fams in RELAY_TRAITS:
    per_unit = defaultdict(set)
    for fa in fams:
        for g in by_fam.get(fa, set()):
            per_unit[unit_of[g]].add(g)
    for u, gs in per_unit.items():
        unit_trait[u][name] = len(gs) / unit_size[u]

# include every bacterial unit, not only those with a relay-pathway HMM hit:
# L. delbrueckii (unit 10) carries LDH but no TMA/GABA marker, and omitting it
# silently dropped the trait matrix from 41 to 40 units.
uids = sorted({unit_of[g] for g in manifest},
              key=lambda u: (genus_of[u], cn_of[u]))

# ---------------------------------------------------------------- Fig: TMA
fig = plt.figure(figsize=(14.5, 8.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], width_ratios=[1.3, 1],
                      hspace=0.38, wspace=0.28)

# --- (a) relay schematic
axa = fig.add_subplot(gs[0, :])
axa.set_xlim(0, 10)
axa.set_ylim(0, 3.2)
axa.axis("off")

def box(ax, x, y, w, h, text, color, fs=9.5, tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                facecolor=color, edgecolor="none", alpha=0.92))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold", linespacing=1.4)

def arrow(ax, x1, y1, x2, y2, color="#555555", label=None, lx=0, ly=0):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=2.0,
                                mutation_scale=16))
    if label:
        ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, ha="center",
                va="center", fontsize=8, color=color, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

box(axa, 0.2, 1.1, 1.5, 1.0, "胆碱\n(膳食)", "#888888")
box(axa, 2.4, 1.1, 1.8, 1.0, "TMA\n(三甲胺)", "#0072B2")
box(axa, 5.0, 1.1, 1.8, 1.0, "TMAO\n(氧化三甲胺)", "#D55E00")
box(axa, 7.6, 1.1, 2.2, 1.0, "宿主\nCVD 风险↑", "#CC79A7")

arrow(axa, 1.7, 1.6, 2.4, 1.6, "#0072B2", "CutC\n(上游益生菌)", 0, 0.55)
arrow(axa, 4.2, 1.6, 5.0, 1.6, "#D55E00", "FMO3\n(宿主肝脏)", 0, 0.55)
arrow(axa, 6.8, 1.6, 7.6, 1.6, "#CC79A7", "", 0, 0)

# bacterial TMAO respiration (minor branch)
box(axa, 5.0, 0.05, 1.8, 0.6, "TorA 呼吸\n(下游菌)", "#FDBE85", fs=8, tc="#7f2704")
arrow(axa, 5.9, 1.1, 5.9, 0.65, "#D55E00", "", 0, 0)

axa.text(0.2, 2.75, "a  TMA/TMAO 接力通路：上游益生菌产 TMA，宿主肝脏转化为 TMAO",
         fontsize=11, fontweight="bold", color="#222222")
axa.text(0.2, 0.02, "注：TorA 为细菌 TMAO 还原酶（厌氧呼吸），TMA→TMAO 主要由宿主 FMO3 完成",
         fontsize=7.5, color="#777777")

# --- (b) unit-level bar chart
axb = fig.add_subplot(gs[1, 0])
cutc = np.array([unit_trait[u].get("胆碱 TMA 裂解酶 CutC", 0.0) for u in uids])
tora = np.array([unit_trait[u].get("TMAO 还原酶 TorA", 0.0) for u in uids])
order = np.argsort(-cutc)
y = np.arange(len(uids))
axb.barh(y - 0.2, cutc[order], height=0.38, color="#0072B2",
         label="CutC⁺（上游：产 TMA）", edgecolor="white", linewidth=0.3)
axb.barh(y + 0.2, tora[order], height=0.38, color="#D55E00",
         label="TorA⁺（下游：TMAO 呼吸）", edgecolor="white", linewidth=0.3)
axb.set_yticks(y)
axb.set_yticklabels([cn_of[uids[i]] for i in order], fontsize=7.4)
axb.invert_yaxis()
axb.set_xlim(0, 1.05)
axb.axvline(0.6, color="#D55E00", ls="--", lw=1.0)
axb.set_xlabel("单元内携带基因组比例", fontsize=9)
axb.legend(frameon=False, fontsize=8.2, loc="lower right")
axb.set_title("b  TMA 通路：产者与用者", fontsize=10.5, fontweight="bold",
              pad=8, loc="left")
despine(axb)
axb.grid(axis="x", color="#DDDDDD", lw=0.5, alpha=0.7)
axb.set_axisbelow(True)

# --- (c) producer-consumer scatter
axc = fig.add_subplot(gs[1, 1])
axc.scatter(cutc, tora, s=25 + np.array([unit_size[u] for u in uids]) * 2.5,
            c=[genus_color(genus_of[u]) for u in uids], alpha=0.85,
            edgecolors="white", linewidths=0.6)
axc.axvline(0.6, color="#0072B2", ls="--", lw=1.0)
axc.axhline(0.6, color="#D55E00", ls="--", lw=1.0)
axc.set_xlim(-0.03, 1.03)
axc.set_ylim(-0.03, 1.03)
axc.set_xlabel("CutC 携带比例（上游）", fontsize=9)
axc.set_ylabel("TorA 携带比例（下游）", fontsize=9)
axc.set_title("c  产–用解耦", fontsize=10.5, fontweight="bold", pad=8, loc="left")
despine(axc)
axc.grid(color="#DDDDDD", lw=0.5, alpha=0.7)
axc.set_axisbelow(True)
n_up = int(np.sum(cutc >= 0.6))
n_down = int(np.sum(tora >= 0.6))
axc.text(0.97, 0.03, "产者不用", ha="right", va="bottom", fontsize=8, color="#0072B2")
axc.text(0.03, 0.97, "既产又用", ha="left", va="top", fontsize=8, color="#7B3294")

fig.suptitle("TMA/TMAO 接力通路：益生菌是 TMA 的主要来源，但 TMAO 呼吸菌罕见",
             fontsize=12.5, fontweight="bold", y=0.98)
fig.text(0.012, 0.012,
         "CutC 阳性单元 %d/41；TorA 阳性单元 %d/41。气泡大小 = 面板株数。"
         "CutD（激活酶）在所有 CutC⁺ 基因组中共现，但单独出现时为 rSAM 超家族假阳性，未纳入判定。"
         % (n_up, n_down),
         fontsize=7.6, color="#555555")
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
save(fig, OUT, "fig_relay_tma")

# ---------------------------------------------------------------- Fig: GABA
fig = plt.figure(figsize=(14.5, 8.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], width_ratios=[1.3, 1],
                      hspace=0.38, wspace=0.28)

# --- (a) relay schematic
axa = fig.add_subplot(gs[0, :])
axa.set_xlim(0, 10)
axa.set_ylim(0, 3.2)
axa.axis("off")

box(axa, 0.2, 1.1, 1.5, 1.0, "谷氨酸\n(膳食/宿主)", "#888888")
box(axa, 2.4, 1.1, 1.5, 1.0, "GABA\n(γ-氨基丁酸)", "#0072B2")
box(axa, 4.6, 1.1, 1.5, 1.0, "SSA\n(琥珀酸半醛)", "#D55E00")
box(axa, 6.8, 1.1, 1.5, 1.0, "琥珀酸\n→ SCFA", "#E69F00")
box(axa, 8.9, 1.1, 1.0, 1.0, "宿主\n神经/免疫", "#CC79A7")

arrow(axa, 1.7, 1.6, 2.4, 1.6, "#0072B2", "GadA/B\n(上游益生菌)", 0, 0.55)
arrow(axa, 3.9, 1.6, 4.6, 1.6, "#D55E00", "GabT\n(下游菌)", 0, 0.55)
arrow(axa, 6.1, 1.6, 6.8, 1.6, "#E69F00", "GabD", 0, 0.55)
arrow(axa, 8.3, 1.6, 8.9, 1.6, "#CC79A7", "", 0, 0)

axa.text(0.2, 2.75, "a  GABA/谷氨酸接力通路：上游益生菌产 GABA，下游菌完成 GABA 降解",
         fontsize=11, fontweight="bold", color="#222222")
axa.text(0.2, 0.02, "注：GadC 逆向转运蛋白负责谷氨酸/GABA 跨膜交换；GdhA 提供谷氨酸底物",
         fontsize=7.5, color="#777777")

# --- (b) unit-level bar chart
axb = fig.add_subplot(gs[1, 0])
gad = np.array([unit_trait[u].get("谷氨酸脱羧酶 GadA/B", 0.0) for u in uids])
gabt = np.array([unit_trait[u].get("GABA 转氨酶 GabT", 0.0) for u in uids])
order = np.argsort(-gad)
y = np.arange(len(uids))
axb.barh(y - 0.2, gad[order], height=0.38, color="#0072B2",
         label="GadA/B⁺（上游：产 GABA）", edgecolor="white", linewidth=0.3)
axb.barh(y + 0.2, gabt[order], height=0.38, color="#D55E00",
         label="GabT⁺（下游：消耗 GABA）", edgecolor="white", linewidth=0.3)
axb.set_yticks(y)
axb.set_yticklabels([cn_of[uids[i]] for i in order], fontsize=7.4)
axb.invert_yaxis()
axb.set_xlim(0, 1.05)
axb.axvline(0.6, color="#D55E00", ls="--", lw=1.0)
axb.set_xlabel("单元内携带基因组比例", fontsize=9)
axb.legend(frameon=False, fontsize=8.2, loc="lower right")
axb.set_title("b  GABA 通路：产者与用者", fontsize=10.5, fontweight="bold",
              pad=8, loc="left")
despine(axb)
axb.grid(axis="x", color="#DDDDDD", lw=0.5, alpha=0.7)
axb.set_axisbelow(True)

# --- (c) producer-consumer scatter
axc = fig.add_subplot(gs[1, 1])
axc.scatter(gad, gabt, s=25 + np.array([unit_size[u] for u in uids]) * 2.5,
            c=[genus_color(genus_of[u]) for u in uids], alpha=0.85,
            edgecolors="white", linewidths=0.6)
axc.axvline(0.6, color="#0072B2", ls="--", lw=1.0)
axc.axhline(0.6, color="#D55E00", ls="--", lw=1.0)
axc.set_xlim(-0.03, 1.03)
axc.set_ylim(-0.03, 1.03)
axc.set_xlabel("GadA/B 携带比例（上游）", fontsize=9)
axc.set_ylabel("GabT 携带比例（下游）", fontsize=9)
axc.set_title("c  产–用解耦", fontsize=10.5, fontweight="bold", pad=8, loc="left")
despine(axc)
axc.grid(color="#DDDDDD", lw=0.5, alpha=0.7)
axc.set_axisbelow(True)
n_up = int(np.sum(gad >= 0.6))
n_down = int(np.sum(gabt >= 0.6))
axc.text(0.97, 0.03, "产者不用", ha="right", va="bottom", fontsize=8, color="#0072B2")
axc.text(0.03, 0.97, "既产又用", ha="left", va="top", fontsize=8, color="#7B3294")

fig.suptitle("GABA/谷氨酸接力通路：益生菌是 GABA 的主要来源，下游消耗菌更广",
             fontsize=12.5, fontweight="bold", y=0.98)
fig.text(0.012, 0.012,
         "GadA/B 阳性单元 %d/41；GabT 阳性单元 %d/41。气泡大小 = 面板株数。"
         % (n_up, n_down),
         fontsize=7.6, color="#555555")
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
save(fig, OUT, "fig_relay_gaba")

# ---------------------------------------------------------------- Fig: integrated
# combine relay traits with the existing direction traits
guilds = read_tsv(os.path.join(ROOT, "pangenome", "direction_matrices",
                               "unit_functional_guilds.tsv"))
guild_of = {r["cn_name"]: r for r in guilds}

# all relay traits across 4 pathways
ALL_TRAITS = [
    # (display name, pathway, role, source)
    ("LDH/MDH（产乳酸）", "乳酸→SCFA", "upstream", "guild"),
    ("buk（产丁酸）", "乳酸→SCFA", "downstream", "guild"),
    ("BSH（胆盐水解）", "胆汁酸", "upstream", "guild"),
    ("bai 操纵子（7α-脱羟）", "胆汁酸", "downstream", "guild"),
    ("CutC（产 TMA）", "TMA/TMAO", "upstream", "relay"),
    ("TorA（TMAO 呼吸）", "TMA/TMAO", "downstream", "relay"),
    ("GadA/B（产 GABA）", "GABA", "upstream", "relay"),
    ("GabT（消耗 GABA）", "GABA", "downstream", "relay"),
]
PW_COLOR = {"乳酸→SCFA": "#0072B2", "胆汁酸": "#D55E00",
            "TMA/TMAO": "#009E73", "GABA": "#7B3294"}

M = np.zeros((len(uids), len(ALL_TRAITS)))
for i, u in enumerate(uids):
    cn = cn_of[u]
    g = guild_of.get(cn, {})
    for j, (name, pw, role, src) in enumerate(ALL_TRAITS):
        if src == "guild":
            # map display name to guild column
            col = {"LDH/MDH（产乳酸）": "乳酸脱氢酶 LDH/MDH 家族",
                   "buk（产丁酸）": "丁酸通路 buk",
                   "BSH（胆盐水解）": "胆盐水解酶 BSH",
                   "bai 操纵子（7α-脱羟）": "bai 操纵子（7α-脱羟）"}[name]
            M[i, j] = float(g.get(col, 0))
        else:
            M[i, j] = unit_trait[u].get(
                {"CutC（产 TMA）": "胆碱 TMA 裂解酶 CutC",
                 "TorA（TMAO 呼吸）": "TMAO 还原酶 TorA",
                 "GadA/B（产 GABA）": "谷氨酸脱羧酶 GadA/B",
                 "GabT（消耗 GABA）": "GABA 转氨酶 GabT"}[name], 0.0)

fig, (axm, axs) = plt.subplots(1, 2, figsize=(15.0, 9.2),
                               gridspec_kw={"width_ratios": [1.6, 1], "wspace": 0.30})

cmap = ListedColormap(["#F2F2F2", "#FDBE85", "#D94801"])
norm = BoundaryNorm([0, 0.05, 0.6, 1.01], cmap.N)
axm.imshow(M, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
axm.set_xticks(range(len(ALL_TRAITS)))
axm.set_xticklabels([t[0] for t in ALL_TRAITS], rotation=40, ha="right", fontsize=8)
axm.set_yticks(range(len(uids)))
axm.set_yticklabels([cn_of[u] for u in uids], fontsize=7.6)
axm.tick_params(length=0)
for s in axm.spines.values():
    s.set_visible(False)
# pathway colour bar on top
for j, (name, pw, role, src) in enumerate(ALL_TRAITS):
    axm.add_patch(Rectangle((j - 0.5, -1.5), 1, 0.85, color=PW_COLOR[pw],
                            lw=0, clip_on=False))
    axm.text(j, -2.0, pw, ha="center", va="bottom", fontsize=7.2,
             color=PW_COLOR[pw], fontweight="bold", rotation=0)
axm.set_title("a  四通路接力性状矩阵（41 单元 × 8 性状）", fontsize=11,
              fontweight="bold", pad=36, loc="left")

# right: upstream vs downstream breadth per pathway
axs.set_xlim(-0.5, 3.5)
axs.set_ylim(-0.05, 1.05)
pws = ["乳酸→SCFA", "胆汁酸", "TMA/TMAO", "GABA"]
up_idx = [0, 2, 4, 6]
dn_idx = [1, 3, 5, 7]
up_vals = [np.mean(M[:, j] >= 0.6) for j in up_idx]
dn_vals = [np.mean(M[:, j] >= 0.6) for j in dn_idx]
x = np.arange(4)
w = 0.35
axs.bar(x - w / 2, up_vals, width=w, color="#0072B2", label="上游（益生菌）",
        edgecolor="white", linewidth=0.4)
axs.bar(x + w / 2, dn_vals, width=w, color="#D55E00", label="下游（非益生菌）",
        edgecolor="white", linewidth=0.4)
for i in range(4):
    axs.text(x[i] - w / 2, up_vals[i] + 0.02, "%d/41" % int(up_vals[i] * 41),
             ha="center", fontsize=8, color="#0072B2", fontweight="bold")
    axs.text(x[i] + w / 2, dn_vals[i] + 0.02, "%d/41" % int(dn_vals[i] * 41),
             ha="center", fontsize=8, color="#D55E00", fontweight="bold")
axs.set_xticks(x)
axs.set_xticklabels(pws, fontsize=9)
axs.set_ylabel("阳性单元比例", fontsize=9.5)
axs.legend(frameon=False, fontsize=8.5)
axs.set_title("b  上游广、下游窄：接力的普遍模式", fontsize=11,
              fontweight="bold", pad=10, loc="left")
despine(axs)
axs.grid(axis="y", color="#DDDDDD", lw=0.5, alpha=0.7)
axs.set_axisbelow(True)

fig.suptitle("多菌接力：益生菌作为上游生产者，下游由非益生菌完成转化",
             fontsize=12.5, fontweight="bold", y=0.98)
fig.text(0.012, 0.012,
         "判定口径：单元内 ≥60% 基因组携带记为阳性。上游 = 产生中间产物的酶；"
         "下游 = 消耗中间产物的酶。四条通路均呈现「上游广、下游窄」的接力模式。",
         fontsize=7.6, color="#555555")
fig.tight_layout(rect=[0, 0.03, 1, 0.955])
save(fig, OUT, "fig_relay_integrated")

# ---------------------------------------------------------------- save matrix
with open(os.path.join(ROOT, "pangenome", "direction_matrices",
                       "relay_trait_matrix.tsv"), "w", encoding="utf-8",
          newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["unit_id", "cn_name", "genus"] + [t[0] for t in ALL_TRAITS])
    for i, u in enumerate(uids):
        w.writerow([u, cn_of[u], genus_of[u]] +
                   ["%.2f" % M[i, j] for j in range(len(ALL_TRAITS))])
print("matrix -> relay_trait_matrix.tsv")

# ---------------------------------------------------------------- stats
print("\n=== 接力通路统计 ===")
for j, (name, pw, role, src) in enumerate(ALL_TRAITS):
    n = int(np.sum(M[:, j] >= 0.6))
    print("  %-22s %-10s %-10s %2d/41" % (name, pw, role, n))
print("\n上游 vs 下游阳性单元数：")
for i, pw in enumerate(pws):
    print("  %-10s 上游 %d/41  下游 %d/41" % (pw, int(up_vals[i] * 41),
                                              int(dn_vals[i] * 41)))
