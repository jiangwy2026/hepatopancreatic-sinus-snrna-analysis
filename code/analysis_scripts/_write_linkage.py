# -*- coding: utf-8 -*-
import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
doc = r"""
# 前期稿件（Cell Reports 投稿 S-26-05162）与本项目单细胞数据的对接方案

生成时间：2026-10-16 ｜ 稿件：*Shrimp hepatopancreatic sinus is a bone marrow-like hematopoietic niche*
（Wenyong Jiang, Chuncheng Song, Fan Wang；汕头大学；投 Cell Reports，Manuscript Draft 59 页）

## 一、前期稿件已确立的结论

| # | 结论 | 关键证据（图中） |
|---|---|---|
| 1 | 抽血（hemospasia）后 **只有肝胰腺窦（Hs）** 出现显著增殖，其它 9 种组织无一致变化 | Fig 1D-E：BrdU + FACS 定量 10 种组织 |
| 2 | Hs 内**多核合胞体（syncytium）→ 单核细胞**的转变发生在 **6 hph**；核面积显著下降（p<0.0001） | Fig 1F：共聚焦 + Hoechst/BrdU + 核面积统计 |
| 3 | THC 在 3–6 h 下降 >30%，6 h 起回升 | Fig 1B |
| 4 | Hs 血浆（Hsp）与循环血浆（Cp）完全不同：Hsp 富含虾青素（红色）、244 个上调代谢物、230 个蛋白 | Fig 2 |
| 5 | **甜菜碱（betaine）** 被体内实验验证为促进血细胞恢复的关键代谢物（候选还有 L-carnitine、niacinamide、limaprost、p-coumaric acid、captopril、oleic acid） | Fig 3 + Key Resources Table |
| 6 | 蛋白呈两阶段调控：6 h 上调细胞周期启动因子（eIF4b、ARHGAP12、dynein light chain）；随后持续上调蛋白合成机器（eEF1α、RPS3A、RPS9）；黏附相关 filamin-A、beta-parvin 持续表达 | Fig 4 |
| 7 | Hs 与淋巴器官都是**晚期、髓系偏向**的造血区室；与人的骨髓共享 14 / 20 个保守同源基因，全部为髓系基因；检测不到经典 HSC marker | Fig 5 + Table 1 |

**已登记数据**：代谢组 MTBLS1478（MetaboLights）、蛋白组 PXD079762（ProteomeXchange）、转录组 PRJNA1478458（SRA，Hs 细胞 bulk RNA-seq，对照 vs 6 h，n=20/组）。
**他们转录组用的参考基因组**：GCF_042767895.1（2025 版本）——与本项目（GCF_003789085.1，2018 科海1号）**不同**。

### 稿件明确写出的技术缺口（对我们最有利的一段）

> "Despite repeated technical attempts, high-quality single-cell transcriptomic data could not be acquired from isolated Hs tissue, most likely due to the extreme fragility of Hs-resident syncytia during sample collection."（Results 第 2 段）

> "Despite extensive attempts to recapitulate the syncytium-to-single-cell transition in vitro, isolated Hs syncytia rapidly lost viability under culture conditions, precluding the reconstitution of their native proliferative microenvironment."（Discussion 局限段）
> "we observed obvious lymphoid organ enlargement upon hemospasia... However, our proliferation assays failed to produce consistent functional outcomes... These progenitors could migrate toward the hepatopancreatic sinus, where they undergo syncytium formation and subsequent cellularization..."（Discussion）

→ 也就是说：**他们没有单细胞数据；合胞体来源、迁移假说、细胞类型身份全部是推测。** 我们的 9,658 细胞 / 9 个亚群数据正是这块空白。

## 二、已核实的接口（本项目数据 × 前期稿件）

### 接口 1：增殖与时间点完全对应
- 我们的 treated 样本 = **抽血后 6 h**，与他们成像/组学的 6 hph 一致
- 我们：PH（原血细胞）在 6 h 扩增 **2.30 倍**（χ²=652.7, P=1.1e−135）；PH 的 G2M 期细胞占 **52.2%**，其余三群 29–33%
- 我们：早期祖细胞中 **ccnd2 +2.02**、**血细胞转谷氨酰胺酶 +2.10** 上调
→ 结论 1（Hs 是增殖中心）得到**独立方法学（单细胞转录组）的验证**

### 接口 2：合胞体假说 ↔ 我们的 PH 标记基因
我们 PH 群最特异的 marker 中有两个与"细胞融合/基质交联"直接相关：
- **LOC113823945（血细胞转谷氨酰胺酶）**：PH 58.7% vs 其他 9.9%，log2FC 4.50（Δ48.8，PH 第 3 特异）
- **SSPO（hemocytin）**：PH 71.3% vs 29.1%
转谷氨酰胺酶是已知的细胞融合/组织交联酶（胎盘合胞体、破骨细胞融合均依赖它）。
→ **可检验假设**：PH 代表合胞体解离后的单核细胞，其"祖细胞程序"（转谷氨酰胺酶 + hemocytin）就是维持合胞体状态的分子基础。

### 接口 3：跨物种保守基因（Table 1）在单细胞上的定位（已做初步映射）
用基因符号在 26,778 个基因中检索他们的虾基因名：**34 个里对上 17 个**（其余因基因组版本不同未命中，需做 ortholog 映射）。已对上者的分布：

| 他们的虾基因 | 人同源基因（他们的注释） | 我们数据中的主要表达群 |
|---|---|---|
| Syne1 | SYNE1（Late Myeloid） | **PH 63%** > GH 40% > MH2 30%（*他们淋巴器官里 SYNE1 落在 "Progenitor cells" 簇 → 互相印证*） |
| Mmel1（neprilysin 家族） | ECE1（Neutrophil） | **PH 61%** ≈ GH 58% > MH2 52% |
| Nlrp3 | NLRP3（Monocyte） | **MH2 21%**（最高） |
| LGMN | LGMN（Macrophages） | **MH2 12%**（最高） |
| Rgs2 | RGS2（Neutrophil） | **MH2 50%**（最高） |
| Tctp | TPT1（Macrophages） | MH1 48%（血细胞中最高） |
| KCNQ5 | KCNQ5（Early Myeloid Progenitor） | **NEU 21%，血细胞几乎不表达** ← 与他们的推断**冲突**，需修正 |
| Act5C | ACTB（Neutrophil） | 全群高（看家基因，符合预期） |

→ 我们的数据可以把他们"组织水平"的保守性结论升级为**细胞类型水平的髓系身份判定**（MH2 = 巨噬/单核样承接 NLRP3/LGMN/RGS2，PH 承接 Syne1），同时纠正个别基因（KCNQ5）的错误归属。

### 接口 4：多组学（蛋白/代谢）的细胞类型归属
他们 6 h 上调的蛋白（eIF4b、ARHGAP12、dynein light chain、eEF1α、RPS3A、RPS9、filamin-A、beta-parvin）与代谢通路（甜菜碱、虾青素、抗氧化）都能在我们的 9 个亚群里定位：哪个亚群表达这套"翻译/周期/黏附"机器，就说明是它驱动了 6 h 的增殖。

### 接口 5：他们的 bulk 转录组可作为我们结论的独立验证
PRJNA1478458（对照 vs 6 h，各 20 只）→ 用我们的单细胞作参考做**反卷积**，看 6 h 时哪个细胞类型占比上升。若 PH 上升，就是对我们 2.30 倍扩增的跨平台复现。

## 三、建议的后续分析（按性价比排序）

### A. 立即可做（数据都在手上，1–3 天）

1. **合胞体假说的直接检验**（最高价值，直接回答他们的核心机制）
   - a) 回查**未过滤矩阵**（ADRL26080060_raw，19,795 细胞）：被当作双细胞/高 RNA 丢弃的条码里，是否存在"同型 PH 双细胞"（= 多核合胞体）。我们的 QC 与厂商 QC 都可能把合胞体当双细胞删掉——**必须先排除这个可能**
   - b) 单条码内 **S 期与 G2M 期评分同时偏高**的细胞（单核细胞中互斥，多核合胞体里可共存）
   - c) 融合机器基因（转谷氨酰胺酶家族、hemocytin、tetraspanin/CD9/CD81、DC-STAMP/Myomaker 同源）随拟时的表达动态
2. **Table 1 保守基因的单细胞定位图**：先做基因组版本间的 ortholog 映射（用他们 Table 1 的 XP 蛋白号 BLAST 到我们参考的蛋白），再做 35 基因 × 9 亚群的气泡图 → 直接接续他们的 Fig 5
3. **6 h 蛋白/代谢通路的细胞归属**：翻译机器、细胞周期启动、黏附（filamin-A/beta-parvin）、抗氧化（cMnSOD3/ferritin/GS）、甜菜碱通路（SLC6A12/BGT1、CHDH、BHMT、NFAT5）
4. **用我们的单细胞反卷积他们的 bulk RNA-seq**（PRJNA1478458）→ 独立验证 PH 在 6 h 扩增
5. **与淋巴器官 scRNA-seq（Sun et al. 2025，李富花组）的标签转移**：检验"淋巴器官前体迁移到窦内"的假说；他们的 Table 1 淋巴器官 marker 序列已公开可用

### B. 需要新实验（写进下一篇或修回补充）

6. **时间序列单细胞**（0/3/6/12/24 h）：现在只有 6 h 单点，加上时间序列才能把"合胞体解离 → PH → 成熟"的顺序钉死（与他们的 THC 曲线一一对应）
7. **切片原位杂交/免疫荧光**：用我们最特异的 PH marker（血细胞转谷氨酰胺酶 LOC113823945、Grik2、PPAF2）标记肝胰腺窦切片——**如果能在他们的合胞体上看到信号，两篇研究就完全融合**（这是最漂亮的一张图，成本也低）
8. **BrdU 脉冲-追踪 + 单细胞**：把"谁增殖、谁分化"从推断变成追踪
9. **甜菜碱处理 + 单细胞**：他们已证明 betaine 促进恢复；我们能回答"betaine 作用在哪群细胞"（若 PH 扩增 ⇒ 机制闭环）
10. **合胞体的单细胞捕获方法学**：温和解离/不消化直接捕获（或用 10x 的 fixed-RNA 方案），争取直接测到多核合胞体

## 四、必须处理的风险与冲突

| 风险 | 说明 | 处理 |
|---|---|---|
| **参考基因组版本不同** | 他们 GCF_042767895.1（2025），我们 GCF_003789085.1（2018）；按基因名只能对上 17/34 | 用他们 Table 1 的 XP 蛋白号做 BLAST/ortholog 映射（OrthoFinder 或 BLAST+ best hit） |
| **合胞体可能在单细胞数据里"消失"** | 我们的 PH 呈单核外观，可能因为处理中解离、或被双细胞过滤删除 | 必做 A-1 的 a/b 项；若确实被过滤，需重新用未过滤数据建库分析 |
| **KCNQ5 冲突** | 他们的"早期髓系祖细胞" marker 在我们数据里属神经内分泌群 | 复核 ortholog 映射；必要时在文中说明差异 |
| **样本定义差异** | 他们测"Hs 细胞"（消化混合），我们测"窦内流出细胞" | 在方法里写清；反卷积时注意组成差异 |
| **发表策略** | 同期投稿需避免自我重复审稿 | 建议定位为 **companion（单细胞解析）**，明确引用他们的投稿/预印本 |

## 五、可以立刻命名的一张"对接图"（建议 Fig 6 或补充图）

**"单细胞分辨率下的肝血窦造血"** 四联图：
- A：9 群 UMAP（PH 高亮）+ PH 的 G2M 比例柱状图（对应他们 Fig 1E）
- B：PH 特异的融合/祖细胞基因（转谷氨酰胺酶、hemocytin、ccnd2、Syne1）的表达与拟时动态（对应他们 Fig 1F 的合胞体→细胞）
- C：Table 1 保守基因在 9 群的气泡图（对应他们 Fig 5 与 Table 1）
- D：6 h 上调蛋白/代谢通路的细胞归属（对应他们 Fig 4）

## 六、附：本次读取产生的中间文件

- 稿件正文文本：`文章所需\cell report\_extract_MANUSCRIPTS_0702.txt`、`_extract_Manuscript_0702.txt`
- 覆盖信文本：`_extract_Cover letter 0702.txt`、`_extract_Cover letter 0703.txt`
- Table 1 的 34 个虾基因在单细胞数据中的检索结果：`analysis_output\manuscript_table1_genes_in_scrna.csv`
- 检索脚本：`analysis_output\scripts\122_table1_genes.py`
"""
open(r"D:\桌面\single cell\analysis_output\CELL_REPORTS_LINKAGE_PLAN.md", "w", encoding="utf8").write(doc)
print("saved:", r"D:\桌面\single cell\analysis_output\CELL_REPORTS_LINKAGE_PLAN.md", len(doc), "chars")
