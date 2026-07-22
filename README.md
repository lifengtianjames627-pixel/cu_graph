# Cu(II) second-sphere VOC sensing — ORCA 6.1.1 workflow

## 科学目标（一句话）

比较三种 Cu(II) 盐的**二级界氢键环境**（BF₄⁻ < Cl⁻ < TFA⁻），看哪一个能让 **丙酮** 比水 / 甲醇 / 乙酸乙酯 / 甲苯 **结合得更牢**，并在光谱上给出可测的 d–d / LMCT 位移——即“丙酮选择性识别”。

热力学模型是**加合**（guest 挤进二级界口袋，抗衡离子留下）：

```
dE_bind(G)      = E(host·G) − E(host) − E(G)                         [kcal/mol]
dE_bind_BBCP(G) = E(host·G) − E(host@guest-ghosts) − E(guest@host-ghosts)
```

更负 = 结合更强。选择性看 `DD_vs_acetone = dE_bind(G) − dE_bind(acetone)`：正值表示丙酮比该竞争者更强。

---

## 方法 / 基组 / 配套物理解释（全局）

引擎：**ORCA 6.1.1**（关键词按官方手册，见 [Counterpoise](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/counterpoise.html)）。

| 设置 | 值 | 为什么 |
|------|-----|--------|
| 泛函 | **PBE0** | 杂化泛函（25% HF）。纯 GGA 对开壳层 Cu(II) 的 d–d 分裂常塌缩；PBE0 对过渡金属配合物几何与激发能更稳。 |
| 色散 | **D3BJ** | 二级界识别靠氢键 + 色散；不用 D3 会系统性低估弱结合（尤其甲苯）。 |
| 波函数 | **UKS** | Cu(II) 是 d⁹，必须开壳层。闭壳层 RKS 会强制成对、完全错。 |
| 多重度 | **MULT=2** | 双重态 (S=1/2)。**绝不要设成 3**——那是三重态，不是 Cu(II) 基态。 |
| Cu 基组 | **def2-TZVP** | Cu 的 d 轨道需要三 ζ；太小会扭曲 JT 几何与金属–配体键长。 |
| 轻原子基组 | **def2-SVP** | H/C/O/F/Cl/B 的筛选级；便宜，opt 能跑完。定稿可升到 def2-TZVP。 |
| gCP | **GCP(DFT/SVP)** | ORCA 6.1.1 官方几何 CP；小基组 BSSE 大，gCP 自动加进 `FINAL SINGLE POINT ENERGY`（手册 Table 2.50）。 |
| Boys–Bernardi CP | ghost atoms (`:`) | ORCA 6 **不再用** `%cp nfrag`；官方写法是 ghost / `GhostFrags`。脚本对每个 complex 再跑 A/B 两个 ghost 单点。 |
| Opt+Freq | `TightOpt Freq` + `Calc_Hess true` | 确认无虚频（真极小）；起始 Hessian 帮助逃出高对称 JT 鞍点。 |
| TD-DFT | `nroots 60` | Cu(II) 的 d–d 弱、LMCT 更高；根数少会漏掉传感相关的激发。`dosoc false`：手册写明 UKS 参考下 SOC 不可用。 |

---

## 文件布局与电荷 / 多重度 `* xyz q m`

ORCA 坐标块写法：`* xyz <charge> <multiplicity>`。下面每个体系都按化学电荷守恒设定。

### Hosts（二级界阴离子不同）

| 文件 | 组成 | 原子数 | `* xyz` | 电荷来源 |
|------|------|--------|---------|----------|
| `host_tfa.xyz` | [Cu(H₂O)₆](CF₃COO)₂ | 33 | **`0 2`** | Cu²⁺ + 2×TFA⁻ = 0；双重态 |
| `host_bf4.xyz` | [Cu(H₂O)₆](BF₄)₂ | 29 | **`0 2`** | Cu²⁺ + 2×BF₄⁻ = 0 |
| `host_cl.xyz`  | Cu(H₂O)₂Cl₂ | 9 | **`0 2`** | Cu²⁺ + 2×Cl⁻ = 0（Cl 内界配位） |

要点：
- **水全是内界**（配位在 Cu 上），没有“外层水”。
- TFA 必须是 **CF₃COO⁻（去质子）**，不能是 CF₃COOH，否则电荷不是 0。
- 原子顺序建议：Cu → 内界水 → 两个阴离子（CP/ghost 按“host 在前”切分）。

### Guests（共享，三个 host 只算一遍）

| 文件 | 分子 | 原子数 | `* xyz` | 为什么 |
|------|------|--------|---------|--------|
| `guests/acetone.xyz` | CH₃COCH₃ | 10 | **`0 1`** | 中性闭壳层；**目标分析物**（纯 H 键受体） |
| `guests/water.xyz` | H₂O | 3 | **`0 1`** | 最强竞争者 / 湿度干扰 |
| `guests/methanol.xyz` | CH₃OH | 6 | **`0 1`** | 极性邻居（给体+受体） |
| `guests/ethyl_acetate.xyz` | CH₃COOCH₂CH₃ | 14 | **`0 1`** | 中等尺寸酯受体 |
| `guests/toluene.xyz` | C₆H₅CH₃ | 15 | **`0 1`** | **阴性对照**（无氢键，只靠色散） |

Guest 用多重度 **1**（单重态），不是 2——它们没有未成对电子。

### Complexes（host·guest 加合物）

| 文件 | 内容 | `* xyz` |
|------|------|---------|
| `complex_<h>_<g>.xyz` | **host 原子在前**，guest 原子在后（二级界口袋里摆好朝向） | **`0 2`** |

电荷：`CHARGE_HOST + CHARGE_GUEST = 0+0 = 0`。  
多重度：**2**（Cu 仍在，仍是双重态）。

---

## 每一步算什么、为什么需要

脚本默认 `HOSTS=(tfa)`。改成 `(tfa bf4 cl)` 后每个 host 各跑一遍下面流水线。

### Step 0 — `prep`：生成 `.inp`

```bash
bash run_workflow.sh prep
```

不跑 ORCA，只把 xyz → 输入文件。检查生成的 `*.inp` 里 `* xyz q m` 是否正确。

---

### Step 1 — Host 几何优化 + 频率  
**对象：** `host_<h>`（裸二级界 Cu 盐）  
**输入：** `01_host_<h>_optfreq.inp` → `* xyz 0 2`  
**方法：** PBE0-D3BJ / UKS / TightOpt+Freq / Cu:def2-TZVP, light:def2-SVP / GCP(DFT/SVP)

**算什么：**
- 找到 host 的平衡结构（含 JT 畸变：轴向 Cu–O 拉长）。
- 频率确认 **无虚频**（真极小，不是鞍点）。
- 给出 `E(host)` —— 结合能公式的基准。

**为什么需要：**  
没有可靠的 host 极小，后面所有 `dE_bind` 都是噪声。JT 几何还会调制二级界口袋形状，影响 guest 怎么坐进去。

---

### Step 2 — Guest 几何优化 + 频率  
**对象：** 五个孤立 VOC（共享）  
**输入：** `02_guest_<g>_optfreq.inp` → `* xyz 0 1`  
**方法：** 同上（闭壳层，ORCA 对 mult=1 走 RKS/UKS 均可）

**算什么：**
- 每个 guest 自己的气相极小与 `E(G)`。
- 无虚频检查。

**为什么需要：**  
结合能是相对**自由 guest** 定义的。若用未优化的骨架坐标，丙酮 vs 甲醇的“强弱”会被几何误差淹没。五个 guest 只算一次，三个 host 共用。

---

### Step 3 — Complex（host·guest）优化 + 频率  
**对象：** guest 已挤进二级界的加合物  
**输入：** `03_complex_<h>_<g>_optfreq.inp` → `* xyz 0 2`  
**方法：** 同 host（开壳层双重态）

**算什么：**
- 结合构象的局部极小与 `E(host·G)`。
- 无虚频；可看 guest 是否仍停在二级界（而不是挤进内界取代水）。

**为什么需要：**  
这是识别真正发生的结构。`dE_bind` 的分子项。每个 guest 建议先手动摆 2–3 个朝向，取能量最低的一个再进正式表。

---

### Step 4 — Boys–Bernardi 基组重叠误差（BSSE）校正  
**对象：** 固定在 complex 几何上的两个 ghost 单点  
**输入（ORCA 6.1.1 官方 ghost 写法）：**
- `04_bbcp_<h>_<g>_A.inp`：host 真原子 + guest **ghost**（元素后加 `:`）→ `* xyz 0 2`
- `04_bbcp_<h>_<g>_B.inp`：host **ghost** + guest 真原子 → `* xyz 0 1`

**算什么：**
```
dE_bind_BBCP = E(complex) − E(A) − E(B)
```
即 dimer 基组下的单体能量，去掉“guest 基函数帮 host 降能（反之亦然）”的假结合。

**为什么需要：**  
你用的是 def2-SVP 级轻原子基组，BSSE 可达数 kcal/mol，足以翻转丙酮 vs 甲醇的排序。gCP（`USE_GCP=1`）已经便宜地修了一部分；ghost-CP（`USE_GHOST_CP=1`）是手册推荐的正式 Boys–Bernardi，用来给最终选择性表定心。

> 旧版 `%cp nfrag` **不是** ORCA 6.1.1 文档写法，已删除。

---

### Step 5 — TD-DFT（host 与 complex）  
**对象：**  
- `05_host_<h>_tddft.inp` → `* xyz 0 2`  
- `05_complex_<h>_<g>_tddft.inp` → `* xyz 0 2`  

**方法：** PBE0-D3BJ UKS + `%tddft nroots 60`（基组仍走 `%basis` 混合，不在 `!` 行堆两个基组）

**算什么：**
- 最低一批激发根（d–d 弱、LMCT 更强）。
- **相对位移**：complex 相对 host 的根能量变化 = 传感信号假说。

**为什么需要：**  
实验上 Cu(II) 传感常靠颜色 / 吸收变化。只有结合能没有光谱，无法对接“可检测”。注意：TD-DFT 的 d–d 绝对值不准，**相对位移**才有用；同时看基态前线轨道分裂作交叉检查。

---

### Step 6 — `analyze`：汇总选择性

```bash
bash run_workflow.sh analyze
```

打印：
- 各体系能量、`<S**2>`（自旋污染；理想 ~0.75）、虚频？
- 每个 guest 的 `dE_bind`、`dE_BBCP`、`DD_vs_acetone`
- TD-DFT 最低几根

**判据：**  
某 host 下若丙酮的 `dE_bind` 最负，且对水/甲醇的 `DD_vs_acetone` 明显为正 → 该阴离子环境支持丙酮选择性。再比较 tfa / bf4 / cl 谁最干净。

---

## 怎么跑

```bash
# 1. 填好 host_tfa.xyz 与 complex_tfa_<guest>.xyz（原子顺序：host 前、guest 后）
# 2. 确认顶部 HOSTS=(tfa)、ORCA 路径、NPROCS

tmux new -s cu
bash run_workflow.sh run 2>&1 | tee logs/run_$(date +%Y%m%d_%H%M).log
# Ctrl+B 再 D 脱离

# 回来
tmux attach -t cu
bash run_workflow.sh analyze
```

- 已完成的 job 会 **自动跳过**（断点续跑）。
- 强制重跑：`FORCE=1 bash run_workflow.sh run`
- 关掉贵的 ghost-CP（只留 gCP）：脚本顶部 `USE_GHOST_CP=0`
- 关掉 gCP：`USE_GCP=0`

---

## 建几何时的硬规则（避免白算）

1. **TFA⁻ 去质子**；整体 `charge=0`。
2. **不要外层水**；二级界是阴离子 + 内界水向外的 H。
3. Host 起点破对称：轴向 Cu–O ≈ 2.3 Å，赤道 ≈ 2.0 Å（避开 JT 鞍点）。
4. TFA⁻ / BF₄⁻ 放在内界水 H 的氢键距离内（O···O ≈ 2.7–2.9 Å），不要飘在真空里。
5. Complex = host 坐标 + guest 坐进二级界口袋；每个 guest 试 2–3 朝向取最低。
6. 甲苯是阴性对照：若它也出现很负的 `dE_bind`，说明口袋太宽或色散主导，选择性叙事要改。

---

## 配置开关（脚本顶部）

| 变量 | 含义 |
|------|------|
| `HOSTS=(tfa)` | 当前跑哪些 host；以后改 `(tfa bf4 cl)` |
| `TARGET=acetone` | 选择性对比的目标 guest |
| `USE_GCP=1` | 能量是否含 gCP |
| `USE_GHOST_CP=1` | 是否跑 Boys–Bernardi ghost 对 |
| `NPROCS` / `ORCA` | 并行核数 / ORCA 可执行文件路径 |

---

## Cu(II) 备注

- `<S**2>` 接近 0.75；远大于 0.8 要警惕自旋污染。
- 虚频 → 不是极小：改起点或加 `Calc_Hess`（已开）。
- 传感结论看 **2D 图**：`dE_bind`（或 BBCP）vs d–d/LMCT 位移，不要只看能量。
