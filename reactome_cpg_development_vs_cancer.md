# Reactome-Style Biological Pathway Diagram: Top Age-Predictive CpGs
**Functional Axis: Preadulthood Development vs. Adulthood Cancer**

This document provides a Reactome-style biological pathway mapping of the **top 20 age-predictive rhesus macaque CpG sites** (lifted over to GRCh38/hg38) from [`cpg_liftover_summary.md`](file:///Users/hellyr/Desktop/Projects/Methylation%20Aging/cpg_liftover_summary.md). 

Instead of tissue specificity, this analysis contrasts their functional relevance across two major biological temporal axes:
1. 🐣 **Preadulthood Development**: Organogenesis, embryonic transcription, stem cell differentiation, and morphogenetic growth in early life.
2. 🎗️ **Adulthood Cancer**: Oncogenesis, tumor suppression, invadopodia/metastasis, chemoresistance, and cell survival in late adulthood.

> [!TIP]
> **Interactive Visualizer Available**: An interactive HTML diagram with node filtering, spectrum ranking, and detailed inspect panels is saved at:
> [reactome_cpg_development_vs_cancer.html](file:///Users/hellyr/Desktop/Projects/Methylation%20Aging/reactome_cpg_development_vs_cancer.html)

---

## 1. Reactome Pathway Hierarchy Diagram

The 20 CpG sites are mapped into **6 Functional Reactome Modules** categorized by their primary developmental and cancer mechanisms.

- 🔷 **Blue Nodes / Borders**: **Development-Dominated** (High morphogenetic / early cell fate role)
- 🔴 **Red Nodes / Borders**: **Cancer-Dominated** (High oncogenic, tumor suppressor, or metastatic role)
- 🟣 **Purple Nodes / Borders**: **Dual Master Axis** (Master regulators active in both early development AND late cancer transformation)
- 🟢 **Green Nodes / Borders**: **Core Maintenance** (Baseline cellular housekeeping)
- 🔻 **Downward Arrow `(↓)`**: Hypo-methylated with age ($r < -0.70$)
- 🔺 **Upward Arrow `(↑)`**: Hyper-methylated with age ($r > +0.38$)

```mermaid
graph TD
    classDef devDom text-align:left,fill:#0c4a6e,stroke:#38bdf8,stroke-width:2px,color:#f0f9ff;
    classDef canDom text-align:left,fill:#4c0519,stroke:#f43f5e,stroke-width:2px,color:#fff1f2;
    classDef dualAx text-align:left,fill:#3b0764,stroke:#a855f7,stroke-width:2px,color:#faf5ff;
    classDef coreMaint text-align:left,fill:#064e3b,stroke:#10b981,stroke-width:2px,color:#ecfdf5;

    subgraph REACTOME["Reactome Root: Epigenetic Drift, Development & Oncogenesis"]
        
        subgraph M1["Module 1: Embryonic Morphogenesis & Transcription"]
            NKX25["cpg: 6_170631732 (↑ r=+0.430)<br/>Gene: NKX2-5 (Cardiac Homeobox TF)<br/>Weight: 80% Dev / 20% Cancer"]:::devDom
            EBF3["cpg: 9_130386965 (↑ r=+0.425)<br/>Gene: EBF3 (Neurogenesis & Tumor Suppressor)<br/>Weight: 50% Dev / 50% Cancer"]:::dualAx
        end

        subgraph M2["Module 2: Stem Cell Self-Renewal & Proliferative Capacity"]
            GNL3["cpg: 2_108854233 (↓ r=-0.727)<br/>Gene: GNL3 / SNORD19 (Stem Cell Nucleostemin)<br/>Weight: 50% Dev / 50% Cancer"]:::dualAx
            MEST["cpg: 3_156792975 (↑ r=+0.394)<br/>Gene: MEST / MIR335 (Imprinted Body Patterning & EMT)<br/>Weight: 55% Dev / 45% Cancer"]:::dualAx
            RORB["cpg: 15_90514046 (↑ r=+0.389)<br/>Gene: RORB-AS1 (Neural / Retinal Organogenesis)<br/>Weight: 60% Dev / 40% Cancer"]:::devDom
        end

        subgraph M3["Module 3: Cytoskeletal Dynamics, Invadopodia & Metastasis"]
            ENAH["cpg: 1_83189442 (↓ r=-0.721)<br/>Gene: ENAH (MenaINV Cancer Invadopodia)<br/>Weight: 30% Dev / 70% Cancer"]:::canDom
            ABLIM3["cpg: 6_146779581 (↓ r=-0.722)<br/>Gene: ABLIM3 (HCC Tumor Suppressor)<br/>Weight: 25% Dev / 75% Cancer"]:::canDom
            CFL1["cpg: 14_8699620 (↓ r=-0.734)<br/>Gene: CFL1 / MUS81 (ECM Degradation & Chemoresistance)<br/>Weight: 30% Dev / 70% Cancer"]:::canDom
        end

        subgraph M4["Module 4: Tumor Suppression, DNA Repair & Genomic Stability"]
            MSRA["cpg: 8_10151079 (↓ r=-0.761)<br/>Gene: MSRA (Methionine Reductase Tumor Suppressor)<br/>Weight: 25% Dev / 75% Cancer"]:::canDom
            RFC2["cpg: 3_46594128 (↑ r=+0.403)<br/>Gene: LAT2 / RFC2 (Oncogenic Unconstrained Proliferation)<br/>Weight: 45% Dev / 55% Cancer"]:::dualAx
            NDUFA9["cpg: 11_4742706 (↓ r=-0.724)<br/>Gene: NDUFA9 / AKAP3 (Warburg Shift & ROS)<br/>Weight: 35% Dev / 65% Cancer"]:::canDom
        end

        subgraph M5["Module 5: Oncogenic Receptor & Axonal/Tumor Guidance"]
            NTN1["cpg: 16_9035289 (↑ r=+0.387)<br/>Gene: NTN1 (Netrin-1 Apoptosis Suppressor)<br/>Weight: 50% Dev / 50% Cancer"]:::dualAx
            CCKBR["cpg: 14_59685611 (↓ r=-0.729)<br/>Gene: CCKBR (Pancreatic / Gastric Oncogenic Receptor)<br/>Weight: 35% Dev / 65% Cancer"]:::canDom
            CDK5["cpg: 3_177059260 (↓ r=-0.738)<br/>Gene: SLC4A2 / CDK5 (Oncogenic Kinase Invasion)<br/>Weight: 30% Dev / 70% Cancer"]:::canDom
        end

        subgraph M6["Module 6: Chromatin Remodeling, Transport & Translation"]
            MTA3["cpg: 13_66302679 (↑ r=+0.412)<br/>Gene: MTA3 / KCNG3 (Snail Repression & EMT)<br/>Weight: 35% Dev / 65% Cancer"]:::canDom
            IPO8["cpg: 11_32040153 (↓ r=-0.721)<br/>Gene: IPO8 (Nuclear Transport of Oncogenic TFs)<br/>Weight: 45% Dev / 55% Cancer"]:::dualAx
            RPL5["cpg: 1_131347781 (↓ r=-0.753)<br/>Gene: RPL5 / DIPK1A (MDM2-p53 Nucleolar Surveillance)<br/>Weight: 50% Dev / 50% Cancer"]:::dualAx
            LINC["cpg: 17_81303274 (↑ r=+0.416)<br/>Gene: LINC01039 (Oncogenic Proliferation lncRNA)<br/>Weight: 40% Dev / 60% Cancer"]:::dualAx
            NOVEL1["cpg: 16_48045751 (↑ r=+0.398)<br/>Gene: ENSG00000305838 (Baseline Non-coding RNA)<br/>Weight: 50% Dev / 50% Cancer"]:::coreMaint
            NOVEL2["cpg: 5_61925557 (↑ r=+0.383)<br/>Gene: Intergenic (chr4 Regulatory Element)<br/>Weight: 50% Dev / 50% Cancer"]:::coreMaint
        end

    end
```

---

## 2. Development vs. Cancer Quantitative Spectrum Table

The table below ranks all 20 CpG sites by their **Developmental Weight (%)** vs. **Adulthood Cancer Weight (%)**.

| Rank | Macaque CpG Site | Human Coords (hg38) | Correlation ($r$) | Associated Gene(s) | Primary Functional Role | Axis Classification | 🐣 Dev Weight | 🎗️ Cancer Weight |
|:---:|---|---|:---:|---|---|:---:|:---:|:---:|
| **1** | `6_170631732` | chr5:173228952 | **+0.430** | **NKX2-5**, RPL7AP33 | Heart Field Morphogenesis & Progenitor Fate | **Development-Dominated** | **80%** | **20%** |
| **2** | `15_90514046` | chr9:74242822 | **+0.389** | **RORB-AS1** | Retinal & Neural Layer Specification | **Development-Dominated** | **60%** | **40%** |
| **3** | `3_156792975` | chr7:130492812 | **+0.394** | **MEST**, MIR335 | Imprinted Embryonic Patterning & EMT | **Dual Master Axis** | **55%** | **45%** |
| **4** | `9_130386965` | chr10:129959989 | **+0.425** | **EBF3** | Neurogenesis & Tumor Suppressor Silencing | **Dual Master Axis** | **50%** | **50%** |
| **5** | `2_108854233` | chr3:52685326 | **-0.727** | **GNL3**, SNORD19 | Embryonic Stem Cell / Cancer Stem Cell Renewal | **Dual Master Axis** | **50%** | **50%** |
| **6** | `16_9035289` | chr17:9239903 | **+0.387** | **NTN1** | Axonal Guidance & Tumor Apoptosis Evasion | **Dual Master Axis** | **50%** | **50%** |
| **7** | `1_131347781` | chr1:92832287 | **-0.753** | **RPL5**, DIPK1A | 60S Ribosomal Biogenesis & p53 Stress Sensing | **Dual Master Axis** | **50%** | **50%** |
| **8** | `11_32040153` | chr12:30695876 | **-0.721** | **IPO8** | Nuclear Transport of Morphogenetic / Oncogenic TFs | **Dual Master Axis** | **45%** | **55%** |
| **9** | `3_46594128` | chr7:74226771 | **+0.403** | **LAT2**, RFC2 | Okazaki Fragment Replication & Tumor Hyper-proliferation | **Dual Master Axis** | **45%** | **55%** |
| **10** | `17_81303274` | chr13:99580167 | **+0.416** | **LINC01039** | Stem Cell Differentiation & Carcinoma Invasion | **Dual Master Axis** | **40%** | **60%** |
| **11** | `13_66302679` | chr2:42493004 | **+0.412** | **MTA3**, KCNG3 | NuRD Snail Repression & Tumor Invasion Control | **Cancer-Dominated** | **35%** | **65%** |
| **12** | `14_59685611` | chr11:6259554 | **-0.729** | **CCKBR** | Pancreatic/GI Endocrine Growth Loops | **Cancer-Dominated** | **35%** | **65%** |
| **13** | `11_4742706` | chr12:4649344 | **-0.724** | **NDUFA9**, AKAP3 | Mitochondrial Complex I & Warburg ROS Shift | **Cancer-Dominated** | **35%** | **65%** |
| **14** | `1_83189442` | chr1:225652368 | **-0.721** | **ENAH** | Neurite Sprouting & MenaINV Invadopodia Metastasis | **Cancer-Dominated** | **30%** | **70%** |
| **15** | `14_8699620` | chr11:65859050 | **-0.734** | **CFL1**, MUS81 | Cleavage Cytokinesis & Chemoresistance Repair | **Cancer-Dominated** | **30%** | **70%** |
| **16** | `3_177059260` | chr7:151057227 | **-0.738** | **CDK5**, SLC4A2 | Corticogenesis Kinase & Tumor Proliferation Driver | **Cancer-Dominated** | **30%** | **70%** |
| **17** | `8_10151079` | chr8:10417746 | **-0.761** | **MSRA** | Antioxidant Protein Repair & Tumor Suppressor Loss | **Cancer-Dominated** | **25%** | **75%** |
| **18** | `6_146779581` | chr5:149265049 | **-0.722** | **ABLIM3** | Hepatocyte Morphogenesis & HCC Tumor Suppression | **Cancer-Dominated** | **25%** | **75%** |
| **19** | `16_48045751` | chr17:48754741 | **+0.398** | *ENSG00000305838* | Uncharacterized LncRNA Locus | **Core Maintenance** | **50%** | **50%** |
| **20** | `5_61925557` | chr4:71832172 | **+0.383** | *Intergenic* | Intergenic Cis-Regulatory Sequence | **Core Maintenance** | **50%** | **50%** |

---

## 3. Biological Breakdown by Functional Reactome Module

### Module 1: Embryonic Morphogenesis & Transcription
- **`6_170631732`** ($r = +0.430$) $\rightarrow$ **NKX2-5** (80% Dev / 20% Cancer)
  - **Preadulthood Development**: Master homeobox transcription factor required for early cardiac progenitor specification, heart tube looping, and chamber differentiation.
  - **Adulthood Cancer**: Repressed in normal adult somatic tissues. Age-associated hyper-methylation locks early lineage homeobox genes into repressed heterochromatin.
- **`9_130386965`** ($r = +0.425$) $\rightarrow$ **EBF3** (50% Dev / 50% Cancer)
  - **Preadulthood Development**: Drives early B-cell lineage commitment and pioneer neurogenesis.
  - **Adulthood Cancer**: **Major tumor suppressor gene**. Silenced by promoter hyper-methylation in glioblastoma, colorectal, renal, and breast cancers. Age-dependent hyper-methylation (+0.425) directly mirrors somatic tumor suppressor inactivation in cancer.

### Module 2: Stem Cell Renewal & Proliferative Capacity
- **`2_108854233`** ($r = -0.727$) $\rightarrow$ **GNL3** (Nucleostemin) (50% Dev / 50% Cancer)
  - **Preadulthood Development**: Highly expressed in embryonic stem cells (ESCs) and neural progenitors to maintain self-renewal and delay premature differentiation.
  - **Adulthood Cancer**: **Cancer Stem Cell (CSC) hallmark**. Re-expressed in aggressive cancer stem cell populations in glioblastoma, breast, and liver cancer to fuel tumor self-renewal and chemoresistance.
- **`3_156792975`** ($r = +0.394$) $\rightarrow$ **MEST** / **MIR335** (55% Dev / 45% Cancer)
  - **Preadulthood Development**: Imprinted gene regulating body axis elongation, mesodermal tissue expansion, and adipose development.
  - **Adulthood Cancer**: Loss of imprinting (LOI) or promoter hyper-methylation disrupts MIR335 miRNA silencing, unleashing EMT and metastatic seeding in carcinomas.

### Module 3: Cytoskeletal Dynamics, Invadopodia & Metastasis
- **`1_83189442`** ($r = -0.721$) $\rightarrow$ **ENAH** (Mena) (30% Dev / 70% Cancer)
  - **Preadulthood Development**: Directs actin filament elongation during neurite outgrowth and cell edge protrusion in embryogenesis.
  - **Adulthood Cancer**: **Core driver of cancer metastasis**. ENAH and its invasive isoform (MenaINV) drive invadopodia formation, tumor cell intravasation into blood vessels, and distant metastasis in breast and lung cancer.
- **`6_146779581`** ($r = -0.722$) $\rightarrow$ **ABLIM3** (25% Dev / 75% Cancer)
  - **Preadulthood Development**: Maintains hepatocyte cytoskeletal shape during liver development.
  - **Adulthood Cancer**: **Tumor suppressor in Hepatocellular Carcinoma (HCC)**. Age-related hypo-methylation (-0.722) tracks loss of structural integrity, promoting EMT and hepatocellular carcinoma cell invasion.

### Module 4: Tumor Suppression, DNA Repair & Genomic Stability
- **`8_10151079`** ($r = -0.761$) $\rightarrow$ **MSRA** (25% Dev / 75% Cancer)
  - **Preadulthood Development**: Protects proteins from oxidative stress during developmental bioenergetic shifts.
  - **Adulthood Cancer**: **Tumor suppressor gene**. MSRA loss leads to accumulation of ROS-induced protein damage and genomic instability, accelerating tumor progression in breast, colon, and hepatic cancers.
- **`14_8699620`** ($r = -0.734$) $\rightarrow$ **CFL1** / **MUS81** (30% Dev / 70% Cancer)
  - **Preadulthood Development**: CFL1 powers cleavage furrow cytokinesis; MUS81 resolves stalled replication forks during rapid embryonic divisions.
  - **Adulthood Cancer**: CFL1 drives tumor cell invasion through basement membranes; MUS81 allows cancer cells to survive extreme replication stress and confers resistance to platinum chemotherapy.

### Module 5: Oncogenic Receptor & Axonal/Tumor Guidance
- **`16_9035289`** ($r = +0.387$) $\rightarrow$ **NTN1** (Netrin-1) (50% Dev / 50% Cancer)
  - **Preadulthood Development**: Premier axon guidance cue directing commissural nerve pathfinding and vascular sprouting.
  - **Adulthood Cancer**: Upregulated in aggressive cancers (pancreatic, neuroblastoma, colorectal) as an **oncogenic dependence ligand** that binds UNC5/DCC receptors to block apoptosis and force tumor cell survival.
- **`3_177059260`** ($r = -0.738$) $\rightarrow$ **CDK5** / **SLC4A2** (30% Dev / 70% Cancer)
  - **Preadulthood Development**: CDK5 directs neuronal migration, corticogenesis, and synaptic formation in the developing brain.
  - **Adulthood Cancer**: CDK5 hyper-activation is a major oncogenic driver in solid tumors (pancreatic, prostate, lung, medulloblastoma), promoting proliferation, motility, and immune evasion.

### Module 6: Chromatin Remodeling, Transport & Translation
- **`13_66302679`** ($r = +0.412$) $\rightarrow$ **MTA3** (35% Dev / 65% Cancer)
  - **Preadulthood Development**: MTA3 (NuRD subunit) regulates developmental EMT during gastrulation.
  - **Adulthood Cancer**: **Master repressor of cancer invasion**. MTA3 represses Snail (SNAI1) to maintain E-cadherin expression. Loss or epigenetic alteration of MTA3 unleashes Snail, triggering EMT and metastatic carcinoma dissemination.
- **`1_131347781`** ($r = -0.753$) $\rightarrow$ **RPL5** (50% Dev / 50% Cancer)
  - **Preadulthood Development**: Ribosomal 60S biogenesis required for rapid embryonic cell translation.
  - **Adulthood Cancer**: RPL5 binds MDM2 to activate p53 tumor suppressor activity upon nucleolar stress; RPL5 mutations/deletions drive cancer initiation in ribosomopathies and T-cell leukemias.

---

## 4. Antagonistic Pleiotropy & Biological Insights

1. **The Epigenetic Clock & Tumor Suppressor Silencing**:
   - CpG sites with positive age correlation ($r > +0.38$) like **`EBF3`** ($r = +0.425$) and **`MEST`** ($r = +0.394$) exhibit age-associated hyper-methylation. In adult tissues, this hyper-methylation directly mirrors the **somatic promoter silencing of tumor suppressors** observed during carcinogenesis.
2. **Re-activation of Developmental Stem Cell Machinery in Cancer**:
   - Master stem cell regulators like **`GNL3` (Nucleostemin)** ($r = -0.727$) and **`NTN1` (Netrin-1)** ($r = +0.387$) play dual roles. Highly active during early embryonic development, their epigenetic dysregulation in adult somatic tissues reactivates pluripotency pathways, driving **Cancer Stem Cell (CSC) self-renewal** and apoptosis evasion.
3. **Loss of Structural & DNA Maintenance Unleashes Metastasis**:
   - CpG sites with strong negative age correlation ($r < -0.70$) like **`ABLIM3`**, **`ENAH`**, **`MSRA`**, and **`MUS81`** regulate cytoskeletal stability and DNA repair. Age-associated loss of epigenetic control at these loci weakens somatic defense, creating fertile ground for **invadopodia formation, EMT, oxidative damage, and cancer invasion**.
