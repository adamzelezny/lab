# PIK3CA Variant Differential Expression & Pathway Analysis Report

This report presents the differential expression (PyDESeq2) and pathway enrichment (GSEA Pre-ranked) results for purified **Tumor Cells** and **Neurons** across the PIK3CA mutant genotypes. The variant **R88Q** is treated as the control baseline group (no activity) to highlight the activated properties of **H1047R** and **C420R**.

## 1. Malignant Tumor Cells Analysis

Tumor cells were defined as cells from Cluster 5 expressing glial/tumor markers (`Gfap` or `Gpc3`) but lacking neuronal marker (`Map2`). Comparing these cells isolates the signals that each tumor type sends to the brain microenvironment.

### Contrast: H1047R vs R88Q control

#### Top 10 Upregulated Genes (Significance padj < 0.05):

| index | baseMean          | log2FoldChange      | padj                  |
| ----- | ----------------- | ------------------- | --------------------- |
| Ubc   | 632.3777675855774 | -1.2006255273606747 | 0.0033282725094746324 |
| Ifng  | 11.43556268234031 | -3.0216256895243383 | 0.0033282725094746324 |

#### Top 10 Downregulated Genes (Significance padj < 0.05):

| index | baseMean          | log2FoldChange      | padj                  |
| ----- | ----------------- | ------------------- | --------------------- |
| Ifng  | 11.43556268234031 | -3.0216256895243383 | 0.0033282725094746324 |
| Ubc   | 632.3777675855774 | -1.2006255273606747 | 0.0033282725094746324 |

#### Top Enriched Biological Process Pathways (NES sorted):

| Term                                                                   | NES                 | NOM p-val | FDR q-val             |
| ---------------------------------------------------------------------- | ------------------- | --------- | --------------------- |
| Mitochondrial Gene Expression (GO:0140053)                             | -2.607677987129238  | 0.0       | 0.0                   |
| Mitochondrial Translation (GO:0032543)                                 | -2.5622343820692906 | 0.0       | 0.0                   |
| Cell-Cell Adhesion Via Plasma-Membrane Adhesion Molecules (GO:0098742) | 2.311697959290873   | 0.0       | 0.011763164183246357  |
| Modulation Of Chemical Synaptic Transmission (GO:0050804)              | 2.2996242271285703  | 0.0       | 0.0058815820916231784 |
| tRNA Modification (GO:0006400)                                         | -2.2807921493731684 | 0.0       | 0.0025575013822039965 |
| Regulation Of Mitochondrial Translation (GO:0070129)                   | -2.229033552584034  | 0.0       | 0.005754378109958992  |
| Positive Regulation Of Protein Tyrosine Kinase Activity (GO:0061098)   | 2.208941913765955   | 0.0       | 0.003921054727748786  |
| Positive Regulation Of Wound Healing (GO:0090303)                      | 2.202531270966968   | 0.0       | 0.0029407910458115892 |
| tRNA Processing (GO:0008033)                                           | -2.19594271143587   | 0.0       | 0.004603502487967193  |
| Synapse Organization (GO:0050808)                                      | 2.1911431084248862  | 0.0       | 0.0023526328366492715 |
| Positive Regulation Of Cell Migration (GO:0030335)                     | 2.147332984198391   | 0.0       | 0.001960527363874393  |
| Translation (GO:0006412)                                               | -2.1452822850954902 | 0.0       | 0.00639375345550999   |
| Regulation Of Chromosome Segregation (GO:0051983)                      | 2.1428764919212546  | 0.0       | 0.0033609040523561017 |
| Fatty Acid Oxidation (GO:0019395)                                      | -2.1201027695109    | 0.0       | 0.010960720209445699  |
| Ribosome Assembly (GO:0042255)                                         | -2.1029450129907716 | 0.0       | 0.009590630183264986  |

### Contrast: C420R vs R88Q control

#### Top 10 Upregulated Genes (Significance padj < 0.05):

*No significant genes found.*

#### Top 10 Downregulated Genes (Significance padj < 0.05):

*No significant genes found.*

#### Top Enriched Biological Process Pathways (NES sorted):

| Term                                                    | NES                 | NOM p-val | FDR q-val             |
| ------------------------------------------------------- | ------------------- | --------- | --------------------- |
| Striated Muscle Contraction (GO:0006941)                | 2.812924332173684   | 0.0       | 0.0                   |
| Skeletal Muscle Contraction (GO:0003009)                | 2.585091261097288   | 0.0       | 0.0                   |
| Response To Type II Interferon (GO:0034341)             | -2.509335781173912  | 0.0       | 0.0                   |
| Muscle Contraction (GO:0006936)                         | 2.427665976782573   | 0.0       | 0.0                   |
| Myofibril Assembly (GO:0030239)                         | 2.3802598276233318  | 0.0       | 0.0                   |
| Regulation Of Muscle Contraction (GO:0006937)           | 2.344998908070245   | 0.0       | 0.0                   |
| Ribonucleoprotein Complex Biogenesis (GO:0022613)       | -2.3411270421228236 | 0.0       | 0.0                   |
| Ribosome Biogenesis (GO:0042254)                        | -2.3265223833595026 | 0.0       | 0.0                   |
| T Cell Receptor Signaling Pathway (GO:0050852)          | -2.3158418252531017 | 0.0       | 0.0                   |
| Skeletal Muscle Organ Development (GO:0060538)          | 2.27975201843418    | 0.0       | 0.0                   |
| Regulation Of Mitotic Spindle Organization (GO:0060236) | -2.227778822004742  | 0.0       | 0.002519893899204244  |
| Muscle Organ Development (GO:0007517)                   | 2.2197959507602025  | 0.0       | 0.0                   |
| Maturation Of SSU-rRNA (GO:0030490)                     | -2.218133895132849  | 0.0       | 0.0020999115826702036 |
| Ribosomal Small Subunit Biogenesis (GO:0042274)         | -2.203717337833127  | 0.0       | 0.0017999242137173171 |
| Cellular Response To Type II Interferon (GO:0071346)    | -2.203629886031054  | 0.0       | 0.0015749336870026524 |

### Contrast: H1047R vs C420R direct

#### Top 10 Upregulated Genes (Significance padj < 0.05):

| index | baseMean          | log2FoldChange     | padj                  |
| ----- | ----------------- | ------------------ | --------------------- |
| Ubc   | 632.3777675855774 | -1.643854449919328 | 6.792980674593296e-08 |

#### Top 10 Downregulated Genes (Significance padj < 0.05):

| index | baseMean          | log2FoldChange     | padj                  |
| ----- | ----------------- | ------------------ | --------------------- |
| Ubc   | 632.3777675855774 | -1.643854449919328 | 6.792980674593296e-08 |

#### Top Enriched Biological Process Pathways (NES sorted):

| Term                                                                | NES                 | NOM p-val | FDR q-val             |
| ------------------------------------------------------------------- | ------------------- | --------- | --------------------- |
| Striated Muscle Contraction (GO:0006941)                            | -2.762759849649181  | 0.0       | 0.0                   |
| Skeletal Muscle Contraction (GO:0003009)                            | -2.498387355321755  | 0.0       | 0.0                   |
| Muscle Contraction (GO:0006936)                                     | -2.420035582993681  | 0.0       | 0.0                   |
| Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775) | -2.363257431235904  | 0.0       | 0.0                   |
| Sarcomere Organization (GO:0045214)                                 | -2.335739258482514  | 0.0       | 0.0                   |
| Granulocyte Chemotaxis (GO:0071621)                                 | 2.3322545097555945  | 0.0       | 0.0                   |
| Cellular Respiration (GO:0045333)                                   | -2.330149293103737  | 0.0       | 0.0                   |
| Lymphocyte Chemotaxis (GO:0048247)                                  | 2.320361878192411   | 0.0       | 0.0                   |
| Regulation Of Chemokine Production (GO:0032642)                     | 2.32035307910058    | 0.0       | 0.0                   |
| Neutrophil Chemotaxis (GO:0030593)                                  | 2.312415169559415   | 0.0       | 0.0                   |
| Aerobic Respiration (GO:0009060)                                    | -2.26796272934089   | 0.0       | 0.0007478103179128614 |
| Positive Regulation Of Cell Cycle Process (GO:0090068)              | 2.267025388647985   | 0.0       | 0.0                   |
| Heart Contraction (GO:0060047)                                      | -2.2591920876668943 | 0.0       | 0.0006543340281737538 |
| Neutrophil Migration (GO:1990266)                                   | 2.241250235313755   | 0.0       | 0.0023687349027579235 |
| Myofibril Assembly (GO:0030239)                                     | -2.2320694841043407 | 0.0       | 0.0011632604945311178 |

## 2. Microenvironmental Neurons Analysis

Neurons were defined as cells expressing neuronal markers (`Map2` or `Dlg4`) but lacking glial marker (`Gfap`). Comparing these cells reveals the transcription changes occurring in host neurons in response to the surrounding tumor variants.

### Contrast: H1047R vs R88Q control

#### Top 10 Upregulated Genes (Significance padj < 0.05):

| index   | baseMean           | log2FoldChange      | padj                  |
| ------- | ------------------ | ------------------- | --------------------- |
| Acaa1b  | 34.95799443471443  | 3.781005491640278   | 4.274436354639093e-08 |
| Slc15a2 | 126.77929114216309 | 1.7577146957298062  | 0.0028559309639679684 |
| Il1a    | 229.8244873887872  | 1.5852794105904031  | 0.0028559309639679684 |
| Eci1    | 863.5607299217365  | -0.7927994329148137 | 0.002208106989004102  |
| Ubc     | 4343.982206894171  | -0.9726148586603368 | 0.0014478955804355816 |
| Chil1   | 111.08280935961152 | -2.912802907280347  | 0.040073495655230766  |
| Gm45351 | 11.675857782974365 | -3.052662372612423  | 0.02483358558931049   |

#### Top 10 Downregulated Genes (Significance padj < 0.05):

| index   | baseMean           | log2FoldChange      | padj                  |
| ------- | ------------------ | ------------------- | --------------------- |
| Gm45351 | 11.675857782974365 | -3.052662372612423  | 0.02483358558931049   |
| Chil1   | 111.08280935961152 | -2.912802907280347  | 0.040073495655230766  |
| Ubc     | 4343.982206894171  | -0.9726148586603368 | 0.0014478955804355816 |
| Eci1    | 863.5607299217365  | -0.7927994329148137 | 0.002208106989004102  |
| Il1a    | 229.8244873887872  | 1.5852794105904031  | 0.0028559309639679684 |
| Slc15a2 | 126.77929114216309 | 1.7577146957298062  | 0.0028559309639679684 |
| Acaa1b  | 34.95799443471443  | 3.781005491640278   | 4.274436354639093e-08 |

#### Top Enriched Biological Process Pathways (NES sorted):

| Term                                                            | NES                 | NOM p-val | FDR q-val |
| --------------------------------------------------------------- | ------------------- | --------- | --------- |
| Mitochondrial Translation (GO:0032543)                          | -2.9591966538367798 | 0.0       | 0.0       |
| Ribosome Biogenesis (GO:0042254)                                | -2.9120737697732744 | 0.0       | 0.0       |
| Mitochondrial Gene Expression (GO:0140053)                      | -2.8605003601753722 | 0.0       | 0.0       |
| Ribonucleoprotein Complex Biogenesis (GO:0022613)               | -2.672916124910662  | 0.0       | 0.0       |
| Translation (GO:0006412)                                        | -2.6582787046962295 | 0.0       | 0.0       |
| Protein Targeting To Mitochondrion (GO:0006626)                 | -2.636213876180383  | 0.0       | 0.0       |
| rRNA Processing (GO:0006364)                                    | -2.6348508052073596 | 0.0       | 0.0       |
| Ribosomal Small Subunit Biogenesis (GO:0042274)                 | -2.620661822333236  | 0.0       | 0.0       |
| tRNA Modification (GO:0006400)                                  | -2.5964171544229964 | 0.0       | 0.0       |
| Aerobic Respiration (GO:0009060)                                | -2.5906263472640507 | 0.0       | 0.0       |
| Positive Regulation Of MAPK Cascade (GO:0043410)                | 2.5409478893668633  | 0.0       | 0.0       |
| Hemopoiesis (GO:0030097)                                        | 2.5283714995850026  | 0.0       | 0.0       |
| NADH Dehydrogenase Complex Assembly (GO:0010257)                | -2.52118434259322   | 0.0       | 0.0       |
| Mitochondrial Respiratory Chain Complex I Assembly (GO:0032981) | -2.52118434259322   | 0.0       | 0.0       |
| Cytoplasmic Translational Initiation (GO:0002183)               | -2.5058446869317517 | 0.0       | 0.0       |

### Contrast: C420R vs R88Q control

#### Top 10 Upregulated Genes (Significance padj < 0.05):

| index   | baseMean           | log2FoldChange     | padj                  |
| ------- | ------------------ | ------------------ | --------------------- |
| Acaa1b  | 34.95799443471443  | 3.93395535000288   | 3.136055964308353e-08 |
| Slc15a2 | 126.77929114216309 | 1.684724954625178  | 0.04553889452409356   |
| Cwc22   | 782.5745558776416  | 1.4322482278916153 | 0.038631191442261084  |

#### Top 10 Downregulated Genes (Significance padj < 0.05):

| index   | baseMean           | log2FoldChange     | padj                  |
| ------- | ------------------ | ------------------ | --------------------- |
| Cwc22   | 782.5745558776416  | 1.4322482278916153 | 0.038631191442261084  |
| Slc15a2 | 126.77929114216309 | 1.684724954625178  | 0.04553889452409356   |
| Acaa1b  | 34.95799443471443  | 3.93395535000288   | 3.136055964308353e-08 |

#### Top Enriched Biological Process Pathways (NES sorted):

| Term                                                                      | NES                 | NOM p-val | FDR q-val |
| ------------------------------------------------------------------------- | ------------------- | --------- | --------- |
| Muscle Contraction (GO:0006936)                                           | 2.82451192011456    | 0.0       | 0.0       |
| Myofibril Assembly (GO:0030239)                                           | 2.7403267649657588  | 0.0       | 0.0       |
| Striated Muscle Contraction (GO:0006941)                                  | 2.709227286939802   | 0.0       | 0.0       |
| Sarcomere Organization (GO:0045214)                                       | 2.5707621709165505  | 0.0       | 0.0       |
| Sister Chromatid Segregation (GO:0000819)                                 | -2.492418426008752  | 0.0       | 0.0       |
| Skeletal Muscle Tissue Development (GO:0007519)                           | 2.4704346003582778  | 0.0       | 0.0       |
| Actomyosin Structure Organization (GO:0031032)                            | 2.450412433244521   | 0.0       | 0.0       |
| Negative Regulation Of Mitotic Metaphase/Anaphase Transition (GO:0045841) | -2.4468489069806747 | 0.0       | 0.0       |
| Heart Contraction (GO:0060047)                                            | 2.405880474433092   | 0.0       | 0.0       |
| Mitotic Nuclear Division (GO:0140014)                                     | -2.3560671681281358 | 0.0       | 0.0       |
| Ribosome Biogenesis (GO:0042254)                                          | -2.3556942518009976 | 0.0       | 0.0       |
| Mitotic Spindle Checkpoint Signaling (GO:0071174)                         | -2.3464012299344477 | 0.0       | 0.0       |
| Spindle Assembly Checkpoint Signaling (GO:0071173)                        | -2.3464012299344477 | 0.0       | 0.0       |
| Mitotic Spindle Assembly Checkpoint Signaling (GO:0007094)                | -2.3464012299344477 | 0.0       | 0.0       |
| Regulation Of DNA-templated DNA Replication (GO:0090329)                  | -2.32452569681761   | 0.0       | 0.0       |

### Contrast: H1047R vs C420R direct

#### Top 10 Upregulated Genes (Significance padj < 0.05):

| index  | baseMean           | log2FoldChange      | padj                   |
| ------ | ------------------ | ------------------- | ---------------------- |
| Slc7a2 | 85.37912230835742  | 1.8129111503130066  | 0.006678581382664082   |
| Ubc    | 4343.982206894171  | -1.4174951256831643 | 3.9444606203465833e-10 |
| Sbsn   | 49.113008931778765 | -2.4819635180378326 | 0.0036720178501482804  |

#### Top 10 Downregulated Genes (Significance padj < 0.05):

| index  | baseMean           | log2FoldChange      | padj                   |
| ------ | ------------------ | ------------------- | ---------------------- |
| Sbsn   | 49.113008931778765 | -2.4819635180378326 | 0.0036720178501482804  |
| Ubc    | 4343.982206894171  | -1.4174951256831643 | 3.9444606203465833e-10 |
| Slc7a2 | 85.37912230835742  | 1.8129111503130066  | 0.006678581382664082   |

#### Top Enriched Biological Process Pathways (NES sorted):

| Term                                                                | NES                 | NOM p-val | FDR q-val |
| ------------------------------------------------------------------- | ------------------- | --------- | --------- |
| Aerobic Respiration (GO:0009060)                                    | -3.0130375931839906 | 0.0       | 0.0       |
| Oxidative Phosphorylation (GO:0006119)                              | -2.9749987498810877 | 0.0       | 0.0       |
| Aerobic Electron Transport Chain (GO:0019646)                       | -2.900314004065589  | 0.0       | 0.0       |
| Proton Motive Force-Driven Mitochondrial ATP Synthesis (GO:0042776) | -2.822390676099021  | 0.0       | 0.0       |
| Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775) | -2.8168898886751124 | 0.0       | 0.0       |
| Cellular Respiration (GO:0045333)                                   | -2.7856729348499494 | 0.0       | 0.0       |
| Mitochondrial Gene Expression (GO:0140053)                          | -2.7602283066501827 | 0.0       | 0.0       |
| Proton Motive Force-Driven ATP Synthesis (GO:0015986)               | -2.7381951424180637 | 0.0       | 0.0       |
| Mitochondrial Respiratory Chain Complex I Assembly (GO:0032981)     | -2.723723550650148  | 0.0       | 0.0       |
| NADH Dehydrogenase Complex Assembly (GO:0010257)                    | -2.723723550650148  | 0.0       | 0.0       |
| Mitochondrial Respiratory Chain Complex Assembly (GO:0033108)       | -2.7222680657651788 | 0.0       | 0.0       |
| Mitochondrial Electron Transport, NADH To Ubiquinone (GO:0006120)   | -2.7083679178999267 | 0.0       | 0.0       |
| Myofibril Assembly (GO:0030239)                                     | -2.6505639349681314 | 0.0       | 0.0       |
| Notch Signaling Pathway (GO:0007219)                                | 2.6308965699644307  | 0.0       | 0.0       |
| Mitochondrial Translation (GO:0032543)                              | -2.6126218386750524 | 0.0       | 0.0       |

