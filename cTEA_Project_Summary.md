# cTEA 项目总结

**版本**: 0.1.0  
**日期**: 2026-04-28  
**目标**: 融合MEGAnE的速度与xTEA的准确性，实现快速、准确的转座子插入（TEI）检测  

---

## 一、项目概述

### 1.1 背景
cTEA（CRAM-optimized Transposable Element Analyzer）是一个用于检测基因组中转座子插入（Transposable Element Insertion，TEI）的工具。它结合了：
- **MEGAnE的速度**：C++实现，单次扫描
- **xTEA的准确性**：严格的过滤逻辑和后处理链  

### 1.2 核心优化vs xTEA
| 维度 | xTEA | cTEA |
|------|------|------|
| **扫描策略** | 分染色体，多次扫描BAM | **单次扫描**，内存存储 |
| **处理时间** | ~4小时（30x WGS） | **目标20-30分钟** |
| **中间文件** | 有（temp files） | **无**（全内存） |
| **CRAM支持** | 部分 | **原生支持**（htslib） |
| **多线程** | 有限 | **OpenMP + ThreadPool** |

### 1.3 当前状态
- ✅ **Phase 1（C++核心）**：完成
- ⚠️ **Phase 2（Python包装和过滤）**：进行中
- ⚠️ **Phase 3（CRAM支持和批量运行）**：待开始

---

## 二、已实现的功能

### 2.1 C++核心模块 ✅

#### **文件列表**：
| 文件 | 功能 | 状态 |
|------|------|------|
| **bam_processor.h/cpp** | BAM/CRAM处理器（htslib） | ✅ 完成 |
| **candidate_detector.h/cpp** | 候选检测器（单次扫描） | ✅ 完成 + 聚集算法 |
| **bindings.cpp** | pybind11 Python绑定 | ✅ 完成 |
| **extract_reads.cpp** | 读取提取器（合并MEGAnE） | ✅ 完成 |
| **kmer_filter.cpp** | k-mer过滤器（C++优化） | ✅ 完成 |
| **Makefile** | 编译配置（支持pybind11） | ✅ 完成 |

#### **关键功能**：
1. ✅ **单次扫描架构**（KEY OPTIMIZATION）：
   ```cpp
   // 单次扫描，内存存储证据
   std::unordered_map<std::string, 
                     std::unordered_map<uint32_t, CandidateEvidence>> evidence;
   bam_processor_->scan_all_evidence(evidence, callback);
   ```

2. ✅ **聚集算法**（Priority 1改进）：
   - 参考xTEA的`chain_regions`函数
   - 100bp窗口（PEAK_WINDOW）
   - 聚类级别证据求和
   - 检测双峰（two-side events）

3. ✅ **阈值调整**（Priority 2优化）：
   - `min_clip_len_`：2 → 1（适应低信号）
   - `clip_cluster_diff_cutoff_`：300 → 100bp（更严格）
   - `sva_clip_cluster_diff_cutoff_`：200 → 100bp

4. ✅ **pybind11绑定**：
   - 暴露`BamProcessor`, `CandidateDetector`类到Python
   - 暴露`BamHeaderInfo`, `CandidateEvidence`, `ScanStatistics`结构体
   - 支持Python回调函数（进度监控）

### 2.2 Python层 ⚠️

| 文件 | 功能 | 状态 |
|------|------|------|
| **ctea_core.py** | C++核心包装器 | ⚠️ 占位符（需更新） |
| **ctea_main.py** | 主流程脚本 | ⚠️ 框架完成，部分占位符 |
| **test_pybind11.py** | pybind11测试脚本 | ✅ 完成 |
| **filters/** | xTEA过滤器移植 | ❌ 未开始 |

### 2.3 测试和文档 ✅

| 文件 | 用途 | 状态 |
|------|------|------|
| **cTEA_Test_Guide.md** | 服务器集群测试指南 | ✅ 完成 |
| **cTEA_Push_Guide.md** | SSH推送指南 | ✅ 完成 |
| **cTEA_Development_Plan.md** | 开发计划 | ✅ 完成 |
| **MEGA-xTEA_Comprehensive_Gap_Analysis.md** | 差距分析 | ✅ 完成 |
| **evidence_clustering_comparison.md** | 聚集算法比较 | ✅ 完成 |

---

## 三、代码结构

### 3.1 目录结构
```
cTEA/
├── cpp/                          # C++核心
│   ├── bam_processor.h/cpp       ✅ BAM/CRAM处理
│   ├── candidate_detector.h/cpp  ✅ 候选检测（+聚集算法）
│   ├── bindings.cpp             ✅ pybind11绑定
│   ├── extract_reads.cpp       ✅ 读取提取
│   ├── kmer_filter.cpp         ✅ k-mer过滤
│   ├── Makefile                ✅ 编译配置
│   ├── ThreadPool.h           ✅ 线程池
│   ├── complementary_seq.hpp    ✅ 辅助函数
│   ├── dna_to_2bit.hpp       ✅ 辅助函数
│   ├── libcTEA.a              ✅ 静态库
│   ├── libcTEA.so             ✅ 共享库
│   └── libcTEA_pybind11.so    ✅ pybind11模块
│
├── python/                       # Python层
│   ├── ctea_core.py           ⚠️ 需更新（使用pybind11）
│   ├── ctea_main.py           ⚠️ 需更新
│   ├── test_pybind11.py       ✅ 测试脚本
│   ├── filters/               ❌ 待移植xTEA过滤器
│   ├── genotyping/            ❌ 待实现
│   └── utils/                 ❌ 待实现
│
├── tests/                       # 测试数据
│   └── generate_test_data.py  ✅ 测试数据生成器
│
├── docker/                      # Docker配置
│   └── .gitkeep
│
├── wdl/                         # WDL脚本
│   └── .gitkeep
│
└── 文档/
    ├── cTEA_Development_Plan.md              ✅
    ├── cTEA_Test_Guide.md                 ✅
    ├── cTEA_Push_Guide.md                 ✅
    ├── MEGA-xTEA_Comprehensive_Gap_Analysis.md  ✅
    └── evidence_clustering_comparison.md       ✅
```

### 3.2 核心类关系
```
cTEA::BamProcessor
    ├── initialize()      // 初始化BAM/CRAM
    ├── scan_all_evidence()  // 单次扫描，收集证据
    ├── get_header_info()   // 获取文件头信息
    ├── get_statistics()   // 获取扫描统计
    ├── is_cram()         // 检查是否为CRAM
    ├── extract_clipped_reads()
    └── extract_discordant_pairs()

cTEA::CandidateDetector
    ├── initialize()              // 初始化检测器
    ├── detect_candidates()       // 主检测函数
    ├── build_candidates()        // 从证据构建候选位点
    ├── cluster_nearby_candidates()  // ✅ NEW: 聚集附近位点
    ├── determine_support_type()   // 判定支持类型
    ├── filter_candidates()       // 应用过滤器
    ├── process_cluster()         // 处理单个聚类
    └── get_statistics()         // 获取检测统计
```

---

## 四、技术特点

### 4.1 性能优化
| 优化点 | xTEA | cTEA |
|---------|------|------|
| **扫描次数** | N次（N=染色体数） | **1次**（单次扫描） |
| **内存使用** | 低（单染色体） | 高（全基因组证据） |
| **中间文件** | 有（CLIP_TMP） | **无** |
| **数据结构** | 文件-based | **内存哈希表** |
| **预计速度** | 基准（4小时） | **快3-5倍**（目标20-30分钟） |

### 4.2 CRAM支持
- ✅ 使用htslib原生支持CRAM格式
- ✅ 需要参考FASTA文件（通过`--reference`参数传递）
- ✅ 自动检测文件类型（`is_cram()`）

### 4.3 Python集成
- ✅ pybind11绑定（直接调用C++，无subprocess开销）
- ✅ ctypes备用方案（标准库）
- ⚠️ 当前`ctea_core.py`还是占位符，未使用pybind11

---

## 五、编译和使用

### 5.1 系统依赖
```bash
# Ubuntu/Debian
sudo apt-get install g++ make libhts-dev zlib1g-dev libbz2-dev liblzma-dev python3 python3-pip

# macOS (Homebrew)
brew install htslib python3
```

### 5.2 Python依赖
```bash
python3 -m venv cTEA_venv
source cTEA_venv/bin/activate
pip install pybind11 numpy scipy scikit-learn
```

### 5.3 编译
```bash
cd cpp/

# 编译所有目标
make clean && make all

# 编译pybind11模块
make libcTEA_pybind11.so

# 创建符号链接
ln -sf libcTEA_pybind11.so libcTEA.so
```

### 5.4 测试
```bash
cd /path/to/cTEA

# 基础测试
python3 test_pybind11.py

# BAM扫描测试
python3 test_cTEA_scan.py /data/test.bam 4

# CRAM扫描测试（需要参考基因组）
python3 test_cTEA_scan.py /data/test.cram /data/genome.fa 4

# 候选检测测试
python3 test_cTEA_detect.py /data/test.cram output.bed /data/genome.fa 4
```

---

## 六、当前限制

### 6.1 未实现的功能
| 功能 | 状态 | 优先级 |
|------|------|----------|
| **mate位置证据收集** | ❌ 未实现 | Priority 3 |
| **consensus比对** | ❌ 未实现 | Priority 4 |
| **xTEA后过滤链移植** | ❌ 未开始 | Phase 2核心 |
| **转导检测** | ❌ 未实现 | Phase 2 |
| **TSD验证** | ❌ 未实现 | Phase 2 |
| **局部组装验证** | ❌ 未实现 | Phase 2 |
| **ML特征提取修复** | ⚠️ 近似值 | Phase 2 |

### 6.2 已知问题
1. ⚠️ **CandidateDetector::detect_candidates** 在pybind11绑定中是占位符（现已连通C++实现）
2. ⚠️ **未使用私有字段警告**：`min_af_cutoff_`, `total_reads_processed_`, `candidates_found_`
3. ⚠️ **过滤链不完整**：当前只有基本过滤，缺少xTEA的1247行`x_post_filter.py`

---

## 七、下一步计划

### 7.1 Phase 2: Python包装和过滤（进行中）
1. **更新ctea_core.py**：使用pybind11而非占位符
2. **移植xTEA过滤器**：
   - `x_post_filter.py` → `python/filters/post_filter.py`
   - `x_transduction.py` → `python/filters/transduction.py`
   - `x_TSD.py` → `python/filters/tsd_validator.py`
   - `x_local_assembly.py` → `python/filters/local_assembly.py`

### 7.2 Phase 3: CRAM支持和批量运行
1. **Docker镜像**：`docker/Dockerfile`
2. **WDL脚本**：`wdl/ctea.wdl`
3. **CRAM性能测试**：不同压缩级别对比

### 7.3 Phase 4: 测试和性能调优
1. **单元测试**：`tests/test_cpp.py`, `tests/test_filters.py`
2. **端到端测试**：30x WGS CRAM全流程
3. **性能基准**：vs xTEA和MEGAnE

---

## 八、项目统计

### 8.1 代码统计
| 语言 | 文件数 | 代码行数（估计） |
|------|--------|------------------|
| **C++** | 6个核心文件 | ~3000行 |
| **Python** | 2个脚本 | ~500行（不含待移植） |
| **Markdown文档** | 6个文档 | ~1500行 |

### 8.2 开发进度
```
Phase 1 (C++核心): ████████████ 100% ✅
Phase 2 (Python过滤): ████░░░░░░░   15% ⚠️
Phase 3 (CRAM/批量):  █░░░░░░░░░    0% ❌
Phase 4 (测试调优):  █░░░░░░░░░    0% ❌
```

### 8.3 GitHub仓库
- **URL**: https://github.com/Mirror-fish/cTEA
- **最新提交**: "Phase 1: Implement pybind11 bindings and test infrastructure"
- **分支**: main
- **SSH**: `git@github.com:Mirror-fish/cTEA.git`

---

## 九、关键成果

### 9.1 技术成果
1. ✅ **单次扫描架构**：比xTEA快3-5倍（理论）
2. ✅ **聚集算法**：参考xTEA的chain_regions，100bp窗口
3. ✅ **pybind11绑定**：Python可以直接调用C++核心
4. ✅ **CRAM原生支持**：使用htslib，自动检测文件类型
5. ✅ **测试基础设施**：测试脚本、操作指南、推送指南

### 9.2 文档成果
1. ✅ **开发计划**：cTEA_Development_Plan.md
2. ✅ **差距分析**：MEGA-xTEA_Comprehensive_Gap_Analysis.md
3. ✅ **聚集比较**：evidence_clustering_comparison.md
4. ✅ **测试指南**：cTEA_Test_Guide.md
5. ✅ **推送指南**：cTEA_Push_Guide.md
6. ✅ **项目总结**：cTEA_Project_Summary.md（本文件）

---

## 十、致谢和参考

### 10.1 参考项目
- **xTEA**: https://github.com/pegasuswang/xTea
- **MEGAnE**: https://github.com/ponnhidev/MEGAnE
- **pybind11**: https://github.com/pybind/pybind11

### 10.2 核心算法来源
- xTEA的`x_TEI_locator.py`：分染色体扫描逻辑
- xTEA的`x_post_filter.py`：后过滤链（1247行）
- MEGAnE的`extract_discordant.cpp`：C++读取提取
- MEGAnE的`1_indiv_call_genotype.py`：主流程框架

---

## 结束

**cTEA Phase 1 完成！**  
**下一步**：继续Phase 2（移植xTEA过滤器）  

**项目仓库**: https://github.com/Mirror-fish/cTEA  
**问题反馈**：通过GitHub Issues或Pull Request  

---
**文档结束** | 生成时间：2026-04-28 18:10 (Europe/Madrid)