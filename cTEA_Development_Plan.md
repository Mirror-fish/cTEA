# cTEA (CRAM-optimized Transposable Element Analyzer)
## 开发设计计划

**目标**: 将30x WGS样本运行时间从4小时降至30分钟，支持CRAM，提供Docker/WDL批量运行

---

## 1. 核心架构设计

### 1.1 整体架构：C++核心 + Python包装

```
cTEA/
├── cpp/                    # C++核心（参考MEGAnE的cpp/）
│   ├── extract_reads.cpp       # 合并extract_discordant + extract_unmapped
│   ├── kmer_filter.cpp        # k-mer过滤（从MEGAnE移植优化）
│   ├── candidate_detector.cpp  # 候选检测（整合xTEA的clip+disc逻辑）
│   ├── bam_processor.h/cpp   # CRAM/BAM处理（htslib）
│   └── ThreadPool.h           # 线程池（从MEGAnE复制）
├── python/
│   ├── ctea_main.py         # 主入口（替代xTEA的x_TEA_main.py）
│   ├── filters/              # xTEA过滤逻辑的高效实现
│   │   ├── sva_filter.py    # SVA专用过滤（从MEGA-xTEA移植优化）
│   │   ├── fp_filter.py     # 假阳性过滤
│   │   └── polyA_detector.py
│   ├── genotyping/           # 基因分型
│   │   ├── ml_genotype.py  # ML基因分型（支持CRAM直接读取特征）
│   │   └── gaussian.py     # Gaussian fallback
│   └── utils/
│       ├── cram_utils.py    # CRAM支持工具
│       └── config.py        # 统一配置（合并xTEA和MEGAnE参数）
├── wdl/
│   ├── ctea.wdl             # WDL工作流
│   └── ctea_inputs.json     # 输入模板
├── docker/
│   └── Dockerfile           # Docker镜像定义
└── tests/                  # 测试套件
```

### 1.2 关键优化策略

| 瓶颈 | xTEA实现 | cTEA优化方案 | 预期加速 |
|------|----------|-----------|----------|
| BAM/CRAM读取 | Python pysam逐文件处理 | C++ htslib多线程读取 + 内存缓存 | 5-10x |
| 分染色体操作 | 按染色体拆分→处理→合并（大量I/O） | 单次扫描+内存数据结构（避免中间文件） | 3-5x |
| 候选检测 | clip和disc分两步，重复扫描 | 单次扫描同时收集clip+disc证据 | 2x |
| k-mer过滤 | 无 | 从MEGAnE移植C++实现 | 10x（vs无过滤） |
| 基因分型 | 需要重新扫描BAM | 首次扫描时缓存特征到内存/临时文件 | 2x |
| 并行化 | Python multiprocessing（进程间序列化开销） | C++线程池 + OpenMP | 2-4x |

---

## 2. 详细开发步骤

### Phase 1: C++核心开发（第1-2周）

#### Step 1.1: BAM/CRAM处理器（3天）
- **文件**: `cpp/bam_processor.h`, `cpp/bam_processor.cpp`
- **功能**: 
  - 使用htslib支持BAM/CRAM（通过`--reference`参数支持CRAM解密）
  - 多线程读取（`OMP_NUM_THREADS`环境变量）
  - 内存缓存常用区域（参考MEGAnE的`extract_discordant.cpp`）
- **参考**: MEGAnE的`cpp/extract_discordant.cpp` lines 1-100

#### Step 1.2: 合并read提取器（3天）
- **文件**: `cpp/extract_reads.cpp`
- **功能**: 
  - 合并MEGAnE的`extract_discordant.cpp`和`extract_unmapped.cpp`
  - 单次扫描同时提取：discordant pairs、unmapped reads、clipped reads
  - 输出：标准化的中间格式（避免xTEA的多临时文件问题）
- **关键优化**: 避免xTEA的`x_TEI_locator.py` lines 112-168的分染色体合并瓶颈

#### Step 1.3: k-mer过滤器（2天）
- **文件**: `cpp/kmer_filter.cpp`
- **功能**: 从MEGAnE移植`0_build_kmer_set.py`的C++版本
- **参考**: MEGAnE的`cpp/convert_rep_to_2bit.cpp`

#### Step 1.4: 候选检测器（4天）
- **文件**: `cpp/candidate_detector.cpp`
- **功能**: 
  - 整合xTEA的clip检测逻辑（`x_TEI_locator.py`的`call_TEI_candidate_sites_from_clip_reads_v2`）
  - 整合discordant pair检测（避免xTEA的两步法）
  - 输出：候选位点+证据计数（BED格式）
- **关键改进**: 不像xTEA那样分染色体处理，而是单次扫描完成

### Phase 2: Python包装和过滤逻辑（第3-4周）

#### Step 2.1: 主入口和管道（3天）
- **文件**: `python/ctea_main.py`
- **功能**: 
  - 类似MEGAnE的`1_indiv_call_genotype.py`的简洁流程
  - 但避免subprocess调用，直接import C++模块（使用pybind11或ctypes）
  - 支持CRAM输入（传递给C++时带`--reference`参数）

#### Step 2.2: xTEA过滤逻辑移植（5天）
- **文件**: `python/filters/sva_filter.py`, `fp_filter.py`
- **功能**: 
  - 从MEGA-xTEA移植并优化（修复之前分析的问题）
  - 使用C++输出的标准化格式，避免重复解析
  - 实现xTEA的`x_post_filter.py`核心逻辑（但用Python高效实现）

#### Step 2.3: PolyA检测优化（2天）
- **文件**: `python/filters/polyA_detector.py`
- **功能**: 从xTEA移植`x_polyA.py`，但支持从C++传递的序列

#### Step 2.4: ML基因分型（3天）
- **文件**: `python/genotyping/ml_genotype.py`
- **功能**: 
  - 从MEGA-xTEA移植，但修复特征提取问题
  - **关键**: 从C++首次扫描时缓存特征，避免重新扫描BAM
  - 支持sklearn和Deep Forest模型

### Phase 3: CRAM支持和批量运行（第5周）

#### Step 3.1: CRAM完整支持（2天）
- **文件**: `python/utils/cram_utils.py`
- **功能**: 
  - 检查CRAM依赖（htslib with CRAM support）
  - 自动传递`--reference`给htslib
  - 测试CRAM vs BAM性能差异

#### Step 3.2: Docker镜像（2天）
- **文件**: `docker/Dockerfile`
- **内容**:
  ```dockerfile
  FROM biocontainers/htslib:v1.9-2-ubuntu
  RUN apt-get update && apt-get install -y python3 python3-pip
  COPY . /ctea/
  WORKDIR /ctea
  RUN cd cpp && make && pip3 install -r python/requirements.txt
  ENTRYPOINT ["python3", "/ctea/python/ctea_main.py"]
  ```

#### Step 3.3: WDL脚本（3天）
- **文件**: `wdl/ctea.wdl`
- **功能**: 
  - 支持单样本和批量运行
  - 输入：CRAM/BAM列表 + reference FASTA
  - 输出：VCF文件（合并MEGAnE和xTEA的输出格式）
  - 集成到现有WDL管道（参考MEGA-xTEA的`wdl/mega_xtea_call.wdl`）

### Phase 4: 测试和性能调优（第6周）

#### Step 4.1: 单元测试（2天）
- **文件**: `tests/test_cpp.py`, `tests/test_filters.py`
- **内容**: 使用合成数据测试每个模块

#### Step 4.2: 端到端测试（2天）
- **测试数据**: 30x WGS CRAM文件（~25GB）
- **对比**: xTEA (4小时) vs cTEA (目标30分钟)

#### Step 4.3: 性能分析和调优（3天）
- **工具**: perf, valgrind
- **目标**: 识别剩余瓶颈，优化内存使用

---

## 3. 关键技术决策

### 3.1 避免xTEA的分染色体问题
**xTEA的问题** (`x_TEI_locator.py` lines 112-168):
```python
# 问题：每个染色体都要读取整个BAM，然后合并
for chrm in m_chrms:
    for i in range(cnt):  # 每个BAM文件
        # 读取临时文件，按染色体过滤
        if tmp_chrm != chrm:
            continue
```

**cTEA解决方案**:
- C++单次扫描BAM，在内存中按染色体组织数据
- 使用哈希表（unordered_map）存储候选位点，避免中间文件
- 输出时再按染色体分组写入（如果需要）

### 3.2 CRAM支持
- 使用htslib的`sam_read1()` API，自动检测BAM/CRAM
- CRAM需要reference FASTA：通过环境变量`CRAM_REF`或命令行参数传递
- 参考MEGAnE的`cpp/parse_fai.hpp`实现FASTA索引读取

### 3.3 整合xTEA的过滤逻辑
**保留xTEA的优势**:
- SVA专用过滤（`x_post_filter.py`的post_processing_SVA）
- AF冲突检查（4-ratio test）
- 低分歧参考TE拷贝过滤

**但用更高效方式实现**:
- 不从BAM重新扫描，而是从首次扫描缓存证据
- 使用numpy进行向量化计算（vs xTEA的纯Python循环）

---

## 4. 开发时间表

| 周次 | 任务 | 里程碑 |
|------|------|--------|
| 第1周 | C++ BAM/CRAM处理器 | 能读取CRAM并提取基本统计 |
| 第2周 | C++候选检测器 | 完成clip+disc检测，输出候选BED |
| 第3周 | Python主入口 + 过滤逻辑 | 能运行完整管道（无CRAM） |
| 第4周 | 过滤逻辑优化 + ML基因分型 | 移植xTEA核心过滤器 |
| 第5周 | CRAM支持 + Docker + WDL | 支持批量CRAM运行 |
| 第6周 | 测试 + 性能调优 | 30x样本 < 30分钟 |

---

## 5. 预期性能对比

| 步骤 | xTEA (4小时) | MEGAnE (30分钟) | cTEA (目标) |
|------|----------------|-------------------|--------------|
| BAM/CRAM读取 | ~60分钟 (pysam) | ~5分钟 (C++) | ~3分钟 (优化C++) |
| 候选检测 | ~120分钟 (两步法) | ~15分钟 (k-mer + C++) | ~10分钟 (合并扫描) |
| 过滤 | ~30分钟 (Python) | ~5分钟 (C++ + Python) | ~5分钟 (优化Python) |
| 基因分型 | ~30分钟 (重新扫描) | ~5分钟 (缓存) | ~2分钟 (首次缓存) |
| **总计** | **4小时** | **30分钟** | **20-30分钟** |

---

## 6. GitHub仓库设置

**仓库名**: cTEA  
**初始化**:
```bash
mkdir cTEA
cd cTEA
git init
# 添加此设计文档
git add cTEA_Development_Plan.md
git commit -m "Initial commit: cTEA development plan"
git remote add origin git@github.com:Mirror-fish/cTEA.git
git push -u origin main
```

---

## 7. 风险缓解

| 风险 | 缓解措施 |
|------|----------|
| C++/Python集成复杂 | 使用pybind11（vs ctypes），提供详细文档 |
| CRAM性能不如BAM | 测试不同压缩级别，考虑缓存解压后的BAM |
| xTEA过滤逻辑复杂 | 分阶段移植，先实现核心过滤器，再添加高级功能 |
| 多线程竞争条件 | 使用线程安全的数据结构，参考MEGAnE的ThreadPool.h |

---

## 8. 后续扩展

- **Phase 2**: 添加体细胞突变检测（从xTEA移植`x_somatic_calling.py`）
- **Phase 3**: 支持单细胞数据（10X Chromium）
- **Phase 4**: 整合长读长（Nanopore/PacBio）支持

---

**请审阅此计划。批准后我将：**
1. 创建GitHub仓库 `cTEA`
2. 开始Phase 1的C++核心开发
3. 每周提供进度更新