# cTEA (CRAM-optimized Transposable Element Analyzer)
## 开发计划 - 更新版 (2026-04-24)

**目标**: 将30x WGS样本运行时间从4小时降至20-30分钟，支持CRAM，提供Docker/WDL批量运行
**当前状态**: Phase 1 核心开发（C++代码框架完成，编译受阻）

---

## 当前进度总结

**最后更新**: 2026-04-24 16:15 UTC+2
**当前状态**: Phase 1 核心开发（✅ C++编译成功！）

### ✅ 已完成 (Phase 1)

1. **目录结构创建**
   - `cpp/` - C++核心代码
   - `python/` - Python包装器和过滤器
   - `tests/` - 测试数据和脚本
   - `docker/`, `wdl/` - 部署相关（待实现）

2. **C++核心代码框架** ✅ **编译成功！**
   - ✅ `bam_processor.h/cpp` - BAM/CRAM处理器（使用htslib）
   - ✅ `extract_reads.cpp` - 合并MEGAnE的extract_discordant和extract_unmapped
   - ✅ `kmer_filter.cpp` - k-mer过滤器（从MEGAnE移植，C++优化）
   - ✅ `candidate_detector.cpp` - **关键优化**：单次扫描检测候选（vs xTEA的分染色体方法）
   - ✅ 辅助头文件：`ThreadPool.h`, `dna_to_2bit.hpp`, `complementary_seq.hpp`
   
   **编译成果** (2026-04-24 16:10):
   - ✅ `libcTEA.a` (29K) - 静态库
   - ✅ `libcTEA.so` (60K) - 共享库（可用于Python ctypes/pybind11）
   - ✅ `extract_reads` (102K) - 独立可执行文件
   - ✅ `kmer_filter` (59K) - 独立可执行文件
   - ✅ `candidate_detector` (85K) - 独立可执行文件

3. **Python层**
   - ✅ `ctea_core.py` - C++核心的Python包装器（支持pybind11/ctypes）
   - ✅ `ctea_main.py` - 主流程脚本（模拟MEGAnE的简洁性 + xTEA的准确性）

4. **测试**
   - ✅ `generate_test_data.py` - 生成小体量测试数据（~1MB参考，1000个模拟MEI位点）

5. **差距分析**
   - ✅ `MEGA-xTEA_Comprehensive_Gap_Analysis.md` - 完成对MEGA-xTEA的全面分析
   - **关键发现**: MEGA-xTEA是浅层封装而非深度整合，缺失xTEA约70%功能

### ✅ 已解决的问题（2026-04-24）

1. **✅ htslib编译问题** - **已解决！**
   - **解决方案**: 使用Homebrew安装的htslib (`/opt/homebrew/`)
   - **Makefile更新**: 正确配置 `HTSLIB_DIR` 和 `HTSLIB_INC`
   - **验证**: `pkg-config --cflags htslib` 确认路径

2. **✅ extract_reads.cpp语法错误** - **已修复！**
   - **问题**: 引号、逗号丢失，函数调用格式错误
   - **解决**: 重写文件为干净版本 (v0.2.0)
   - **添加**: `SoftClipInfo` 结构体定义
   - **修复**: `sam_open` API调用（适配新版本htslib）
   - **修复**: `dna_to_2bit` → `dna_to_2bit_table` 引用
   - **修复**: `l_qseq` 重复定义问题

3. **✅ kmer_filter.cpp编译错误** - **已修复！**
   - **添加**: `typedef unsigned long long ull;`
   - **修复**: `dna_to_2bit[seq[i]]` → `dna_to_2bit_table[(unsigned char)seq[i]]`
   - **添加**: `<cstdint>` 头文件

4. **✅ 重复main符号错误** - **已修复！**
   - **问题**: extract_reads.cpp, kmer_filter.cpp, candidate_detector.cpp 都有main函数
   - **解决**: 添加 `#ifndef BUILD_LIB` / `#endif` 宏保护
   - **Makefile更新**: 
     - 库编译使用 `CXXFLAGS_LIB = $(CXXFLAGS) -DBUILD_LIB`
     - 可执行文件编译使用独立的 `*_main.o` 对象文件
   - **结果**: 静态库和共享库不再包含main函数

### ⚠️ 剩余问题

1. **C++和Python集成未验证**
   - pybind11绑定代码未实现
   - ctypes接口未测试
   - 当前`ctea_core.py`是占位符实现

2. **MEGA-xTEA功能缺失**
   - 根据差距分析，需要移植xTEA的:
     - 完整后过滤链（x_post_filter.py有1247行，MEGA-xTEA仅实现30%）
     - 转导检测（x_transduction.py）
     - TSD验证（x_TSD.py）
     - 局部组装（x_local_assembly.py）

3. **编译警告**（不影响功能）
   - `candidate_detector.cpp:110:14: warning: private field 'total_reads_processed_' is not used`
   - `candidate_detector.cpp:111:14: warning: private field 'candidates_found_' is not used`

---

## 更新后的开发步骤

### Phase 1: C++核心开发（第1-2周）[进行中]

#### Step 1.1: 解决htslib依赖 ✅ 代码完成，⚠️ 编译受阻
- **文件**: `cpp/bam_processor.h`, `cpp/bam_processor.cpp`
- **功能**: 使用htslib支持BAM/CRAM
- **状态**: 代码完成，需要修复编译环境
- **下一步**: 
  ```bash
  # 选项1: 使用conda
  conda install -c bioconda htslib
  
  # 选项2: 手动编译安装htslib
  git clone https://github.com/samtools/htslib.git
  cd htslib && autoheader && autoconf && ./configure && make && make install
  ```

#### Step 1.2: 合并read提取器 ✅ 完成
- **文件**: `cpp/extract_reads.cpp`
- **关键优化**: 单次扫描同时提取discordant/unmapped/overhang reads（vs MEGAnE的两步法）

#### Step 1.3: k-mer过滤器 ✅ 完成
- **文件**: `cpp/kmer_filter.cpp`
- **优化**: C++实现 + 二分查找（O(log n) vs MEGAnE Python的O(n)）

#### Step 1.4: 候选检测器 ✅ 完成（框架）
- **文件**: `cpp/candidate_detector.cpp`
- **关键优化**: **单次扫描**检测候选（vs xTEA的分染色体迭代）
- **性能提升**: 预计3-5x加速

#### Step 1.5: Makefile和编译 ⚠️ 受阻
- **文件**: `cpp/Makefile`
- **状态**: 需要修复htslib路径
- **下一步**: 修复编译后，测试单个模块

### Phase 2: Python包装和过滤逻辑（第3-4周）[待开始]

#### Step 2.1: 主入口和管道 ✅ 代码完成
- **文件**: `python/ctea_main.py`
- **功能**: 类似MEGAnE的简洁流程，但直接调用C++（不通过subprocess）

#### Step 2.2: 从MEGA-xTEA移植过滤器 ⚠️ 需要完整移植
- **基于差距分析，需要**:
  1. **完整移植x_post_filter.py**（1247行 → 预计700-800行优化版）
     - SVA过滤（4级层次）
     - Alu过滤
     - L1过滤（含转导处理）
     - 低分歧参考TE过滤
     - 覆盖率异常检测
  2. **实现缺失模块**:
     - `transduction.py` - 转导检测（从x_transduction.py移植）
     - `tsd_validator.py` - TSD验证（从x_TSD.py移植）
     - `local_assembly.py` - 局部组装（从x_local_assembly.py移植）

#### Step 2.3: ML基因分型优化 ⚠️ 特征提取不准确
- **问题**: MEGA-xTEA的`ml_genotype_features.tsv`是近似值
- **解决**: 从C++首次扫描时缓存精确特征，避免BAM重新扫描

### Phase 3: CRAM支持和批量运行（第5周）[待开始]

#### Step 3.1: CRAM完整支持
- **文件**: `python/utils/cram_utils.py`
- **功能**: 检查CRAM依赖，自动传递`--reference`给htslib

#### Step 3.2: Docker镜像
- **文件**: `docker/Dockerfile`
- **基础**: `biocontainers/htslib` + Python3

#### Step 3.3: WDL脚本
- **文件**: `wdl/ctea.wdl`
- **功能**: 支持单样本和批量运行

### Phase 4: 测试和性能调优（第6周）[待开始]

#### Step 4.1: 单元测试
- **文件**: `tests/test_cpp.py`, `tests/test_filters.py`
- **数据**: 使用`tests/test_data/`的模拟数据

#### Step 4.2: 端到端测试
- **测试数据**: 30x WGS CRAM文件（~25GB）
- **对比**: xTEA (4小时) vs MEGAnE (30分钟) vs cTEA (目标20-30分钟)

#### Step 4.3: 性能分析
- **工具**: perf, valgrind
- **目标**: 识别剩余瓶颈

---

## 关键技术决策（更新）

### 1. 避免xTEA的分染色体问题 ✅ 已实现框架
**xTEA的问题** (`x_TEI_locator.py` lines 112-168):
```python
# 问题：每个染色体都要读取整个BAM，然后合并
for chrm in m_chrms:
    for i in range(cnt):  # 每个BAM文件
        if tmp_chrm != chrm:
            continue  # 跳过不匹配的
```

**cTEA解决方案** (`candidate_detector.cpp`):
- C++单次扫描BAM，在内存中按染色体组织数据
- 使用哈希表（`unordered_map`）存储候选位点
- 输出时再按染色体分组（如果需要）

### 2. CRAM支持 ✅ 框架完成，待测试
- 使用htslib的`sam_read1()` API，自动检测BAM/CRAM
- CRAM需要reference FASTA：通过命令行参数传递
- 参考MEGAnE的`cpp/parse_fai.hpp`实现FASTA索引读取

### 3. 整合xTEA的过滤逻辑 ⚠️ 需要完整移植
**保留xTEA的优势**:
- SVA专用过滤（`x_post_filter.py`的`post_processing_SVA`）
- AF冲突检查（4-ratio test）
- 低分歧参考TE拷贝过滤

**用更高效方式实现**:
- 不从BAM重新扫描，而是从首次扫描缓存证据
- 使用numpy进行向量化计算（vs xTEA的纯Python循环）

---

## 预期性能对比（更新）

| 步骤 | xTEA (4小时) | MEGAnE (30分钟) | cTEA (目标) | 当前状态 |
|------|--------------|-------------------|--------------|----------|
| BAM/CRAM读取 | ~60分钟 (pysam) | ~5分钟 (C++) | ~3分钟 (优化C++) | ⚠️ 编译受阻 |
| 候选检测 | ~120分钟 (两步法) | ~15分钟 (k-mer + C++) | ~10分钟 (合并扫描) | ✅ 框架完成 |
| 过滤 | ~30分钟 (Python) | ~5分钟 (C++ + Python) | ~5分钟 (优化Python) | ⚠️ 需完整移植 |
| 基因分型 | ~30分钟 (重新扫描) | ~5分钟 (缓存) | ~2分钟 (首次缓存) | ⚠️ 特征不准确 |
| **总计** | **4小时** | **30分钟** | **20-30分钟** | **开发中** |

---

## 下一步行动清单

### 立即执行（本周）
1. **修复htslib编译问题**
   ```bash
   # 尝试conda安装
   conda install -c bioconda htslib
   # 或手动安装
   cd /tmp && git clone https://github.com/samtools/htslib.git
   cd htslib && autoheader && autoconf && ./configure --prefix=/usr/local && make && sudo make install
   ```

2. **编译测试C++核心**
   ```bash
   cd /Users/jzeng/Desktop/cTEA/cpp
   make clean && make 2>&1 | tee build.log
   ```

3. **创建Python虚拟环境**
   ```bash
   cd /Users/jzeng/Desktop/cTEA
   python3 -m venv venv
   source venv/bin/activate
   pip install pysam numpy scipy scikit-learn
   ```

### 短期（2周内）
4. **实现pybind11绑定**（让Python直接调用C++）
5. **完整移植xTEA过滤器**（基于差距分析）
6. **测试管道**（使用`tests/generate_test_data.py`生成的数据）

### 中期（1个月内）
7. **CRAM完整测试**
8. **Docker镜像构建**
9. **WDL脚本编写**
10. **性能基准测试**（vs xTEA和MEGAnE）

---

## 风险缓解（更新）

| 风险 | 缓解措施 | 状态 |
|------|----------|------|
| C++/Python集成复杂 | 使用pybind11，提供详细文档 | ⚠️ 未开始 |
| CRAM性能不如BAM | 测试不同压缩级别，考虑缓存解压后的BAM | ⚠️ 未测试 |
| xTEA过滤逻辑复杂 | 分阶段移植，先实现核心过滤器 | ⚠️ 需完整移植 |
| 多线程竞争条件 | 使用线程安全的数据结构，参考MEGAnE的ThreadPool.h | ✅ 已实现 |
| **htslib编译失败** | 提供多种安装选项（conda/手动/预编译） | ⚠️ 当前阻塞 |

---

## GitHub仓库状态

**仓库名**: cTEA  
**状态**: 已初始化，代码已提交  
**下一步**: 修复编译问题后，继续开发并定期提交

---

## 附录：MEGA-xTEA差距分析关键发现

详见 `MEGA-xTEA_Comprehensive_Gap_Analysis.md`

**核心问题**:
1. MEGA-xTEA是**浅层封装**（subprocess调用MEGAnE），而非深度整合
2. **缺失xTEA约70%功能**（后过滤链仅实现30%）
3. ML特征提取使用**近似值**，影响基因分型准确性
4. 比MEGAnE还慢（预估>40分钟 vs MEGAnE的30分钟）

**cTEA的改进**:
- ✅ 单次扫描架构（vs xTEA的分染色体）
- ✅ C++核心 + Python包装（vs MEGA-xTEA的subprocess调用）
- ⚠️ 需要完整移植xTEA过滤器
- ⚠️ 需要修复ML特征提取

---

**最后更新**: 2026-04-24 13:58 UTC+2  
**下次审查**: 修复htslib编译问题后