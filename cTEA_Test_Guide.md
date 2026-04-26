# cTEA 测试操作指南
**版本**: 0.1.0  
**日期**: 2026-04-26  
**目标**: 在服务器集群上测试cTEA的BAM/CRAM处理功能

---

## 一、依赖条件

### 1.1 系统依赖（必需）
| 依赖 | 版本要求 | 安装命令（Ubuntu/Debian） | 验证命令 |
|------|----------|------------------------|----------|
| **C++编译器** | g++ ≥ 7.0 | `sudo apt-get install g++` | `g++ --version` |
| **make** | ≥ 4.0 | `sudo apt-get install make` | `make --version` |
| **htslib** | ≥ 1.9 | `sudo apt-get install libhts-dev` 或 从源码编译 | `pkg-config --modversion htslib` |
| **zlib** | ≥ 1.2 | `sudo apt-get install zlib1g-dev` | - |
| **bzip2** | ≥ 1.0 | `sudo apt-get install libbz2-dev` | - |
| **lzma** | ≥ 5.2 | `sudo apt-get install liblzma-dev` | - |
| **Python** | ≥ 3.7 | `sudo apt-get install python3 python3-pip` | `python3 --version` |

### 1.2 Python依赖（必需）
```bash
# 创建虚拟环境（推荐）
python3 -m venv cTEA_venv
source cTEA_venv/bin/activate

# 安装依赖
pip install pybind11 numpy scipy scikit-learn

# 验证安装
python3 -c "import pybind11; print('pybind11 version:', pybind11.__version__)"
```

### 1.3 可选依赖（用于生成测试数据）
| 依赖 | 用途 | 安装命令 |
|------|------|----------|
| **ART** | 模拟RNA-seq/基因组测序数据 | `pip install ART` 或从源码编译 |
| **dwgsim** | 模拟基因组测序数据 | 从https://github.com/dhully/dwgsim 编译 |

---

## 二、项目编译

### 2.1 获取源代码
```bash
# 克隆仓库
git clone https://github.com/Mirror-fish/cTEA.git
cd cTEA

# 或使用SSH（如果已配置）
git clone git@github.com:Mirror-fish/cTEA.git
cd cTEA
```

### 2.2 编译C++核心
```bash
cd cpp

# 清理旧编译
make clean

# 编译静态库、共享库和可执行文件
make all 2>&1 | tee build.log

# 编译pybind11模块（Python绑定）
# 注意：需要确保pybind11已安装
make libcTEA_pybind11.so 2>&1 | tee pybind11_build.log

# 验证编译结果
ls -lh libcTEA.* libcTEA_pybind11.so
# 应该看到：
# libcTEA.a        (静态库)
# libcTEA.so       (共享库)
# libcTEA_pybind11.so (pybind11模块)
```

### 2.3 创建Python模块链接
```bash
cd cpp

# 创建符号链接，让Python可以导入libcTEA
ln -sf libcTEA_pybind11.so libcTEA.so

# 验证
ls -la libcTEA.so
```

---

## 三、测试数据准备

### 3.1 测试数据要求
- **BAM或CRAM文件**：包含基因组测序数据
- **参考FASTA文件**：对应参考基因组（CRAM必需）
- **推荐测试数据**：
  - 30x WGS（全基因组测序）数据
  - 染色体21（chr21）数据（较小，便于测试）
  - 已知包含TEI位点的样本（如NA19129）

### 3.2 数据示例路径
假设您的测试数据位于：
```
/data/projects/test/
├── NA19129_chr21_unmapped.cram    # 测试CRAM文件
├── genome.fa                              # 参考基因组
└── genome.fa.fai                          # FASTA索引（samtools faidx生成）
```

---

## 四、测试步骤

### 4.1 基础功能测试（Python层）

创建测试脚本 `test_cTEA_basic.py`：
```python
#!/usr/bin/env python3
"""
cTEA基础功能测试
测试BAM/CRAM处理、候选检测等核心功能
"""

import sys
import os

# 添加cpp目录到Python路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'cpp'))

try:
    import libcTEA
    print("✅ 成功导入libcTEA模块")
    print(f"   版本: {libcTEA.__version__}")
except ImportError as e:
    print(f"❌ 导入libcTEA失败: {e}")
    sys.exit(1)

# 测试1: 创建BamProcessor对象
print("\n=== 测试1: BamProcessor类 ===")
try:
    bp = libcTEA.BamProcessor()
    print("✅ 创建BamProcessor对象成功")
except Exception as e:
    print(f"❌ 创建BamProcessor失败: {e}")

# 测试2: 创建CandidateDetector对象
print("\n=== 测试2: CandidateDetector类 ===")
try:
    detector = libcTEA.CandidateDetector()
    print("✅ 创建CandidateDetector对象成功")
except Exception as e:
    print(f"❌ 创建CandidateDetector失败: {e}")

print("\n✅ 基础功能测试完成！")
```

运行测试：
```bash
cd /path/to/cTEA
source cTEA_venv/bin/activate  # 如果在虚拟环境中
python3 test_cTEA_basic.py
```

### 4.2 BAM/CRAM扫描测试

创建测试脚本 `test_cTEA_scan.py`：
```python
#!/usr/bin/env python3
"""
cTEA BAM/CRAM扫描测试
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'cpp'))

import libcTEA

def test_bam_scan(bam_path, ref_path="", n_threads=4):
    """测试BAM/CRAM扫描功能"""
    print(f"\n=== 测试BAM/CRAM扫描 ===")
    print(f"输入文件: {bam_path}")
    print(f"参考文件: {ref_path if ref_path else '无（BAM不需要）'}")
    print(f"线程数: {n_threads}")
    
    try:
        # 创建BamProcessor
        bp = libcTEA.BamProcessor()
        
        # 初始化
        print("\n1. 初始化BamProcessor...")
        result = bp.initialize(bam_path, ref_path, n_threads)
        if not result:
            print("❌ 初始化失败")
            return False
        print("✅ 初始化成功")
        
        # 检查是否为CRAM
        is_cram = bp.is_cram()
        print(f"   文件类型: {'CRAM' if is_cram else 'BAM'}")
        
        # 获取头信息
        print("\n2. 获取文件头信息...")
        header = bp.get_header_info()
        print(f"   染色体数量: {header.n_chromosomes}")
        print(f"   染色体列表: {header.chromosome_names[:5]}...")  # 只显示前5个
        
        # 扫描证据（带进度回调）
        print("\n3. 扫描证据...")
        def progress_callback(current, total):
            if total == 0:
                print(f"\r   已处理 {current} 条reads...", end="", flush=True)
            else:
                percent = int(current * 100 / total)
                print(f"\r   进度: {current}/{total} ({percent}%)...", end="", flush=True)
        
        evidence = bp.scan_all_evidence(progress_callback)
        print("\n✅ 扫描完成")
        
        # 统计证据
        total_sites = sum(len(pos_map) for pos_map in evidence.values())
        print(f"   收集到 {total_sites} 个候选位点的证据")
        
        # 获取统计信息
        print("\n4. 获取统计信息...")
        stats = bp.get_statistics()
        print(f"   总reads数: {stats.total_reads}")
        print(f"   已映射reads: {stats.mapped_reads}")
        print(f"   软剪切reads: {stats.clipped_reads}")
        print(f"   异常配对reads: {stats.discordant_pairs}")
        print(f"   未映射reads: {stats.unmapped_reads}")
        print(f"   平均读长: {stats.avg_read_length:.1f}")
        print(f"   估计覆盖度: {stats.estimated_coverage:.1f}x")
        
        return True
        
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python3 test_cTEA_scan.py <input.bam/cram> [reference.fa] [n_threads]")
        print("示例: python3 test_cTEA_scan.py /data/NA19129_chr21.cram /data/genome.fa 4")
        sys.exit(1)
    
    bam_path = sys.argv[1]
    ref_path = sys.argv[2] if len(sys.argv) > 2 else ""
    n_threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    
    if not os.path.isfile(bam_path):
        print(f"❌ 输入文件不存在: {bam_path}")
        sys.exit(1)
    
    if ref_path and not os.path.isfile(ref_path):
        print(f"❌ 参考文件不存在: {ref_path}")
        sys.exit(1)
    
    test_bam_scan(bam_path, ref_path, n_threads)
```

运行扫描测试：
```bash
cd /path/to/cTEA

# 测试BAM文件（不需要参考基因组）
python3 test_cTEA_scan.py /data/NA19129_chr21.bam 4

# 测试CRAM文件（需要参考基因组）
python3 test_cTEA_scan.py /data/NA19129_chr21.cram /data/genome.fa 4
```

### 4.3 候选检测测试

创建测试脚本 `test_cTEA_detect.py`：
```python
#!/usr/bin/env python3
"""
cTEA候选检测测试
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'cpp'))

import libcTEA

def test_candidate_detection(bam_path, output_bed, ref_path="", n_threads=4):
    """测试候选检测功能"""
    print(f"\n=== 测试候选检测 ===")
    print(f"输入文件: {bam_path}")
    print(f"输出BED: {output_bed}")
    print(f"参考文件: {ref_path if ref_path else '无'}")
    
    try:
        # 创建CandidateDetector
        detector = libcTEA.CandidateDetector()
        
        # 初始化
        print("\n1. 初始化CandidateDetector...")
        result = detector.initialize(bam_path, ref_path, n_threads)
        if not result:
            print("❌ 初始化失败")
            return False
        print("✅ 初始化成功")
        
        # 执行检测
        print("\n2. 执行候选检测...")
        # 注意：目前detect_candidates是占位符，可能返回False
        result = detector.detect_candidates(output_bed, "", "")
        if result:
            print("✅ 检测成功")
            
            # 获取统计
            stats = detector.get_statistics()
            print(f"   总reads数: {stats['total_reads']}")
            print(f"   候选位点: {stats['candidates_found']}")
            
            # 检查输出文件
            if os.path.isfile(output_bed):
                with open(output_bed) as f:
                    lines = f.readlines()
                    print(f"   输出文件: {len(lines)} 行（含表头）")
        else:
            print("⚠️  检测返回False（可能功能未完成）")
        
        return True
        
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python3 test_cTEA_detect.py <input.bam/cram> <output.bed> [reference.fa] [n_threads]")
        sys.exit(1)
    
    bam_path = sys.argv[1]
    output_bed = sys.argv[2]
    ref_path = sys.argv[3] if len(sys.argv) > 3 else ""
    n_threads = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    
    test_candidate_detection(bam_path, output_bed, ref_path, n_threads)
```

运行检测测试：
```bash
cd /path/to/cTEA

# 测试候选检测
python3 test_cTEA_detect.py \
    /data/NA19129_chr21.cram \
    /data/output/candidates.bed \
    /data/genome.fa 4
```

---

## 五、预期测试结果

### 5.1 成功指标
| 测试项 | 预期结果 | 备注 |
|--------|----------|------|
| 模块导入 | ✅ 成功导入libcTEA | 版本应为0.1.0 |
| BamProcessor创建 | ✅ 对象创建成功 | - |
| CandidateDetector创建 | ✅ 对象创建成功 | - |
| BAM初始化 | ✅ 返回True | - |
| CRAM初始化 | ✅ 返回True | 需要参考FASTA |
| 扫描证据 | ✅ 返回证据字典 | 可能大量位点 |
| 候选检测 | ⚠️  可能返回False | 目前是占位符实现 |

### 5.2 当前限制
1. ⚠️ **CandidateDetector.detect_candidates** 目前是**占位符实现**，可能不产生实际输出
2. ⚠️ **过滤链** 尚未移植（post-filter、转导检测等）
3. ⚠️ **Python包装器** `ctea_core.py` 还是占位符，未使用pybind11

### 5.3 可测试的功能
✅ **可以测试**：
- BAM/CRAM文件读取
- 文件头信息获取
- 证据扫描收集（scan_all_evidence）
- 扫描统计信息获取
- CRAM格式支持验证

❌ **暂不可测试**：
- 完整候选检测流程
- 后过滤链
- 基因分型
- 最终BED/VCF输出

---

## 六、故障排除

### 6.1 常见问题

**问题1**: `ImportError: No module named 'libcTEA'`
```bash
# 解决：检查符号链接
cd /path/to/cTEA/cpp
ls -la libcTEA.so
# 应该指向libcTEA_pybind11.so

# 重新创建链接
ln -sf libcTEA_pybind11.so libcTEA.so
```

**问题2**: `Undefined symbols` 或 `symbol(s) not found`
```bash
# 解决：检查htslib安装
pkg-config --modversion htslib

# 如果没有，安装htslib
sudo apt-get install libhts-dev
# 或从源码编译
```

**问题3**: `CRAM读取失败`
```bash
# 确保提供参考FASTA文件
# CRAM需要参考基因组进行解压缩
python3 test_cTEA_scan.py input.cram reference.fa 4
```

**问题4**: `Python.h` 找不到（编译时）
```bash
# 安装Python开发头文件
sudo apt-get install python3-dev
```

---

## 七、性能基准测试（可选）

### 7.1 运行时间对比
```bash
# 记录开始时间
START=$(date +%s)

# 运行cTEA扫描
python3 test_cTEA_scan.py /data/30x_sample.cram /data/genome.fa 8

# 记录结束时间
END=$(date +%s)
echo "总运行时间: $((END-START)) 秒"
```

### 7.2 内存使用监控
```bash
# 使用time命令监控资源使用
/usr/bin/time -v python3 test_cTEA_scan.py input.cram output.bed ref.fa 4 2>&1 | grep -E "(Maximum resident set size|Elapsed time)"
```

---

## 八、测试完成检查清单

- [ ] 依赖软件已安装（g++, make, htslib, Python）
- [ ] pybind11已安装（pip install pybind11）
- [ ] 项目已成功编译（make all && make libcTEA_pybind11.so）
- [ ] 符号链接已创建（libcTEA.so -> libcTEA_pybind11.so）
- [ ] BAM扫描测试通过
- [ ] CRAM扫描测试通过（如有CRAM文件）
- [ ] 候选检测测试运行（即使返回False，检查错误信息）
- [ ] 性能数据已记录（运行时间、内存使用）

---

## 九、报告测试结果

测试完成后，请提供以下信息：

1. **系统信息**：
   ```bash
   uname -a
   g++ --version
   python3 --version
   pkg-config --modversion htslib
   ```

2. **编译日志**：`cpp/build.log` 和 `cpp/pybind11_build.log`

3. **测试输出**：运行测试脚本的完整输出

4. **问题报告**：任何错误信息或异常行为

---

## 附录：快速测试命令汇总

```bash
# 1. 编译
cd /path/to/cTEA/cpp
make clean && make all && make libcTEA_pybind11.so

# 2. 创建链接
cd /path/to/cTEA/cpp
ln -sf libcTEA_pybind11.so libcTEA.so

# 3. 基础测试
cd /path/to/cTEA
python3 test_cTEA_basic.py

# 4. 扫描测试（BAM）
python3 test_cTEA_scan.py /data/test.bam 4

# 5. 扫描测试（CRAM）
python3 test_cTEA_scan.py /data/test.cram /data/genome.fa 4

# 6. 检测测试
python3 test_cTEA_detect.py /data/test.cram /data/output.bed /data/genome.fa 4
```

---

**文档结束** | 如有问题，请参考仓库：https://github.com/Mirror-fish/cTEA