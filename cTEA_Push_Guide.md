# cTEA SSH Push 指南
**日期**: 2026-04-26  
**目标**: 将最新代码通过SSH推送到GitHub仓库  

---

## 一、检查SSH密钥配置

### 1.1 检查现有SSH密钥
```bash
# 查看现有密钥
ls -la ~/.ssh/

# 应该看到类似：
# id_ed25519       (私钥)
# id_ed25519.pub   (公钥)
# 或者
# id_rsa
# id_rsa.pub
```

### 1.2 如果没有密钥，生成新的
```bash
# 生成ed25519密钥（推荐）
ssh-keygen -t ed25519 -C "your_email@example.com"

# 或者生成RSA密钥
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"

# 按提示操作：
# - 输入密钥保存路径（默认~/.ssh/即可）
# - 输入密码（可选，但推荐设置）
```

### 1.3 添加SSH密钥到ssh-agent
```bash
# 启动ssh-agent
eval "$(ssh-agent -s)"

# 添加私钥
ssh-add ~/.ssh/id_ed25519
# 或者
ssh-add ~/.ssh/id_rsa
```

### 1.4 添加公钥到GitHub
```bash
# 复制公钥内容
cat ~/.ssh/id_ed25519.pub
# 或者
cat ~/.ssh/id_rsa.pub

# 复制输出的全部内容，然后：
# 1. 登录GitHub
# 2. 点击右上角头像 → Settings
# 3. 左侧菜单选择 "SSH and GPG keys"
# 4. 点击 "New SSH key"
# 5. 粘贴公钥内容，设置标题（如：My Laptop）
# 6. 点击 "Add SSH key"
```

### 1.5 测试SSH连接
```bash
# 测试GitHub SSH连接
ssh -T git@github.com

# 成功会显示：
# Hi Mirror-fish! You've successfully authenticated, but GitHub does not provide shell access.
```

---

## 二、检查仓库配置

### 2.1 查看当前仓库状态
```bash
cd /Users/jzeng/Desktop/cTEA

# 查看仓库状态
git status

# 查看远程仓库配置
git remote -v

# 应该显示：
# origin  git@github.com:Mirror-fish/cTEA.git (fetch)
# origin  git@github.com:Mirror-fish/cTEA.git (push)
```

### 2.2 如果远程仓库不是SSH，修改它
```bash
# 如果不是SSH格式，修改为SSH
git remote set-url origin git@github.com:Mirror-fish/cTEA.git

# 再次验证
git remote -v
```

---

## 三、准备提交的新文件

### 3.1 查看修改和新增的文件
```bash
cd /Users/jzeng/Desktop/cTEA

# 查看仓库状态
git status

# 应该看到类似：
# modified:   cpp/Makefile
# modified:   cpp/candidate_detector.cpp
# modified:   cpp/candidate_detector.h
# new file:   cpp/bindings.cpp
# new file:   python/test_pybind11.py
# new file:   cTEA_Test_Guide.md
# new file:   cTEA_Push_Guide.md
```

### 3.2 添加文件到暂存区
```bash
cd /Users/jzeng/Desktop/cTEA

# 添加所有修改和新增文件
git add cpp/bindings.cpp
git add cpp/candidate_detector.h
git add cpp/candidate_detector.cpp
git add cpp/Makefile
git add python/test_pybind11.py
git add cTEA_Test_Guide.md
git add cTEA_Push_Guide.md

# 或者添加所有更改（谨慎使用）
# git add -A

# 验证暂存区
git status
# 应该显示 "Changes to be committed"
```

---

## 四、提交更改

### 4.1 创建提交
```bash
cd /Users/jzeng/Desktop/cTEA

# 提交更改
git commit -m "Phase 1: Implement pybind11 bindings and test infrastructure

- Add cpp/bindings.cpp: pybind11 bindings for C++ core
  - Expose BamProcessor, CandidateDetector classes to Python
  - Add BamHeaderInfo, CandidateEvidence, ScanStatistics structs
  - Implement scan_all_evidence with Python callback
  
- Update cpp/candidate_detector.h:
  - Fix namespace consistency (cTEA namespace)
  - Separate declaration and implementation
  - Add proper struct and class definitions
  
- Update cpp/candidate_detector.cpp:
  - Clean up implementation (remove duplicate class definition)
  - Fix constructor and destructor
  - Implement all methods declared in header
  
- Update cpp/Makefile:
  - Add pybind11 compilation support
  - Add Python framework linking for macOS
  - Create libcTEA_pybind11.so target
  
- Add python/test_pybind11.py:
  - Test script for pybind11 module
  - Verify import, class creation, basic functionality
  
- Add cTEA_Test_Guide.md:
  - Comprehensive testing guide for server cluster
  - Include dependencies, compilation, testing steps
  - Document expected results and troubleshooting"

# 验证提交
git log -1
```

---

## 五、推送到GitHub

### 5.1 推送前先拉取最新代码（避免冲突）
```bash
cd /Users/jzeng/Desktop/cTEA

# 拉取远程最新更改
git pull origin main

# 如果有冲突，解决后再次提交
# git add -A
# git commit -m "Resolve merge conflicts"
```

### 5.2 推送到GitHub
```bash
cd /Users/jzeng/Desktop/cTEA

# 推送更改到main分支
git push origin main

# 如果是第一次推送，可能需要：
# git push -u origin main
```

### 5.3 输入SSH密钥密码
```bash
# 推送时会提示输入SSH密钥密码（如果设置了）
# Enter passphrase for key '/Users/jzeng/.ssh/id_ed25519':

# 输入密码后，应该看到：
# Enumerating objects: XX, done.
# Counting objects: XX, done.
# Compressing objects: XX% (XX/XX), done.
# Writing objects: XX% (XX/XX), XX KiB | XX MiB/s, done.
# Total XX (delta XX), reused XX (delta XX), pack-reused XX
# To github.com:Mirror-fish/cTEA.git
#    xxxxxxx..yyyyyyy  main -> main
```

---

## 六、验证推送成功

### 6.1 查看GitHub仓库
```bash
# 在浏览器打开
open https://github.com/Mirror-fish/cTEA

# 应该看到最新的提交：
# "Phase 1: Implement pybind11 bindings and test infrastructure"
# 以及所有新增和修改的文件
```

### 6.2 本地验证
```bash
cd /Users/jzeng/Desktop/cTEA

# 查看最新提交
git log -1

# 查看文件列表
ls -la cpp/bindings.cpp
ls -la cTEA_Test_Guide.md
ls -la python/test_pybind11.py
```

---

## 七、完整的推送命令汇总

```bash
# ============ 完整流程 ============

# 1. 检查SSH连接
cd /Users/jzeng/Desktop/cTEA
ssh -T git@github.com

# 2. 检查仓库状态
git status
git remote -v

# 3. 添加文件到暂存区
git add cpp/bindings.cpp
git add cpp/candidate_detector.h
git add cpp/candidate_detector.cpp
git add cpp/Makefile
git add python/test_pybind11.py
git add cTEA_Test_Guide.md
git add cTEA_Push_Guide.md

# 4. 提交更改
git commit -m "Phase 1: Implement pybind11 bindings and test infrastructure"

# 5. 拉取最新代码（避免冲突）
git pull origin main

# 6. 推送到GitHub
git push origin main

# ============ 完成 ============
```

---

## 八、常见问题解决

### 问题1: `Permission denied (publickey)`
```bash
# 解决：确保SSH密钥已添加到GitHub
ssh -T git@github.com
# 如果还是失败，重新添加密钥到agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

### 问题2: `fatal: Authentication failed`
```bash
# 解决：检查远程URL是否为SSH格式
git remote -v
# 如果不是，修改：
git remote set-url origin git@github.com:Mirror-fish/cTEA.git
```

### 问题3: `! [rejected] main -> main (fetch first)`
```bash
# 解决：先拉取远程更改
git pull origin main
# 解决可能的冲突后，再次推送
git push origin main
```

### 问题4: 忘记SSH密钥密码
```bash
# 如果忘记密码，需要重新生成密钥
ssh-keygen -t ed25519 -C "your_email@example.com"
# 然后重新添加到GitHub（参考步骤1.4）
```

---

## 九、推送内容清单

✅ **新增文件**：
- `cpp/bindings.cpp` - pybind11绑定代码
- `python/test_pybind11.py` - pybind11测试脚本
- `cTEA_Test_Guide.md` - 测试操作指南
- `cTEA_Push_Guide.md` - 本Push指南

✅ **修改文件**：
- `cpp/candidate_detector.h` - 修复头文件
- `cpp/candidate_detector.cpp` - 清理实现
- `cpp/Makefile` - 添加pybind11支持

✅ **编译产物**（不需要push，已在.gitignore中）：
- `cpp/libcTEA.a`
- `cpp/libcTEA.so`
- `cpp/libcTEA_pybind11.so`
- `cpp/*.o`

---

## 十、推送后下一步

推送成功后，您可以在服务器集群上：

```bash
# 1. 克隆最新代码
git clone git@github.com:Mirror-fish/cTEA.git
cd cTEA

# 2. 编译项目
cd cpp
make clean && make all && make libcTEA_pybind11.so

# 3. 创建符号链接
ln -sf libcTEA_pybind11.so libcTEA.so

# 4. 运行测试
cd ..
python3 test_pybind11.py
python3 test_cTEA_basic.py
python3 test_cTEA_scan.py /path/to/input.cram /path/to/reference.fa 4
```

---

**指南结束** | 如有问题，请参考：https://docs.github.com/en/authentication/connecting-to-github-with-ssh