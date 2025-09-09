# RadonPy 中文使用文档

## 目录
1. [简介](#简介)
2. [安装指南](#安装指南)
3. [快速开始](#快速开始)
4. [核心功能](#核心功能)
5. [详细教程](#详细教程)
6. [API参考](#api参考)
7. [实例代码](#实例代码)
8. [常见问题](#常见问题)
9. [参考文献](#参考文献)

---

## 简介

RadonPy 是第一个用于聚合物材料物性全自动计算的开源 Python 库，使用全原子经典分子动力学（MD）模拟。对于给定的聚合物重复单元化学结构，整个 MD 模拟过程可以完全自动执行，包括分子建模、平衡态和非平衡态 MD 模拟、平衡状态的自动判定、失败时的自动重启调度以及后处理中的物性计算。

### 主要特点

- 🚀 **全自动化流程**：从 SMILES 输入到物性输出的一站式解决方案
- 📊 **15种物性计算**：密度、热容、体积模量、热导率等
- 🔬 **多种聚合物类型**：均聚物、共聚物（交替、无规、嵌段）
- 🛠 **开源软件栈**：基于 LAMMPS、Psi4、RDKit 等开源工具
- 📈 **数据科学工具**：力场描述符、聚合物信息学工具

### 系统架构

```
输入 (SMILES) → 预处理 → MD模拟 → 后处理 → 物性输出
                  ↓         ↓         ↓
              构象搜索   LAMMPS    性质计算
              力场分配    Psi4     数据分析
              链生成
```

---

## 安装指南

### 系统要求

- Python 3.7-3.12
- LAMMPS >= 3Mar20
- RDKit >= 2020.03
- Psi4 >= 1.5
- 8GB+ RAM（推荐 16GB）
- Linux/macOS（Windows 需要 WSL）

### 方法一：使用 Conda（推荐）

#### Psi4 >= 1.8 版本
```bash
# 1. 创建 conda 环境
conda create -n radonpy python=3.11
conda activate radonpy

# 2. 安装依赖包
conda install -c conda-forge/label/libint_dev -c conda-forge -c psi4 \
    rdkit psi4 resp mdtraj matplotlib scipy pandas

# 3. 安装 LAMMPS
conda install -c conda-forge lammps

# 4. 安装 RadonPy
pip install radonpy-pypi
```

#### Psi4 <= 1.7 版本
```bash
# 1. 创建 conda 环境
conda create -n radonpy python=3.9
conda activate radonpy

# 2. 安装依赖包
conda install -c psi4 -c conda-forge \
    rdkit psi4 resp mdtraj matplotlib scipy pandas

# 3. 安装 LAMMPS
conda install -c conda-forge lammps

# 4. 安装 RadonPy
pip install radonpy-pypi
```

### 方法二：使用 pip

```bash
# 最小安装（不含 MD 模拟功能）
pip install radonpy-pypi

# 包含 LAMMPS（可进行 MD 模拟，但无 DFT 计算）
pip install radonpy-pypi[lammps]
```

### 方法三：从源码安装

```bash
# 克隆仓库
git clone https://github.com/RadonPy/RadonPy.git
cd RadonPy

# 安装
pip install -e .
```

### 环境变量设置

如果手动编译 LAMMPS：
```bash
export LAMMPS_EXEC=/path/to/lammps/binary
```

---

## 快速开始

### 最简单的例子：计算聚苯乙烯的性质

```python
from radonpy.core import utils, poly
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.sim.preset import eq

# 1. 定义聚苯乙烯的重复单元
smiles = '*CC(*)c1ccccc1'  # 聚苯乙烯 SMILES
mol = utils.mol_from_smiles(smiles)

# 2. 生成聚合物链（约1000个原子）
ter = utils.mol_from_smiles('*C')  # 端基
n = poly.calc_n_from_num_atoms(mol, 1000, terminal1=ter)
polymer = poly.polymerize_rw(mol, n)
polymer = poly.terminate_rw(polymer, ter)

# 3. 分配力场参数
ff = GAFF2_mod()
ff.ff_assign(polymer)

# 4. 创建模拟单元
cell = poly.amorphous_cell(polymer, 10, density=0.05)

# 5. 运行平衡态MD模拟
eqmd = eq.EQ21step(cell, work_dir='./PS_simulation')
cell = eqmd.exec(temp=300, press=1.0)

# 6. 分析结果
analy = eqmd.analyze()
properties = analy.get_all_prop(temp=300, press=1.0, save=True)

# 7. 打印结果
print(f"密度: {properties['density']:.3f} g/cm³")
print(f"比热容: {properties['Cp']:.1f} J/kg·K")
print(f"体积模量: {properties['bulk_modulus']/1e9:.2f} GPa")
```

---

## 核心功能

### 1. 聚合物链生成（radonpy.core.poly）

#### 均聚物
```python
# 生成聚丙烯
smiles = '*CC(*)(C)'
mol = utils.mol_from_smiles(smiles)
polymer = poly.polymerize_rw(mol, n=50, tacticity='isotactic')
```

#### 共聚物
```python
# 交替共聚物
mol1 = utils.mol_from_smiles('*CC(*)')
mol2 = utils.mol_from_smiles('*CC(*)C')
copolymer = poly.copolymerize_rw([mol1, mol2], n=50)

# 无规共聚物（50:50比例）
random_copoly = poly.random_copolymerize_rw([mol1, mol2], n=50, ratio=[0.5, 0.5])

# 嵌段共聚物
block_copoly = poly.block_copolymerize_rw([mol1, mol2], [25, 25])
```

### 2. 构象搜索（radonpy.sim.qm）

自动搜索聚合物重复单元的稳定构象：

```python
from radonpy.sim import qm

# 构象搜索流程：
# 1. ETKDG生成1000个初始构象
# 2. GAFF2_mod力场优化
# 3. Butina聚类
# 4. DFT优化前4个构象
mol, energy = qm.conformation_search(
    mol, 
    ff=ff,
    nconf=1000,        # 生成构象数
    dft_nconf=4,       # DFT优化构象数
    opt_method='wb97m-d3bj',
    opt_basis='6-31G(d,p)'
)
```

### 3. 电子性质计算

```python
# RESP电荷
qm.assign_charges(mol, charge='RESP')

# HOMO/LUMO和偶极矩
properties = qm.sp_prop(mol)
print(f"HOMO: {properties['qm_homo']:.2f} eV")
print(f"LUMO: {properties['qm_lumo']:.2f} eV")

# 极化率
polar = qm.polarizability(mol)
print(f"极化率: {polar['qm_polarizability']:.2f} Å³")
```

### 4. 力场分配（radonpy.ff）

支持多种力场：

```python
# GAFF2原始参数
from radonpy.ff.gaff2 import GAFF2
ff = GAFF2()

# GAFF2修正参数（针对含氟聚合物优化）
from radonpy.ff.gaff2_mod import GAFF2_mod
ff = GAFF2_mod()

# GAFF原始参数
from radonpy.ff.gaff import GAFF
ff = GAFF()

# 分配力场
success = ff.ff_assign(polymer)
```

### 5. 模拟单元生成

```python
# 无定形单组分
cell = poly.amorphous_cell(polymer, n=10, density=0.05)

# 无定形混合物
cell = poly.amorphous_mixture_cell(
    [polymer1, polymer2], 
    [5, 5],  # 各5条链
    density=0.1
)

# 取向结构（向列相）
cell = poly.nematic_cell(polymer, n=10, director=[1,0,0])

# 结晶结构
cell = poly.crystal_cell(polymer, a=10, b=10, c=20)
```

### 6. MD模拟预设（radonpy.sim.preset）

#### 平衡态模拟
```python
from radonpy.sim.preset import eq

# 21步压缩/解压缩协议
eqmd = eq.EQ21step(cell, work_dir='./equilibration')
cell = eqmd.exec(temp=300, press=1.0, mpi=4, omp=2)

# 额外平衡（5ns NPT）
eqmd = eq.Additional(cell, work_dir='./additional_eq')
cell = eqmd.exec(temp=300, press=1.0)

# 检查平衡状态
analy = eqmd.analyze()
is_equilibrated = analy.check_eq()
```

#### 非平衡态模拟（热导率）
```python
from radonpy.sim.preset import tc

# Müller-Plathe方法
nemd = tc.NEMD_MP(cell, work_dir='./thermal', axis='x')
cell = nemd.exec(temp=300, decomp=True, mpi=10)

# 分析热导率
analy = nemd.analyze()
thermal_cond = analy.calc_tc(decomp=True, save=True)
print(f"热导率: {thermal_cond:.3f} W/m·K")

# 分解分析
decomp = analy.TCdecomp_data
print(f"键贡献: {decomp['bond']:.3f}")
print(f"角贡献: {decomp['angle']:.3f}")
```

---

## 详细教程

### 教程1：完整的自动化工作流

```python
import os
from radonpy.core import utils, poly
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.sim import qm
from radonpy.sim.preset import eq, tc

def calculate_polymer_properties(smiles, name, temp=300, press=1.0):
    """
    完整的聚合物性质计算流程
    
    参数:
        smiles: 聚合物重复单元的SMILES字符串
        name: 聚合物名称（用于文件夹命名）
        temp: 温度(K)
        press: 压力(atm)
    """
    work_dir = f'./{name}_calculation'
    os.makedirs(work_dir, exist_ok=True)
    
    # 初始化力场
    ff = GAFF2_mod()
    
    # 步骤1：构象搜索
    print("步骤1：构象搜索...")
    mol = utils.mol_from_smiles(smiles)
    mol, energy = qm.conformation_search(
        mol, ff=ff, 
        work_dir=work_dir,
        log_name=f'{name}_conf'
    )
    
    # 步骤2：电子性质计算
    print("步骤2：计算电子性质...")
    qm.assign_charges(mol, charge='RESP', work_dir=work_dir)
    elec_props = qm.sp_prop(mol, work_dir=work_dir)
    polar = qm.polarizability(mol, work_dir=work_dir)
    
    # 步骤3：生成聚合物链
    print("步骤3：生成聚合物链...")
    ter = utils.mol_from_smiles('*C')
    qm.assign_charges(ter, charge='RESP', work_dir=work_dir)
    
    n = poly.calc_n_from_num_atoms(mol, 1000, terminal1=ter)
    polymer = poly.polymerize_rw(mol, n, tacticity='atactic')
    polymer = poly.terminate_rw(polymer, ter)
    
    # 步骤4：力场分配
    print("步骤4：分配力场参数...")
    if not ff.ff_assign(polymer):
        raise ValueError("力场分配失败")
    
    # 步骤5：创建模拟单元
    print("步骤5：创建模拟单元...")
    cell = poly.amorphous_cell(polymer, 10, density=0.05)
    
    # 步骤6：平衡态MD
    print("步骤6：运行平衡态MD...")
    eqmd = eq.EQ21step(cell, work_dir=work_dir)
    cell = eqmd.exec(temp=temp, press=press, mpi=4)
    
    # 检查平衡
    analy = eqmd.analyze()
    props = analy.get_all_prop(temp=temp, press=press, save=True)
    is_eq = analy.check_eq()
    
    # 额外平衡（如需要）
    max_additional = 4
    for i in range(max_additional):
        if is_eq:
            break
        print(f"  额外平衡 {i+1}/{max_additional}...")
        eqmd = eq.Additional(cell, work_dir=work_dir)
        cell = eqmd.exec(temp=temp, press=press, mpi=4)
        analy = eqmd.analyze()
        props = analy.get_all_prop(temp=temp, press=press, save=True)
        is_eq = analy.check_eq()
    
    if not is_eq:
        print("警告：未达到完全平衡状态")
    
    # 步骤7：非平衡态MD（热导率）
    print("步骤7：计算热导率...")
    nemd = tc.NEMD_MP(cell, work_dir=work_dir)
    cell = nemd.exec(temp=temp, decomp=True, mpi=10)
    nemd_analy = nemd.analyze()
    tc_value = nemd_analy.calc_tc(decomp=True, save=True)
    
    # 汇总结果
    results = {
        '聚合物': name,
        'SMILES': smiles,
        '温度': temp,
        '压力': press,
        **props,
        '热导率': tc_value,
        'HOMO': elec_props['qm_homo'],
        'LUMO': elec_props['qm_lumo'],
        '带隙': elec_props['qm_lumo'] - elec_props['qm_homo'],
        '极化率': polar['qm_polarizability']
    }
    
    return results

# 使用示例
if __name__ == '__main__':
    # 计算聚乙烯的性质
    pe_results = calculate_polymer_properties(
        smiles='*CC(*)',
        name='polyethylene'
    )
    
    # 打印结果
    for key, value in pe_results.items():
        if isinstance(value, float):
            print(f"{key}: {value:.3f}")
        else:
            print(f"{key}: {value}")
```

### 教程2：批量计算多种聚合物

```python
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# 聚合物数据库
polymer_database = [
    {'name': 'PE', 'smiles': '*CC(*)', 'cn_name': '聚乙烯'},
    {'name': 'PP', 'smiles': '*CC(*)(C)', 'cn_name': '聚丙烯'},
    {'name': 'PS', 'smiles': '*CC(*)c1ccccc1', 'cn_name': '聚苯乙烯'},
    {'name': 'PMMA', 'smiles': '*CC(*)(C)C(=O)OC', 'cn_name': '聚甲基丙烯酸甲酯'},
    {'name': 'PVC', 'smiles': '*CC(*)Cl', 'cn_name': '聚氯乙烯'},
    {'name': 'PTFE', 'smiles': '*C(F)C(*)F', 'cn_name': '聚四氟乙烯'},
]

def batch_calculate(polymer_list, max_workers=2):
    """批量计算聚合物性质"""
    results = []
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for polymer in polymer_list:
            future = executor.submit(
                calculate_polymer_properties,
                polymer['smiles'],
                polymer['name']
            )
            futures.append((polymer, future))
        
        for polymer, future in futures:
            try:
                result = future.result(timeout=3600)  # 1小时超时
                result['中文名'] = polymer['cn_name']
                results.append(result)
                print(f"✓ {polymer['cn_name']} 计算完成")
            except Exception as e:
                print(f"✗ {polymer['cn_name']} 计算失败: {e}")
    
    # 保存结果
    df = pd.DataFrame(results)
    df.to_csv('polymer_properties.csv', index=False)
    df.to_excel('polymer_properties.xlsx', index=False)
    
    return df

# 运行批量计算
results_df = batch_calculate(polymer_database)
print(results_df)
```

### 教程3：自定义MD模拟

```python
from radonpy.sim import md, lammps

class CustomMD:
    """自定义MD模拟类"""
    
    def __init__(self, mol, work_dir='./custom_md'):
        self.mol = mol
        self.work_dir = work_dir
        self.lmp = lammps.LAMMPS(mol, work_dir=work_dir)
    
    def run_nvt(self, temp, steps, dt=1.0):
        """运行NVT模拟"""
        # 创建LAMMPS输入文件
        self.lmp.make_dat()
        
        # 设置NVT参数
        md_params = md.MdInput()
        md_params.temp = temp
        md_params.dt = dt
        md_params.nstep = steps
        md_params.ensemble = 'nvt'
        md_params.dump_freq = 1000
        md_params.thermo_freq = 100
        
        # 生成输入脚本
        self.lmp.make_input(md_params)
        
        # 运行模拟
        self.lmp.run(mpi=4, omp=2)
        
        # 分析结果
        return self.analyze_trajectory()
    
    def run_npt(self, temp, press, steps, dt=1.0):
        """运行NPT模拟"""
        md_params = md.MdInput()
        md_params.temp = temp
        md_params.press = press
        md_params.dt = dt
        md_params.nstep = steps
        md_params.ensemble = 'npt'
        md_params.barostat = 'iso'
        
        self.lmp.make_input(md_params)
        self.lmp.run(mpi=4, omp=2)
        
        return self.analyze_trajectory()
    
    def run_annealing(self, temp_start, temp_end, steps):
        """退火模拟"""
        md_params = md.MdInput()
        md_params.temp = temp_start
        md_params.temp_end = temp_end
        md_params.nstep = steps
        md_params.ensemble = 'nvt'
        md_params.temp_damp = 100
        
        self.lmp.make_input(md_params)
        self.lmp.run(mpi=4, omp=2)
        
        return self.analyze_trajectory()
    
    def analyze_trajectory(self):
        """分析轨迹"""
        import mdtraj
        
        # 读取轨迹
        traj = mdtraj.load(f'{self.work_dir}/traj.dcd', 
                          top=f'{self.work_dir}/data.pdb')
        
        # 计算性质
        rg = mdtraj.compute_rg(traj)
        
        results = {
            'n_frames': len(traj),
            'rg_mean': rg.mean(),
            'rg_std': rg.std(),
            'volume': traj.unitcell_volumes.mean()
        }
        
        return results

# 使用自定义MD
custom = CustomMD(cell, work_dir='./my_simulation')

# 退火
annealing_results = custom.run_annealing(600, 300, 100000)

# NVT平衡
nvt_results = custom.run_nvt(300, 500000)

# NPT产生
npt_results = custom.run_npt(300, 1.0, 1000000)
```

### 教程4：力场描述符计算

```python
from radonpy.ff.descriptor import ForceFieldDescriptor

def calculate_ff_descriptors(polymer):
    """计算力场描述符用于机器学习"""
    
    # 创建描述符计算器
    ffd = ForceFieldDescriptor()
    
    # 计算各种描述符
    descriptors = {}
    
    # 1. 原子类型统计
    atom_types = ffd.get_atom_types(polymer)
    descriptors['n_atom_types'] = len(atom_types)
    
    # 2. 键参数统计
    bond_params = ffd.get_bond_parameters(polymer)
    descriptors['mean_bond_k'] = bond_params['k'].mean()
    descriptors['std_bond_k'] = bond_params['k'].std()
    descriptors['mean_bond_r0'] = bond_params['r0'].mean()
    
    # 3. 角参数统计
    angle_params = ffd.get_angle_parameters(polymer)
    descriptors['mean_angle_k'] = angle_params['k'].mean()
    descriptors['mean_angle_theta0'] = angle_params['theta0'].mean()
    
    # 4. 二面角参数统计
    dihedral_params = ffd.get_dihedral_parameters(polymer)
    descriptors['n_dihedrals'] = len(dihedral_params)
    
    # 5. 非键相互作用
    vdw_params = ffd.get_vdw_parameters(polymer)
    descriptors['mean_epsilon'] = vdw_params['epsilon'].mean()
    descriptors['mean_sigma'] = vdw_params['sigma'].mean()
    
    # 6. 电荷分布
    charges = ffd.get_charges(polymer)
    descriptors['total_charge'] = charges.sum()
    descriptors['charge_std'] = charges.std()
    descriptors['max_charge'] = charges.max()
    descriptors['min_charge'] = charges.min()
    
    return descriptors

# 计算描述符
ff_descriptors = calculate_ff_descriptors(polymer)

# 用于机器学习
import numpy as np
from sklearn.ensemble import RandomForestRegressor

# 准备数据（示例）
X = []  # 描述符矩阵
y = []  # 目标性质

for polymer in polymer_list:
    desc = calculate_ff_descriptors(polymer)
    X.append(list(desc.values()))
    y.append(polymer_properties)  # 如玻璃化温度

X = np.array(X)
y = np.array(y)

# 训练模型
model = RandomForestRegressor(n_estimators=100)
model.fit(X, y)

# 预测新聚合物性质
new_polymer_desc = calculate_ff_descriptors(new_polymer)
prediction = model.predict([list(new_polymer_desc.values())])
```

---

## API参考

### radonpy.core.utils

| 函数 | 描述 | 参数 | 返回值 |
|-----|------|------|--------|
| `mol_from_smiles(smiles)` | 从SMILES创建分子 | smiles: str | RDKit Mol |
| `mol_to_pdb(mol, filename)` | 保存为PDB文件 | mol: Mol, filename: str | None |
| `mol_to_pkl(mol, filename)` | 保存为pickle文件 | mol: Mol, filename: str | None |
| `mol_from_pkl(filename)` | 从pickle读取 | filename: str | RDKit Mol |

### radonpy.core.poly

| 函数 | 描述 | 参数 | 返回值 |
|-----|------|------|--------|
| `polymerize_rw(mol, n, **kwargs)` | 生成均聚物 | mol: Mol, n: int | Mol |
| `copolymerize_rw(mols, n)` | 生成交替共聚物 | mols: list, n: int | Mol |
| `random_copolymerize_rw(mols, n, ratio)` | 生成无规共聚物 | mols: list, n: int, ratio: list | Mol |
| `block_copolymerize_rw(mols, n_list)` | 生成嵌段共聚物 | mols: list, n_list: list | Mol |
| `terminate_rw(mol, ter)` | 添加端基 | mol: Mol, ter: Mol | Mol |
| `amorphous_cell(mol, n, density)` | 创建无定形单元 | mol: Mol, n: int, density: float | Mol |

### radonpy.sim.qm

| 函数 | 描述 | 参数 | 返回值 |
|-----|------|------|--------|
| `conformation_search(mol, **kwargs)` | 构象搜索 | mol: Mol | (Mol, energy) |
| `assign_charges(mol, charge='RESP')` | 计算原子电荷 | mol: Mol, charge: str | bool |
| `sp_prop(mol, **kwargs)` | 计算电子性质 | mol: Mol | dict |
| `polarizability(mol, **kwargs)` | 计算极化率 | mol: Mol | dict |

### radonpy.sim.preset.eq

| 类 | 描述 | 主要方法 |
|----|------|----------|
| `EQ21step` | 21步平衡协议 | exec(), analyze() |
| `Additional` | 额外NPT平衡 | exec(), analyze() |

### radonpy.sim.preset.tc

| 类 | 描述 | 主要方法 |
|----|------|----------|
| `NEMD_MP` | Müller-Plathe热导率 | exec(), analyze() |

---

## 实例代码

### 实例1：含氟聚合物的特殊处理

```python
# 聚四氟乙烯(PTFE)需要使用修正的力场参数
from radonpy.ff.gaff2_mod import GAFF2_mod

smiles = '*C(F)C(*)F'  # PTFE
mol = utils.mol_from_smiles(smiles)

# 使用GAFF2_mod（针对氟优化）
ff = GAFF2_mod()

# 生成聚合物
polymer = poly.polymerize_rw(mol, 50, tacticity='syndiotactic')
polymer = poly.terminate_rw(polymer, utils.mol_from_smiles('*F'))

# 力场分配
success = ff.ff_assign(polymer)

# 由于PTFE密度较高，调整初始密度
cell = poly.amorphous_cell(polymer, 10, density=0.1)  # 更高的初始密度
```

### 实例2：计算玻璃化转变温度

```python
def calculate_tg(polymer, temp_range=(200, 500), n_points=10):
    """通过密度-温度曲线计算Tg"""
    
    temperatures = np.linspace(temp_range[0], temp_range[1], n_points)
    densities = []
    
    for T in temperatures:
        # 在每个温度下平衡
        eqmd = eq.Additional(polymer, work_dir=f'./Tg_T{T}')
        cell = eqmd.exec(temp=T, press=1.0)
        
        # 获取密度
        analy = eqmd.analyze()
        props = analy.get_all_prop(temp=T, press=1.0)
        densities.append(props['density'])
    
    # 拟合并找到转折点
    from scipy.optimize import curve_fit
    
    def bilinear(T, Tg, m1, m2, b):
        """双线性模型"""
        return np.where(T < Tg, m1*T + b, m2*T + b + (m1-m2)*Tg)
    
    popt, _ = curve_fit(bilinear, temperatures, densities, p0=[350, -0.001, -0.0005, 1.0])
    Tg = popt[0]
    
    return Tg, temperatures, densities
```

### 实例3：聚合物共混物模拟

```python
def simulate_polymer_blend(polymer1, polymer2, ratio=(0.5, 0.5), compatibility=True):
    """模拟聚合物共混物"""
    
    # 计算链数
    total_chains = 10
    n1 = int(total_chains * ratio[0])
    n2 = int(total_chains * ratio[1])
    
    # 创建共混物
    blend = poly.amorphous_mixture_cell(
        [polymer1, polymer2],
        [n1, n2],
        density=0.05
    )
    
    if compatibility:
        # 相容体系：正常平衡
        eqmd = eq.EQ21step(blend, work_dir='./blend_compatible')
    else:
        # 不相容体系：可能需要更长时间
        eqmd = eq.EQ21step(blend, work_dir='./blend_incompatible')
        # 修改平衡参数
        eqmd.eq_nstep = 10000000  # 更长的平衡时间
    
    blend = eqmd.exec(temp=300, press=1.0)
    
    # 分析相分离
    analy = eqmd.analyze()
    
    # 计算径向分布函数判断混合程度
    rdf = analy.calc_rdf(type1='polymer1', type2='polymer2')
    
    return blend, rdf
```

### 实例4：结晶聚合物模拟

```python
def create_crystalline_polymer(mol, crystal_params):
    """创建结晶聚合物"""
    
    # 生成伸直链
    polymer = poly.polymerize_rw(
        mol, 
        n=20,
        tacticity='isotactic',  # 等规立构有利于结晶
        extended=True  # 伸直链构象
    )
    
    # 创建晶胞
    crystal = poly.crystal_cell(
        polymer,
        a=crystal_params['a'],
        b=crystal_params['b'], 
        c=crystal_params['c'],
        alpha=crystal_params.get('alpha', 90),
        beta=crystal_params.get('beta', 90),
        gamma=crystal_params.get('gamma', 90),
        n_cells=(2, 2, 1)  # 2x2x1超胞
    )
    
    # 优化晶体结构
    from radonpy.sim import md
    
    # 能量最小化
    md_params = md.MdInput()
    md_params.minimize = True
    md_params.min_style = 'cg'
    md_params.min_tol = 1e-4
    
    lmp = lammps.LAMMPS(crystal, work_dir='./crystal_opt')
    lmp.make_input(md_params)
    lmp.run()
    
    return crystal
```

---

## 常见问题

### Q1: 安装时出现依赖冲突怎么办？

**A:** 建议创建新的conda环境，严格按照版本要求安装：
```bash
conda create -n radonpy_clean python=3.11
conda activate radonpy_clean
# 然后按照安装指南重新安装
```

### Q2: LAMMPS找不到怎么办？

**A:** 设置环境变量：
```bash
# 查找LAMMPS位置
which lmp_serial
# 或
which lmp_mpi

# 设置环境变量
export LAMMPS_EXEC=/path/to/lammps
```

### Q3: 内存不足错误

**A:** 调整模拟参数：
- 减少聚合物链数：`poly.amorphous_cell(polymer, 5)`  # 从10减到5
- 减少原子数：使用更小的n值
- 使用更多的MPI进程分配内存

### Q4: 力场分配失败

**A:** 检查以下几点：
1. SMILES是否正确（使用`*`标记连接点）
2. 是否包含不支持的元素（GAFF2仅支持H,C,N,O,F,P,S,Cl,Br,I）
3. 尝试使用不同的力场（GAFF、GAFF2、GAFF2_mod）

### Q5: 平衡状态一直达不到

**A:** 可能的解决方案：
1. 增加额外平衡次数
2. 提高模拟温度加速平衡
3. 检查初始结构是否合理
4. 对于高Tg聚合物，在Tg以上温度平衡

### Q6: 热导率计算结果异常

**A:** 检查：
1. 系统是否完全平衡
2. 温度梯度是否线性（`Tgrad_check`）
3. 模拟时间是否足够长
4. 系统尺寸是否足够大

### Q7: 如何提高计算速度？

**A:** 优化策略：
```python
# 1. 使用并行计算
exec(mpi=8, omp=2, gpu=1)  # 使用8个MPI进程，2个OpenMP线程，1个GPU

# 2. 减少DFT计算
qm.conformation_search(mol, dft_nconf=2)  # 只优化2个构象

# 3. 使用更大的时间步长（小心使用）
md_params.dt = 2.0  # fs

# 4. 减少输出频率
md_params.dump_freq = 10000  # 减少轨迹输出
```

### Q8: 如何处理特殊聚合物？

**特殊处理示例：**

```python
# 导电聚合物（如聚苯胺）
# 需要考虑不同氧化态
pani_reduced = '*Nc1ccc(cc1)N*'
pani_oxidized = '*N=C1C=CC(=N*)C=C1'

# 生物可降解聚合物（如PLA）
# 需要正确的立体化学
pla_L = '*OC(C)C(=O)*'  # L-乳酸
pla_D = '*OC([C@@H](C))C(=O)*'  # D-乳酸

# 离子聚合物
# 需要添加抗衡离子
polymer_ion = poly.add_counterions(polymer, charge=-1, ion='Na+')
```

---

## 参考文献

1. **RadonPy主要论文**
   - Y. Hayashi, J. Shiomi, J. Morikawa, R. Yoshida, "RadonPy: Automated Physical Property Calculation using All-atom Classical Molecular Dynamics Simulations for Polymer Informatics," npj Comput. Mater., 8:222 (2022)

2. **力场相关**
   - J. Wang et al., "Development and testing of a general amber force field," J. Comput. Chem., 25, 1157-1174 (2004)
   - J. Trag, D. Zahn, "Improved GAFF2 parameters for fluorinated alkanes," J. Mol. Model., 25, 39 (2019)

3. **MD协议**
   - G.S. Larsen et al., "Molecular Simulations of PIM-1-like Polymers," Macromolecules, 44, 6944-6951 (2011) [21步协议]

4. **相关项目**
   - [XenonPy](https://github.com/yoshida-lab/XenonPy) - 材料信息学机器学习工具
   - [SMiPoly](https://github.com/PEJpOhno/SMiPoly) - 虚拟聚合物生成器

---

## 贡献与支持

### 报告问题
如果遇到bug或有功能建议，请在GitHub提交issue：
https://github.com/RadonPy/RadonPy/issues

### 贡献代码
欢迎提交Pull Request！请确保：
1. 代码符合PEP 8规范
2. 添加适当的测试
3. 更新相关文档

### 引用
如果在研究中使用了RadonPy，请引用：
```bibtex
@article{hayashi2022radonpy,
  title={RadonPy: Automated Physical Property Calculation using All-atom Classical Molecular Dynamics Simulations for Polymer Informatics},
  author={Hayashi, Yoshihiro and Shiomi, Junichiro and Morikawa, Junko and Yoshida, Ryo},
  journal={npj Computational Materials},
  volume={8},
  pages={222},
  year={2022},
  publisher={Nature Publishing Group}
}
```

### 许可证
RadonPy采用BSD-3-Clause许可证。

### 联系方式
- 主要开发者：Yoshihiro Hayashi (yhayashi@ism.ac.jp)
- 研究机构：统计数理研究所（The Institute of Statistical Mathematics）

---

## 附录

### A. 支持的物理性质列表

| 性质 | 单位 | 计算方法 | 所需模拟 |
|-----|------|----------|----------|
| 密度 | g/cm³ | 统计平均 | NPT |
| 比热容(Cp) | J/kg·K | 焓涨落 | NPT |
| 比热容(Cv) | J/kg·K | 能量涨落 | NVT |
| 线膨胀系数 | K⁻¹ | 长度-温度相关 | NPT |
| 体积膨胀系数 | K⁻¹ | 体积-温度相关 | NPT |
| 压缩率 | Pa⁻¹ | 体积涨落 | NPT |
| 体积模量 | Pa | 压缩率倒数 | NPT |
| 等熵压缩率 | Pa⁻¹ | 声速相关 | NPT |
| 等熵体积模量 | Pa | 等熵压缩率倒数 | NPT |
| 静态介电常数 | - | 偶极矩涨落 | NVT |
| 折射率 | - | Lorentz-Lorenz | DFT+MD |
| 回转半径 | Å | 质心距离 | Any |
| 端到端距离 | Å | 链端距离 | Any |
| 向列序参数 | - | 取向有序度 | Any |
| 热导率 | W/m·K | Müller-Plathe | NEMD |

### B. SMILES输入规范

**基本规则：**
- 使用`*`表示聚合物连接点（必须有2个）
- 支持手性中心：`@`和`@@`
- 支持顺反异构：`/`和`\`
- 支持芳香环：小写字母

**示例：**
```python
# 简单烷基链
'*CC(*)'  # 聚乙烯
'*CC(*)(C)'  # 聚丙烯

# 含官能团
'*CC(*)C(=O)O'  # 聚丙烯酸
'*CC(*)O'  # 聚乙烯醇

# 芳香聚合物
'*CC(*)c1ccccc1'  # 聚苯乙烯
'*Oc1ccc(cc1)O*'  # 聚对苯二酚

# 杂环聚合物
'*n1cccc1*'  # 聚吡咯
'*c1cncs1*'  # 聚噻唑
```

### C. 计算资源估算

| 聚合物体系 | 原子数 | CPU核数 | 内存 | 计算时间 |
|-----------|--------|---------|------|----------|
| 小体系 | ~1,000 | 4 | 4GB | 2-4小时 |
| 中等体系 | ~5,000 | 8 | 8GB | 6-12小时 |
| 大体系 | ~10,000 | 16 | 16GB | 12-24小时 |
| 超大体系 | ~50,000 | 32+ | 32GB+ | 1-3天 |

### D. 故障排除检查清单

- [ ] Python版本是否正确（3.7-3.12）
- [ ] 所有依赖包是否安装完整
- [ ] LAMMPS路径是否正确设置
- [ ] 工作目录是否有写入权限
- [ ] SMILES格式是否正确
- [ ] 初始密度是否合理（0.01-0.1 g/cm³）
- [ ] MPI进程数是否匹配CPU核数
- [ ] 是否有足够的磁盘空间（至少10GB）

---

**文档版本：** 1.0.0  
**最后更新：** 2025年1月  
**RadonPy版本：** 0.2.10

本文档持续更新中，最新版本请访问：https://github.com/RadonPy/RadonPy