#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP HTML可视化验证脚本

按照 .claude/skills/fep-visualization.md 的要求全面检查HTML输出。

Usage:
    python fep-visualization-check.py <mapping.html>

Examples:
    python fep-visualization-check.py oplsaa/GMX_PROLIG_FEP/common/hybrid/mapping.html
"""

import sys
import json
import re
from pathlib import Path
from collections import Counter
from typing import Dict, Tuple, Any


def load_html(html_path: str) -> str:
    """Load HTML file content."""
    path = Path(html_path)
    if not path.exists():
        print(f"❌ HTML file not found: {html_path}")
        sys.exit(1)
    return path.read_text()


def extract_data(html: str) -> Dict[str, Any]:
    """Extract key data from HTML."""
    data = {}

    # Extract atoms and bonds
    for mol_id in ["A", "B"]:
        # Atoms
        match = re.search(rf"const ATOMS_{mol_id} = (\[.*?\]);", html, re.DOTALL)
        if match:
            data[f"atoms_{mol_id.lower()}"] = json.loads(match.group(1))

        # Bonds
        match = re.search(rf"const BONDS_{mol_id} = (\[.*?\]);", html, re.DOTALL)
        if match:
            data[f"bonds_{mol_id.lower()}"] = json.loads(match.group(1))

    # Extract total charges
    match_a = re.search(r"molecule-label molecule-label-a.*?Total Charge: ([+-]?[\d.]+)", html, re.DOTALL)
    match_b = re.search(r"molecule-label molecule-label-b.*?Total Charge: ([+-]?[\d.]+)", html, re.DOTALL)
    if match_a:
        data["charge_a"] = float(match_a.group(1))
    if match_b:
        data["charge_b"] = float(match_b.group(1))

    # Extract statistics
    match = re.search(r"Common: (\d+)", html)
    if match:
        data["common_count_html"] = int(match.group(1))

    # Extract correspondence map
    match = re.search(r"const CORRESPONDENCE = (\{.*?\});", html, re.DOTALL)
    if match:
        data["correspondence"] = json.loads(match.group(1))

    return data


def check_gray_atoms(html: str) -> Tuple[bool, str]:
    """a. 灰色原子检测"""
    gray_count = html.count("rgb(200, 200, 200)")

    if gray_count == 0:
        return True, f"✓ 灰色原子: {gray_count} (应为0)"
    else:
        return False, f"✗ 灰色原子: {gray_count} (应为0)"


def check_unknown_atoms(data: Dict[str, Any]) -> Tuple[bool, str]:
    """a. Unknown atoms检测"""
    unknown_a = sum(1 for a in data.get("atoms_a", []) if a.get("classification") == "unknown")
    unknown_b = sum(1 for a in data.get("atoms_b", []) if a.get("classification") == "unknown")

    if unknown_a == 0 and unknown_b == 0:
        return True, f"✓ Unknown classification: A={unknown_a}, B={unknown_b} (均为0)"
    else:
        return False, f"✗ Unknown classification: A={unknown_a}, B={unknown_b} (应为0)"


def check_atom_placeholders(data: Dict[str, Any]) -> Tuple[bool, str]:
    """b. 原子标签检查"""
    placeholders_a = sum(1 for a in data.get("atoms_a", []) if a["name"].startswith("Atom"))
    placeholders_b = sum(1 for a in data.get("atoms_b", []) if a["name"].startswith("Atom"))

    # Some valid input formats legitimately use a small number of single-letter
    # atom names (for example S or N in MOL2). Treat this as a problem only if it
    # becomes a widespread placeholder pattern.
    only_element_a = sum(1 for a in data.get("atoms_a", []) if len(a["name"]) == 1 and a["name"].isupper())
    only_element_b = sum(1 for a in data.get("atoms_b", []) if len(a["name"]) == 1 and a["name"].isupper())

    if placeholders_a == 0 and placeholders_b == 0:
        result = f"✓ AtomXX占位符: A={placeholders_a}, B={placeholders_b} (均为0)"
    else:
        result = f"✗ AtomXX占位符: A={placeholders_a}, B={placeholders_b} (应为0)"

    element_name_limit = 3
    if only_element_a <= element_name_limit and only_element_b <= element_name_limit:
        result += f"\n✓ Single-letter labels accepted where chemically valid: A={only_element_a}, B={only_element_b}"
    else:
        result += f"\n✗ Too many single-letter labels: A={only_element_a}, B={only_element_b}"

    success = (
        placeholders_a == 0
        and placeholders_b == 0
        and only_element_a <= element_name_limit
        and only_element_b <= element_name_limit
    )
    return success, result


def check_total_charge(data: Dict[str, Any]) -> Tuple[bool, str]:
    """b. 总电荷显示（映射后）"""
    charge_a = data.get("charge_a", 0.0)
    charge_b = data.get("charge_b", 0.0)

    # 总电荷必须 ≈ 0（允许浮点误差）
    tolerance = 0.01
    ok_a = abs(charge_a) < tolerance
    ok_b = abs(charge_b) < tolerance

    if ok_a and ok_b:
        return True, f"✓ 总电荷: A={charge_a:.6f}, B={charge_b:.6f} (均≈0)"
    else:
        return False, f"✗ 总电荷: A={charge_a:.6f}, B={charge_b:.6f} (应≈0)"


def check_common_count_consistency(data: Dict[str, Any]) -> Tuple[bool, str]:
    """c. Common数量一致性"""
    common_a = sum(1 for a in data.get("atoms_a", []) if a["classification"] == "common")
    common_b = sum(1 for a in data.get("atoms_b", []) if a["classification"] == "common")

    if common_a == common_b:
        return True, f"✓ Common数量: A={common_a}, B={common_b} (相等)"
    else:
        return False, f"✗ Common数量: A={common_a}, B={common_b} (不相等)"


def check_common_atoms_mapping(data: Dict[str, Any]) -> Tuple[bool, str]:
    """Check that common atoms have correspondence entries.

    Atom names are not required to be identical across states. This is
    especially important for OPLS/LigParGen-style inputs where the native atom
    names can legitimately differ while the mapping is still correct.
    """
    atoms_a = data.get("atoms_a", [])
    atoms_b = data.get("atoms_b", [])
    correspondence = data.get("correspondence", {})

    common_a = [(i, a) for i, a in enumerate(atoms_a) if a["classification"] == "common"]
    common_b = [(i, a) for i, a in enumerate(atoms_b) if a["classification"] == "common"]

    missing_a = [a["name"] for i, a in common_a if f"a_{i}" not in correspondence]
    missing_b = [a["name"] for i, a in common_b if f"b_{i}" not in correspondence]

    if not missing_a and not missing_b and len(common_a) == len(common_b):
        return True, f"✓ Common原子有一一对应: A={len(common_a)}, B={len(common_b)}"
    return False, f"✗ Common原子对应缺失: A缺失={missing_a}, B缺失={missing_b}"


def check_bond_orders(data: Dict[str, Any]) -> Tuple[bool, str]:
    """d. 键级渲染检查"""
    results = []

    for mol_name, bond_key in [("Ligand A", "bonds_a"), ("Ligand B", "bonds_b")]:
        bonds = data.get(bond_key, [])
        atoms = data.get(f"atoms_{mol_name.split()[-1].lower()}", [])

        if not bonds:
            results.append(f"✗ {mol_name}: 无键数据")
            continue

        # Count bond types
        cnt = Counter(b.get("type", 1) for b in bonds)

        # Check if we have aromatic and double bonds
        has_aromatic = cnt.get(12, 0) > 0
        has_double = cnt.get(2, 0) > 0
        has_triple = cnt.get(3, 0) > 0

        # Check bond count (should be ≈ atom count - 1)
        total_bonds = len(bonds)
        total_atoms = len(atoms)
        bond_count_ok = abs(total_bonds - (total_atoms - 1)) <= 3

        status = "✓" if (has_aromatic or has_double or has_triple) and bond_count_ok else "✗"
        results.append(
            f"{status} {mol_name}: {total_bonds} bonds, "
            f"SINGLE={cnt[1]}, DOUBLE={cnt.get(2,0)}, TRIPLE={cnt.get(3,0)}, AROMATIC={cnt.get(12,0)}"
        )

    if all("✓" in r for r in results):
        return True, "\n".join(results)
    else:
        return False, "\n".join(results)


def check_atom_properties(data: Dict[str, Any]) -> Tuple[bool, str]:
    """检查原子属性完整性"""
    issues = []

    for mol_name, atom_key in [("Ligand A", "atoms_a"), ("Ligand B", "atoms_b")]:
        atoms = data.get(atom_key, [])

        for i, atom in enumerate(atoms[:10]):  # Check first 10 atoms
            required_fields = ["name", "element", "charge", "classification"]
            missing = [f for f in required_fields if f not in atom or atom[f] is None]

            if missing:
                issues.append(f"✗ {mol_name} atom[{i}] 缺少字段: {missing}")

        # Check for zero charges (unless truly neutral atom)
        zero_charge_a = sum(1 for a in atoms if a.get("charge", 0) == 0.0)
        zero_charge_b = sum(1 for a in atoms if a.get("charge", 0) == 0.0)

        if zero_charge_a > len(atoms) * 0.5:  # More than half have zero charge
            issues.append(f"⚠ {mol_name}: {zero_charge_a}/{len(atoms)} 原子电荷为0")

    if not issues:
        return True, "✓ 原子属性完整"
    else:
        return False, "\n".join(issues)


def validate_html(html_path: str) -> bool:
    """主验证函数"""
    print(f"=== FEP HTML可视化验证 ===")
    print(f"文件: {html_path}\n")

    html = load_html(html_path)
    data = extract_data(html)

    checks = [
        ("灰色原子检测", lambda: check_gray_atoms(html)),
        ("Unknown atoms", lambda: check_unknown_atoms(data)),
        ("原子标签", lambda: check_atom_placeholders(data)),
        ("总电荷", lambda: check_total_charge(data)),
        ("Common数量一致性", lambda: check_common_count_consistency(data)),
        ("Common原子对应", lambda: check_common_atoms_mapping(data)),
        ("键级渲染", lambda: check_bond_orders(data)),
        ("原子属性", lambda: check_atom_properties(data)),
    ]

    all_pass = True
    for name, check_func in checks:
        try:
            passed, message = check_func()
            print(f"[{name}]")
            print(f"  {message}\n")
            if not passed:
                all_pass = False
        except Exception as e:
            print(f"[{name}]")
            print(f"  ✗ 检查失败: {e}\n")
            all_pass = False

    print("=" * 70)
    if all_pass:
        print("✅ 所有检查通过！")
        return True
    else:
        print("✗ 部分检查失败，需要修复")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python fep-visualization-check.py <mapping.html>")
        print("\nExample:")
        print("  python fep-visualization-check.py oplsaa/GMX_PROLIG_FEP/common/hybrid/mapping.html")
        sys.exit(1)

    html_path = sys.argv[1]
    success = validate_html(html_path)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
