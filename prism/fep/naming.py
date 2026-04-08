#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM FEP System Naming Standardization

This module provides standardized naming conventions for FEP system output directories.
Ensures consistent directory naming across different force field combinations.
"""

import re
from typing import Optional, Dict
from pathlib import Path


class FEPSystemNamer:
    """FEP 系统目录命名标准化工具

    标准格式: <protein_ff>-mut_<ligand_ff>
    示例: amber14sb_ol15-mut_gaff2
    """

    # 力场名称标准化映射
    FORCEFIELD_NORMALIZATION = {
        # AMBER force fields
        "amber14sb_OL15": "amber14sb_ol15",
        "amber14sb": "amber14sb",
        "amber99sb-ildn": "amber99sb_ildn",
        "amber19sb": "amber19sb",
        "amber03": "amber03",
        # CHARMM force fields
        "charmm36_jul2022": "charmm36_jul2022",
        "charmm36m": "charmm36m",
        "charmm27": "charmm27",
        # OPLS force fields
        "oplsaa": "oplsaa",
        "oplsaam": "oplsaam",
        # Ligand force fields
        "gaff2": "gaff2",
        "gaff": "gaff",
        "openff": "openff",
        "cgenff": "cgenff",
        "opls": "opls",
    }

    @classmethod
    def generate_name(cls, protein_ff: str, ligand_ff: str) -> str:
        """生成标准 FEP 系统目录名

        Args:
            protein_ff: 蛋白质力场名称 (如 amber14sb_OL15)
            ligand_ff: 配体力场名称 (如 gaff2, opls)

        Returns:
            标准化目录名，格式: <protein_ff>-mut_<ligand_ff>
            例如: amber14sb_ol15-mut_gaff2

        Examples:
            >>> FEPSystemNamer.generate_name("amber14sb_OL15", "gaff2")
            'amber14sb_ol15-mut_gaff2'
            >>> FEPSystemNamer.generate_name("charmm36_jul2022", "cgenff")
            'charmm36_jul2022-mut_cgenff'
        """
        # 标准化力场名称
        protein = cls._normalize_forcefield_name(protein_ff)
        ligand = cls._normalize_forcefield_name(ligand_ff)

        # 生成标准目录名
        standard_name = f"{protein}-mut_{ligand}"

        return standard_name

    @classmethod
    def validate_name(cls, name: str) -> bool:
        """验证目录名是否符合 FEP 系统命名规范

        规范: <protein_ff>-mut_<ligand_ff>
        - 只能包含小写字母、数字、下划线、连字符
        - 必须包含 "-mut_" 分隔符

        Args:
            name: 目录名称

        Returns:
            是否符合规范

        Examples:
            >>> FEPSystemNamer.validate_name("amber14sb_ol15-mut_gaff2")
            True
            >>> FEPSystemNamer.validate_name("oplsaa")
            False
            >>> FEPSystemNamer.validate_name("amber14sb_ol15_gaff2")
            False
        """
        # 检查格式: <protein_ff>-mut_<ligand_ff>
        pattern = r"^[a-z0-9_]+-mut_[a-z0-9_]+$"
        return bool(re.match(pattern, name))

    @classmethod
    def suggest_name(cls, protein_ff: str, ligand_ff: str) -> str:
        """建议目录名（带提示信息）

        Args:
            protein_ff: 蛋白质力场名称
            ligand_ff: 配体力场名称

        Returns:
            建议的标准目录名
        """
        standard = cls.generate_name(protein_ff, ligand_ff)
        return standard

    @classmethod
    def parse_name(cls, name: str) -> Optional[Dict[str, str]]:
        """解析标准目录名，提取力场信息

        Args:
            name: 标准目录名

        Returns:
            包含 protein_ff 和 ligand_ff 的字典，如果解析失败返回 None

        Examples:
            >>> FEPSystemNamer.parse_name("amber14sb_ol15-mut_gaff2")
            {'protein_ff': 'amber14sb_ol15', 'ligand_ff': 'gaff2'}
        """
        if not cls.validate_name(name):
            return None

        parts = name.split("-mut_")
        if len(parts) == 2:
            return {"protein_ff": parts[0], "ligand_ff": parts[1]}
        return None

    @classmethod
    def _normalize_forcefield_name(cls, ff_name: str) -> str:
        """标准化力场名称

        转换规则：
        1. 转小写
        2. 移除多余分隔符（- 和 _）
        3. 使用映射表进行标准化

        Args:
            ff_name: 原始力场名称

        Returns:
            标准化后的力场名称
        """
        # Check mapping case-insensitively
        normalized_lower = ff_name.lower()
        for key, value in cls.FORCEFIELD_NORMALIZATION.items():
            if key.lower() == normalized_lower:
                return value

        # 标准化处理：转小写，移除特殊字符
        normalized = ff_name.lower()
        # 移除所有连字符和下划线
        normalized = re.sub(r"[-_]+", "", normalized)

        return normalized

    @classmethod
    def get_default_output_dir(cls, protein_ff: str, ligand_ff: str, base_dir: str = ".") -> Path:
        """获取默认输出目录路径

        Args:
            protein_ff: 蛋白质力场名称
            ligand_ff: 配体力场名称
            base_dir: 基础目录（默认为当前目录）

        Returns:
            完整的输出目录路径
        """
        standard_name = cls.generate_name(protein_ff, ligand_ff)
        return Path(base_dir) / standard_name


def generate_fep_system_name(protein_ff: str, ligand_ff: str) -> str:
    """生成 FEP 系统目录名的便捷函数

    Args:
        protein_ff: 蛋白质力场名称
        ligand_ff: 配体力场名称

    Returns:
        标准化目录名

    Examples:
        >>> generate_fep_system_name("amber14sb_OL15", "gaff2")
        'amber14sb_ol15-mut_gaff2'
    """
    return FEPSystemNamer.generate_name(protein_ff, ligand_ff)


def validate_fep_system_name(name: str) -> bool:
    """验证 FEP 系统目录名的便捷函数

    Args:
        name: 目录名称

    Returns:
        是否符合规范
    """
    return FEPSystemNamer.validate_name(name)


def parse_fep_system_name(name: str) -> Optional[Dict[str, str]]:
    """解析 FEP 系统目录名的便捷函数

    Args:
        name: 标准目录名

    Returns:
        包含力场信息的字典，解析失败返回 None
    """
    return FEPSystemNamer.parse_name(name)
