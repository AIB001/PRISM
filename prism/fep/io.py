#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compatibility shim for legacy ``prism.fep.io`` imports."""

from prism.fep.common import io as _impl

globals().update({name: getattr(_impl, name) for name in dir(_impl) if not name.startswith("__")})

__all__ = getattr(_impl, "__all__", [name for name in globals() if not name.startswith("__")])
