# from .base import ForceFieldGeneratorBase
# from .gaff import GAFFForceFieldGenerator
# from .openff import OpenFFForceFieldGenerator

# __all__ = [
#     "ForceFieldGeneratorBase",
#     "GAFFForceFieldGenerator",
#     "OpenFFForceFieldGenerator",
# ]

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Force field generators for PRISM
"""

try:
    from .base import ForceFieldGeneratorBase
    from .gaff import GAFFForceFieldGenerator
    from .openff import OpenFFForceFieldGenerator
    
    __all__ = [
        "ForceFieldGeneratorBase",
        "GAFFForceFieldGenerator", 
        "OpenFFForceFieldGenerator",
    ]
except ImportError as e:
    # Handle missing dependencies gracefully
    print(f"Warning: Some force field generators may not be available: {e}")
    
    # Try to import what we can
    try:
        from .base import ForceFieldGeneratorBase
        __all__ = ["ForceFieldGeneratorBase"]
    except ImportError:
        __all__ = []
    
    try:
        from .gaff import GAFFForceFieldGenerator
        __all__.append("GAFFForceFieldGenerator")
    except ImportError:
        pass
        
    try:
        from .openff import OpenFFForceFieldGenerator
        __all__.append("OpenFFForceFieldGenerator")
    except ImportError:
        pass
    