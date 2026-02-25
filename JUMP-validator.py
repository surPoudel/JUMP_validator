#!/usr/bin/env python
"""Legacy entrypoint wrapper.

This keeps the historical filename (`JUMP-validator.py`) working, while routing
execution to the new cross-platform CLI implementation.

Usage (same as before):

    python JUMP-validator.py path/to/jump_validator.params

"""

from jump_validator.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
