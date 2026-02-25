#!/usr/bin/env python
"""GUI entrypoint wrapper.

This provides a convenient, legacy-style filename for launching the GUI:

    python JUMP-validator-gui.py

The GUI simply asks for a `.params` file and then runs the same backend
pipeline as the CLI (`jump_validator.validator.run`).
"""

from jump_validator.gui import main


if __name__ == "__main__":
    raise SystemExit(main())
