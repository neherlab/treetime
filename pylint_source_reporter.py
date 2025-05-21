"""
Adds source code snippets to pylint output

Usage:

pip3 install pylint colorama pygments

PYTHONPATH=. pylint treetime/ --output-format=pylint_source_reporter.SourceCodeReporter
"""

from colorama import Fore, Style, init
from pygments import highlight
from pygments.formatters import Terminal256Formatter
from pygments.lexers import PythonLexer
from pylint.reporters.text import TextReporter

init(autoreset=True)


class SourceCodeReporter(TextReporter):
    def handle_message(self, msg):
        # Color and bold diagnostic message
        level_color = {
            "fatal": Fore.RED,
            "error": Fore.RED,
            "warning": Fore.YELLOW,
            "refactor": Fore.MAGENTA,
            "convention": Fore.CYAN,
            "info": Fore.GREEN,
        }.get(msg.category, "")

        self.writeln(
            f"{Style.BRIGHT}{level_color}"
            f"{msg.msg_id}: {msg.msg} ({msg.symbol}) @ {msg.path}:{msg.line}:{msg.column}"
            f"{Style.RESET_ALL}"
        )

        try:
            with open(msg.path, encoding="utf-8") as f:
                lines = f.readlines()
            start = max(msg.line - 5, 0)
            end = min(msg.line + 4, len(lines))
            raw_block = "".join(lines[start:end])

            highlighted_block = highlight(raw_block, PythonLexer(), Terminal256Formatter(style="material"))
            highlighted_lines = highlighted_block.splitlines()

            for i, rendered in enumerate(highlighted_lines, start=start + 1):
                prefix = ">>" if i == msg.line else "  "
                style = Style.BRIGHT + Fore.RED if i == msg.line else Style.DIM
                self.writeln(f"{style}{prefix} {i:4}:{Style.RESET_ALL} {rendered}")
        except Exception:
            pass

        self.writeln("")
        self.writeln("")
