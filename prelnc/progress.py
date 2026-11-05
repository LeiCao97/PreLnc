from __future__ import absolute_import

"""Progress reporting and end-of-run summaries.

Progress lines go to stderr so that stdout stays usable for piping. On a
terminal a single line is refreshed in place; when stderr is redirected the
same information is printed as one line per 5% so logs stay readable.
"""

import shutil
import sys
import time
import unicodedata

PREFIX = "[PreLnc]"

REFRESH_INTERVAL = 0.1
BUCKET_COUNT = 20


def _display_width(text):
    width = 0
    for char in text:
        width += 2 if unicodedata.east_asian_width(char) in ("W", "F") else 1
    return width


def pad(text, width):
    return text + " " * max(0, width - _display_width(text))


def truncate(text, width):
    if _display_width(text) <= width:
        return text
    result = []
    used = 0
    for char in text:
        char_width = 2 if unicodedata.east_asian_width(char) in ("W", "F") else 1
        if used + char_width > width:
            break
        result.append(char)
        used += char_width
    return "".join(result)


def _write(stream, text):
    try:
        stream.write(text)
    except UnicodeEncodeError:
        encoding = getattr(stream, "encoding", None) or "ascii"
        stream.write(text.encode(encoding, "replace").decode(encoding, "replace"))
    stream.flush()


def format_duration(seconds):
    seconds = max(0.0, float(seconds))
    if seconds < 60:
        return "%.0fs" % seconds
    minutes, remainder = divmod(int(seconds), 60)
    if minutes < 60:
        return "%dm%02ds" % (minutes, remainder)
    hours, minutes = divmod(minutes, 60)
    return "%dh%02dm" % (hours, minutes)


def info(message, stream=None):
    _write(stream or sys.stderr, "%s %s\n" % (PREFIX, message))


def warn(message, stream=None):
    _write(stream or sys.stderr, "%s WARNING: %s\n" % (PREFIX, message))


class Progress(object):
    """A counter with percentage, elapsed time and a rough ETA."""

    def __init__(self, total, label, stream=None, enabled=True, unit="items"):
        self.label = label
        self.total = max(int(total), 0)
        self.unit = unit
        self.stream = stream or sys.stderr
        self.enabled = bool(enabled) and self.total > 0
        self.count = 0
        self.started = time.time()
        self.finished = False
        self._last_render = 0.0
        self._last_reported = 0
        self._interactive = bool(getattr(self.stream, "isatty", lambda: False)())
        self._bucket = max(1, self.total // BUCKET_COUNT) if self.total > BUCKET_COUNT else 1
        if self.enabled and self._interactive:
            self._render(force=True)

    def _line(self):
        elapsed = time.time() - self.started
        if self.count > 0:
            remaining = elapsed * (self.total - self.count) / float(self.count)
        else:
            remaining = 0.0
        percent = 100.0 * self.count / self.total
        return "%s %s [%3d%%] %d/%d %s · elapsed %s · eta %s" % (
            PREFIX,
            self.label,
            int(percent),
            self.count,
            self.total,
            self.unit,
            format_duration(elapsed),
            format_duration(remaining),
        )

    def _render(self, force=False):
        if not self.enabled:
            return
        now = time.time()
        if not force and now - self._last_render < REFRESH_INTERVAL:
            return
        self._last_render = now
        line = self._line()
        columns = shutil.get_terminal_size((80, 20)).columns - 1
        _write(self.stream, "\r" + truncate(line, max(columns, 20)) + "\x1b[K")

    def advance(self, step=1):
        if not self.enabled:
            return
        self.count = min(self.count + step, self.total)
        if self._interactive:
            self._render()
        elif self.count - self._last_reported >= self._bucket or self.count >= self.total:
            self._last_reported = self.count
            _write(self.stream, self._line() + "\n")

    def finish(self):
        if not self.enabled or self.finished:
            return
        self.finished = True
        self.count = self.total
        line = self._line()
        if self._interactive:
            _write(self.stream, "\r" + truncate(line, max(shutil.get_terminal_size((80, 20)).columns - 1, 20)) + "\x1b[K\n")
        elif self._last_reported != self.total:
            _write(self.stream, line + "\n")


def print_summary(rows, title=None, stream=None):
    """Print an aligned key/value block, e.g. ``[PreLnc] Output   : path``."""
    stream = stream or sys.stderr
    if title:
        _write(stream, "%s %s\n" % (PREFIX, title))
    width = max(_display_width(label) for label, _ in rows) if rows else 0
    for label, value in rows:
        _write(stream, "%s %s : %s\n" % (PREFIX, pad(label, width), value))
