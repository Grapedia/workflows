#!/usr/bin/env python3
"""Minimal inline-SVG chart helpers for the TITAN HTML reports (no dependencies).

Every function returns an SVG string that takes its colours from CSS custom properties
(--c1 ... --c6, --ink, --ink2, --grid, --surface), so a single stylesheet gives the light
and the dark theme.  Hover tooltips are <title> elements; a data table can be attached by
the caller (see `data_table`).

Categorical colour order (validated with the dataviz palette validator, light and dark):
blue, orange, aqua, yellow, violet, red  ->  --c1 ... --c6.
"""
import html

W = 760  # viewBox width; SVGs scale to the container


def esc(text):
    return html.escape(str(text), quote=True)


def fmt_int(n):
    return "{:,}".format(int(round(n))).replace(",", " ")


def dec_fr(x, digits=1):
    """Number with a French decimal comma."""
    return (("%." + str(digits) + "f") % x).replace(".", ",")


def fmt_pct(x, digits=1):
    return dec_fr(x, digits) + "\u202f%"


def _svg(height, body, label):
    return ('<svg class="chart" viewBox="0 0 %d %d" role="img" aria-label="%s" '
            'preserveAspectRatio="xMinYMin meet">%s</svg>' % (W, height, esc(label), body))


def _text(x, y, s, cls="t2", anchor="start", extra=""):
    return '<text x="%.1f" y="%.1f" class="%s" text-anchor="%s" %s>%s</text>' % (x, y, cls, anchor, extra, esc(s))


def hbars(items, label, label_w=250, bar_h=22, gap=9, value_fmt=None, max_val=None, right_pad=90):
    """items: list of dict(label, value, color=1..6, text=None, tip=None)."""
    max_val = max_val or max(i["value"] for i in items) or 1
    plot_w = W - label_w - right_pad
    h = len(items) * (bar_h + gap) + gap
    out = []
    for k, it in enumerate(items):
        y = gap + k * (bar_h + gap)
        w = max(1.0, plot_w * it["value"] / max_val)
        txt = it.get("text") or (value_fmt(it["value"]) if value_fmt else fmt_int(it["value"]))
        tip = it.get("tip") or "%s : %s" % (it["label"], txt)
        out.append(_text(label_w - 10, y + bar_h * 0.7, it["label"], "t1", "end"))
        out.append('<g><title>%s</title><rect x="%d" y="%.1f" width="%.1f" height="%d" rx="3" '
                   'class="f%d"/></g>' % (esc(tip), label_w, y, w, bar_h, it.get("color", 1)))
        out.append(_text(label_w + w + 8, y + bar_h * 0.7, txt, "t1"))
    return _svg(h, "".join(out), label)


def stacked_hbars(rows, segments, label, label_w=170, bar_h=26, gap=12, pct=True, min_label_w=46,
                  legend=True):
    """rows: list of (row_label, [values...]); segments: list of (name, color_index).
    Each bar is normalised to 100 %."""
    plot_w = W - label_w - 20
    legend_h = 26 if legend else 0
    h = legend_h + len(rows) * (bar_h + gap) + gap
    out = []
    if legend:
        x = label_w
        for name, c in segments:
            out.append('<rect x="%d" y="6" width="12" height="12" rx="2" class="f%d"/>' % (x, c))
            out.append(_text(x + 17, 16, name, "t2"))
            x += 17 + 7.2 * len(name) + 16
    for r, (rl, vals) in enumerate(rows):
        y = legend_h + gap + r * (bar_h + gap)
        tot = float(sum(vals)) or 1.0
        out.append(_text(label_w - 10, y + bar_h * 0.68, rl, "t1", "end"))
        x = label_w
        for (name, c), v in zip(segments, vals):
            w = plot_w * v / tot
            if w <= 0:
                continue
            tip = "%s | %s : %s (%s)" % (rl, name, fmt_int(v), fmt_pct(100 * v / tot, 1))
            out.append('<g><title>%s</title><rect x="%.1f" y="%d" width="%.1f" height="%d" class="f%d" '
                       'stroke="var(--surface)" stroke-width="2"/></g>' % (esc(tip), x, y, w, bar_h, c))
            if w >= min_label_w:
                s = fmt_pct(100 * v / tot, 1) if pct else fmt_int(v)
                out.append(_text(x + w / 2, y + bar_h * 0.68, s, "td" if c in (3, 4) else ("tr%d" % (c - 6) if c >= 7 else "tw"), "middle"))
            x += w
    return _svg(h, "".join(out), label)


def grouped_hbars(categories, series, label, label_w=250, bar_h=15, gap=14, unit="%", max_val=None,
                  right_pad=80):
    """categories: list of labels; series: list of (name, [values per category], color_index)."""
    max_val = max_val or max(max(v) for _, v, _ in series) or 1
    plot_w = W - label_w - right_pad
    legend_h = 26
    block = bar_h * len(series) + 4
    h = legend_h + len(categories) * (block + gap) + gap
    out = []
    x = label_w
    for name, _, c in series:
        out.append('<rect x="%d" y="6" width="12" height="12" rx="2" class="f%d"/>' % (x, c))
        out.append(_text(x + 17, 16, name, "t2"))
        x += 17 + 7.2 * len(name) + 16
    for k, cat in enumerate(categories):
        y0 = legend_h + gap + k * (block + gap)
        out.append(_text(label_w - 10, y0 + block / 2 + 4, cat, "t1", "end"))
        for s, (name, vals, c) in enumerate(series):
            y = y0 + s * (bar_h + 2)
            w = max(1.0, plot_w * vals[k] / max_val) if vals[k] > 0 else 0
            tip = "%s | %s : %s%s" % (cat, name, dec_fr(vals[k]), unit)
            out.append('<g><title>%s</title><rect x="%d" y="%.1f" width="%.1f" height="%d" rx="3" class="f%d"/></g>'
                       % (esc(tip), label_w, y, w, bar_h, c))
            out.append(_text(label_w + w + 6, y + bar_h * 0.78, dec_fr(vals[k]) + unit.replace("%", "\u202f%"),
                             "t2"))
    return _svg(h, "".join(out), label)


def line_chart(series, label, xlabel, ylabel, height=300, xmax=None, ymax=None, bands=None,
               xticks=None, yticks=None, yfmt=None, left=64, xlog=False):
    """series: list of (name, xs, ys, color_index); bands: list of (xs, lo, hi, color_index)."""
    import math
    right, top, bottom = 24, 36, 46
    pw, ph = W - left - right, height - top - bottom
    xmax = xmax or max(max(xs) for _, xs, _, _ in series)
    xmin = min(min(xs) for _, xs, _, _ in series)
    ymax = ymax or max(max(ys) for _, _, ys, _ in series)
    yfmt = yfmt or (lambda v: fmt_int(v))

    def X(v):
        if xlog:
            return left + pw * (math.log10(v) - math.log10(xmin)) / (math.log10(xmax) - math.log10(xmin))
        return left + pw * (v - xmin) / float(xmax - xmin)

    def Y(v):
        return top + ph * (1 - v / float(ymax))

    out = []
    x = left
    for name, _, _, c in series:
        out.append('<rect x="%d" y="8" width="12" height="12" rx="2" class="f%d"/>' % (x, c))
        out.append(_text(x + 17, 18, name, "t2"))
        x += 17 + 7.2 * len(name) + 16
    for v in (yticks or [ymax * i / 4.0 for i in range(5)]):
        out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" class="grid"/>' % (left, left + pw, Y(v), Y(v)))
        out.append(_text(left - 8, Y(v) + 4, yfmt(v), "t2", "end"))
    for v in (xticks or [xmin + (xmax - xmin) * i / 5.0 for i in range(6)]):
        out.append(_text(X(v), top + ph + 18, ("%g" % v).replace(".", ","), "t2", "middle"))
    out.append(_text(left + pw / 2, height - 6, xlabel, "t2", "middle"))
    out.append('<text transform="translate(14,%d) rotate(-90)" class="t2" text-anchor="middle">%s</text>'
               % (top + ph / 2, esc(ylabel)))
    for xs, lo, hi, c in (bands or []):
        pts = ["%.1f,%.1f" % (X(a), Y(b)) for a, b in zip(xs, hi)] + \
              ["%.1f,%.1f" % (X(a), Y(b)) for a, b in reversed(list(zip(xs, lo)))]
        out.append('<polygon points="%s" class="f%d" opacity="0.18"/>' % (" ".join(pts), c))
    for name, xs, ys, c in series:
        pts = " ".join("%.1f,%.1f" % (X(a), Y(b)) for a, b in zip(xs, ys))
        out.append('<polyline points="%s" fill="none" class="s%d" stroke-width="2" stroke-linejoin="round" '
                   'stroke-linecap="round"/>' % (pts, c))
        # hover points
        for a, b in list(zip(xs, ys))[::max(1, len(xs) // 24)]:
            out.append('<g><title>%s | x=%g : %s</title><circle cx="%.1f" cy="%.1f" r="6" fill="transparent"/>'
                       '<circle cx="%.1f" cy="%.1f" r="2.5" class="f%d"/></g>'
                       % (esc(name), a, esc(yfmt(b)), X(a), Y(b), X(a), Y(b), c))
    return _svg(height, "".join(out), label)


def heatmap(row_labels, col_labels, values, texts, label, label_w=210, cell_w=78, cell_h=30, head_h=120):
    """values in [0,1] drive the blue sequential ramp; texts are the printed cell strings."""
    ramp = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    w = label_w + cell_w * len(col_labels) + 110  # room for the rotated column headers
    h = head_h + cell_h * len(row_labels) + 8
    out = []
    for j, cl in enumerate(col_labels):
        cx = label_w + j * cell_w + cell_w / 2
        out.append('<text transform="translate(%.1f,%d) rotate(-30)" class="t2">%s</text>' % (cx - 8, head_h - 8, esc(cl)))
    for i, rl in enumerate(row_labels):
        y = head_h + i * cell_h
        out.append(_text(label_w - 10, y + cell_h * 0.65, rl, "t1", "end"))
        for j, v in enumerate(values[i]):
            idx = 0 if v is None else min(len(ramp) - 1, int(v * len(ramp)))
            dark_text = idx >= 3
            out.append('<g><title>%s | %s : %s</title><rect x="%d" y="%d" width="%d" height="%d" fill="%s" '
                       'stroke="var(--surface)" stroke-width="2"/></g>'
                       % (esc(rl), esc(col_labels[j]), esc(texts[i][j]), label_w + j * cell_w, y, cell_w, cell_h, ramp[idx]))
            out.append('<text x="%d" y="%d" text-anchor="middle" class="%s">%s</text>'
                       % (label_w + j * cell_w + cell_w / 2, y + cell_h * 0.65, "hm-w" if dark_text else "hm-d",
                          esc(texts[i][j])))
    svg = ('<svg class="chart" viewBox="0 0 %d %d" role="img" aria-label="%s" preserveAspectRatio="xMinYMin meet">%s</svg>'
           % (max(w, W), h, esc(label), "".join(out)))
    return svg


def data_table(headers, rows, caption=""):
    """Plain HTML table used as the accessible table view of a chart."""
    th = "".join("<th>%s</th>" % esc(h) for h in headers)
    body = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % esc(c) for c in r) for r in rows)
    cap = "<caption>%s</caption>" % esc(caption) if caption else ""
    return '<div class="tablewrap"><table>%s<thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>' % (cap, th, body)


def figure(svg, title, caption, table_html=""):
    det = ""
    if table_html:
        det = '<details class="tv"><summary>Voir les données</summary>%s</details>' % table_html
    return ('<figure><figcaption class="ft">%s</figcaption>%s<p class="fc">%s</p>%s</figure>'
            % (title, svg, caption, det))
