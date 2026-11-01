#!/usr/bin/env node
/* Markdown -> submission-ready .docx for the rs4354668 manuscript set.
 *
 * Handles: headings, bold/italic/code, unicode super/subscripts, <sup> tags,
 * pipe tables (with dual widths), fenced code blocks, bullet lists, and
 * horizontal rules. Emits double-spaced body text with continuous line
 * numbers and page numbers, as journals require.
 */
const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  PageNumber, Footer, LineNumberRestartFormat, convertInchesToTwip,
} = require('docx');

const FONT = 'Times New Roman';
const PAGE_W = 11906, PAGE_H = 16838;          // A4 in DXA
const MARGIN = convertInchesToTwip(1);
const TEXT_W = PAGE_W - 2 * MARGIN;

// ---------------------------------------------------------------- inline
// Recursive so that nested emphasis -- **bold with *italic* inside** -- is
// parsed rather than emitted as literal asterisks.
function inline(md, base = {}) {
  const runs = [];
  const re = /(\*\*[\s\S]+?\*\*|\*[^*\n]+?\*|`[^`]+?`|<sup>[^<]*<\/sup>|<sub>[^<]*<\/sub>)/;
  let rest = md;
  while (rest.length) {
    const m = re.exec(rest);
    if (!m) { pushPlain(runs, rest, base); break; }
    pushPlain(runs, rest.slice(0, m.index), base);
    const tok = m[0];
    if (tok.startsWith('**')) {
      runs.push(...inline(tok.slice(2, -2), { ...base, bold: true }));
    } else if (tok.startsWith('<sup>')) {
      runs.push(new TextRun({ text: tok.replace(/<\/?sup>/g, ''), font: FONT,
                              ...base, superScript: true }));
    } else if (tok.startsWith('<sub>')) {
      runs.push(new TextRun({ text: tok.replace(/<\/?sub>/g, ''), font: FONT,
                              ...base, subScript: true }));
    } else if (tok.startsWith('`')) {
      runs.push(new TextRun({ text: tok.slice(1, -1), font: 'Courier New',
                              size: base.size ? base.size - 2 : 18, ...base,
                              font: 'Courier New' }));
    } else {
      runs.push(...inline(tok.slice(1, -1), { ...base, italics: true }));
    }
    rest = rest.slice(m.index + tok.length);
  }
  return runs.length ? runs : [new TextRun({ text: '', font: FONT, ...base })];
}

function pushPlain(runs, text, base) {
  if (text === '') return;
  runs.push(new TextRun({ text: clean(text), font: FONT, ...base }));
}

function clean(s) {
  return s.replace(/\\\*/g, '*').replace(/\\_/g, '_');
}

// ---------------------------------------------------------------- tables
function parseTable(lines) {
  const rows = lines
    .filter(l => l.trim().startsWith('|'))
    .map(l => l.trim().replace(/^\||\|$/g, '').split('|').map(c => c.trim()));
  if (rows.length < 2) return null;
  const sep = rows[1].every(c => /^:?-{2,}:?$/.test(c));
  const header = rows[0];
  const body = sep ? rows.slice(2) : rows.slice(1);
  const n = header.length;
  // Proportional to the widest cell in each column, clamped so no column is
  // squeezed small enough to break words. Equal widths wrap "Schizophrenia".
  const raw = Array.from({ length: n }, (_, i) => {
    const cells = [header[i] || '', ...body.map(r => r[i] || '')];
    const longest = Math.max(...cells.map(c => {
      const words = c.replace(/[*`]/g, '').split(/\s+/);
      return Math.max(c.replace(/[*`]/g, '').length * 0.55,
                      ...words.map(w => w.length));
    }));
    return Math.max(6, longest);
  });
  const tot = raw.reduce((a, b) => a + b, 0);
  const widths = raw.map(r => Math.floor(TEXT_W * r / tot));
  widths[n - 1] = TEXT_W - widths.slice(0, -1).reduce((a, b) => a + b, 0);

  const mkCell = (txt, isHeader) => new TableCell({
    width: { size: widths[0], type: WidthType.DXA },
    shading: isHeader
      ? { type: ShadingType.CLEAR, fill: 'EFEFEF' } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({
      spacing: { line: 240, before: 0, after: 0 },
      children: inline(txt, { size: 18, bold: isHeader || undefined }),
    })],
  });

  const mkRow = (cells, isHeader) => new TableRow({
    tableHeader: isHeader,
    children: cells.slice(0, n).concat(Array(Math.max(0, n - cells.length)).fill(''))
      .map((c, i) => {
        const cell = mkCell(c, isHeader);
        cell.root.find(x => x && x.rootKey === 'w:tcPr');
        return new TableCell({
          width: { size: widths[i], type: WidthType.DXA },
          shading: isHeader ? { type: ShadingType.CLEAR, fill: 'EFEFEF' } : undefined,
          margins: { top: 60, bottom: 60, left: 100, right: 100 },
          children: [new Paragraph({
            spacing: { line: 240, before: 0, after: 0 },
            children: inline(c, { size: 18, bold: isHeader || undefined }),
          })],
        });
      }),
  });

  return new Table({
    columnWidths: widths,
    width: { size: TEXT_W, type: WidthType.DXA },
    rows: [mkRow(header, true), ...body.map(r => mkRow(r, false))],
  });
}

// ---------------------------------------------------------------- convert
function convert(md, opts = {}) {
  const spacing = opts.double
    ? { line: 480, before: 0, after: 0 }
    : { line: 276, before: 0, after: 120 };
  const out = [];
  const lines = md.split('\n');
  let i = 0;

  while (i < lines.length) {
    const L = lines[i];

    // fenced code block
    if (L.trim().startsWith('```')) {
      i++;
      const buf = [];
      while (i < lines.length && !lines[i].trim().startsWith('```')) buf.push(lines[i++]);
      i++;
      buf.forEach(b => out.push(new Paragraph({
        spacing: { line: 240, before: 0, after: 0 },
        children: [new TextRun({ text: b || ' ', font: 'Courier New', size: 16 })],
      })));
      out.push(new Paragraph({ text: '', spacing: { after: 120 } }));
      continue;
    }

    // table
    if (L.trim().startsWith('|')) {
      const buf = [];
      while (i < lines.length && lines[i].trim().startsWith('|')) buf.push(lines[i++]);
      const t = parseTable(buf);
      if (t) { out.push(t); out.push(new Paragraph({ text: '', spacing: { after: 160 } })); }
      continue;
    }

    // horizontal rule
    if (/^-{3,}$/.test(L.trim())) {
      out.push(new Paragraph({
        spacing: { before: 120, after: 120 },
        border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: 'BBBBBB' } },
        children: [new TextRun({ text: '', font: FONT })],
      }));
      i++; continue;
    }

    // headings
    const h = L.match(/^(#{1,4})\s+(.*)$/);
    if (h) {
      const lvl = h[1].length;
      const map = { 1: HeadingLevel.HEADING_1, 2: HeadingLevel.HEADING_2,
                    3: HeadingLevel.HEADING_3, 4: HeadingLevel.HEADING_4 };
      out.push(new Paragraph({
        heading: map[lvl],
        spacing: { before: lvl <= 2 ? 300 : 220, after: 120, line: 276 },
        keepNext: true,
        children: inline(h[2], { bold: true, size: lvl === 1 ? 28 : lvl === 2 ? 24 : 22 }),
      }));
      i++; continue;
    }

    // numbered entry (reference list) -- one paragraph each, hanging indent
    if (/^\d+\.\s+/.test(L.trim())) {
      const buf2 = [L.trim()];
      i++;
      while (i < lines.length && lines[i].trim() !== ''
             && !/^\d+\.\s/.test(lines[i].trim())
             && !/^[#>|`-]/.test(lines[i].trim())) buf2.push(lines[i++].trim());
      out.push(new Paragraph({
        spacing: { line: opts.double ? 360 : 240, before: 0, after: 80 },
        indent: { left: 420, hanging: 420 },
        children: inline(buf2.join(' ')),
      }));
      continue;
    }

    // bullet list
    if (/^[-*]\s+/.test(L.trim())) {
      out.push(new Paragraph({
        bullet: { level: 0 },
        spacing,
        children: inline(L.trim().replace(/^[-*]\s+/, '')),
      }));
      i++; continue;
    }

    // blockquote
    if (L.trim().startsWith('>')) {
      out.push(new Paragraph({
        spacing, indent: { left: 420 },
        children: inline(L.trim().replace(/^>\s?/, ''), { italics: true }),
      }));
      i++; continue;
    }

    // blank
    if (L.trim() === '') { i++; continue; }

    // paragraph (join wrapped lines)
    const buf = [L];
    i++;
    while (i < lines.length && lines[i].trim() !== ''
           && !/^[#>|`-]/.test(lines[i].trim())
           && !/^\d+\.\s/.test(lines[i].trim())
           && !/^[-*]\s/.test(lines[i].trim())) buf.push(lines[i++]);
    out.push(new Paragraph({
      spacing, alignment: AlignmentType.LEFT,
      children: inline(buf.join(' ')),
    }));
  }
  return out;
}

// ---------------------------------------------------------------- build
function build(children, opts) {
  return new Document({
    creator: 'Alper Karagol; Taner Karagol',
    title: opts.title,
    styles: {
      default: {
        document: { run: { font: FONT, size: 22 } },
        heading1: { run: { font: FONT, size: 28, bold: true, color: '000000' } },
        heading2: { run: { font: FONT, size: 24, bold: true, color: '000000' } },
        heading3: { run: { font: FONT, size: 22, bold: true, color: '000000' } },
        heading4: { run: { font: FONT, size: 22, bold: true, italics: true, color: '000000' } },
      },
    },
    sections: [{
      properties: {
        page: {
          size: { width: PAGE_W, height: PAGE_H },
          margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN },
        },
        ...(opts.lineNumbers ? {
          lineNumbers: { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS, distance: 360 },
        } : {}),
      },
      footers: {
        default: new Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 18 })],
          })],
        }),
      },
      children,
    }],
  });
}

module.exports = { convert, build, inline, Paragraph, TextRun, FONT, AlignmentType };

if (require.main === module) {
  const [, , inFile, outFile, mode] = process.argv;
  const md = fs.readFileSync(inFile, 'utf8');
  const doc = build(convert(md, { double: mode !== 'supp' }), {
    title: inFile, lineNumbers: mode !== 'supp',
  });
  Packer.toBuffer(doc).then(b => { fs.writeFileSync(outFile, b); console.log('wrote', outFile); });
}
