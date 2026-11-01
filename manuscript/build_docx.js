#!/usr/bin/env node
/* Build the submission .docx files. */
const fs = require('fs');
const { Packer, Paragraph, TextRun, AlignmentType, PageBreak, ImageRun } = require('docx');
const path = require('path');
const ROOT = path.resolve(__dirname, '..');
const FIGDIR = process.env.FIGDIR || path.join(ROOT, 'results');
const OUTDIR = process.env.OUTDIR || path.join(ROOT, 'results');
const M = require('./md2docx.js');
const FONT = M.FONT;

const p = (text, o = {}) => new Paragraph({
  spacing: { line: o.line || 276, before: o.before || 0, after: o.after ?? 120 },
  alignment: o.align || AlignmentType.LEFT,
  children: o.runs || [new TextRun({ text, font: FONT, size: o.size || 22,
                                     bold: o.bold, italics: o.italics })],
});

// ===================================================== MANUSCRIPT
let md = fs.readFileSync(path.join(__dirname, 'Manuscript_submission.md'), 'utf8');
md = md.split('# ⚠ PRE-SUBMISSION NOTES')[0];        // internal notes: never ship

const titleMd = md.match(/^#\s+(.*)$/m)[1];
const title = titleMd.replace(/\*/g, '');
const body = md.slice(md.indexOf('## Abstract'));

const titlePage = [
  p(null, { after: 300, line: 360, runs: M.inline(titleMd, { bold: true, size: 28 }) }),
  p('Running title: Anatomy of a post-GWAS candidate-gene finding',
    { italics: true, after: 300 }),

  new Paragraph({
    spacing: { after: 200, line: 276 },
    children: [
      new TextRun({ text: 'Alper Karagöl', font: FONT, size: 22 }),
      new TextRun({ text: '1,2,3,†,*', font: FONT, size: 22, superScript: true }),
      new TextRun({ text: ' and Taner Karagöl', font: FONT, size: 22 }),
      new TextRun({ text: '1,4,†,*', font: FONT, size: 22, superScript: true }),
    ],
  }),

  ...[
    ['1', 'Istanbul University, Istanbul Medical Faculty, Istanbul, Turkey'],
    ['2', 'Present address: Broad Institute of MIT and Harvard, Cambridge, MA, USA'],
    ['3', 'Present address: Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA'],
    ['4', 'Present address: McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA'],
  ].map(([n, t]) => new Paragraph({
    spacing: { after: 40, line: 240 },
    children: [
      new TextRun({ text: n, font: FONT, size: 20, superScript: true }),
      new TextRun({ text: ' ' + t, font: FONT, size: 20 }),
    ],
  })),

  new Paragraph({
    spacing: { before: 160, after: 40, line: 240 },
    children: [
      new TextRun({ text: '†', font: FONT, size: 20, superScript: true }),
      new TextRun({ text: ' These authors contributed equally.', font: FONT, size: 20 }),
    ],
  }),
  new Paragraph({
    spacing: { after: 240, line: 240 },
    children: [
      new TextRun({ text: '*', font: FONT, size: 20, superScript: true }),
      new TextRun({ text: ' Joint corresponding authors.', font: FONT, size: 20 }),
    ],
  }),

  p('Correspondence', { bold: true, after: 60, size: 20 }),
  p('Alper Karagöl — akaragol@mit.edu — ORCID 0009-0001-7864-0732',
    { size: 20, after: 40, line: 240 }),
  p('Taner Karagöl — karagolt@mit.edu — ORCID 0009-0005-1011-7661',
    { size: 20, after: 280, line: 240 }),

  p('Word count: Abstract 248; main text ~4,400', { size: 20, after: 40, line: 240 }),
  p('Display items: 4 tables, 3 figures', { size: 20, after: 40, line: 240 }),
  p('Supplementary: Tables S1–S12, Supplementary Methods S1–S2, Supplementary Software',
    { size: 20, after: 40, line: 240 }),
  p('Keywords: SLC1A2; EAAT2; major depressive disorder; candidate gene; excess significance; meta-research; astrocyte',
    { size: 20, after: 0, line: 240 }),

  new Paragraph({ children: [new PageBreak()] }),
];

const figPages = [];
[['Figure 1', path.join(FIGDIR, 'Figure1.png'), 560],
 ['Figure 2', path.join(FIGDIR, 'Figure2.png'), 560],
 ['Figure 3', path.join(FIGDIR, 'Figure3.png'), 580]].forEach(([label, file, w], i) => {
  if (!fs.existsSync(file)) return;
  const dim = require('child_process')
    .execSync(`python3 -c "from PIL import Image;im=Image.open('${file}');print(im.width,im.height)"`)
    .toString().trim().split(' ').map(Number);
  figPages.push(new Paragraph({ children: [new PageBreak()] }));
  figPages.push(p(label, { bold: true, after: 160 }));
  figPages.push(new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 120 },
    children: [new ImageRun({
      type: 'png', data: fs.readFileSync(file),
      transformation: { width: w, height: Math.round(w * dim[1] / dim[0]) },
    })],
  }));
});

const doc = M.build([...titlePage, ...M.convert(body, { double: true }), ...figPages],
                    { title, lineNumbers: true });
Packer.toBuffer(doc).then(b => {
  fs.writeFileSync(path.join(OUTDIR, 'Manuscript_Karagol_SLC1A2.docx'), b);
  console.log('wrote', path.join(OUTDIR, 'Manuscript_Karagol_SLC1A2.docx'));
});

// ===================================================== SUPPLEMENTARY
const smd = fs.readFileSync(path.join(__dirname, 'Supplementary_Information.md'), 'utf8');
const sHead = [
  p('Supplementary Information', { size: 28, bold: true, after: 200, line: 360 }),
  p(null, { after: 160, line: 300, runs: M.inline(titleMd, { bold: true, size: 22 }) }),
  p('Alper Karagöl and Taner Karagöl', { after: 300 }),
];
const sBody = smd.slice(smd.indexOf('## Contents'));
const sdoc = M.build([...sHead, ...M.convert(sBody, { double: false })],
                     { title: 'Supplementary Information', lineNumbers: false });
Packer.toBuffer(sdoc).then(b => {
  fs.writeFileSync(path.join(OUTDIR, 'Supplementary_Karagol_SLC1A2.docx'), b);
  console.log('wrote', path.join(OUTDIR, 'Supplementary_Karagol_SLC1A2.docx'));
});
