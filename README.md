# C6-Multiplatform

**Cross-platform toolkit for synthetic biology and molecular biology operations**

[![npm version](https://img.shields.io/npm/v/c6-sim.svg?style=flat)](https://www.npmjs.com/package/c6-sim)  
[![jsDelivr hits](https://data.jsdelivr.com/v1/package/npm/c6-sim/badge)](https://www.jsdelivr.com/package/npm/c6-sim)

---

## 🚀 Overview

C6 is a comprehensive library for parsing, simulating, and analyzing DNA construction workflows. It supports a range of synthetic biology operations and automates reasoning about transcriptional structure and expression.

---

## 📦 Installation

```bash
npm install c6-sim
```

Or load via CDN:

```html
<script src="https://cdn.jsdelivr.net/npm/c6-sim/dist/c6-sim.min.js"></script>
```

---

## 🧪 Quickstart Example

```javascript
import C6 from 'c6-sim';

const cfText = `
PCR Fwd Rev on Template, Product
oligo Fwd ATGCGT
oligo Rev TACGCA
plasmid Template ATGCGTACGCATAGC
`;

const cf = C6.parseCF(cfText);
const result = C6.simCF(cf);
console.log(result);
```

---

## 🔧 Core Capabilities

### Construction File Simulation

C6 can simulate PCR, Gibson Assembly, Golden Gate, ligation, and transformation:

```js
const cfText = `
GoldenGate oF oR on Vector, Digest
Gibson Digest Insert, Product
Transform Product into Mach1 with Amp, Colony
`;

const cf = C6.parseCF(cfText);
const result = C6.simCF(cf);
```

---

### Feature Annotation

Automatically detect annotated elements in DNA sequences:

```js
const features = C6.annotateSequenceSmart("ATGAGTGAAGAGAGGAGAAATACTAGATGGCGTCT...");
console.log(features); // List of features with type, label, and position
```

---

### Transcriptional Unit Inference

Predict transcriptional structure from annotated features:

```js
const tus = C6.inferTranscriptionalUnits(features);
console.log(tus);
```

---

### Protein Expression Inference

Identify which CDSs are likely expressed:

```js
const expressed = C6.inferExpressedProteins(tus);
console.log(expressed);
```

---

### Utilities

Handy tools for formatting and data extraction:

```js
C6.merge("5'", "ATG", "CTG", "3'", ".");
C6.makeJSON([["name", "G00101"], ["seq", "ATGC"]]);
C6.field('{"name":"G00101"}', "name"); // → "G00101"
```

---

## 🧪 Running Tests

This project uses [Vitest](https://vitest.dev):

```bash
npm test
```

Test coverage includes core simulation, sequence utilities, oligo design, and gene annotation.

---

## 📚 Repository

[GitHub – UCB BioE Anderson Lab](https://github.com/UCB-BioE-Anderson-Lab/C6-Multiplatform)

---

## 🧬 License

MIT License © University of California, Berkeley
